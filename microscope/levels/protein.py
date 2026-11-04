"""
Nivel 1 (proteínas): la proteína plegada, átomo a átomo (unidad: nanómetro).

Igual que la doble hélice: cada átomo se dibuja como su nube de electrones y la forma de la
proteína sale de la suma de las nubes. Además:

  Modo didáctico
    - núcleos marcados y enlaces covalentes con el color de cada elemento
    - la "traza" del esqueleto: una línea que une los carbonos α (uno por aminoácido),
      coloreada por el tipo de aminoácido (hidrófobo, polar, con carga + o -...): así se
      ve cómo se pliega la cadena
  Modo realista
    - sin líneas ni colores inventados; si la estructura es una predicción, las nubes se
      tiñen con su confianza (pLDDT, los colores de AlphaFold DB)

Los aminoácidos resaltados (mutaciones, al comparar con otra proteína) llevan una esfera
magenta en su carbono α en los dos modos.

Al acercarse a un átomo se entra en el mismo nivel de átomo que en el ADN (atom.py): la
molécula de proteína tiene la misma forma que la de ADN (dnainvent/scene/protein.py).
"""
from typing import Callable, Optional

import numpy as np

from atomicworld.scene import orbitals
from dnainvent.constants import AMINO_ACID_CLASSES
from dnainvent.protein import confidence
from dnainvent.scene.elements import atom_type
from microscope.levels.base import Child, Focus, Level, Target, aimed
from microscope.render.points import LineSet, PointCloud, PointStyle

REALISTIC_COLOR = np.array([0.55, 0.70, 1.00])
TARGET_POINTS = 350_000
ENTER_FRACTION = 0.45
VIB_AMP_NM = 0.006
OVERVIEW_DISTANCE = 3.0           # nm: más lejos que esto, el HUD describe la proteína entera
CA_TRACE_MAX = 0.45               # nm: dos Cα consecutivos de la cadena están a ~0,38 nm
TUBE_STEPS = 6                    # puntos por tramo del "tubo" del esqueleto
TUBE_SIZE = 0.07                  # nm
HIGHLIGHT_COLOR = (1.0, 0.25, 0.85)
CLASS_COLORS = {
    "hydrophobic": (0.95, 0.78, 0.30),
    "special": (0.80, 0.80, 0.80),
    "polar": (0.35, 0.90, 0.50),
    "positive": (0.35, 0.55, 1.00),
    "negative": (1.00, 0.35, 0.35),
}


def _hex(color: str) -> tuple:
    return tuple(int(color[i:i + 2], 16) / 255 for i in (1, 3, 5))


def residue_letter(mol, i: int) -> str:
    return mol.sequence[int(mol.pair_index[i]) - mol.start]


def protein_atom_context(tr, mol, i: int) -> str:
    """'esqueleto · Glutamato (Glu) 7'"""
    return tr("microscope.protein.where", part=tr(f"microscope.protein.parts.{mol.part[i]}"),
              residue=tr.amino_acid(residue_letter(mol, i)), number=int(mol.pair_index[i]) + 1)


class ProteinLevel(Level):
    depth = 0
    unit_m = 1e-9

    def __init__(self, molecule, tr, label: Callable[[], str] | str = "",
                 highlight: dict[int, str] | None = None, tint=None):
        """label: nombre para el HUD ("Hemoglobin subunit beta · P68871").
        highlight: aminoácido (índice en la cadena) -> texto ("Glu7Val").
        tint: color fijo para todas las nubes (al superponer dos proteínas)."""
        super().__init__(tr)
        self.mol = molecule
        self.label = label if callable(label) else (lambda text=label: text)
        self.highlight = highlight or {}
        self.tint = tint
        self.vdw = molecule.radius_vdw()
        self.clouds, self.lines = [], []
        self.tube = None
        self.ca = self._alpha_carbons()

    @property
    def name(self):
        return self.tr("microscope.levels.protein")

    def _alpha_carbons(self) -> dict[int, int]:
        """Aminoácido -> índice del átomo de su Cα."""
        return {int(self.mol.pair_index[i]): i for i, n in enumerate(self.mol.names) if n == "CA"}

    # ------------------------------------------------------------ construcción
    def _cloud_colors(self, didactic: bool) -> np.ndarray:
        mol = self.mol
        if self.tint is not None:
            return np.tile(self.tint, (mol.n_atoms, 1))
        if didactic:
            return np.array([atom_type(e).color for e in mol.elements])
        if mol.confidence is not None:
            return np.array([_hex(confidence.color(c)) for c in mol.confidence])
        return np.tile(REALISTIC_COLOR, (mol.n_atoms, 1))

    def build(self, renderer, didactic: bool):
        mol = self.mol
        rng = np.random.default_rng(1)
        types = {s: atom_type(s) for s in set(mol.elements)}
        total_e = sum(types[e].electrons for e in mol.elements)
        ppe = int(np.clip(TARGET_POINTS / total_e, 4, 60))
        cache = {s: orbitals.sample_atom(t.number, t.shells, ppe * 40, didactic, rng)[0] / 1000.0
                 for s, t in types.items()}
        colors = self._cloud_colors(didactic)
        alpha = 0.10 if self.tint is not None else 0.16

        pos, col, size, seed = [], [], [], []
        for i, sym in enumerate(mol.elements):
            t, src = types[sym], cache[sym]
            n = t.electrons * ppe
            pos.append(mol.positions[i] + src[rng.integers(0, len(src), n)])
            col.append(np.tile([*colors[i], alpha], (n, 1)))
            size.append(np.full(n, 0.010))
            seed.append(np.full(n, float(i)))
        self.clouds.append(PointCloud(renderer, np.vstack(pos), np.vstack(col),
                                      np.concatenate(size), np.concatenate(seed)))

        if didactic and self.tint is None:
            cpk = np.array([[*types[e].color, 0.9] for e in mol.elements])
            self.clouds.append(PointCloud(renderer, mol.positions, cpk, self.vdw * 0.22,
                                          np.arange(mol.n_atoms, dtype="f4")))
            a, b = mol.positions[mol.bonds[:, 0]], mol.positions[mol.bonds[:, 1]]
            ca_ = np.array([[*types[mol.elements[i]].color, 0.45] for i in mol.bonds[:, 0]])
            cb_ = np.array([[*types[mol.elements[i]].color, 0.45] for i in mol.bonds[:, 1]])
            self.lines.append((LineSet(renderer, a, b, ca_, cb_), 0.0))
        if didactic or self.tint is not None:
            self._trace(renderer)
        self._highlights(renderer)

    def _trace(self, renderer):
        """Línea del esqueleto: Cα de cada aminoácido con el siguiente."""
        mol, a, b, ca_, cb_ = self.mol, [], [], [], []
        for k in sorted(self.ca):
            i, j = self.ca[k], self.ca.get(k + 1)
            if j is None or np.linalg.norm(mol.positions[i] - mol.positions[j]) > CA_TRACE_MAX:
                continue
            a.append(mol.positions[i])
            b.append(mol.positions[j])
            ca_.append([*self._trace_color(k), 0.95])
            cb_.append([*self._trace_color(k + 1), 0.95])
        if not a:
            self.tube = None
            return
        a, b, ca_, cb_ = np.array(a), np.array(b), np.array(ca_), np.array(cb_)
        self.lines.append((LineSet(renderer, a, b, ca_, cb_), 0.0))
        # El mismo esqueleto como un tubo de puntos: de lejos las líneas (1 píxel) se pierden
        # entre las nubes; los puntos tienen un tamaño mínimo en pantalla y siempre se ven.
        t = np.linspace(0, 1, TUBE_STEPS, endpoint=False)[None, :, None]
        pos = (a[:, None, :] * (1 - t) + b[:, None, :] * t).reshape(-1, 3)
        col = (ca_[:, None, :] * (1 - t) + cb_[:, None, :] * t).reshape(-1, 4)
        col[:, 3] = 0.85
        self.tube = PointCloud(renderer, pos, col, np.full(len(pos), TUBE_SIZE),
                               np.arange(len(pos), dtype="f4"))

    def _trace_color(self, k: int):
        if self.tint is not None:
            return self.tint
        if k in self.highlight:
            return HIGHLIGHT_COLOR
        letter = self.mol.sequence[k - self.mol.start]
        return CLASS_COLORS[AMINO_ACID_CLASSES.get(letter, "special")]

    def _highlights(self, renderer):
        marks = [self.ca[k] for k in self.highlight if k in self.ca]
        if marks:
            p = self.mol.positions[marks]
            self.clouds.append(PointCloud(renderer, p, np.tile([*HIGHLIGHT_COLOR, 1.0], (len(p), 1)),
                                          np.full(len(p), 0.18), np.arange(len(p), dtype="f4")))
            self.highlight_cloud = len(self.clouds) - 1
        else:
            self.highlight_cloud = None

    def release(self):
        for c in self.clouds + ([self.tube] if self.tube else []):
            c.release()
        for line, _ in self.lines:
            line.release()
        self.clouds, self.lines, self.tube = [], [], None

    # ------------------------------------------------------------ dibujo
    def render(self, renderer, cam, didactic, fade):
        style = PointStyle(vib_amp=VIB_AMP_NM, min_px=1.2, max_px=400, softness=1.0)
        renderer.points(self.clouds[0], cam, style, fade)
        if self.tube is not None:
            renderer.points(self.tube, cam, PointStyle(min_px=2.5, max_px=60, softness=0.45),
                            fade)
        for k, cloud in enumerate(self.clouds[1:], start=1):
            if k == self.highlight_cloud:
                renderer.points(cloud, cam, PointStyle(min_px=6.0, max_px=200, softness=0.5),
                                fade)
            else:
                renderer.points(cloud, cam, PointStyle(vib_amp=VIB_AMP_NM, min_px=2.0,
                                                       max_px=200, softness=0.35), fade * 0.9)
        for lines, dash in self.lines:
            renderer.lines(lines, cam, fade, dash)

    # ------------------------------------------------------------ navegación
    def _dist(self, pos):
        return np.linalg.norm(self.mol.positions - pos, axis=1)

    def surface_distance(self, pos):
        return float(self._dist(pos).min())

    def child_to_enter(self, pos, only=None):
        d = self._dist(pos)
        i = only[1] if only is not None else int(np.argmin(d))
        if d[i] < ENTER_FRACTION * self.vdw[i]:
            from microscope.levels.atom import AtomLevel
            return Child(key=("atom", i), center=self.mol.positions[i].astype(float),
                         scale=1000.0, create=lambda: AtomLevel(self.mol, i, self.tr))
        return None

    def dive_target(self, pos, forward, didactic):
        from microscope.levels.helix import atom_title
        mol = self.mol
        i = aimed(mol.positions, pos, forward, 5)
        if i is None:
            v = mol.positions - pos
            d = np.linalg.norm(v, axis=1)
            cos = (v @ forward) / np.maximum(d, 1e-12)
            heavy = np.array([e != "H" for e in mol.elements])
            i = int(np.argmin(np.where(heavy, d * (2.0 - cos), np.inf)))
        p = mol.positions[i].astype(float)
        tr = self.tr
        return Target(lambda: tr("microscope.autopilot.atom_target_residue",
                                 atom=atom_title(tr, mol, i),
                                 residue=tr.amino_acid(residue_letter(mol, i)),
                                 number=int(mol.pair_index[i]) + 1),
                      lambda: p, key=("atom", i))

    def focus(self, pos, forward, didactic):
        if self._dist(pos).min() > OVERVIEW_DISTANCE:       # de lejos: la proteína entera
            return self.overview(didactic)
        i = aimed(self.mol.positions, pos, forward, 5)
        return None if i is None else self.atom_focus(i, float(self._dist(pos)[i]), didactic)

    def overview(self, didactic) -> Focus:
        mol, tr = self.mol, self.tr
        if mol.start == 0 and len(mol.sequence) == mol.length:
            detail = tr("microscope.protein.residues", residues=len(mol.sequence),
                        atoms=mol.n_atoms)
        else:
            detail = tr("microscope.protein.window", first=mol.start + 1, last=mol.end,
                        length=mol.length, atoms=mol.n_atoms)
        if didactic:
            expl = tr.lines("microscope.protein.explanation")
        elif mol.confidence is not None:
            expl = tr.lines("microscope.protein.explanation_confidence")
        else:
            expl = []
        if self.highlight:
            expl = expl + [tr("microscope.protein.highlighted",
                              list=", ".join(list(self.highlight.values())[:8]))]
        return Focus(tr("microscope.protein.title", name=self.label()), detail, expl)

    def atom_focus(self, i, dist, didactic) -> Focus:
        from microscope.levels.helix import atom_title
        mol, tr = self.mol, self.tr
        k = int(mol.pair_index[i])
        expl = []
        if didactic:
            kind = AMINO_ACID_CLASSES.get(residue_letter(mol, i), "special")
            expl = (tr.lines(f"microscope.protein.explanation_{mol.part[i]}")
                    + [tr(f"microscope.protein.classes.{kind}")])
        if mol.confidence is not None:
            band = tr(confidence.i18n_key(float(mol.confidence[i])))
            expl.append(tr("microscope.protein.atom_confidence",
                           plddt=tr.decimal(float(mol.confidence[i]), 1),
                           band=band[:1].lower() + band[1:]))
        if k in self.highlight:
            expl.append(tr("microscope.protein.atom_highlighted", label=self.highlight[k]))
        return Focus(atom_title(tr, mol, i), protein_atom_context(tr, mol, i), expl,
                     pair_index=k, distance=dist)

    def time_info(self):
        return self.tr("microscope.helix.time")


class OverlayLevel(Level):
    """Las dos proteínas en la misma escena, cada una con su color (ya superpuestas: comparten
    coordenadas). Para entrar en un átomo se elige el más cercano de cualquiera de las dos."""
    depth = 0
    unit_m = 1e-9
    TINTS = ((0.45, 0.85, 1.00), (1.00, 0.62, 0.25))        # natural: cian; la tuya: naranja

    def __init__(self, levels: list[ProteinLevel], tr):
        super().__init__(tr)
        self.levels = levels
        for level, tint in zip(levels, self.TINTS):
            level.tint = tint

    @property
    def name(self):
        return self.tr("microscope.levels.overlay")

    @property
    def mol(self):
        return self.levels[0].mol

    def build(self, renderer, didactic):
        for level in self.levels:
            level.build(renderer, didactic)

    def release(self):
        for level in self.levels:
            level.release()

    def render(self, renderer, cam, didactic, fade):
        for level in self.levels:
            level.render(renderer, cam, didactic, fade)

    def surface_distance(self, pos):
        return min(level.surface_distance(pos) for level in self.levels)

    def _nearest(self, pos) -> int:
        return int(np.argmin([level.surface_distance(pos) for level in self.levels]))

    def child_to_enter(self, pos, only=None):
        if only is not None:
            side, inner = only[2], only[:2]
        else:
            side, inner = self._nearest(pos), None
        child = self.levels[side].child_to_enter(pos, inner)
        if child is None:
            return None
        return Child(key=(*child.key, side), center=child.center, scale=child.scale,
                     create=child.create)

    def dive_target(self, pos, forward, didactic):
        targets = [level.dive_target(pos, forward, didactic) for level in self.levels]
        side = int(np.argmin([np.linalg.norm(t.point() - pos) for t in targets]))
        t = targets[side]
        return Target(t.label, t.point, key=(*t.key, side), arrive_radius=t.arrive_radius)

    def focus(self, pos, forward, didactic):
        if self.surface_distance(pos) > OVERVIEW_DISTANCE:
            return self.overview(didactic)
        best: Optional[tuple] = None
        for side, level in enumerate(self.levels):
            i = aimed(level.mol.positions, pos, forward, 5)
            if i is not None:
                d = float(np.linalg.norm(level.mol.positions[i] - pos))
                if best is None or d < best[0]:
                    best = (d, side, i)
        if best is None:
            return None
        d, side, i = best
        focus = self.levels[side].atom_focus(i, d, didactic)
        who = self.tr(f"microscope.protein.sides.{('natural', 'yours')[side]}")
        focus.detail = f"{who} · {focus.detail}"
        return focus

    def overview(self, didactic) -> Focus:
        tr = self.tr
        lines = tr.lines("microscope.protein.overlay_explanation")
        highlight = self.levels[0].highlight
        if highlight:
            lines = lines + [tr("microscope.protein.highlighted",
                                list=", ".join(list(highlight.values())[:8]))]
        return Focus(tr("microscope.protein.overlay_title", natural=self.levels[0].label(),
                        yours=self.levels[1].label()), "", lines)

    def time_info(self):
        return self.tr("microscope.helix.time")
