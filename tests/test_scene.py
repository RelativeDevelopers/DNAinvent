"""Tests de la molécula 3D (dnainvent.scene). No necesitan tarjeta gráfica."""

import json

import numpy as np
import pytest

pytest.importorskip("numpy")
pytest.importorskip("atomicworld")

from dnainvent import InvalidParameterError, InvalidSequenceError, helix3d
from dnainvent.scene import describe, elements
from dnainvent.scene.limits import MAX_WEB_PAIRS, RECOMMENDED_PAIRS
from dnainvent.scene.molecule import build_molecule, helical_step

EXAMPLE = "ATGCGTACGTTAGC"


@pytest.fixture(scope="module")
def mol():
    return build_molecule(EXAMPLE)


def test_hydrogen_bonds_come_from_geometry(mol):
    # no se escriben a mano: A–T debe dar 2 y G–C 3
    per_pair = np.bincount(mol.pair_index[mol.hbonds[:, 0]], minlength=len(EXAMPLE))
    expected = [2 if b in "AT" else 3 for b in EXAMPLE]
    assert per_pair.tolist() == expected


def test_helical_step_is_the_one_of_helix3d():
    R, t = helical_step()
    twist = np.degrees(np.arccos((np.trace(R) - 1) / 2))
    assert twist == pytest.approx(360 / helix3d.BASE_PAIRS_PER_TURN, abs=1e-6)   # 34,29°
    assert abs(t[2]) == pytest.approx(helix3d.RISE)                               # 3,38 Å


def test_ten_and_a_half_pairs_per_turn():
    # el par 21 está girado 2 vueltas exactas respecto al par 0
    m = build_molecule("A" * 22)
    p = lambda k: m.positions[(m.pair_index == k) & (np.array(m.names) == "P") & (m.strand == 1)][0]
    a, b = p(0), p(21)
    angle = np.degrees(np.arctan2(a[1], a[0]) - np.arctan2(b[1], b[0])) % 360
    assert min(angle, 360 - angle) < 1.0
    assert abs(b[2] - a[2]) == pytest.approx(21 * helix3d.RISE / 10, rel=1e-3)   # en nm


def test_backbone_bonds_between_pairs_have_a_real_length(mol):
    names = np.array(mol.names)
    between = mol.pair_index[mol.bonds[:, 0]] != mol.pair_index[mol.bonds[:, 1]]
    links = mol.bonds[between]
    assert len(links) == 2 * (len(EXAMPLE) - 1)        # un O3'–P por hebra y paso
    assert set(map(frozenset, names[links])) == {frozenset({"O3'", "P"})}
    lengths = np.linalg.norm(mol.positions[links[:, 0]] - mol.positions[links[:, 1]], axis=1)
    assert np.all((lengths > 0.150) & (lengths < 0.170))   # P–O ≈ 1,6 Å


def test_every_covalent_bond_is_short(mol):
    lengths = np.linalg.norm(mol.positions[mol.bonds[:, 0]] - mol.positions[mol.bonds[:, 1]], axis=1)
    assert lengths.max() < 0.175 and lengths.min() > 0.09


def test_molecule_is_centered_and_has_the_right_size(mol):
    assert np.allclose(mol.positions.mean(axis=0), 0, atol=1e-4)
    extent = np.ptp(mol.positions, axis=0)
    assert 1.8 < extent[0] < 2.4 and 1.8 < extent[1] < 2.4        # ~2 nm de diámetro
    assert extent[2] == pytest.approx(13 * 0.338 + 0.6, abs=0.3)


def test_only_dna_elements(mol):
    assert set(mol.elements) == set(elements.DNA_ELEMENTS)


def test_window_keeps_real_pair_numbers():
    dna = "ACGT" * 50
    m = build_molecule(dna, start=120, count=10)
    assert m.sequence == dna[120:130]
    assert (m.start, m.end, m.length) == (120, 130, 200)
    assert m.pair_index.min() == 120 and m.pair_index.max() == 129


def test_window_is_cut_at_the_end():
    m = build_molecule("ACGTACGT", start=6, count=10)
    assert m.sequence == "GT"


def test_window_errors():
    with pytest.raises(InvalidParameterError):
        build_molecule("ACGT", start=4)
    with pytest.raises(InvalidParameterError):
        build_molecule("ACGT", count=0)
    with pytest.raises(InvalidSequenceError):
        build_molecule("ACGU")


def test_lowercase_and_spaces_are_cleaned():
    assert build_molecule("acg t").sequence == "ACGT"


def test_element_data_comes_from_atomicworld():
    c = elements.atom_type("C")
    assert (c.number, c.isotope.a, c.isotope.neutrons) == (6, 12, 6)
    assert elements.atom_type("P").isotope.a == 31
    assert c.vdw_nm == pytest.approx(0.170)
    assert all(0 <= v <= 1 for v in elements.atom_type("N").color)


def test_limits():
    assert RECOMMENDED_PAIRS <= MAX_WEB_PAIRS


def test_describe_structure_is_json_and_consistent():
    d = describe.describe_structure(EXAMPLE)
    text = json.dumps(d)
    assert len(text) < 200_000
    atoms = d["atoms"]
    n = atoms["count"]
    assert len(atoms["positions_nm"]) == 3 * n
    assert len(atoms["elements"]) == len(atoms["strand"]) == len(atoms["part"]) == n
    assert len(atoms["names"]) == len(atoms["pair"]) == n
    assert max(d["bonds"]) < n and len(d["bonds"]) % 2 == 0
    assert len(d["hbonds"]) == 2 * sum(2 if b in "AT" else 3 for b in EXAMPLE)
    assert set(d["elements"]) == set(atoms["elements"])
    carbon = d["elements"]["C"]
    assert carbon["element"]["number"] == 6 and carbon["isotope"]["a"] == 12
    assert len(carbon["nucleus"]["centers_fm"]) == 12


def test_describe_window():
    d = describe.describe_structure("ACGT" * 40, start=100, count=20)
    assert d["start"] == 100 and d["length"] == 160 and len(d["sequence"]) == 20
    assert min(d["atoms"]["pair"]) == 100
