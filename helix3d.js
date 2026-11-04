/* Visor 3D de la doble hélice.
 *
 * Python (dnainvent/helix3d.py) calcula dónde está cada par de bases.
 * Aquí solo se gira la hélice alrededor de su eje (la misma matriz que
 * rotate_z) y se proyecta en un <canvas> 2D, 60 veces por segundo.
 *
 * Cómo se consigue el efecto 3D sin ninguna librería:
 *   - Perspectiva: lo que está cerca se dibuja un poco más grande.
 *   - Orden de pintado: se dibuja de atrás hacia delante (algoritmo del pintor),
 *     así lo que está cerca tapa lo que está detrás.
 *   - Profundidad: lo que está detrás se funde con el color del fondo.
 */
"use strict";

window.HelixViewer = (() => {
  const SPIN_SPEED = 0.7;        // radianes por segundo
  const MAX_SCALE = 14;          // píxeles por ángstrom como máximo (secuencias muy cortas)
  const SIDE_MARGIN = 34;        // espacio para las etiquetas 5' y 3'
  const CAMERA = 12;             // distancia de la cámara, en radios de la hélice
  // Con perspectiva, lo más cercano se ve hasta un 9 % más grande: hay que reservarle sitio.
  const PERSPECTIVE = CAMERA / (CAMERA - 1);

  const hexToRgb = (hex) => {
    const h = hex.replace("#", "");
    return [0, 2, 4].map((i) => parseInt(h.slice(i, i + 2), 16));
  };
  const mix = (a, b, amount) => `rgb(${a.map((v, i) => Math.round(v + (b[i] - v) * amount)).join(",")})`;

  function readColors() {
    const style = getComputedStyle(document.documentElement);
    const color = (name) => hexToRgb(style.getPropertyValue(name).trim());
    return {
      A: color("--a"), C: color("--c"), G: color("--g"), T: color("--t"),
      backbone: color("--muted"), ink: color("--ink"), background: color("--surface"),
    };
  }

  /**
   * Monta el visor dentro de root.
   *   dna    secuencia validada (solo A, C, G, T)
   *   api    función para llamar a la API
   *   t      función de traducción
   *   view   estado que sobrevive a los redibujados: { angle, playing, start }
   * Devuelve { destroy } para parar la animación al salir de la herramienta.
   */
  function mount(root, { dna, api, t, view }) {
    const canvas = root.querySelector("canvas");
    const ctx = canvas.getContext("2d");
    const tip = root.querySelector(".helix-tip");
    const toggle = root.querySelector(".helix-toggle");
    const rangeLabel = root.querySelector(".helix-from");
    const range = root.querySelector(".helix-range");
    const caption = root.querySelector(".helix-window");
    const reduceMotion = matchMedia("(prefers-reduced-motion: reduce)");
    const darkMode = matchMedia("(prefers-color-scheme: dark)");

    if (view.playing === undefined) view.playing = !reduceMotion.matches;
    let colors = readColors();
    let data = null;
    let width = 0, height = 0, scale = 1;
    let rungs = [];            // puntos medios de los peldaños en pantalla, para el puntero
    let dragging = null;
    let visible = true;
    let dirty = true;
    let frame = 0, lastTime = 0, request = 0, debounce = 0;
    let alive = true;

    const windowSize = () => (root.clientWidth < 520 ? 10 : 24);

    // ---------- Datos ----------

    async function load() {
      const count = Math.min(windowSize(), dna.length);
      view.start = Math.max(0, Math.min(view.start ?? 0, dna.length - count));
      const id = ++request;
      try {
        const result = await api("/api/helix", { dna, start: view.start, count });
        if (!alive || id !== request) return;     // llegó tarde: ya se pidió otra ventana
        data = result;
      } catch {
        return;
      }
      const first = data.start + 1, last = data.start + data.pairs.length;
      caption.textContent = t("web.helix.window", { first, last, length: data.length });
      canvas.setAttribute("aria-label", t("web.helix.label", { first, last }));
      rangeLabel.hidden = data.length <= count;
      range.max = String(data.length - count + 1);
      range.value = String(data.start + 1);
      resize();
    }

    // ---------- Tamaño ----------

    function resize() {
      width = canvas.parentElement.clientWidth;
      const count = data ? data.pairs.length : windowSize();
      const rise = data?.rise ?? 3.38, radius = data?.radius ?? 10;
      const length = Math.max(count - 1, 1) * rise;
      scale = Math.min((width - 2 * SIDE_MARGIN) / (length * PERSPECTIVE + 2 * radius * 0.2),
                       MAX_SCALE, 320 / (2.4 * radius));
      height = Math.ceil(2 * radius * scale * PERSPECTIVE * 1.1 + 28);
      const ratio = window.devicePixelRatio || 1;
      canvas.style.height = `${height}px`;
      canvas.width = Math.round(width * ratio);
      canvas.height = Math.round(height * ratio);
      ctx.setTransform(ratio, 0, 0, ratio, 0, 0);
      dirty = true;
    }

    // ---------- Proyección ----------

    function projector() {
      const c = Math.cos(view.angle), s = Math.sin(view.angle);
      const pairs = data.pairs;
      const zMid = (pairs[0].strand1[2] + pairs[pairs.length - 1].strand1[2]) / 2;
      const camera = CAMERA * data.radius;   // perspectiva suave
      return ([x, y, z]) => {
        // Rotación alrededor del eje Z (igual que rotate_z en Python).
        const rx = c * x - s * y;
        const ry = s * x + c * y;
        // El eje Z va de izquierda a derecha; Y apunta hacia quien mira.
        const f = camera / (camera - ry);
        return { x: width / 2 + (z - zMid) * scale * f, y: height / 2 - rx * scale * f, depth: ry, f };
      };
    }

    /** Punto del esqueleto para una posición fraccionaria t (misma fórmula que backbone_point). */
    function backbone(tIndex, offset) {
      const angle = ((tIndex * data.twist + offset) * Math.PI) / 180;
      return [data.radius * Math.cos(angle), data.radius * Math.sin(angle), tIndex * data.rise];
    }

    // ---------- Dibujo ----------

    function draw() {
      ctx.clearRect(0, 0, width, height);
      if (!data || !data.pairs.length) return;
      const project = projector();
      const radius = data.radius;
      const shade = (rgb, depth) => mix(rgb, colors.background, 0.62 * (1 - (depth + radius) / (2 * radius)));
      const items = [];

      // Esqueletos: curvas suaves, 8 tramos entre cada par de bases.
      const first = data.pairs[0].index, last = data.pairs[data.pairs.length - 1].index;
      for (const offset of [0, data.groove]) {
        const steps = Math.max(1, (last - first) * 8);
        let previous = project(backbone(first, offset));
        for (let k = 1; k <= steps; k++) {
          const point = project(backbone(first + ((last - first) * k) / steps, offset));
          const a = previous, depth = (a.depth + point.depth) / 2, f = (a.f + point.f) / 2;
          items.push({ depth, paint: () => line(a, point, shade(colors.backbone, depth), 1.15 * scale * f) });
          previous = point;
        }
      }

      // Peldaños: mitad del color de la base, mitad del de su pareja.
      rungs = [];
      for (const pair of data.pairs) {
        const p1 = project(pair.strand1), p2 = project(pair.strand2);
        const mid = { x: (p1.x + p2.x) / 2, y: (p1.y + p2.y) / 2, depth: (p1.depth + p2.depth) / 2, f: (p1.f + p2.f) / 2 };
        const width1 = 0.7 * scale * mid.f;
        items.push({ depth: (p1.depth + mid.depth) / 2, paint: () => line(p1, mid, shade(colors[pair.base], (p1.depth + mid.depth) / 2), width1) });
        items.push({ depth: (p2.depth + mid.depth) / 2, paint: () => line(mid, p2, shade(colors[pair.partner], (p2.depth + mid.depth) / 2), width1) });
        rungs.push({ ...mid, pair });
      }

      items.sort((a, b) => a.depth - b.depth).forEach((item) => item.paint());
      labelEnds(project);
    }

    function line(a, b, color, lineWidth) {
      ctx.strokeStyle = color;
      ctx.lineWidth = Math.max(1, lineWidth);
      ctx.lineCap = "round";
      ctx.beginPath();
      ctx.moveTo(a.x, a.y);
      ctx.lineTo(b.x, b.y);
      ctx.stroke();
    }

    /** 5' y 3' en los extremos: las dos hebras van en sentidos opuestos. */
    function labelEnds(project) {
      const first = data.pairs[0], last = data.pairs[data.pairs.length - 1];
      ctx.fillStyle = mix(colors.backbone, colors.background, 0);
      ctx.font = `600 13px ${getComputedStyle(document.body).fontFamily}`;
      ctx.textBaseline = "middle";
      const label = (text, point, side) => {
        ctx.textAlign = side < 0 ? "right" : "left";
        ctx.fillText(text, point.x + side * 10, point.y);
      };
      label("5'", project(first.strand1), -1);
      label("3'", project(last.strand1), 1);
      label("3'", project(first.strand2), -1);
      label("5'", project(last.strand2), 1);
    }

    // ---------- Animación ----------

    function tick(time) {
      if (!alive) return;
      const elapsed = lastTime ? Math.min(time - lastTime, 100) / 1000 : 0;
      lastTime = time;
      if (view.playing && !dragging && visible) {
        view.angle += SPIN_SPEED * elapsed;
        dirty = true;
      }
      if (dirty && visible) { draw(); dirty = false; }
      frame = requestAnimationFrame(tick);
    }

    function setPlaying(playing) {
      view.playing = playing;
      toggle.textContent = t(playing ? "web.helix.pause" : "web.helix.play");
    }

    // ---------- Interacción ----------

    function showTip(event) {
      const box = canvas.getBoundingClientRect();
      const x = event.clientX - box.left, y = event.clientY - box.top;
      let best = null;
      for (const rung of rungs) {
        const distance = Math.hypot(rung.x - x, rung.y - y);
        if (distance < 18 && (!best || rung.depth > best.depth)) best = rung;
      }
      if (!best) { tip.hidden = true; return; }
      const { pair } = best;
      tip.textContent = t("web.helix.pair", {
        position: pair.index + 1,
        base: `${pair.base} (${t(`web.bases.${pair.base}`)})`,
        partner: `${pair.partner} (${t(`web.bases.${pair.partner}`)})`,
      });
      tip.style.left = `${best.x}px`;
      tip.style.top = `${best.y}px`;
      tip.hidden = false;
    }

    const onPointerDown = (event) => {
      dragging = { x: event.clientX, angle: view.angle };
      canvas.setPointerCapture(event.pointerId);
      tip.hidden = true;
    };
    const onPointerMove = (event) => {
      if (dragging) {
        view.angle = dragging.angle + (event.clientX - dragging.x) * 0.012;
        dirty = true;
      } else if (event.pointerType === "mouse") {
        showTip(event);
      }
    };
    const onPointerUp = (event) => {
      if (dragging && Math.abs(event.clientX - dragging.x) < 4) showTip(event);   // un toque muestra el par
      dragging = null;
    };
    const onLeave = () => { if (!dragging) tip.hidden = true; };
    const onKey = (event) => {
      if (event.key === "ArrowLeft" || event.key === "ArrowRight") {
        view.angle += event.key === "ArrowLeft" ? -0.2 : 0.2;
        dirty = true;
        event.preventDefault();
      } else if (event.key === " ") {
        setPlaying(!view.playing);
        event.preventDefault();
      }
    };
    const onToggle = () => setPlaying(!view.playing);
    const onRange = () => {
      clearTimeout(debounce);
      debounce = setTimeout(() => { view.start = Number(range.value) - 1; load(); }, 90);
    };
    const onTheme = () => { colors = readColors(); dirty = true; };

    canvas.addEventListener("pointerdown", onPointerDown);
    canvas.addEventListener("pointermove", onPointerMove);
    canvas.addEventListener("pointerup", onPointerUp);
    canvas.addEventListener("pointercancel", onPointerUp);
    canvas.addEventListener("pointerleave", onLeave);
    canvas.addEventListener("keydown", onKey);
    toggle.addEventListener("click", onToggle);
    range.addEventListener("input", onRange);
    darkMode.addEventListener("change", onTheme);

    // Si cambia el ancho, cambia cuántos pares caben: se pide otra ventana.
    let lastWindow = windowSize();
    const resizer = new ResizeObserver(() => {
      if (windowSize() !== lastWindow) { lastWindow = windowSize(); load(); } else resize();
    });
    resizer.observe(root);
    // Fuera de la pantalla no se dibuja: ahorra batería.
    const watcher = new IntersectionObserver(([entry]) => { visible = entry.isIntersecting; dirty = true; });
    watcher.observe(canvas);

    setPlaying(view.playing);
    resize();
    load();
    frame = requestAnimationFrame(tick);

    return {
      destroy() {
        alive = false;
        cancelAnimationFrame(frame);
        clearTimeout(debounce);
        resizer.disconnect();
        watcher.disconnect();
        darkMode.removeEventListener("change", onTheme);
      },
    };
  }

  return { mount };
})();
