/* msa_grid.js - colour-coded multiple sequence alignment grid for NEOHELIX.
 *
 * Standalone: no dependencies, injects its own CSS, so style.css / app.js /
 * index.html do not need to change to try it. It reuses your theme variables
 * (--panel, --line, --text, --muted, --cyan) when they exist.
 *
 * Usage (data is exactly what /api/msa already returns):
 *   MSAGrid.render(containerElement, { rows: data.rows, analysis: data.analysis });
 *   // rows: [{name, sequence}] with gaps as "-"; analysis: {consensus, conservation}
 *   // analysis is optional - it is recomputed here if missing.
 *
 * Three views: "bases" (colour by nucleotide), "heat" (colour by column
 * conservation), "diff" (dots where a row matches the consensus).
 */
(function (global) {
  "use strict";

  const CSS = `
.mg-wrap{--cell:22px;background:var(--panel,#0e1c27);border:1px solid var(--line,#24404e);border-radius:12px;padding:14px;color:var(--text,#e7f1f5);font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace}
.mg-ui{font-family:Inter,system-ui,-apple-system,sans-serif;font-size:13px;color:var(--muted,#8da6b2)}
.mg-toolbar{display:flex;flex-wrap:wrap;gap:16px;align-items:center;margin-bottom:12px}
.mg-seg{display:inline-flex;border:1px solid var(--line,#24404e);border-radius:8px;overflow:hidden}
.mg-seg button{background:transparent;border:0;color:inherit;padding:6px 12px;font:inherit;cursor:pointer}
.mg-seg button.on{background:var(--cyan,#4ce4ed);color:#04141a;font-weight:600}
.mg-zoom{display:flex;align-items:center;gap:8px}
.mg-scroll{overflow-x:auto;padding-bottom:8px}
.mg-grid{width:max-content;min-width:100%}
.mg-row{display:flex;align-items:center;height:calc(var(--cell) + 3px)}
.mg-name{position:sticky;left:0;z-index:2;flex:0 0 120px;width:120px;padding-right:10px;background:var(--panel,#0e1c27);font:12px Inter,system-ui,sans-serif;color:var(--text,#e7f1f5);white-space:nowrap;overflow:hidden;text-overflow:ellipsis}
.mg-col{flex:0 0 var(--cell);width:var(--cell);height:var(--cell);margin-right:1px;display:flex;align-items:center;justify-content:center;font-size:calc(var(--cell)*.56);font-weight:600;border-radius:3px}
.mg-ruler{height:18px}
.mg-ruler .mg-col{height:16px;justify-content:flex-end;overflow:visible;white-space:nowrap;font-size:10px;font-weight:400;color:var(--muted,#8da6b2)}
.mg-bars{height:44px;align-items:flex-end;margin-top:4px}
.mg-bars .mg-col{height:42px;align-items:flex-end;border-radius:0}
.mg-bars i{display:block;width:100%;border-radius:2px 2px 0 0}
.mg-cons{margin-top:2px}
.mg-A{background:rgba(95,227,161,.22);color:#5fe3a1}
.mg-T{background:rgba(255,138,138,.22);color:#ff8a8a}
.mg-G{background:rgba(255,200,87,.22);color:#ffc857}
.mg-C{background:rgba(106,169,255,.22);color:#6aa9ff}
.mg-N{background:rgba(183,196,203,.16);color:#b7c4cb}
.mg-gap{color:#46606c;font-weight:400}
.mg-dot{color:#38515d;font-weight:400}
.mg-hit{box-shadow:inset 0 0 0 1px currentColor}
.mg-legend{display:flex;flex-wrap:wrap;gap:14px;margin-top:10px}
.mg-legend b{display:inline-block;width:10px;height:10px;border-radius:3px;margin-right:6px;vertical-align:-1px}
.mg-tip{position:fixed;z-index:50;display:none;pointer-events:none;background:#04141a;border:1px solid var(--line,#24404e);border-radius:8px;padding:7px 10px;font:12px/1.5 Inter,system-ui,sans-serif;color:#e7f1f5;white-space:nowrap}`;

  const BASE = { A: "A", T: "T", U: "T", G: "G", C: "C", N: "N" };
  const LEGEND = [["A", "#5fe3a1"], ["T / U", "#ff8a8a"], ["G", "#ffc857"], ["C", "#6aa9ff"], ["gap", "#46606c"]];

  const esc = (s) => String(s).replace(/[&<>"']/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c]));

  function ensureStyle() {
    if (document.getElementById("msa-grid-styles")) return;
    const s = document.createElement("style");
    s.id = "msa-grid-styles";
    s.textContent = CSS;
    document.head.appendChild(s);
  }

  // Same definition as backend alignment_analysis: % of ALL rows sharing the top base.
  function analyse(rows) {
    const consensus = [], conservation = [];
    for (let c = 0; c < rows[0].sequence.length; c++) {
      const counts = {}; let best = "-", top = 0;
      for (const r of rows) {
        const ch = r.sequence[c];
        if (ch === "-") continue;
        counts[ch] = (counts[ch] || 0) + 1;
        if (counts[ch] > top) { top = counts[ch]; best = ch; }
      }
      consensus.push(best);
      conservation.push(top ? Math.round((top / rows.length) * 1000) / 10 : 0);
    }
    return { consensus: consensus.join(""), conservation };
  }

  function render(container, data, opts) {
    opts = opts || {};
    ensureStyle();
    const rows = data.rows;
    if (!rows || !rows.length) throw new Error("MSAGrid: no rows supplied");
    const width = rows[0].sequence.length;
    if (rows.some((r) => r.sequence.length !== width)) throw new Error("MSAGrid: all rows must be the same length");

    const an = data.analysis && data.analysis.consensus && data.analysis.conservation ? data.analysis : analyse(rows);
    const cons = an.conservation, consensus = an.consensus;
    const state = { mode: opts.mode || "bases", cell: opts.cellSize || 22 };
    const resPos = rows.map((r) => { let k = 0; return Array.from(r.sequence, (ch) => (ch === "-" ? k : ++k)); });
    const fullyConserved = cons.filter((v) => v === 100).length;

    container.innerHTML =
      `<div class="mg-wrap mg-ui">
        <div class="mg-toolbar">
          <div class="mg-seg" role="group" aria-label="Colour mode">
            <button data-mode="bases">Bases</button><button data-mode="heat">Conservation</button><button data-mode="diff">Differences</button>
          </div>
          <label class="mg-zoom">Cell size <input type="range" min="10" max="34" value="${state.cell}"></label>
          <span>${rows.length} sequences &times; ${width} columns &middot; ${fullyConserved} fully conserved</span>
        </div>
        <div class="mg-scroll"><div class="mg-grid"></div></div>
        <div class="mg-legend">${LEGEND.map(([t, c]) => `<span><b style="background:${c}"></b>${t}</span>`).join("")}<span id="mg-hint"></span></div>
        <div class="mg-tip"></div>
      </div>`;

    const wrap = container.querySelector(".mg-wrap");
    const grid = container.querySelector(".mg-grid");
    const tip = container.querySelector(".mg-tip");
    const hint = container.querySelector("#mg-hint");
    const buttons = container.querySelectorAll(".mg-seg button");
    const HINTS = { bases: "Colour = nucleotide", heat: "Brighter = more conserved column", diff: "Dots match the consensus; coloured bases differ" };

    function cell(r, c) {
      const ch = rows[r].sequence[c];
      const head = `<div class="mg-col `, tail = `" data-r="${r}" data-c="${c}"`;
      if (ch === "-") return `${head}mg-gap${tail}>-</div>`;
      const cls = "mg-" + (BASE[ch] || "N");
      if (state.mode === "heat") {
        const a = 0.1 + 0.85 * (cons[c] / 100);
        return `${head}${tail} style="background:rgba(76,228,237,${a.toFixed(2)});color:${a > 0.55 ? "#04141a" : "#cfeff2"}">${esc(ch)}</div>`;
      }
      if (state.mode === "diff") {
        return ch === consensus[c] ? `${head}mg-dot${tail}>&middot;</div>` : `${head}${cls} mg-hit${tail}>${esc(ch)}</div>`;
      }
      return `${head}${cls}${tail}>${esc(ch)}</div>`;
    }

    function build() {
      wrap.style.setProperty("--cell", state.cell + "px");
      hint.textContent = HINTS[state.mode];
      buttons.forEach((b) => b.classList.toggle("on", b.dataset.mode === state.mode));
      let h = '<div class="mg-row mg-ruler"><div class="mg-name"></div>';
      for (let c = 0; c < width; c++) h += `<div class="mg-col">${c === 0 || (c + 1) % 10 === 0 ? c + 1 : ""}</div>`;
      h += "</div>";
      rows.forEach((r, i) => {
        h += `<div class="mg-row"><div class="mg-name" title="${esc(r.name)}">${esc(r.name)}</div>`;
        for (let c = 0; c < width; c++) h += cell(i, c);
        h += "</div>";
      });
      h += '<div class="mg-row mg-bars"><div class="mg-name">Conservation</div>';
      for (let c = 0; c < width; c++) h += `<div class="mg-col" data-c="${c}"><i style="height:${cons[c]}%;background:rgba(76,228,237,${(0.25 + 0.75 * cons[c] / 100).toFixed(2)})"></i></div>`;
      h += '</div><div class="mg-row mg-cons"><div class="mg-name">Consensus</div>';
      for (let c = 0; c < width; c++) {
        const ch = consensus[c];
        h += `<div class="mg-col ${ch === "-" ? "mg-gap" : "mg-" + (BASE[ch] || "N")}" data-c="${c}">${esc(ch)}</div>`;
      }
      grid.innerHTML = h + "</div>";
    }

    wrap.addEventListener("click", (e) => {
      const b = e.target.closest("[data-mode]");
      if (b) { state.mode = b.dataset.mode; build(); }
    });
    wrap.querySelector("input[type=range]").addEventListener("input", (e) => { state.cell = +e.target.value; wrap.style.setProperty("--cell", state.cell + "px"); });
    wrap.querySelector(".mg-scroll").addEventListener("mousemove", (e) => {
      const t = e.target.closest("[data-c]");
      if (!t) { tip.style.display = "none"; return; }
      const c = +t.dataset.c;
      let html = "";
      if (t.dataset.r !== undefined) {
        const r = +t.dataset.r, ch = rows[r].sequence[c];
        html = `<b>${esc(rows[r].name)}</b><br>${ch === "-" ? "Gap" : "Base " + esc(ch) + " (residue " + resPos[r][c] + ")"} at column ${c + 1}<br>`;
      } else html = `Column ${c + 1}<br>`;
      tip.innerHTML = html + `Conservation ${cons[c]}% &middot; consensus ${esc(consensus[c])}`;
      tip.style.display = "block";
      tip.style.left = e.clientX + 14 + "px";
      tip.style.top = e.clientY + 14 + "px";
    });
    wrap.querySelector(".mg-scroll").addEventListener("mouseleave", () => { tip.style.display = "none"; });

    build();
    return { setMode: (m) => { state.mode = m; build(); }, setCellSize: (n) => { state.cell = n; build(); } };
  }

  global.MSAGrid = { render };
})(window);
