
/* msa_tree.js - UPGMA tree visualization for NEOHELIX */
(function (global) {
  "use strict";

  const CSS = `
  .mt-wrap {
    background: var(--panel, #0e1c27);
    border: 1px solid var(--line, #24404e);
    border-radius: 12px;
    padding: 16px;
    color: var(--text, #e7f1f5);
    overflow-x: auto;
  }
  .mt-wrap h3 {
    margin: 0 0 8px;
    font-family: Inter, system-ui, sans-serif;
  }
  .mt-subtitle {
    color: var(--muted, #8da6b2);
    font: 12px Inter, system-ui, sans-serif;
    margin-bottom: 12px;
  }
  .mt-svg {
    display: block;
    min-width: 420px;
    width: 100%;
    overflow: visible;
  }
  .mt-branch {
    fill: none;
    stroke: var(--cyan, #4ce4ed);
    stroke-width: 2;
  }
  .mt-label {
    fill: var(--text, #e7f1f5);
    font: 12px Inter, system-ui, sans-serif;
  }
  .mt-axis {
    fill: var(--muted, #8da6b2);
    font: 10px Inter, system-ui, sans-serif;
  }
  .mt-error {
    color: #ff8a8a;
    font: 13px Inter, system-ui, sans-serif;
  }`;

  function addStyle() {
    if (document.getElementById("msa-tree-styles")) return;
    const style = document.createElement("style");
    style.id = "msa-tree-styles";
    style.textContent = CSS;
    document.head.appendChild(style);
  }

  // Convert percentage identity into distance.
  function makeDistanceMatrix(identity) {
    return identity.map(row =>
      row.map(value => 1 - Number(value) / 100)
    );
  }

  function buildUPGMA(names, identity) {
    const distances = makeDistanceMatrix(identity);

    let clusters = names.map((name, i) => ({
      id: i,
      name,
      size: 1,
      height: 0,
      left: null,
      right: null,
      leftLength: 0,
      rightLength: 0
    }));

    let nextId = names.length;
    const distanceMap = new Map();

    function key(a, b) {
      return a < b ? `${a},${b}` : `${b},${a}`;
    }

    for (let i = 0; i < clusters.length; i++) {
      for (let j = i + 1; j < clusters.length; j++) {
        distanceMap.set(key(i, j), distances[i][j]);
      }
    }

    while (clusters.length > 1) {
      let bestDistance = Infinity;
      let bestI = -1;
      let bestJ = -1;

      for (let i = 0; i < clusters.length; i++) {
        for (let j = i + 1; j < clusters.length; j++) {
          const d = distanceMap.get(
            key(clusters[i].id, clusters[j].id)
          );

          if (d < bestDistance) {
            bestDistance = d;
            bestI = i;
            bestJ = j;
          }
        }
      }

      const a = clusters[bestI];
      const b = clusters[bestJ];
      const parentHeight = bestDistance / 2;

      const parent = {
        id: nextId++,
        name: "",
        size: a.size + b.size,
        height: parentHeight,
        left: a,
        right: b,
        leftLength: Math.max(0, parentHeight - a.height),
        rightLength: Math.max(0, parentHeight - b.height)
      };

      // UPGMA: weighted average of distances between clusters.
      for (const other of clusters) {
        if (other.id === a.id || other.id === b.id) continue;

        const da = distanceMap.get(key(a.id, other.id));
        const db = distanceMap.get(key(b.id, other.id));

        const combined =
          (da * a.size + db * b.size) / (a.size + b.size);

        distanceMap.set(key(parent.id, other.id), combined);
      }

      clusters = clusters.filter(
        (_, index) => index !== bestI && index !== bestJ
      );
      clusters.push(parent);
    }

    return clusters[0];
  }

  function render(container, data) {
    addStyle();

    const rows = data.rows;
    const identity = data.analysis && data.analysis.identity;

    if (!rows || rows.length < 2) {
      throw new Error("At least two sequences are required for a tree.");
    }

    if (!identity || identity.length !== rows.length ||
        identity.some(row => !row || row.length !== rows.length)) {
      throw new Error("Pairwise identity data is missing or invalid.");
    }

    const names = rows.map((row, i) => row.name || `Sequence ${i + 1}`);
    const root = buildUPGMA(names, identity);

    const width = Math.max(520, names.length * 110 + 180);
    const height = Math.max(180, names.length * 58 + 50);
    const left = 35;
    const right = width - 155;
    const top = 25;
    const bottom = height - 30;
    const maxDistance = root.height || 1;
    const plotWidth = right - left;

    let leafIndex = 0;
    const positions = new Map();

    function assignPositions(node) {
      if (!node.left && !node.right) {
        const y = names.length === 1
          ? (top + bottom) / 2
          : top + leafIndex * ((bottom - top) / (names.length - 1));

        leafIndex++;
        positions.set(node.id, {
          x: right,
          y
        });
        return y;
      }

      const leftY = assignPositions(node.left);
      const rightY = assignPositions(node.right);
      const y = (leftY + rightY) / 2;
      const x = left + ((maxDistance - node.height) / maxDistance) * plotWidth;

      positions.set(node.id, { x, y });
      return y;
    }

    assignPositions(root);

    let branches = "";
    let labels = "";

    function draw(node) {
      const p = positions.get(node.id);

      if (!node.left && !node.right) {
        labels += `<text class="mt-label" x="${p.x + 9}" y="${p.y + 4}">${escapeText(node.name)}</text>`;
        return;
      }

      const l = positions.get(node.left.id);
      const r = positions.get(node.right.id);

      branches += `
        <path class="mt-branch"
          d="M ${p.x} ${l.y} H ${l.x} M ${p.x} ${r.y} H ${r.x} M ${p.x} ${l.y} V ${r.y}" />
      `;

      draw(node.left);
      draw(node.right);
    }

    function escapeText(value) {
      return String(value).replace(/[&<>"']/g, ch => ({
        "&": "&amp;", "<": "&lt;", ">": "&gt;",
        '"': "&quot;", "'": "&#39;"
      })[ch]);
    }

    draw(root);

    // Horizontal distance scale.
    let axis = "";
    for (let i = 0; i <= 4; i++) {
      const x = left + (plotWidth * i / 4);
      const distance = maxDistance * (1 - i / 4);

      axis += `
        <line x1="${x}" y1="${bottom + 8}"
              x2="${x}" y2="${bottom + 13}"
              stroke="currentColor" opacity=".5" />
        <text class="mt-axis" x="${x}" y="${bottom + 27}"
              text-anchor="middle">${distance.toFixed(2)}</text>
      `;
    }

    container.innerHTML = `
      <div class="mt-wrap">
        <h3>UPGMA Guide Tree</h3>
        <div class="mt-subtitle">
          Clustered using pairwise sequence identity. Horizontal scale represents distance.
        </div>
        <svg class="mt-svg"
             viewBox="0 0 ${width} ${height}"
             role="img"
             aria-label="UPGMA sequence similarity tree">
          ${branches}
          ${labels}
          ${axis}
        </svg>
        <div class="mt-subtitle">
          Distance = 1 − pairwise identity. Smaller distances indicate greater similarity.
        </div>
      </div>
    `;

    return root;
  }

  global.MSATree = { render };
})(window);
