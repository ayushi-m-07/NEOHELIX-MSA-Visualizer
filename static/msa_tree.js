
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
  }
  .mt-grid { stroke: var(--line, #24404e); stroke-dasharray: 3 4; }
  .mt-leaf { fill: var(--cyan, #4ce4ed); }
  .mt-node { fill: var(--panel, #0e1c27); stroke: var(--cyan, #4ce4ed); stroke-width: 2; cursor: help; }
  .mt-node:hover { fill: var(--cyan, #4ce4ed); }
  .mt-len { fill: var(--muted, #8da6b2); font: 10px Inter, system-ui, sans-serif; }
  .mt-pct { fill: var(--green, #91f2c0); font: 10px Inter, system-ui, sans-serif; }
  .mt-sub { fill: var(--muted, #8da6b2); font: 10px ui-monospace, SFMono-Regular, monospace; }`;

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
        // Summary of pairwise similarity
    let closest = { v: -1, i: 0, j: 1 };
    let farthest = { v: Infinity, i: 0, j: 1 };
    let total = 0, pairs = 0;
    for (let i = 0; i < names.length; i++) {
      for (let j = i + 1; j < names.length; j++) {
        const v = Number(identity[i][j]);
        total += v; pairs++;
        if (v > closest.v) closest = { v, i, j };
        if (v < farthest.v) farthest = { v, i, j };
      }
    }
    const meanIdentity = pairs ? total / pairs : 0;

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
    let lengths = "";
    let nodes = "";
    let grid = "";
    function leavesOf(node) {
      if (!node.left && !node.right) return [node.name];
      return leavesOf(node.left).concat(leavesOf(node.right));
    }

    function draw(node) {
      const p = positions.get(node.id);

      if (!node.left && !node.right) {
        const info = rows[node.id] && rows[node.id].stats;
        const sub = info ? `${info.length} bp · GC ${info.gc_percent}%` : "";
        nodes += `<circle class="mt-leaf" cx="${p.x}" cy="${p.y}" r="4"><title>${escapeText(node.name)}${sub ? " - " + escapeText(sub) : ""}</title></circle>`;
        labels += `<text class="mt-label" x="${p.x + 12}" y="${p.y + 4}">${escapeText(node.name)}</text>`;
        if (sub) labels += `<text class="mt-sub" x="${p.x + 12}" y="${p.y + 18}">${escapeText(sub)}</text>`;
        return;
      }

      const l = positions.get(node.left.id);
      const r = positions.get(node.right.id);

      branches += `
        <path class="mt-branch"
          d="M ${p.x} ${l.y} H ${l.x} M ${p.x} ${r.y} H ${r.x} M ${p.x} ${l.y} V ${r.y}" />
      `;

      // Branch-length labels (only where there is room)
      [[node.leftLength, l], [node.rightLength, r]].forEach(([len, pos]) => {
        if (len > 0 && pos.x - p.x > 40) {
          lengths += `<text class="mt-len" x="${(p.x + pos.x) / 2}" y="${pos.y - 5}" text-anchor="middle">${len.toFixed(2)}</text>`;
        }
      });

      // Identity at this merge: distance = 1 - identity, and height = distance / 2
      const identityPct = ((1 - 2 * node.height) * 100).toFixed(1);
      if (Math.min(l.x, r.x) - p.x > 56) {
        lengths += `<text class="mt-pct" x="${p.x + 8}" y="${p.y - 7}">${identityPct}%</text>`;
      }

      const members = leavesOf(node);
      nodes += `<circle class="mt-node" cx="${p.x}" cy="${p.y}" r="5"><title>Cluster of ${members.length}: ${escapeText(members.join(", "))}\nAverage identity between its two branches: ${identityPct}%</title></circle>`;

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
      grid += `<line class="mt-grid" x1="${x}" y1="${top - 10}" x2="${x}" y2="${bottom + 8}" />`;

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
         Pairwise distance = 1 − identity. The scale shows height above the leaves (half the distance between the clusters being joined). Numbers on branches are branch lengths. Hover a node to see its cluster members and average identity.
        </div>
        <div class="stats">
          <div class="stat"><b>${closest.v.toFixed(1)}%</b>Closest pair: ${escapeText(names[closest.i])} &amp; ${escapeText(names[closest.j])}</div>
          <div class="stat"><b>${farthest.v.toFixed(1)}%</b>Most divergent: ${escapeText(names[farthest.i])} &amp; ${escapeText(names[farthest.j])}</div>
          <div class="stat"><b>${meanIdentity.toFixed(1)}%</b>Mean pairwise identity</div>
        </div>        
        <svg class="mt-svg"
             viewBox="0 0 ${width} ${height}"
             role="img"
             aria-label="UPGMA sequence similarity tree">
          ${grid}
          ${branches}
          ${lengths}
          ${nodes}
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
