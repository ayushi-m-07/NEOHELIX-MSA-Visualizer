const byId = (id) => document.getElementById(id);

function showError(target, message) {
  target.classList.remove("hidden");
  target.innerHTML = `<p class="error">${escapeHtml(message)}</p>`;
}

function escapeHtml(value) {
  return String(value).replace(/[&<>"']/g, (char) => ({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;","'":"&#39;"}[char]));
}


function showAlignment(target, title, result) {
  target.classList.remove("hidden");

  const a = result.aligned_seq1;
  const b = result.aligned_seq2;

  let markers = "";
  for (let i = 0; i < a.length; i++) {
    markers += a[i] === b[i] ? "|" : " ";
  }

  let html = `
    <h3>${escapeHtml(title)}</h3>

    <div class="stats">
      <div class="stat">
        <b>${result.score}</b>
        Alignment score
      </div>
      <div class="stat">
        <b>${a.length}</b>
        Alignment columns
      </div>
    </div>

    <h3 style="margin-top:18px">Optimal alignment</h3>

    <div class="mono">${escapeHtml(a)}
<span class="matchline">${escapeHtml(markers)}</span>
${escapeHtml(b)}</div>
  `;

  // Display the complete DP matrix only for small sequences.
  const maxLength = Math.max(a.length, b.length);

  if (result.matrix && maxLength <= 50) {
    html += `
      <h3 style="margin-top:18px">
        Dynamic-programming score matrix
      </h3>
      <p>Highlighted cells show the optimal traceback path.</p>
      ${matrixTable(result.matrix, result.path)}
    `;
  } else {
    html += `
      <h3 style="margin-top:18px">
        Dynamic-programming score matrix
      </h3>
      <p>
        The full matrix is hidden for large sequences
        to keep the interface responsive.
        The alignment score and optimal alignment are
        still available above.
      </p>
    `;
  }

  target.innerHTML = html;
}


function matrixTable(matrix, path) {
  const pathKeys = new Set(path.map(([i,j]) => `${i},${j}`));
  let html = '<div style="overflow-x:auto"><table class="data-table"><tbody>';
  matrix.forEach((row, i) => {
    html += "<tr>";
    row.forEach((value, j) => {
      const tag = pathKeys.has(`${i},${j}`) ? ' class="conserved"' : "";
      html += `<td${tag}>${value}</td>`;
    });
    html += "</tr>";
  });
  return html + "</tbody></table></div>";
}

byId("pairButton").addEventListener("click", async () => {
  const button = byId("pairButton");
  const target = byId("pairResult");
  button.disabled = true;
  try {
    const response = await fetch("/api/pairwise", {method:"POST", headers:{"Content-Type":"application/json"}, body:JSON.stringify({seq1:byId("pair1").value,seq2:byId("pair2").value,mode:byId("mode").value,match:byId("match").value,mismatch:byId("mismatch").value,gap:byId("gap").value})});
    const data = await response.json();
    if (!response.ok) throw new Error(data.error);
    showAlignment(target, byId("mode").value === "global" ? "Needleman–Wunsch global alignment" : "Smith–Waterman local alignment", data.result);
    if (byId("mode").value === "local") target.innerHTML += `<p>Local region: sequence 1 [${data.result.start1}:${data.result.end1}], sequence 2 [${data.result.start2}:${data.result.end2}]</p>`;
  } catch (error) { showError(target, error.message); }
  finally { button.disabled = false; }
});

byId("msaButton").addEventListener("click", async () => {
  const button = byId("msaButton");
  const target = byId("msaResult");
  button.disabled = true;
  try {
    const response = await fetch("/api/msa", {method:"POST",headers:{"Content-Type":"application/json"},body:JSON.stringify({sequences:byId("msaInput").value,match:byId("msaMatch").value,mismatch:byId("msaMismatch").value,gap:byId("msaGap").value})});
    const contentType = response.headers.get("content-type") || "";

    if (!contentType.includes("application/json")) {
      const body = await response.text();
      throw new Error(
      `Server returned HTML instead of JSON. HTTP ${response.status}. ` +
      body.slice(0, 250)
  );
}

    const data = await response.json();

    if (!response.ok) {
      throw new Error(data.error || "MSA failed.");
}
    if (!response.ok) throw new Error(data.error);
    let html = `<h3>Progressive star multiple sequence alignment</h3><div class="stats"><div class="stat"><b>${data.rows.length}</b>Sequences</div><div class="stat"><b>${data.rows[0].sequence.length}</b>Alignment columns</div><div class="stat"><b>${escapeHtml(data.rows[data.center_index].name)}</b>Chosen centre sequence</div></div>`;
    data.rows.forEach(row => { html += `<div class="alignment"><div class="seqname">${escapeHtml(row.name)}</div><div class="seq">${escapeHtml(row.sequence)}</div></div>`; });
    html += `<div id="msaGrid" style="margin:14px 0"></div>`;
    html += `<div id="msaTree" style="margin:14px 0"></div>`;
    html += `<h3 style="margin-top:18px">Consensus sequence</h3><div class="mono">${escapeHtml(data.analysis.consensus)}</div>`;
    html += `<h3 style="margin-top:18px">Conservation by position (%)</h3><div class="mono">${data.analysis.conservation.map((v,i)=>`${i+1}:${v}%`).join("   ")}</div>`;
    html += `<h3 style="margin-top:18px">Pairwise identity (%)</h3><table class="data-table"><tbody>`;
    html += `<tr><th>Sequence</th>${data.rows.map(r=>`<th>${escapeHtml(r.name)}</th>`).join("")}</tr>`;
    data.analysis.identity.forEach((row,i)=>{html += `<tr><th>${escapeHtml(data.rows[i].name)}</th>${row.map(v=>`<td>${v}%</td>`).join("")}</tr>`;});
    html += `</tbody></table><button class="download" id="downloadFasta">Download aligned FASTA</button>`;
    target.classList.remove("hidden"); target.innerHTML = html;
    MSAGrid.render(byId("msaGrid"), { rows: data.rows, analysis: data.analysis });
    MSATree.render(byId("msaTree"), {rows: data.rows,analysis: data.analysis});
    byId("downloadFasta").addEventListener("click", () => {
      const fasta = data.rows.map(row => `>${row.name}\n${row.sequence}`).join("\n") + "\n";
      const blob = new Blob([fasta], {type:"text/plain"});
      const url = URL.createObjectURL(blob);
      const link = document.createElement("a"); link.href = url; link.download = "neohelix_alignment.fasta"; link.click(); URL.revokeObjectURL(url);
    });
  } catch (error) { showError(target, error.message); }
  finally { button.disabled = false; }
});

byId("ncbiButton").addEventListener("click", async () => {
  const button = byId("ncbiButton");
  const target = byId("ncbiResult");
  button.disabled = true;
  try {
    const response = await fetch("/api/ncbi", {method:"POST",headers:{"Content-Type":"application/json"},body:JSON.stringify({accession:byId("accession").value,email:byId("email").value,organism:byId("organism").value})});
    const data = await response.json();
    if (!response.ok) throw new Error(data.error);
    target.classList.remove("hidden");
    target.innerHTML = `<h3>${escapeHtml(data.description)}</h3><div class="stats"><div class="stat"><b>${data.length}</b>Base pairs</div><div class="stat"><b>${escapeHtml(data.accession)}</b>Accession</div></div><div class="mono">${escapeHtml(data.sequence)}</div><button class="download" id="addFetched">Add fetched sequence to MSA input</button>`;
   byId("addFetched").addEventListener("click", (event) => {
  const box = byId("msaInput");
  box.value += `\n>${data.accession}\n${data.sequence}`;

  const addButton = event.currentTarget;
  addButton.textContent = "Added to MSA input ✓";
  addButton.disabled = true;

  // Show a styled notice (uses the .notice CSS)
  if (!byId("addedNotice")) {
    const notice = document.createElement("p");
    notice.id = "addedNotice";
    notice.className = "notice";
    notice.textContent = `${data.accession} was added to the MSA input in section 02.`;
    addButton.insertAdjacentElement("afterend", notice);
  }

  // Bring the MSA box into view
  box.scrollIntoView({ behavior: "smooth", block: "center" });
});
 } catch (error) { showError(target, error.message); }
  finally { button.disabled = false; }
});

// Organism type-ahead for the NCBI fetch box (NCBI fetch only)
// Custom dropdown (the browser's native datalist is unreliable)
(() => {
  const input = byId("organism");
  if (!input) return;
  input.removeAttribute("list");
  input.setAttribute("autocomplete", "off");

  const COMMON = [
    ["Homo sapiens","human"],["Mus musculus","house mouse"],["Rattus norvegicus","Norway rat"],
    ["Danio rerio","zebrafish"],["Drosophila melanogaster","fruit fly"],["Caenorhabditis elegans","nematode"],
    ["Saccharomyces cerevisiae","baker's yeast"],["Schizosaccharomyces pombe","fission yeast"],
    ["Escherichia coli","E. coli"],["Bacillus subtilis",""],["Arabidopsis thaliana","thale cress"],
    ["Oryza sativa","rice"],["Zea mays","maize"],["Gallus gallus","chicken"],["Bos taurus","cattle"],
    ["Sus scrofa","pig"],["Canis lupus familiaris","dog"],["Macaca mulatta","rhesus macaque"],
    ["Pan troglodytes","chimpanzee"],["Xenopus laevis","African clawed frog"]
  ].map(([name, common]) => ({ name, common }));

  const wrap = input.parentElement;
  wrap.style.position = "relative";
  const menu = document.createElement("div");
  menu.style.cssText = "display:none;position:absolute;left:0;right:0;top:100%;margin-top:4px;z-index:50;max-height:240px;overflow-y:auto;background:#0e1c27;border:1px solid #4ce4ed;border-radius:7px;box-shadow:0 10px 30px #000a;";
  wrap.appendChild(menu);

  let items = [];
  let active = -1;
  const cache = new Map();
  let timer = null;
  let controller = null;

  function close() { menu.style.display = "none"; active = -1; }
  function choose(name) { input.value = name; close(); }
  function highlight() {
    Array.from(menu.children).forEach((el, k) => {
      el.style.background = k === active ? "#1d3a4a" : "transparent";
    });
    if (active >= 0 && menu.children[active]) menu.children[active].scrollIntoView({ block: "nearest" });
  }

  function render(list) {
    const seen = new Set();
    items = list.filter((o) => !seen.has(o.name) && seen.add(o.name)).slice(0, 12);
    menu.innerHTML = "";
    active = -1;
    if (!items.length) { close(); return; }
    items.forEach((o) => {
      const row = document.createElement("div");
      row.style.cssText = "padding:8px 11px;cursor:pointer;font:13px/1.4 ui-monospace,SFMono-Regular,monospace;color:#e7f1f5;";
      row.textContent = o.name;
      if (o.common) {
        const c = document.createElement("span");
        c.textContent = "  " + o.common;
        c.style.cssText = "color:#8da6b2;font-size:12px;";
        row.appendChild(c);
      }
      row.addEventListener("mousedown", (e) => { e.preventDefault(); choose(o.name); });
      row.addEventListener("mouseenter", () => { active = items.indexOf(o); highlight(); });
      menu.appendChild(row);
    });
    menu.style.display = "block";
  }

  function localMatches(q) {
    const t = q.toLowerCase();
    return COMMON.filter((o) => !t || o.name.toLowerCase().includes(t) || o.common.toLowerCase().includes(t));
  }

  function update() {
    clearTimeout(timer);
    if (controller) controller.abort();
    const q = input.value.trim();
    const local = localMatches(q);
    render(local);
    if (q.length < 3) return;

    const key = q.toLowerCase();
    if (cache.has(key)) { render(local.concat(cache.get(key))); return; }

    timer = setTimeout(async () => {
      controller = new AbortController();
      try {
        const email = encodeURIComponent(byId("email").value.trim());
        const response = await fetch(`/api/organisms?q=${encodeURIComponent(q)}&email=${email}`, { signal: controller.signal });
        const data = await response.json();
        const live = data.organisms || [];
        if (live.length) cache.set(key, live);
        if (input.value.trim() === q && document.activeElement === input) render(local.concat(live));
      } catch (error) { /* suggestions are optional / request cancelled */ }
    }, 250);
  }

  input.addEventListener("input", update);
  input.addEventListener("focus", update);
  input.addEventListener("click", () => { if (menu.style.display === "none") update(); });
  input.addEventListener("blur", close);
  input.addEventListener("keydown", (e) => {
    if (e.key === "ArrowDown") {
      if (menu.style.display === "none") { update(); } else { active = Math.min(active + 1, items.length - 1); highlight(); }
      e.preventDefault();
    } else if (e.key === "ArrowUp") {
      active = Math.max(active - 1, 0); highlight(); e.preventDefault();
    } else if (e.key === "Enter" && active >= 0 && menu.style.display !== "none") {
      choose(items[active].name); e.preventDefault();
    } else if (e.key === "Escape") { close(); }
  });
})();
