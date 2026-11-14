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
  for (let i = 0; i < a.length; i++) markers += a[i] === b[i] ? "|" : " ";
  let html = `<h3>${escapeHtml(title)}</h3><div class="stats"><div class="stat"><b>${result.score}</b>Alignment score</div><div class="stat"><b>${a.length}</b>Alignment columns</div></div>`;
  html += `<div class="mono">${escapeHtml(a)}\n<span class="matchline">${escapeHtml(markers)}</span>\n${escapeHtml(b)}</div>`;
  html += `<h3 style="margin-top:18px">Dynamic-programming score matrix</h3>${matrixTable(result.matrix, result.path)}`;
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
    const data = await response.json();
    if (!response.ok) throw new Error(data.error);
    let html = `<h3>Progressive star multiple sequence alignment</h3><div class="stats"><div class="stat"><b>${data.rows.length}</b>Sequences</div><div class="stat"><b>${data.rows[0].sequence.length}</b>Alignment columns</div><div class="stat"><b>${escapeHtml(data.rows[data.center_index].name)}</b>Chosen centre sequence</div></div>`;
    data.rows.forEach(row => { html += `<div class="alignment"><div class="seqname">${escapeHtml(row.name)}</div><div class="seq">${escapeHtml(row.sequence)}</div></div>`; });
    html += `<h3 style="margin-top:18px">Consensus sequence</h3><div class="mono">${escapeHtml(data.analysis.consensus)}</div>`;
    html += `<h3 style="margin-top:18px">Conservation by position (%)</h3><div class="mono">${data.analysis.conservation.map((v,i)=>`${i+1}:${v}%`).join("   ")}</div>`;
    html += `<h3 style="margin-top:18px">Pairwise identity (%)</h3><table class="data-table"><tbody>`;
    html += `<tr><th>Sequence</th>${data.rows.map(r=>`<th>${escapeHtml(r.name)}</th>`).join("")}</tr>`;
    data.analysis.identity.forEach((row,i)=>{html += `<tr><th>${escapeHtml(data.rows[i].name)}</th>${row.map(v=>`<td>${v}%</td>`).join("")}</tr>`;});
    html += `</tbody></table><button class="download" id="downloadFasta">Download aligned FASTA</button>`;
    target.classList.remove("hidden"); target.innerHTML = html;
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
    const response = await fetch("/api/ncbi", {method:"POST",headers:{"Content-Type":"application/json"},body:JSON.stringify({accession:byId("accession").value,email:byId("email").value})});
    const data = await response.json();
    if (!response.ok) throw new Error(data.error);
    target.classList.remove("hidden");
    target.innerHTML = `<h3>${escapeHtml(data.description)}</h3><div class="stats"><div class="stat"><b>${data.length}</b>Base pairs</div><div class="stat"><b>${escapeHtml(data.accession)}</b>Accession</div></div><div class="mono">${escapeHtml(data.sequence)}</div><button class="download" id="addFetched">Add fetched sequence to MSA input</button>`;
    byId("addFetched").addEventListener("click", () => { const box = byId("msaInput"); box.value += `\n>${data.accession}\n${data.sequence}`; });
  } catch (error) { showError(target, error.message); }
  finally { button.disabled = false; }
});
