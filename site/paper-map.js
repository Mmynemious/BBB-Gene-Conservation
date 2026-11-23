// Paper map for the Literature Review page.
// To add a paper: add a node (x, y in a 900 x 560 box) and its links below.
(function () {
  const groups = {
    expr:  { label: "BBB across species",       color: "#2C6C4C" },
    model: { label: "Human models",             color: "#6B4C9A" },
    pipe:  { label: "Cross-species pipelines",  color: "#1F4FA8" },
    sel:   { label: "Immune gene evolution",    color: "#875110" },
    wider: { label: "Wider reading",            color: "#7E8885" },
    self:  { label: "This project",             color: "#141817" },
  };

  const nodes = [
    { id: "project",    label: "This project",       sub: "CDS conservation",       g: "self",  x: 450, y: 285 },
    { id: "song",       label: "Song 2020",          sub: "Human · mouse, LCM",     g: "expr",  x: 110, y: 70,  href: "papers/song2020.qmd" },
    { id: "miao",       label: "Miao 2026",          sub: "Mouse · human, scRNA",   g: "expr",  x: 300, y: 145, href: "papers/miao2026.qmd" },
    { id: "kumabe",     label: "Kumabe 2025",        sub: "Primate proteomics",     g: "expr",  x: 520, y: 140, href: "papers/kumabe2025.qmd" },
    { id: "obrown",     label: "O'Brown 2018",       sub: "Review across organisms",g: "expr",  x: 420, y: 40,  href: "papers/obrown2018.qmd" },
    { id: "park",       label: "Park 2019",          sub: "Human BBB chip",         g: "model", x: 700, y: 80,  href: "papers/park2019.qmd" },
    { id: "sciam",      label: "Scientific American",sub: "The Healthy Brain",      g: "wider", x: 790, y: 190, href: "reading/sciam-healthy-brain-2026.qmd" },
    { id: "loverso",    label: "LoVerso 2015",       sub: "Pipeline in one graph",  g: "pipe",  x: 110, y: 480, href: "papers/loverso2015.qmd" },
    { id: "sheng",      label: "Sheng 2022",         sub: "CoRMAP, OrthoMCL",       g: "pipe",  x: 110, y: 330, href: "papers/sheng2022.qmd" },
    { id: "beker",      label: "Beker 2025",         sub: "Human BMEC enrichment",  g: "pipe",  x: 330, y: 500, href: "papers/beker2025.qmd" },
    { id: "barr",       label: "Barr 2025",          sub: "CrossFilt",              g: "pipe",  x: 290, y: 395, href: "papers/barr2025.qmd" },
    { id: "nandakumar", label: "Nandakumar 2025",    sub: "Rodent immune genes",    g: "sel",   x: 640, y: 470, href: "papers/nandakumar2025.qmd" },
    { id: "shultz",     label: "Shultz 2019",        sub: "Birds and mammals",      g: "sel",   x: 790, y: 370, href: "papers/shultz2019.qmd" },
    { id: "cury",       label: "Cury 2024",          sub: "Immune genes from bacteria", g: "sel", x: 800, y: 515, href: "papers/cury2024.qmd" },
  ];

  // kind: agrees | method | limit | shared
  const links = [
    ["song", "miao", "agrees", "The same transporter differences in both: Slco1c1 and Slc22a8 mouse-enriched, A2M human-specific. Different tissue, different methods, six years apart."],
    ["song", "park", "shared", "Same lab, Palecek and Shusta. Song is the descriptive groundwork; Park is the engineered model built on it."],
    ["song", "obrown", "limit", "Both point at brain region as an uncontrolled variable."],
    ["obrown", "kumabe", "agrees", "O'Brown cites primate and rodent proteomics from the Ohtsuki group. Kumabe is the current instalment of that work."],
    ["obrown", "park", "limit", "Both undercut the premise from different sides: barrier identity is induced by the tissue around it, and a human chip skips the species jump."],
    ["obrown", "sciam", "agrees", "The leaky barrier thread from the special edition gets a source: permeability differs between regions in mice."],
    ["park", "kumabe", "shared", "How faithful in vitro models are: claudin-5 absent from the immortalised hCMEC/D3 line in Kumabe, present on Park's chip."],
    ["miao", "kumabe", "shared", "The orthologue problem: Miao restricts to one-to-one pairs, Kumabe names that orthologues are not always one to one."],
    ["miao", "barr", "method", "If Miao aligned each species to its own genome, that is close to the approach CrossFilt shows has the most false positives. Ribosomal genes: Miao excludes them, CrossFilt filters their reads."],
    ["miao", "beker", "shared", "Miao compares mouse and human; Beker compares brain endothelium against other endothelium within human only."],
    ["loverso", "barr", "method", "Both put the two species on the same features before counting: LoVerso lifts exons across with axtLift, CrossFilt lifts reads across UCSC chains."],
    ["loverso", "sheng", "method", "Two whole pipelines for cross-species expression: LoVerso aligns to existing genomes, CoRMAP assembles de novo and groups genes with OrthoMCL."],
    ["loverso", "beker", "method", "LoVerso lays out the steps; Beker shows how to handle multiple datasets and quality control."],
    ["nandakumar", "shultz", "agrees", "Immune genes are enriched for positive selection, within rodents and across birds and mammals. Both rely on BUSTED."],
    ["cury", "shultz", "shared", "Both follow immune genes across very long evolutionary distances: Shultz across birds and mammals, Cury back to bacteria."],
    ["miao", "project", "shared", "Same two species and the same question, measured by expression instead of coding sequence. The SLC family diverges most."],
    ["kumabe", "project", "agrees", "Macaque sits closer to human than mouse does, from protein expression instead of sequence."],
    ["song", "project", "limit", "Song's pericyte markers, PDGFRB, CSPG4, KCNJ8, ABCC9, NOTCH3 and VTN, are all on the gene list."],
    ["obrown", "project", "limit", "If barrier identity is induced by the tissue, coding sequence conservation is a weaker proxy for barrier function."],
    ["park", "project", "limit", "Claudin-5 is central to a working barrier but missing from the enrichment-based gene list."],
    ["barr", "project", "method", "A ready-made algorithm for the expression side of the method-driven pipeline, and more support for one-to-one orthologues."],
    ["sheng", "project", "shared", "Runs into all the orthologue problems this project had."],
    ["loverso", "project", "method", "The frame for where each workflow fits in the method-driven pipeline."],
    ["nandakumar", "project", "method", "dN/dS on coding sequence against a random control set, the same comparison as BBB genes against control genes."],
    ["shultz", "project", "method", "dN/dS on coding sequence across a whole genome, with one sequence per species in the species tree set."],
    ["cury", "project", "method", "Uses conservation the other way round: homology search and phylogeny from a known system to find new immune genes in humans."],
  ];

  const kinds = {
    agrees: { label: "Agrees or replicates", color: "#2C6C4C", dash: "" },
    method: { label: "Supplies a method",    color: "#1F4FA8", dash: "" },
    limit:  { label: "Raises a limitation",  color: "#A02D21", dash: "6 4" },
    shared: { label: "Same problem or lab",  color: "#7E8885", dash: "2 4" },
  };

  const root = document.getElementById("paper-map");
  if (!root) return;
  const NS = "http://www.w3.org/2000/svg";
  const byId = Object.fromEntries(nodes.map(n => [n.id, n]));
  const W = 148, H = 46;

  const legend = document.createElement("div");
  legend.className = "pm-legend";
  Object.entries(kinds).forEach(([k, v]) => {
    const b = document.createElement("button");
    b.type = "button"; b.dataset.kind = k; b.setAttribute("aria-pressed", "true");
    b.innerHTML = `<svg width="26" height="8" aria-hidden="true"><line x1="0" y1="4" x2="26" y2="4" stroke="${v.color}" stroke-width="2.5" stroke-dasharray="${v.dash}"/></svg>${v.label}`;
    b.addEventListener("click", () => {
      const on = b.getAttribute("aria-pressed") !== "true";
      b.setAttribute("aria-pressed", String(on));
      svg.querySelectorAll(`.pm-link[data-kind="${k}"]`).forEach(l => l.classList.toggle("pm-off", !on));
    });
    legend.appendChild(b);
  });
  root.appendChild(legend);

  const glegend = document.createElement("div");
  glegend.className = "pm-glegend";
  glegend.innerHTML = Object.values(groups)
    .map(v => `<span><i style="background:${v.color}"></i>${v.label}</span>`).join("");
  root.appendChild(glegend);

  const svg = document.createElementNS(NS, "svg");
  svg.setAttribute("viewBox", "0 0 900 560");
  svg.setAttribute("class", "pm-svg");
  svg.setAttribute("role", "img");
  svg.setAttribute("aria-label", "Map of how the papers in the literature review relate to each other and to this project");
  const scroll = document.createElement("div");
  scroll.className = "pm-scroll";
  scroll.appendChild(svg);
  root.appendChild(scroll);

  const panel = document.createElement("div");
  panel.className = "pm-panel";
  root.appendChild(panel);

  const gLinks = document.createElementNS(NS, "g");
  const gNodes = document.createElementNS(NS, "g");
  svg.append(gLinks, gNodes);

  links.forEach(([a, b, kind, text], i) => {
    const A = byId[a], B = byId[b];
    const l = document.createElementNS(NS, "line");
    l.setAttribute("x1", A.x); l.setAttribute("y1", A.y);
    l.setAttribute("x2", B.x); l.setAttribute("y2", B.y);
    l.setAttribute("stroke", kinds[kind].color);
    l.setAttribute("stroke-dasharray", kinds[kind].dash);
    l.setAttribute("class", "pm-link");
    l.dataset.kind = kind; l.dataset.a = a; l.dataset.b = b; l.dataset.i = i;
    gLinks.appendChild(l);
  });

  nodes.forEach(n => {
    const g = document.createElementNS(NS, "g");
    g.setAttribute("class", "pm-node" + (n.g === "self" ? " pm-self" : ""));
    g.setAttribute("transform", `translate(${n.x - W / 2},${n.y - H / 2})`);
    g.setAttribute("tabindex", "0");
    g.setAttribute("role", "button");
    g.setAttribute("aria-label", `${n.label}: show its connections`);
    g.dataset.id = n.id;
    const r = document.createElementNS(NS, "rect");
    r.setAttribute("width", W); r.setAttribute("height", H); r.setAttribute("rx", 5);
    r.setAttribute("stroke", groups[n.g].color);
    const bar = document.createElementNS(NS, "rect");
    bar.setAttribute("width", 5); bar.setAttribute("height", H); bar.setAttribute("rx", 2);
    bar.setAttribute("fill", groups[n.g].color);
    const t1 = document.createElementNS(NS, "text");
    t1.setAttribute("x", 14); t1.setAttribute("y", 20); t1.setAttribute("class", "pm-t1");
    t1.textContent = n.label;
    const t2 = document.createElementNS(NS, "text");
    t2.setAttribute("x", 14); t2.setAttribute("y", 36); t2.setAttribute("class", "pm-t2");
    t2.textContent = n.sub;
    g.append(r, bar, t1, t2);
    g.addEventListener("click", () => select(n.id));
    g.addEventListener("keydown", e => { if (e.key === "Enter" || e.key === " ") { e.preventDefault(); select(n.id); } });
    gNodes.appendChild(g);
  });

  let current = null;
  function select(id) {
    current = current === id ? null : id;
    svg.classList.toggle("pm-focus", !!current);
    svg.querySelectorAll(".pm-node").forEach(g => {
      const lit = !current || g.dataset.id === current ||
        links.some(([a, b]) => (a === current && b === g.dataset.id) || (b === current && a === g.dataset.id));
      g.classList.toggle("pm-dim", !lit);
      g.classList.toggle("pm-sel", g.dataset.id === current);
    });
    svg.querySelectorAll(".pm-link").forEach(l => {
      l.classList.toggle("pm-hi", !!current && (l.dataset.a === current || l.dataset.b === current));
    });
    render();
  }

  function render() {
    if (!current) {
      panel.innerHTML = `<p class="pm-hint">Click a paper to see how it connects to the others. Click it again to show everything. The legend buttons switch each kind of link on and off.</p>`;
      return;
    }
    const n = byId[current];
    const mine = links.filter(([a, b]) => a === current || b === current);
    const head = n.href
      ? `<a href="${n.href.replace(/\.qmd$/, ".html")}">${n.label}</a>`
      : n.label;
    panel.innerHTML = `<h4>${head} <span>${groups[n.g].label}</span></h4><ul>` +
      mine.map(([a, b, kind, text]) => {
        const o = byId[a === current ? b : a];
        const name = o.href ? `<a href="${o.href.replace(/\.qmd$/, ".html")}">${o.label}</a>` : o.label;
        return `<li><span class="pm-k" style="--k:${kinds[kind].color}">${kinds[kind].label}</span><b>${name}.</b> ${text}</li>`;
      }).join("") + "</ul>";
  }
  render();
})();
