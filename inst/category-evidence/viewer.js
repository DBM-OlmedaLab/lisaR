/* Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
   Científicas (CSIC). Author: David Olmeda Casadomé.
   This file is part of lisaR, free software under the GNU General Public
   License version 3 (GPL-3). See the DESCRIPTION file and
   <https://www.gnu.org/licenses/gpl-3.0.html>. */

(function () {
  'use strict';
  const D = JSON.parse(document.getElementById('evidence-data').textContent);
  const $ = id => document.getElementById(id);
  const arr = x => Array.isArray(x) ? x : (x == null ? [] : [x]);
  const finite = x => typeof x === 'number' && Number.isFinite(x);
  const fmt = x => !finite(x) ? 'NA' : (Math.abs(x) > 0 && Math.abs(x) < .001 ? x.toExponential(2) : Number(x.toPrecision(3)).toString());
  const el = (tag, text, cls) => { const n = document.createElement(tag); if (text != null) n.textContent = String(text); if (cls) n.className = cls; return n; };
  const unique = x => [...new Set(x)];
  const lexical = (a, b) => a < b ? -1 : a > b ? 1 : 0;
  const ascending = (a, b) => finite(a) ? (finite(b) ? a - b : -1) : finite(b) ? 1 : 0;
  const descendingAbs = (a, b) => ascending(finite(a) ? -Math.abs(a) : null, finite(b) ? -Math.abs(b) : null);
  const fdrOrder = (a, b) => ascending(a.gsea_fdr, b.gsea_fdr) || descendingAbs(a.NES, b.NES) || lexical(a.pathway, b.pathway);
  const sortedSets = x => x.slice().sort((a, b) => {
    const mode = $('set-sort').value;
    if (mode === 'name') return lexical(a.pathway, b.pathway);
    if (mode === 'abs_nes') return descendingAbs(a.NES, b.NES) || fdrOrder(a, b);
    if (mode === 'direction') return ({positive: 0, zero: 1, negative: 2, not_evaluable: 3}[a.direction] - {positive: 0, zero: 1, negative: 2, not_evaluable: 3}[b.direction]) || fdrOrder(a, b);
    return fdrOrder(a, b);
  });
  const positive = '#b2182b', negative = '#2166ac', neutral = '#64748b', missing = '#d1d5db';
  const directionColor = x => !finite(x) ? missing : x > 0 ? positive : x < 0 ? negative : neutral;
  const effectColor = (x, max) => {
    if (!finite(x)) return missing;
    const t = Math.min(1, Math.abs(x) / max), rgb = x >= 0 ? [178, 24, 43] : [33, 102, 172];
    return `rgb(${rgb.map(v => Math.round(255 + (v - 255) * t)).join(',')})`;
  };
  const genes = new Map(D.genes.map(g => [g.symbol, g]));
  const setsByCat = new Map(), edgesByCat = new Map();
  for (const s of D.sets) { if (!setsByCat.has(s.category_id)) setsByCat.set(s.category_id, []); setsByCat.get(s.category_id).push(s); }
  for (const e of D.leading_edges) { if (!edgesByCat.has(e.category_id)) edgesByCat.set(e.category_id, []); edgesByCat.get(e.category_id).push(e); }
  let category = '', setPage = 0, genePage = 0, tablePage = 0, selectedSet = '', selectedGene = '', visibleSets = [], visibleGenes = [];
  let currentTableRows = [], currentTableColumns = [], currentCandidateCount = 0, currentGeneCount = 0;
  const recurrenceByCat = new Map();
  for (const c of D.categories) {
    const sig = new Set((setsByCat.get(c.category_id) || []).filter(s => s.significant).map(s => s.pathway));
    const counts = new Map();
    for (const e of edgesByCat.get(c.category_id) || []) if (sig.has(e.pathway)) counts.set(e.symbol, (counts.get(e.symbol) || 0) + 1);
    recurrenceByCat.set(c.category_id, counts);
  }
  const categoryData = () => D.categories.find(c => c.category_id === category);
  const recurrence = symbol => recurrenceByCat.get(category)?.get(symbol) || 0;
  const categoryMean = () => categoryData()?.mean_significant_set_NES;
  const sortedGenes = symbols => symbols.slice().sort((a, b) => {
    const ga = genes.get(a) || {}, gb = genes.get(b) || {}, mode = $('gene-sort').value;
    if (mode === 'name') return lexical(a, b);
    if (mode === 'de_fdr') return ascending(ga.de_fdr, gb.de_fdr) || lexical(a, b);
    if (mode === 'abs_effect') return descendingAbs(ga.log2FC, gb.log2FC) || lexical(a, b);
    if (mode === 'aligned_effect') {
      const mean = categoryMean();
      if (finite(mean) && mean !== 0) return ascending(finite(ga.log2FC) ? -Math.sign(mean) * ga.log2FC : null, finite(gb.log2FC) ? -Math.sign(mean) * gb.log2FC : null) || lexical(a, b);
      return lexical(a, b);
    }
    return recurrence(b) - recurrence(a) || lexical(a, b);
  });
  function settings() {
    return {schema_version: '1.1', analysis_id: D.metadata.analysis_id, collection: D.metadata.collection, tier: D.metadata.tier,
      category_id: category, set_sort: $('set-sort').value, gene_sort: $('gene-sort').value, search: $('search').value,
      include_all: $('include-all').checked, set_page: setPage, gene_page: genePage, max_sets: maxSets, max_genes: maxGenes,
      selected_set: selectedSet, selected_gene: selectedGene, gsea_padj_cutoff: D.metadata.gsea_padj_cutoff,
      de_padj_cutoff: D.metadata.de_padj_cutoff, mean_significant_set_NES: categoryMean() ?? null,
      recurrence_scope: 'all_significant_unique_category_sets', recurrence_denominator: categoryData()?.n_significant_sets || 0,
      shown_sets: visibleSets.map(s => s.pathway), shown_genes: visibleGenes, matching_sets: currentCandidateCount, available_genes: currentGeneCount,
      set_sort_label: $('set-sort').selectedOptions[0].textContent, gene_sort_label: $('gene-sort').selectedOptions[0].textContent};
  }
  function syncShell() {
    if (!window.lisaReportShell) return;
    const c = categoryData();
    window.lisaReportShell.updateContext({selection: c ? c.category_display_name + (selectedGene ? ' / ' + selectedGene : '') : '',
      exact_ids: {analysis_id: D.metadata.analysis_id, collection: D.metadata.collection, tier: D.metadata.tier, category_id: category, gene: selectedGene}});
  }
  document.addEventListener('lisa:shell-ready', syncShell);
  function updateGeneLink() {
    if (!D.metadata.gene_explorer_href) { $('gene-explorer').hidden = true; return; }
    $('gene-explorer').hidden = false; const params = new URLSearchParams(); if (selectedGene) params.set('gene', selectedGene); if (D.metadata.analysis_id) params.set('analysis_id', D.metadata.analysis_id);
    $('gene-explorer').href = D.metadata.gene_explorer_href + (params.size ? '?' + params : '');
    $('gene-explorer').textContent = selectedGene ? `View all evidence for ${selectedGene}` : 'Search all evidence for a gene';
  }
  function updateLocation() {
    const params = new URLSearchParams(location.search); params.set('category', category); if (selectedGene) params.set('gene', selectedGene); else params.delete('gene');
    try { history.replaceState(null, '', '?' + params + '#category-' + encodeURIComponent(category)); } catch (_) { /* File viewers may disable history; interaction still works. */ }
    updateGeneLink(); syncShell();
  }
  const maxSets = D.metadata.max_sets, maxGenes = D.metadata.max_genes, rowsPerPage = 50;
  const edges = () => edgesByCat.get(category) || [];
  const sets = () => setsByCat.get(category) || [];
  const query = () => $('search').value.trim().toLowerCase();
  const match = value => String(value == null ? '' : value).toLowerCase().includes(query());
  function link(label, path) { const a = el('a', label); a.href = path; a.download = ''; return a; }
  function safeCategoryAsset(path) {
    if (typeof path !== 'string' || path.length === 0 ||
        /^[A-Za-z][A-Za-z0-9+.-]*:|^[\\/]|[\\\x00-\x1f\x7f]/.test(path)) return false;
    const parts = path.split('/');
    const safePart = part => {
      if (!part || part === '.' || part === '..') return false;
      let decoded;
      try { decoded = decodeURIComponent(part); } catch (_) { return false; }
      return decoded !== '.' && decoded !== '..' && !/[\/\\\x00-\x1f\x7f]/.test(decoded);
    };
    if (parts.every(safePart)) return true;
    const extensionArtifact = parts.length > 4 &&
      parts.slice(0, 3).every(part => part === '..') &&
      parts[3] === 'artifacts' && parts.slice(4).every(safePart);
    const runArtifact = parts.length > 5 &&
      parts.slice(0, 4).every(part => part === '..') &&
      (parts[4] === 'outputs' || parts[4] === 'artifacts') &&
      parts.slice(5).every(safePart);
    return extensionArtifact || runArtifact;
  }
  // Ported from inst/report_assets/lisa_shell.js (startPending/rescroll): a
  // fragment jump on this page lands in the wrong place because the product
  // and member-chart images are lazy and have zero height until they load, so
  // the target keeps moving after a one-shot scrollIntoView. Re-apply the
  // scroll on the next frame, on window "load" and on every document.body
  // resize, until the first user-initiated input or a bounded deadline.
  let pendingTarget = null, pendingDeadline = 0, pendingObserver = null, pendingPoller = null, interactionBound = false;
  function stopPending() {
    pendingTarget = null; pendingDeadline = 0;
    if (pendingObserver) { pendingObserver.disconnect(); pendingObserver = null; }
    if (pendingPoller) { clearInterval(pendingPoller); pendingPoller = null; }
  }
  function rescroll() {
    if (!pendingTarget || !pendingTarget.isConnected) { stopPending(); return; }
    pendingTarget.scrollIntoView({block: 'start', behavior: 'instant'});
    if (Date.now() >= pendingDeadline) stopPending();
  }
  function startPending(target) {
    stopPending();
    if (!target) return;
    pendingTarget = target; pendingDeadline = Date.now() + 2000;
    if (!interactionBound) {
      interactionBound = true;
      for (const type of ['wheel', 'keydown', 'touchstart'])
        window.addEventListener(type, () => { if (pendingTarget) stopPending(); }, {capture: true, passive: true});
    }
    requestAnimationFrame(rescroll);
    if (document.readyState !== 'complete')
      window.addEventListener('load', function onLoad() { window.removeEventListener('load', onLoad); if (pendingTarget) requestAnimationFrame(rescroll); });
    if (window.ResizeObserver) { pendingObserver = new ResizeObserver(() => { if (pendingTarget) rescroll(); }); pendingObserver.observe(document.body); }
    else pendingPoller = setInterval(() => { if (pendingTarget) rescroll(); else { clearInterval(pendingPoller); pendingPoller = null; } }, 200);
    setTimeout(() => { if (pendingTarget) stopPending(); }, 2000);
  }
  // The report generator injects a static announcement block near the top of
  // the page. Its per-category count is client-side state, so it is written
  // here; the block itself and its href stay plain HTML so the page still
  // announces its extras with JavaScript disabled.
  function updateProductsAnchor(count) {
    const anchor = $('evidence-products-anchor');
    if (!anchor) return;
    const current = anchor.querySelector('[data-lisa-products-current]');
    if (current) current.textContent = count
      ? `This category: ${count} saved product${count === 1 ? '' : 's'}.`
      : 'This category has no saved products.';
    const jump = anchor.querySelector('.evidence-products-jump');
    if (jump) jump.hidden = !count;
  }
  function bindProductsJump() {
    const jump = document.querySelector('.evidence-products-jump');
    if (!jump) return;
    // Deliberately not native fragment navigation: this page's hashchange
    // handler reads the hash as a category selector, so letting the browser
    // set #category-products would clear the selected category.
    jump.addEventListener('click', event => {
      const target = $('category-products');
      if (!target || target.hidden) return;
      event.preventDefault();
      startPending(target);
    });
  }
  // A PNG product without its SVG twin is flagged when the run requested SVG,
  // exactly like figure_card() in the main report; KEGG diagrams are raster
  // backgrounds and are labelled as an explicit exception.
  function svgStatus(target, product, assets) {
    const requested = arr(D.formats).map(f => String(f).toLowerCase()).includes('svg');
    const fmts = assets.map(a => String((a && a.format) || '').toLowerCase());
    if (!requested || !fmts.includes('png') || fmts.includes('svg')) return;
    const pill = document.createElement('span'); pill.className = 'svg-status';
    pill.textContent = /kegg_pathway_map|contrast_kegg_map|painted/i.test(String(product.product || '') + ' ' + assets.map(a => a.href).join(' ')) ? 'SVG not available for this product type' : 'SVG missing';
    target.append(pill);
  }
  function renderCategoryProducts() {
    const root = $('category-products'), cards = $('category-product-cards');
    cards.replaceChildren();
    // Optional attachment payload used by full-report assembly.  Do not
    // recover a category from a filename: an item must carry the active exact
    // ID and every local asset URL must pass the offline-path gate.
    const products = arr(D.category_products).filter(product => product && product.product !== 'kegg' && product.category_id === category &&
      Array.isArray(product.assets) && product.assets.some(asset => asset && safeCategoryAsset(asset.href)));
    root.hidden = !products.length;
    updateProductsAnchor(products.length);
    for (const product of products) {
      const card = el('article', null, 'category-product-card');
      card.append(el('h3', product.label || product.product || 'Saved product'));
      const downloads = el('div', null, 'downloads');
      const assets = product.assets.filter(asset => asset && safeCategoryAsset(asset.href));
      for (const asset of assets) {
        const a = link((/^r$/i.test(asset.format) ? 'R script' : String(asset.format || 'file').toUpperCase()), asset.href);
        if (typeof asset.source_name === 'string' && !/[\\/\x00-\x1f]/.test(asset.source_name)) {
          a.download = asset.source_name; a.title = asset.source_name;
        }
        downloads.append(a);
      }
      card.append(downloads);
      svgStatus(downloads, product, assets);
      const image = assets.find(asset => /^(png|svg)$/i.test(String(asset.format || '')));
      // The thumb container reserves its box through CSS aspect-ratio before
      // the lazy image loads, so the panel does not grow underneath a jump
      // that is already in flight.
      if (image) { const box = el('a', null, 'category-product-thumb'); box.href=image.href; box.target='_blank'; box.rel='noopener'; const img = el('img'); img.src = image.href; img.alt = product.label || product.product || 'Saved category product'; img.loading = 'lazy'; img.decoding = 'async'; box.append(img); card.append(box); }
      if (!image) { const pdf = assets.find(asset => /^pdf$/i.test(String(asset.format || ''))); if (pdf) { const preview = el('object'); preview.data = pdf.href; preview.type = 'application/pdf'; preview.style.width = '100%'; preview.style.height = '360px'; preview.append(link('Open PDF', pdf.href)); card.append(preview); } }
      cards.append(card);
    }
  }
  $('scope').textContent = `Analysis: ${D.metadata.analysis_id || 'not recorded'} · Collection: ${D.metadata.collection || 'not recorded'} · Tier: ${D.metadata.tier}`;
  $('contrast').textContent = `Positive NES / DE contrast: ${D.metadata.positive_contrast} · GSEA FDR ≤ ${D.metadata.gsea_padj_cutoff} · DE FDR ≤ ${D.metadata.de_padj_cutoff}`;
  $('provenance').textContent = `Exact-ID joins; schema ${D.metadata.schema_version}. Universe ledger ${D.metadata.universe_ledger_supplied ? 'supplied' : 'not supplied (mapped-universe completeness cannot be independently checked)'}. ${D.metadata.excluded_unclassified_rows} unclassified annotation rows excluded from category graphics. No cross-tier or longitudinal comparison is performed.`;
  for (const c of D.categories) { const o = el('option', `${c.category_display_name} [${c.category_id}]`); o.value = c.category_id; $('category').append(o); }
  for (const [name, path] of Object.entries(D.downloads)) $('downloads').append(link(name.replaceAll('_', ' ') + ' TSV', path));
  $('downloads').append(link('Figure index TSV', 'tables/figure_index.tsv'));
  $('provenance').append(' Overlap export: ' + D.metadata.overlap_export + '. ' + D.metadata.overlap_scope);
  if (D.recipes) $('downloads').append(link('R script', 'reproduce_category_evidence.R'));
  function selectSet(id) { selectedSet = selectedSet === id ? '' : id; selectedGene = ''; updateLocation(); renderMatrix(); renderTable(); }
  function selectGene(id) { selectedGene = selectedGene === id ? '' : id; selectedSet = ''; updateLocation(); renderMatrix(); renderTable(); }
  function renderCategory() {
    const c = D.categories.find(x => x.category_id === category); if (!c) return;
    $('category').value = category;
    syncShell();
    const info = $('category-info'); info.replaceChildren(el('h2', c.category_display_name));
    if (c.category_description) info.append(el('p', c.category_description));
    if (D.metadata.hommel_support_schema && c.analysis_id === D.metadata.analysis_id && c.collection === D.metadata.collection && !D.metadata.contrast_id) {
      const support = el('section', null, 'hommel-support-summary');
      support.dataset.analysisId = c.analysis_id; support.dataset.collection = c.collection; support.dataset.categoryId = c.category_id;
      support.append(el('h3', 'Category significance and support'), el('p', 'Adjusted category P and minimum supported enrichment', 'note'));
      const help = el('button', 'How to read category significance', 'hommel-help-button');
      help.type = 'button'; help.dataset.hommelHelp = ''; help.setAttribute('aria-haspopup', 'dialog'); help.setAttribute('aria-controls', 'hommel-support-help');
      const statistics = el('dl', null, 'hommel-support-statistics');
      const evaluated = c.n_sets_evaluable > 0;
      const values = [
        ['Minimum enriched / evaluable sets (d/N)', evaluated ? c.support_count : 'NE (N=0)'],
        ['Minimum supported percentage', evaluated ? '≥' + c.minimum_support_pct_text : 'Not evaluable'],
        ['Adjusted category P', finite(c.category_p_adjusted) ? String(c.category_p_adjusted) : 'Not evaluable'],
        ['Category significance', !evaluated ? 'Not evaluable' : c.significant ? 'Significant (adjusted P ≤ 0.05)' : 'Not significant (adjusted P > 0.05)'],
        ['Coverage: raw-P evaluable / total mapped', `${c.n_sets_evaluable} / ${c.n_sets_total}`],
        ['Missing raw P / ineligible or absent', `${c.n_sets_missing_p} / ${c.n_sets_not_eligible_or_absent}`]
      ];
      for (const [label, value] of values) { const cell = el('div'); cell.append(el('dt', label), el('dd', value)); statistics.append(cell); }
      support.append(statistics, el('p', c.support_label));
      const interpretation = !evaluated ? 'No eligible gene set has a usable P value, so this category could not be evaluated.' : c.minimum_enriched_sets > 0 ? `In this analysis, enrichment is supported in at least ${c.minimum_enriched_sets} of the ${c.n_sets_evaluable} evaluated gene sets. This is a minimum, not an exact count: more may be enriched. It does not measure effect size or show that the biological process is activated.` : `The data do not support a positive minimum among the ${c.n_sets_evaluable} evaluated gene sets. This does not mean that none are enriched.`;
      support.append(el('p', interpretation, 'hommel-biological-interpretation'), el('p', 'The adjusted category P tests whether any evaluated gene set is enriched; d/N is the minimum number supported at 95% simultaneous confidence. Their count and NES direction below describe a different part of the result.', 'note'), help);
      info.append(support);
    }
    info.append(el('p', `${c.macrogroup_name || c.macrogroup_id || 'Macrogroup not recorded'} · ${c.support_state.replaceAll('_', ' ')}`, 'note'));
    info.append(el('p', `Category mean NES across all significant unique member sets: ${fmt(c.mean_significant_set_NES)}. ${c.support_state === 'mixed_direction_support' ? 'Mixed NES directions: the mean does not imply uniform support; opposing sets and genes remain visible.' : 'This descriptive mean is not a category-level hypothesis test.'}`, 'note'));
    const counts = [['Mapped sets', c.n_mapped_sets], ['Evaluable / significant sets', `${c.n_evaluable_sets} / ${c.n_significant_sets}`], ['Significant +NES / −NES', `${c.n_significant_positive_sets} / ${c.n_significant_negative_sets}`], ['All evaluable +NES / −NES', `${c.n_positive_sets} / ${c.n_negative_sets}`], ['Unique LE genes (all / significant)', `${c.n_unique_le_genes} / ${c.n_significant_unique_le_genes}`], ['LE assignments (all / significant)', `${c.n_le_assignments} / ${c.n_significant_le_assignments}`]];
    const n = c.n_significant_sets, percent = value => n ? `${fmt(100 * value / n)}%` : 'NA';
    counts[2] = ['Positive / negative significant sets', `${c.n_significant_positive_sets} (${percent(c.n_significant_positive_sets)}) / ${c.n_significant_negative_sets} (${percent(c.n_significant_negative_sets)})`];
    const concordant = c.mean_significant_set_NES > 0 ? c.n_significant_positive_sets : c.mean_significant_set_NES < 0 ? c.n_significant_negative_sets : null;
    counts.push(['Directional concordance / significant sets', concordant === null ? 'Not directional' : `${percent(concordant)} (${concordant}/${n})`]);
    $('counts').replaceChildren(...counts.map(([label, value]) => { const n = el('div', null, 'count'); n.append(el('strong', value), el('span', label)); return n; }));
    renderCategoryProducts();
    const fi = D.figure_index.find(x => x.category_id === category); $('figure-downloads').replaceChildren();
    if (fi && arr(D.formats).length) {
      $('figure-downloads').append(el('span', 'Default capped view:', 'note'));
      for (const f of arr(D.formats)) $('figure-downloads').append(link(f.toUpperCase(), `figures/${fi.stem}.${f}`));
      if (D.source_data) $('figure-downloads').append(link('Plotted source TSV', `figures/${fi.stem}_source.tsv`));
    }
    const member = arr(D.member_figure_index).find(x => x.category_id === category);
    $('member-chart').replaceChildren(); $('member-chart-downloads').replaceChildren();
    if (member) {
      const imageFormat = arr(D.formats).find(x => x === 'png' || x === 'svg');
      if (imageFormat) { const img = el('img'); img.src = `figures/${member.stem}.${imageFormat}`; img.alt = `Gene-set NES and FDR for ${c.category_display_name}`; img.loading = 'lazy'; $('member-chart').append(img); }
      if (!imageFormat && arr(D.formats).includes('pdf')) { const preview = el('object'); preview.data = `figures/${member.stem}.pdf`; preview.type = 'application/pdf'; preview.style.width = '100%'; preview.style.height = '360px'; preview.append(link('Open PDF', preview.data)); $('member-chart').append(preview); }
      for (const f of arr(D.formats)) $('member-chart-downloads').append(link(f.toUpperCase(), `figures/${member.stem}.${f}`));
      if (D.source_data && member.source_tsv) $('member-chart-downloads').append(link('All member states and plotted source (TSV)', member.source_tsv));
      if (D.recipes && member.recipe_r) $('member-chart-downloads').append(link('R script', member.recipe_r));
    } else $('member-chart').append(el('p', 'No static formats were requested. All member states are available in the evidence table.'));
    renderMatrix(); renderTable();
  }
  function membershipMap() { const m = new Map(); for (const e of edges()) { if (!m.has(e.pathway)) m.set(e.pathway, new Set()); m.get(e.pathway).add(e.symbol); } return m; }
  function renderMatrix() {
    const mm = membershipMap();
    const candidates = sortedSets(sets().filter(s => ($('include-all').checked || s.significant) && (!query() || match(s.pathway) || [...(mm.get(s.pathway) || [])].some(match))));
    const pages = Math.max(1, Math.ceil(candidates.length / maxSets)); setPage = Math.min(setPage, pages - 1);
    visibleSets = candidates.slice(setPage * maxSets, (setPage + 1) * maxSets);
    let symbols = sortedGenes(unique(visibleSets.flatMap(s => [...(mm.get(s.pathway) || [])])));
    if (query() && symbols.some(match)) symbols = symbols.filter(match);
    const genePages = Math.max(1, Math.ceil(symbols.length / maxGenes)); genePage = Math.min(genePage, genePages - 1);
    visibleGenes = symbols.slice(genePage * maxGenes, (genePage + 1) * maxGenes);
    currentCandidateCount = candidates.length; currentGeneCount = symbols.length;
    const mean = categoryMean();
    $('active-order').textContent = `Set order: ${$('set-sort').selectedOptions[0].textContent}. Gene order: ${$('gene-sort').selectedOptions[0].textContent}. Recurrence denominator: all ${categoryData()?.n_significant_sets || 0} significant category sets, independent of this page. ` + ($('gene-sort').value === 'aligned_effect' ? (finite(mean) && mean !== 0 ? `Alignment uses sign of category mean NES (${fmt(mean)}); opposite effects are retained.` : 'Category mean NES is zero or unavailable: alignment is undefined; using exact symbol order, not a fabricated direction.') : '');
    updateGeneLink();
    $('set-page').textContent = `Sets ${candidates.length ? setPage * maxSets + 1 : 0}–${Math.min((setPage + 1) * maxSets, candidates.length)} / ${candidates.length}`;
    $('gene-page').textContent = `Genes ${symbols.length ? genePage * maxGenes + 1 : 0}–${Math.min((genePage + 1) * maxGenes, symbols.length)} / ${symbols.length}`;
    $('set-prev').disabled = setPage === 0; $('set-next').disabled = setPage >= pages - 1; $('gene-prev').disabled = genePage === 0; $('gene-next').disabled = genePage >= genePages - 1;
    $('matrix-state').textContent = candidates.length ? `Showing ${visibleSets.length} of ${candidates.length} ${$('include-all').checked ? 'mapped' : 'significant'} matching sets and ${visibleGenes.length} of ${symbols.length} leading-edge genes in this set page.${symbols.length ? '' : ' No leading-edge genes are recorded for these sets.'}` : `No ${$('include-all').checked ? 'mapped' : 'significant'} sets match this view. All set states remain in the tables below.`;
    $('selection').textContent = selectedSet ? `Selected set: ${selectedSet}; ${mm.get(selectedSet)?.size || 0} category LE genes` : selectedGene ? `Selected gene: ${selectedGene}; ${edges().filter(e => e.symbol === selectedGene).length} set connections` : 'Select a set or gene to highlight its connections.';
    const table = el('table'), thead = el('thead'), head = el('tr'); head.append(el('th', 'Exact set ID', 'set-label'));
    head.append(el('th', 'NES', 'stat'), el('th', 'GSEA FDR / state', 'stat'), el('th', 'LE shown / full', 'stat'));
    for (const symbol of visibleGenes) { const th = el('th', null, 'gene-label' + (selectedGene === symbol ? ' picked' : '')); const b = el('button', symbol); b.title = `Select gene ${symbol}`; b.onclick = () => selectGene(symbol); th.append(b); head.append(th); }
    thead.append(head);
    const de = visibleGenes.map(g => genes.get(g) || {}), maxFC = Math.max(1, ...de.filter(g => finite(g.log2FC)).map(g => Math.abs(g.log2FC)));
    for (const [label, field] of [['Gene DE log2FC', 'log2FC'], ['Gene DE FDR', 'de_fdr'], [`LE recurrence / ${categoryData()?.n_significant_sets || 0} significant sets`, 'recurrence']]) {
      const tr = el('tr'); const th = el('th', label, 'set-label track-label'); th.colSpan = 4; tr.append(th);
      for (let i = 0; i < visibleGenes.length; i++) { const g = de[i], td = el('td', field === 'recurrence' ? String(recurrence(visibleGenes[i])) : fmt(g[field]), 'track'); td.title = `${visibleGenes[i]} · gene_id: ${g.gene_id || 'not recorded'} · DE log2FC ${fmt(g.log2FC)} · FDR ${fmt(g.de_fdr)} · ${g.de_state || 'missing'}`;
        if (field === 'log2FC') { if (!finite(g.log2FC)) td.style.background = '#d1d5db'; else { const t = Math.abs(g.log2FC) / maxFC; td.style.background = effectColor(g.log2FC, maxFC); if (t > .65) td.style.color = 'white'; } } tr.append(td); }
      thead.append(tr);
    }
    table.append(thead); const body = el('tbody');
    for (const s of visibleSets) { const tr = el('tr'), title = el('th', null, 'set-label' + (selectedSet === s.pathway ? ' picked' : '')); const b = el('button', s.pathway); b.title = `${s.pathway} · ${s.state} · LE ${s.leading_edge_state}`; b.onclick = () => selectSet(s.pathway); title.append(b); tr.append(title);
      const nes = el('td', fmt(s.NES), 'stat'); nes.style.color = directionColor(s.NES); tr.append(nes, el('td', `${fmt(s.gsea_fdr)} · ${s.state.replaceAll('_', ' ')}`, 'stat')); const shownLE = visibleGenes.filter(g => mm.get(s.pathway)?.has(g)).length; tr.append(el('td', `${shownLE} / ${s.n_leading_edge_genes}`, 'stat'));
      for (const symbol of visibleGenes) { const member = mm.get(s.pathway)?.has(symbol); const td = el('td'); td.title = `${s.pathway} × ${symbol}: ${member ? 'GSEA leading-edge member' : 'not a recorded LE member'}`;
        if (member) { const dot = el('span', null, `member ${s.NES > 0 ? 'nes-pos' : s.NES < 0 ? 'nes-neg' : 'zero'}`); dot.setAttribute('aria-label', 'Leading-edge member'); if ((selectedSet && selectedSet !== s.pathway) || (selectedGene && selectedGene !== symbol)) dot.classList.add('faded'); td.append(dot); }
        if (selectedSet === s.pathway || selectedGene === symbol) td.classList.add('picked'); tr.append(td); }
      body.append(tr);
    }
    table.append(body); $('matrix').replaceChildren(table);
  }
  function overlapRows() {
    const mm = membershipMap(), pool = selectedSet ? sets() : visibleSets, pairs = [];
    const chosen = selectedSet ? pool.find(s => s.pathway === selectedSet) : null;
    const candidates = chosen ? pool.filter(s => s.pathway !== chosen.pathway).map(s => [chosen, s]) : pool.flatMap((s, i) => pool.slice(i + 1).map(t => [s, t]));
    for (const [aa, bb] of candidates) {
      const ga = mm.get(aa.pathway) || new Set(), gb = mm.get(bb.pathway) || new Set(); const intersection = [...ga].filter(g => gb.has(g)).length, union = ga.size + gb.size - intersection;
      pairs.push({pathway_a: aa.pathway, pathway_b: bb.pathway, direction_a: aa.direction, direction_b: bb.direction, direction_relation: !['positive', 'negative'].includes(aa.direction) || !['positive', 'negative'].includes(bb.direction) ? 'not_directional' : aa.direction === bb.direction ? 'same' : 'opposite', n_a: ga.size, n_b: gb.size, intersection_n: intersection, union_n: union, jaccard: ga.size && gb.size ? intersection / union : null, overlap_state: ga.size && gb.size ? 'defined' : 'empty_or_unrecorded_membership'});
    } return pairs;
  }
  function renderTable() {
    const kind = $('table-kind').value; let rows, columns;
    $('table-note').textContent = '';
    if (kind === 'sets') { const matchingPathways = new Set(query() ? edges().filter(e => match(e.symbol)).map(e => e.pathway) : []); rows = sortedSets(sets()).filter(s => !query() || match(s.pathway) || matchingPathways.has(s.pathway)); columns = ['pathway', 'NES', 'gsea_fdr', 'state', 'non_evaluable_reason', 'direction', 'n_leading_edge_genes', 'leading_edge_state', 'source_family', 'display_order']; }
    if (kind === 'genes') { rows = sortedGenes(unique(edges().map(e => e.symbol))).map(s => ({...(genes.get(s) || {symbol: s}), n_significant_set_connections: recurrence(s), n_significant_sets_denominator: categoryData()?.n_significant_sets || 0})).filter(g => !query() || match(g.symbol) || match(g.gene_id)); columns = ['symbol', 'gene_id', 'log2FC', 'de_fdr', 'de_state', 'n_significant_set_connections', 'n_significant_sets_denominator']; }
    if (kind === 'assignments') { rows = edges().filter(e => !query() || match(e.pathway) || match(e.symbol)); columns = ['pathway', 'symbol']; }
    if (kind === 'overlaps') { rows = overlapRows(); columns = ['pathway_a', 'pathway_b', 'direction_a', 'direction_b', 'direction_relation', 'n_a', 'n_b', 'intersection_n', 'union_n', 'jaccard', 'overlap_state']; $('table-note').textContent = selectedSet ? 'Scope: selected set versus every other mapped set in this category. Download matching table to retain every computed pair (no pagination cap).' : 'Scope: currently visible set page only. Select one set for its overlap against every mapped set; download matching rows here. The default exported pair file covers only the default displayed set selection. Empty/unrecorded membership gives NA.'; }
    currentTableRows = rows; currentTableColumns = columns;
    const pages = Math.max(1, Math.ceil(rows.length / rowsPerPage)); tablePage = Math.min(tablePage, pages - 1); const page = rows.slice(tablePage * rowsPerPage, (tablePage + 1) * rowsPerPage);
    const table = el('table'), head = el('tr'); for (const c of columns) head.append(el('th', c)); const thead = el('thead'); thead.append(head); table.append(thead); const tbody = el('tbody');
    for (const r of page) { const tr = el('tr'); for (const c of columns) { const td = el('td'); if ((c === 'pathway' || c === 'symbol') && r[c]) { const b = el('button', r[c]); b.onclick = () => c === 'pathway' ? selectSet(r[c]) : selectGene(r[c]); td.append(b); } else td.textContent = r[c] == null ? 'NA' : typeof r[c] === 'number' ? fmt(r[c]) : String(r[c]); tr.append(td); } tbody.append(tr); }
    table.append(tbody); $('detail-table').replaceChildren(table); $('table-count').textContent = `${rows.length} matching rows`; $('table-page').textContent = `Page ${tablePage + 1} / ${pages}`; $('table-prev').disabled = tablePage === 0; $('table-next').disabled = tablePage >= pages - 1;
  }
  $('table-download').onclick = () => {
    const metadata = {analysis_id: D.metadata.analysis_id, collection: D.metadata.collection, tier: D.metadata.tier, category_id: category,
      gsea_padj_cutoff: D.metadata.gsea_padj_cutoff, pair_scope: $('table-kind').value === 'overlaps' ? (selectedSet ? 'selected_set_vs_all_mapped_sets' : 'current_visible_set_page') : 'all_matching_category_rows'};
    const columns = [...Object.keys(metadata), ...currentTableColumns.filter(c => !Object.hasOwn(metadata, c))];
    const value = x => x == null ? 'NA' : /[\t\r\n"]/.test(String(x)) ? '"' + String(x).replaceAll('"', '""') + '"' : String(x);
    const lines = [columns.join('\t'), ...currentTableRows.map(row => columns.map(c => value(Object.hasOwn(metadata, c) ? metadata[c] : row[c])).join('\t'))];
    const url = URL.createObjectURL(new Blob([lines.join('\n') + '\n'], {type: 'text/tab-separated-values;charset=utf-8'}));
    const a = link('download', url); a.download = category.replace(/[^A-Za-z0-9._-]/g, '_') + '_' + $('table-kind').value + '.tsv'; document.body.append(a); a.click(); a.remove(); setTimeout(() => URL.revokeObjectURL(url), 1000);
  };
  function download(name, content, type) {
    const url = URL.createObjectURL(new Blob([content], {type})); const a = link('download', url); a.download = name; document.body.append(a); a.click(); a.remove(); setTimeout(() => URL.revokeObjectURL(url), 1000);
  }
  const fileStem = () => category.replace(/[^A-Za-z0-9._-]/g, '_') + '_current_view';
  const tsvValue = x => x == null ? 'NA' : /[\t\r\n"]/.test(String(x)) ? '"' + String(x).replaceAll('"', '""') + '"' : String(x);
  function sourceRows() {
    const c = categoryData(), state = settings(), mm = membershipMap();
    const base = {figure_id: 'category_evidence__' + category, figure_type: 'category_evidence', row_type: 'metadata',
      category_id: category, category_display_name: c.category_display_name, pathway: '', symbol: '', set_order: '', gene_order: '', NES: '', gsea_fdr: '', log2FC: '', de_fdr: '', de_state: '', leading_edge_state: '', n_leading_edge_genes: '', n_displayed_leading_edge_genes: '',
      analysis_id: D.metadata.analysis_id, collection: D.metadata.collection, tier: D.metadata.tier, positive_contrast: D.metadata.positive_contrast,
      gsea_padj_cutoff: D.metadata.gsea_padj_cutoff, de_padj_cutoff: D.metadata.de_padj_cutoff, n_mapped_sets: c.n_mapped_sets,
      n_evaluable_sets: c.n_evaluable_sets, n_significant_sets: c.n_significant_sets, n_plotted_sets: visibleSets.length,
      n_available_plot_genes: currentGeneCount, n_plotted_genes: visibleGenes.length, n_matching_sets: currentCandidateCount,
      n_significant_set_connections: '', n_significant_sets_denominator: c.n_significant_sets,
      mean_significant_set_NES: c.mean_significant_set_NES, set_sort_label: state.set_sort_label, gene_sort_label: state.gene_sort_label,
      selection_mode: state.include_all ? 'all mapped matching sets' : 'significant matching sets',
      selected_gene: selectedGene, selected_set: selectedSet};
    const rows = [base];
    visibleSets.forEach((s, i) => rows.push({...base, row_type: 'set', pathway: s.pathway, set_order: i + 1, NES: s.NES, gsea_fdr: s.gsea_fdr, leading_edge_state: s.leading_edge_state, n_leading_edge_genes: s.n_leading_edge_genes, n_displayed_leading_edge_genes: visibleGenes.filter(g => mm.get(s.pathway)?.has(g)).length}));
    visibleGenes.forEach((symbol, i) => { const g = genes.get(symbol) || {}; rows.push({...base, row_type: 'gene', symbol, gene_order: i + 1, log2FC: g.log2FC ?? null, de_fdr: g.de_fdr ?? null, de_state: g.de_state || 'not_in_standardized_DE', n_significant_set_connections: recurrence(symbol)}); });
    visibleSets.forEach((s, i) => visibleGenes.forEach((symbol, j) => { if (mm.get(s.pathway)?.has(symbol)) rows.push({...base, row_type: 'membership', pathway: s.pathway, symbol, set_order: i + 1, gene_order: j + 1}); }));
    return rows.map((r, i) => ({...r, source_row_order: i + 1}));
  }
  for (const [id, enabled] of [['view-svg', arr(D.formats).includes('svg')], ['view-source', D.source_data !== false], ['view-state', D.recipes !== false]]) { $(id).hidden = !enabled; $(id).disabled = !enabled; }
  $('view-state').onclick = () => download(fileStem() + '_settings.json', JSON.stringify(settings(), null, 2) + '\n', 'application/json');
  $('view-source').onclick = () => { const rows = sourceRows(), columns = Object.keys(rows[0]); download(fileStem() + '_source.tsv', [columns.join('\t'), ...rows.map(r => columns.map(k => tsvValue(r[k])).join('\t'))].join('\n') + '\n', 'text/tab-separated-values;charset=utf-8'); };
  $('view-svg').onclick = () => {
    const ns = 'http://www.w3.org/2000/svg', svg = document.createElementNS(ns, 'svg'), state = settings(), c = categoryData(), mm = membershipMap();
    const left = 430, cell = 34, startY = 370, rowH = 28, matrixW = Math.max(1, visibleGenes.length) * cell, annotationX = left + matrixW + 28;
    const width = Math.max(1100, annotationX + 340), height = Math.max(650, startY + Math.max(1, visibleSets.length) * rowH + 115);
    svg.setAttribute('xmlns', ns); svg.setAttribute('viewBox', `0 0 ${width} ${height}`); svg.setAttribute('width', width); svg.setAttribute('height', height);
    function shape(tag, attrs, text) { const n = document.createElementNS(ns, tag); for (const [k, v] of Object.entries(attrs || {})) n.setAttribute(k, v); if (text != null) n.textContent = text; svg.append(n); return n; }
    const text = (x, y, value, attrs = {}) => shape('text', {x, y, 'font-family': 'Arial,sans-serif', 'font-size': 12, fill: '#172033', ...attrs}, value);
    shape('rect', {width, height, fill: 'white'});
    shape('metadata', {}, JSON.stringify({settings: state, source_rows: sourceRows()}));
    text(24, 34, `${c.category_display_name} [${category}]`, {'font-size': 20, 'font-weight': 'bold'});
    text(24, 58, `${D.metadata.analysis_id} · ${D.metadata.collection} · ${D.metadata.tier}`);
    text(24, 80, `Positive contrast: ${D.metadata.positive_contrast} · GSEA FDR ≤ ${D.metadata.gsea_padj_cutoff}`);
    text(24, 102, `Mapped / evaluable / significant sets: ${c.n_mapped_sets} / ${c.n_evaluable_sets} / ${c.n_significant_sets}; category mean NES: ${fmt(c.mean_significant_set_NES)}`);
    text(24, 124, `Current view: ${visibleSets.length} / ${currentCandidateCount} matching sets; ${visibleGenes.length} / ${currentGeneCount} genes in this set page.`);
    text(24, 146, `Sets: ${state.set_sort_label}. Genes: ${state.gene_sort_label}.`);
    text(24, 168, c.support_state === 'mixed_direction_support' ? 'Mixed NES directions: opposing sets and genes are retained; the category mean is descriptive.' : 'Category mean and recurrence are descriptive, not additional significance tests.');
    const maxFC = Math.max(1, ...visibleGenes.map(symbol => genes.get(symbol)?.log2FC).filter(finite).map(Math.abs));
    text(left - 16, 202, 'Gene DE log2FC', {'text-anchor': 'end'}); text(left - 16, 230, 'Gene DE FDR', {'text-anchor': 'end'}); text(left - 16, 258, `LE recurrence / ${c.n_significant_sets} significant sets`, {'text-anchor': 'end'});
    visibleGenes.forEach((symbol, j) => { const g = genes.get(symbol) || {}, x = left + j * cell;
      shape('rect', {x, y: 185, width: cell - 1, height: 24, fill: effectColor(g.log2FC, maxFC)});
      text(x + cell / 2, 201, fmt(g.log2FC), {'text-anchor': 'middle', 'font-size': 10, fill: finite(g.log2FC) && Math.abs(g.log2FC) > maxFC * .65 ? 'white' : '#172033'});
      text(x + cell / 2, 230, fmt(g.de_fdr), {'text-anchor': 'middle', 'font-size': 10}); text(x + cell / 2, 258, recurrence(symbol), {'text-anchor': 'middle', 'font-size': 10});
      text(x + cell / 2, startY - 14, symbol, {transform: `rotate(-60 ${x + cell / 2} ${startY - 14})`, 'text-anchor': 'start', 'font-size': 11});
    });
    text(annotationX, startY - 36, 'Set NES', {'font-weight': 'bold'}); text(annotationX + 90, startY - 36, 'GSEA FDR', {'font-weight': 'bold'}); text(annotationX + 210, startY - 48, 'LE genes', {'font-weight': 'bold'}); text(annotationX + 210, startY - 30, 'shown / full', {'font-weight': 'bold'});
    visibleSets.forEach((s, i) => { const y = startY + i * rowH;
      shape('rect', {x: left, y, width: matrixW, height: rowH - 1, fill: i % 2 ? 'white' : '#f1f5f9'});
      const label = text(left - 14, y + 19, s.pathway.length > 56 ? s.pathway.slice(0, 53) + '...' : s.pathway, {'text-anchor': 'end', 'font-size': 11}); const title = document.createElementNS(ns, 'title'); title.textContent = s.pathway; label.append(title);
      text(annotationX, y + 19, fmt(s.NES), {fill: directionColor(s.NES)}); text(annotationX + 90, y + 19, fmt(s.gsea_fdr)); text(annotationX + 210, y + 19, `${visibleGenes.filter(g => mm.get(s.pathway)?.has(g)).length} / ${s.n_leading_edge_genes}`);
      visibleGenes.forEach((symbol, j) => { const x = left + j * cell;
        if (selectedSet === s.pathway || selectedGene === symbol) shape('rect', {x, y, width: cell, height: rowH - 1, fill: '#fef3c7'});
        if (mm.get(s.pathway)?.has(symbol)) shape('rect', {x: x + 11, y: y + 8, width: 12, height: 12, fill: directionColor(s.NES), opacity: (selectedSet && selectedSet !== s.pathway) || (selectedGene && selectedGene !== symbol) ? '.18' : '1'});
      });
    });
    const footer = height - 68; text(24, footer, 'Red = positive; blue = negative; grey = missing. Dots: set NES direction. Gene DE log2FC and FDR remain distinct.');
    text(24, footer + 22, 'Recurrence denominator: ALL significant category sets, independent of this page. A blank cell is not a recorded LE connection.');
    text(24, footer + 44, 'All exact identifiers, displayed memberships and active view settings are embedded as SVG metadata and available as TSV / JSON.');
    download(fileStem() + '.svg', '<?xml version="1.0" encoding="UTF-8"?>\n' + new XMLSerializer().serializeToString(svg), 'image/svg+xml;charset=utf-8');
  };

  $('category').onchange = () => { category = $('category').value; setPage = genePage = tablePage = 0; selectedSet = selectedGene = ''; updateLocation(); renderCategory(); };
  $('search').oninput = () => { setPage = genePage = tablePage = 0; renderMatrix(); renderTable(); };
  $('include-all').onchange = () => { setPage = genePage = tablePage = 0; renderMatrix(); renderTable(); };
  $('clear-selection').onclick = () => { selectedSet = selectedGene = ''; updateLocation(); renderMatrix(); renderTable(); };
  for (const id of ['set-sort', 'gene-sort']) $(id).onchange = () => { setPage = genePage = tablePage = 0; renderMatrix(); renderTable(); };
  $('reset-view').onclick = () => { $('set-sort').value = 'fdr'; $('gene-sort').value = 'recurrence'; $('search').value = ''; $('include-all').checked = false; selectedSet = selectedGene = ''; setPage = genePage = tablePage = 0; updateLocation(); renderCategory(); };
  $('set-prev').onclick = () => { setPage--; genePage = 0; renderMatrix(); renderTable(); }; $('set-next').onclick = () => { setPage++; genePage = 0; renderMatrix(); renderTable(); };
  $('gene-prev').onclick = () => { genePage--; renderMatrix(); }; $('gene-next').onclick = () => { genePage++; renderMatrix(); };
  $('table-kind').onchange = () => { tablePage = 0; renderTable(); }; $('table-prev').onclick = () => { tablePage--; renderTable(); }; $('table-next').onclick = () => { tablePage++; renderTable(); };
  function loadHash() {
    const params = new URLSearchParams(location.search);
    for (const key of ['analysis_id', 'collection', 'tier']) if (params.has(key) && params.get(key) !== String(D.metadata[key])) {
      $('scope-error').hidden = false; $('scope-error').className = 'scope-warning'; $('scope-error').textContent = `Evidence scope mismatch: requested ${key} does not match this page. No category is displayed under the wrong scope.`;
      $('matrix').replaceChildren(); return;
    }
    let requested = params.get('category') || '';
    if (!requested) try { requested = decodeURIComponent(location.hash.replace(/^#category-/, '')); } catch (_) {}
    if (requested && !D.categories.some(c => c.category_id === requested)) {
      $('scope-error').hidden = false; $('scope-error').className = 'scope-warning'; $('scope-error').textContent = 'The exact requested category does not occur in this analysis, collection and tier. No fallback category was substituted.'; return;
    }
    $('scope-error').hidden = true; category = requested || D.categories[0]?.category_id || ''; selectedSet = ''; selectedGene = params.get('gene') || ''; setPage = genePage = tablePage = 0;
    if (selectedGene) {
      const mm = membershipMap(), candidates = sortedSets(sets().filter(s => s.significant));
      const at = candidates.findIndex(s => mm.get(s.pathway)?.has(selectedGene));
      if (at >= 0) { setPage = Math.floor(at / maxSets); const page = candidates.slice(setPage * maxSets, (setPage + 1) * maxSets); const symbols = sortedGenes(unique(page.flatMap(s => [...(mm.get(s.pathway) || [])]))); genePage = Math.max(0, Math.floor(symbols.indexOf(selectedGene) / maxGenes)); }
    }
    renderCategory();
  }
  window.addEventListener('hashchange', () => { const params = new URLSearchParams(location.search); params.delete('category'); try { history.replaceState(null, '', '?' + params + location.hash); } catch (_) {} loadHash(); });
  bindProductsJump();
  if (!D.categories.length) { $('category-info').append(el('h3', 'No classified category assignments are available.')); $('matrix-state').textContent = 'No category graphics were generated. Complete source states remain in the downloadable tables.'; }
  else loadHash();
})();
