/* Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
   Científicas (CSIC). Author: David Olmeda Casadomé.
   This file is part of lisaR, free software under the GNU General Public
   License version 3 (GPL-3). See the DESCRIPTION file and
   <https://www.gnu.org/licenses/gpl-3.0.html>. */

(function () {
  'use strict';
  const data = JSON.parse(document.getElementById('gene-evidence-data').textContent);
  const el = id => document.getElementById(id);
  const columns = (names, values) => Object.fromEntries(names.map((name, i) => [name, values[i]]));
  const de = columns(data.de_columns, data.de);
  const sets = columns(data.set_columns, data.sets);
  const symbols = data.symbols;
  const symbolIndex = new Map(symbols.map((symbol, i) => [symbol, i]));
  const idIndex = new Map();
  const deIndex = new Map();
  const setIndex = new Map();
  const analyses = new Map(data.analyses.map(row => [row.analysis_id, row.positive_contrast]));
  for (let i = 0; i < de.symbol.length; i++) {
    const gene = de.symbol[i];
    if (!deIndex.has(gene)) deIndex.set(gene, []);
    deIndex.get(gene).push(i);
    const id = de.gene_id[i];
    if (id && !['rowname', 'rownames', '.rownames'].includes((de.gene_id_source[i] || '').toLowerCase())) {
      if (!idIndex.has(id)) idIndex.set(id, new Set());
      idIndex.get(id).add(gene);
    }
  }
  for (let i = 0; i < sets.pathway.length; i++) {
    const pathway = sets.pathway[i];
    if (!setIndex.has(pathway)) setIndex.set(pathway, []);
    setIndex.get(pathway).push(i);
  }
  const state = { gene: null, rows: [], deRows: [], page: 0 };
  const pageSize = 100;
  const finite = x => typeof x === 'number' && Number.isFinite(x);
  const format = x => finite(x) ? (x === 0 ? '0' : (Math.abs(x) < 0.001 ? x.toExponential(3) : Number(x.toPrecision(5)).toString())) : 'unavailable';
  const signedClass = x => !finite(x) ? 'missing' : x > 0 ? 'positive' : x < 0 ? 'negative' : 'neutral';
  const create = (tag, text, className) => {
    const node = document.createElement(tag);
    if (text !== undefined) node.textContent = text;
    if (className) node.className = className;
    return node;
  };
  const scalarText = value => value === null || value === undefined ? 'NA' : String(value);
  const download = (name, rows, keys) => {
    const escape = value => scalarText(value).replace(/[\t\r\n]/g, ' ');
    const text = [keys.join('\t'), ...rows.map(row => keys.map(key => escape(row[key])).join('\t'))].join('\n') + '\n';
    const url = URL.createObjectURL(new Blob([text], { type: 'text/tab-separated-values;charset=utf-8' }));
    const a = create('a'); a.href = url; a.download = name; a.click(); setTimeout(() => URL.revokeObjectURL(url), 1000);
  };
  const baseTable = (container, headings) => {
    container.replaceChildren();
    const table = create('table'); const head = create('thead'); const row = create('tr');
    headings.forEach(label => { const th = create('th', label); th.scope = 'col'; row.append(th); });
    head.append(row); table.append(head); const body = create('tbody'); table.append(body); container.append(table); return body;
  };
  el('scope').textContent = `Dictionary tier: ${data.metadata.tier} · GSEA FDR cutoff: ${data.metadata.gsea_padj_cutoff} · DE FDR cutoff: ${data.metadata.de_padj_cutoff}`;
  [...analyses.keys()].forEach(id => { const option = create('option', id); option.value = id; el('analysis-filter').append(option); });
  [...new Set(data.scopes.map(row => row.collection))].sort().forEach(id => { const option = create('option', id); option.value = id; el('collection-filter').append(option); });

  // URL scope helpers: only fixed local report routes are constructed here.
  function resolveGeneScope(query) {
    const analysis = query.get('analysis_id') || '', collection = query.get('collection') || '';
    if (query.has('tier') && query.get('tier') !== String(data.metadata.tier))
      return {error: 'This link requests a different dictionary tier. No gene was selected.'};
    if (analysis && !data.analyses.some(row => row.analysis_id === analysis))
      return {error: 'The requested analysis is not included in this gene report. No gene was selected.'};
    if (collection && !data.scopes.some(row => row.collection === collection))
      return {error: 'The requested collection is not included in this gene report. No gene was selected.'};
    if (analysis && collection && !data.scopes.some(row => row.analysis_id === analysis && row.collection === collection))
      return {error: 'The requested analysis and collection are not an available scope. No gene was selected.'};
    return {analysis, collection, includeAll: ['1', 'true'].includes(query.get('include_all')),
      leadingOnly: ['1', 'true'].includes(query.get('leading_only')), error: ''};
  }
  function resolveContrastOrigin(query) {
    // Older category/gene links without a return route remain valid.
    if (!['return_contrast_scope', 'return_analysis_a', 'return_analysis_b', 'return_side',
      'return_analysis_id', 'return_collection', 'return_tier'].some(key => query.has(key)))
      return {origin: null, error: ''};
    // Legacy links share the initial view scope. Once filters change, pin the
    // return scope independently so the original A/B side survives a reload.
    const pinned = ['return_analysis_id', 'return_collection', 'return_tier'].some(key => query.has(key));
    const origin = {scope: query.get('return_contrast_scope'), contrast_id: query.get('return_contrast_id'),
      analysis_a: query.get('return_analysis_a'), analysis_b: query.get('return_analysis_b'),
      category_id: query.get('return_category_id'), side: query.get('return_side'),
      analysis_id: query.get(pinned ? 'return_analysis_id' : 'analysis_id'),
      collection: query.get(pinned ? 'return_collection' : 'collection'),
      tier: query.get(pinned ? 'return_tier' : 'tier')};
    const segment = value => typeof value === 'string' && value.length > 0 &&
      !['.', '..'].includes(value) && !/[\/\\\u0000-\u001f\u007f]/.test(value);
    const valid = segment(origin.scope) && segment(origin.collection) && origin.contrast_id && origin.category_id &&
      ['a', 'b'].includes(origin.side) && origin.tier === String(data.metadata.tier) &&
      origin.analysis_a !== origin.analysis_b && origin.analysis_id === origin['analysis_' + origin.side] &&
      ['a', 'b'].every(side => data.scopes.some(row => row.analysis_id === origin['analysis_' + side] && row.collection === origin.collection));
    return valid ? {origin, error: ''} : {origin: null,
      error: 'The originating contrast scope is incomplete or does not match this report. Clear it to search this report.'};
  }
  function contrastReturnHref(origin, gene) {
    const query = new URLSearchParams({contrast_id: origin.contrast_id, analysis_a: origin.analysis_a,
      analysis_b: origin.analysis_b, collection: origin.collection, tier: origin.tier, category_id: origin.category_id});
    if (gene) query.set('gene', gene);
    return `../contrast_evidence/${encodeURIComponent(origin.scope)}/${encodeURIComponent(origin.collection)}/index.html#${query}`;
  }
  function queryForGeneView(query, view, origin) {
    const next = new URLSearchParams(query);
    for (const [key, value] of Object.entries({gene: view.gene, analysis_id: view.analysis,
      collection: view.collection, tier: String(data.metadata.tier),
      include_all: view.includeAll ? '1' : '', leading_only: view.leadingOnly ? '1' : ''})) {
      if (value) next.set(key, value); else next.delete(key);
    }
    if (origin) {
      for (const [key, value] of Object.entries({return_contrast_scope: origin.scope,
        return_contrast_id: origin.contrast_id, return_category_id: origin.category_id,
        return_analysis_a: origin.analysis_a, return_analysis_b: origin.analysis_b, return_side: origin.side,
        return_analysis_id: origin.analysis_id, return_collection: origin.collection, return_tier: origin.tier}))
        next.set(key, value);
    }
    return next;
  }
  // End URL scope helpers.
  const incoming = new URLSearchParams(location.search);
  const initialScope = resolveGeneScope(incoming), returnScope = resolveContrastOrigin(incoming);
  const scopeError = initialScope.error || returnScope.error;
  function syncViewQuery() {
    if (scopeError) return;
    const query = queryForGeneView(location.search, {gene: state.gene === null ? '' : symbols[state.gene],
      analysis: el('analysis-filter').value, collection: el('collection-filter').value,
      includeAll: el('include-all').checked, leadingOnly: el('leading-only').checked}, returnScope.origin);
    try { history.replaceState(null, '', `${location.pathname}?${query}${location.hash || ''}`); }
    catch (_) { /* file hosts may disallow history; selection still works */ }
  }
  function syncContrastReturn() {
    const container = el('contrast-return'), link = el('contrast-return-link');
    const origin = returnScope.origin;
    container.hidden = !!scopeError || !origin || state.gene === null;
    if (container.hidden) { link.removeAttribute('href'); link.textContent = ''; return; }
    link.href = contrastReturnHref(origin, symbols[state.gene]);
    link.textContent = `Return to contrast category ${origin.category_id} (from ${origin.side.toUpperCase()}: ${origin['analysis_' + origin.side]})`;
  }

  function leadingStatus(index) {
    if (sets.leading_edge_state[index] === 'unavailable') return 'unavailable';
    return data.leading_edges[index].includes(state.gene) ? 'yes' : 'no';
  }
  function setRow(index) {
    const row = Object.fromEntries(data.set_columns.map(key => [key, sets[key][index]]));
    row.pathway = data.pathways[row.pathway]; row.symbol = symbols[state.gene]; row.tier = data.metadata.tier;
    row.leading_edge = leadingStatus(index); return row;
  }
  function showSets() {
    const body = baseTable(el('set-table'), ['Analysis', 'Collection', 'Gene set', 'NES', 'GSEA FDR', 'Significant', 'Leading edge', 'LISA categories / evidence', 'Evaluation state']);
    const visible = state.rows.slice(state.page * pageSize, (state.page + 1) * pageSize);
    visible.forEach(index => {
      const row = setRow(index); const tr = create('tr');
      tr.append(create('td', row.analysis_id), create('td', row.collection), create('td', row.pathway, 'identifier'));
      tr.append(create('td', format(row.NES), `numeric ${signedClass(row.NES)}`), create('td', format(row.gsea_fdr), 'numeric'));
      tr.append(create('td', row.evaluable ? (row.significant ? 'yes' : 'no') : 'unavailable'), create('td', row.leading_edge));
      const category = create('td');
      if (row.category_ids) row.category_ids.split(';').forEach(id => {
        const link = create('a', id, 'category-link');
        const query = new URLSearchParams({ category: id, gene: row.symbol, tier: data.metadata.tier, analysis_id: row.analysis_id, collection: row.collection });
        link.href = `../evidence/${encodeURIComponent(row.analysis_id)}/${encodeURIComponent(row.collection)}/index.html?${query}`;
        category.append(link);
      }); else category.textContent = 'unclassified (table only)';
      tr.append(category, create('td', row.evaluable ? 'evaluated' : row.eligibility_status)); body.append(tr);
    });
    if (!visible.length) { const row = create('tr'); const cell = create('td', 'No member gene sets match these filters. This does not imply no biological effect.', 'empty'); cell.colSpan = 9; row.append(cell); body.append(row); }
    const total = state.rows.length; const pages = Math.max(1, Math.ceil(total / pageSize));
    el('page-state').textContent = `Page ${state.page + 1} of ${pages} · ${visible.length} of ${total} matching rows shown`;
    el('previous').disabled = state.page === 0; el('next').disabled = state.page + 1 >= pages;
    el('set-count').textContent = `${total} analysis × collection × member-set records. Order: significant first, then leading-edge yes / no / unavailable, GSEA FDR, absolute NES and exact IDs. Downloads include all matching rows.`;
  }
  function syncShell() {
    syncContrastReturn();
    if (!window.lisaReportShell) return;
    const id = el('analysis-filter').value, collection = el('collection-filter').value;
    const record = data.analyses.find(x => x.analysis_id === id);
    window.lisaReportShell.updateContext({selection: state.gene === null ? '' : symbols[state.gene],
      analysis: id ? (record?.display_title || record?.label || id) : 'All individual analyses', collection,
      direction: id ? analyses.get(id) : '',
      exact_ids: {analysis_id: id, collection, tier: data.metadata.tier, gene: state.gene === null ? '' : symbols[state.gene]}});
  }
  document.addEventListener('lisa:shell-ready', syncShell);
  function filterSets() {
    syncViewQuery();
    if (state.gene === null) { syncShell(); return; }
    const member = data.member_sets[state.gene] || [];
    const all = member.flatMap(index => setIndex.get(index) || []);
    const analysis = el('analysis-filter').value; const collection = el('collection-filter').value;
    state.rows = all.filter(index => (!analysis || sets.analysis_id[index] === analysis) && (!collection || sets.collection[index] === collection) &&
      (el('include-all').checked || sets.significant[index]) && (!el('leading-only').checked || leadingStatus(index) === 'yes'));
    const leadOrder = { yes: 0, no: 1, unavailable: 2 };
    state.rows.sort((a, b) => Number(sets.significant[b]) - Number(sets.significant[a]) || leadOrder[leadingStatus(a)] - leadOrder[leadingStatus(b)] ||
      (finite(sets.gsea_fdr[a]) ? sets.gsea_fdr[a] : Infinity) - (finite(sets.gsea_fdr[b]) ? sets.gsea_fdr[b] : Infinity) ||
      Math.abs(sets.NES[b] || 0) - Math.abs(sets.NES[a] || 0) || sets.analysis_id[a].localeCompare(sets.analysis_id[b]) ||
      sets.collection[a].localeCompare(sets.collection[b]) || data.pathways[sets.pathway[a]].localeCompare(data.pathways[sets.pathway[b]]));
    state.page = 0; showSets(); syncShell();
  }
  function selectGene(index) {
    state.gene = index; el('evidence').hidden = false; el('matches').replaceChildren();
    const symbol = symbols[index]; el('gene-title').textContent = symbol; el('gene-query').value = symbol;
    state.deRows = (deIndex.get(index) || []).map(i => ({...Object.fromEntries(data.de_columns.map(key => [key, key === 'symbol' ? symbol : de[key][i]])), positive_contrast: analyses.get(de.analysis_id[i])}));
    const identifiers = [...new Set(state.deRows.filter(row => !['rowname','rownames','.rownames'].includes((row.gene_id_source || '').toLowerCase())).map(row => row.gene_id).filter(Boolean))];
    el('identifier-summary').textContent = `Exact symbol: ${symbol} · Recorded biological identifiers: ${identifiers.length ? identifiers.join(', ') : 'none available beyond symbol'}. Input row indexes are labelled separately and are not searchable gene aliases. All distinct DE rows are retained.`;
    el('search-state').textContent = `Selected exact symbol ${symbol}.`;
    const body = baseTable(el('de-table'), ['Analysis', 'Input identifier / row ID', 'Positive contrast means', 'log2 fold change', 'Rank value', 'Test statistic', 'DE P value', 'DE FDR', 'DE significance']);
    [...analyses.keys()].forEach(analysis => {
      const rows = state.deRows.filter(row => row.analysis_id === analysis);
      if (!rows.length) rows.push({ analysis_id: analysis, gene_id: '', missing: true });
      rows.forEach(row => {
        const rowId = ['rowname','rownames','.rownames'].includes((row.gene_id_source || '').toLowerCase());
        const idText = row.gene_id ? `${rowId ? 'row index: ' : ''}${row.gene_id} (source: ${row.gene_id_source || 'not recorded'})` : 'unavailable';
        const tr = create('tr'); tr.append(create('td', analysis), create('td', idText), create('td', analyses.get(analysis)));
        for (const key of ['log2FC', 'rank_value', 'statistic', 'pvalue', 'de_fdr']) tr.append(create('td', format(row[key]), `numeric ${key === 'log2FC' ? signedClass(row[key]) : ''}`));
        tr.append(create('td', finite(row.de_fdr) && row.de_fdr >= 0 && row.de_fdr <= 1 ? (row.de_fdr <= data.metadata.de_padj_cutoff ? 'yes' : 'no') : 'unavailable')); body.append(tr);
      });
    });
    el('de-state').textContent = 'Red = positive gene effect; blue = negative; grey = unavailable. A ranking is shown as a test statistic only when the recorded source-column mapping explicitly identifies it as such. Rows missing from standardised DE remain unavailable even when the gene belongs to a significant set.';
    filterSets();
  }
  function search(query) {
    query = query.trim(); el('matches').replaceChildren(); el('evidence').hidden = true; state.gene = null; syncShell();
    if (scopeError) { el('search-state').textContent = 'Clear the incoming scope above before searching this report.'; return; }
    if (!query) { el('search-state').textContent = 'Enter a recorded symbol or identifier.'; return; }
    const exact = new Set(idIndex.get(query) || []); if (symbolIndex.has(query)) exact.add(symbolIndex.get(query));
    if (exact.size === 1) { selectGene([...exact][0]); return; }
    let candidates = [...exact];
    if (exact.size) el('search-state').textContent = `Ambiguous exact query: ${exact.size} symbols match. Choose the intended symbol; evidence is not merged.`;
    else {
      const lower = query.toLowerCase(); const suggestions = new Set();
      symbols.forEach((symbol, i) => { if (symbol.toLowerCase().includes(lower)) suggestions.add(i); });
      idIndex.forEach((indices, id) => { if (id.toLowerCase().includes(lower)) indices.forEach(index => suggestions.add(index)); });
      candidates = [...suggestions];
      el('search-state').textContent = candidates.length ? `No exact match. ${candidates.length} case-insensitive / partial suggestions; choose explicitly (up to 50 displayed).` : 'No recorded exact symbol, identifier or partial match. No aliases have been inferred.';
    }
    candidates.sort((a, b) => symbols[a].localeCompare(symbols[b])).slice(0, 50).forEach(index => {
      const button = create('button', symbols[index], 'match'); button.type = 'button'; button.addEventListener('click', () => selectGene(index)); el('matches').append(button);
    });
  }
  el('gene-search').addEventListener('submit', event => { event.preventDefault(); search(el('gene-query').value); });
  ['analysis-filter', 'collection-filter', 'include-all', 'leading-only'].forEach(id => el(id).addEventListener('change', () => {
    const analysis = el('analysis-filter').value, collection = el('collection-filter').value;
    // A newly chosen scope takes precedence over an incompatible prior filter.
    if (analysis && collection && !data.scopes.some(row => row.analysis_id === analysis && row.collection === collection)) {
      if (id === 'analysis-filter') el('collection-filter').value = '';
      if (id === 'collection-filter') el('analysis-filter').value = '';
    }
    filterSets();
  }));
  el('previous').addEventListener('click', () => { if (state.page > 0) { state.page--; showSets(); } });
  el('next').addEventListener('click', () => { if ((state.page + 1) * pageSize < state.rows.length) { state.page++; showSets(); } });
  el('de-download').addEventListener('click', () => download(`${symbols[state.gene]}_DE.tsv`, state.deRows, [...data.de_columns, 'positive_contrast']));
  el('sets-download').addEventListener('click', () => download(`${symbols[state.gene]}_member_sets.tsv`, state.rows.map(setRow), ['symbol', 'analysis_id', 'collection', 'tier', 'pathway', 'NES', 'gsea_pvalue', 'gsea_fdr', 'evaluable', 'significant', 'leading_edge', 'category_ids', 'classification_status', 'eligibility_status']));
  if (scopeError) {
    el('scope-error').textContent = scopeError; el('scope-error').hidden = false; el('reset-scope').hidden = false;
    el('evidence').hidden = true; el('gene-query').value = ''; el('matches').replaceChildren();
  } else {
    el('analysis-filter').value = initialScope.analysis;
    el('collection-filter').value = initialScope.collection;
    el('include-all').checked = initialScope.includeAll;
    el('leading-only').checked = initialScope.leadingOnly;
    const query = incoming.get('gene'); if (query) { el('gene-query').value = query; search(query); }
  }
  syncShell();
})();
