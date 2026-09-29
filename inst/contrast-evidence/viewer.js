/* Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
   Científicas (CSIC). Author: David Olmeda Casadomé.
   This file is part of lisaR, free software under the GNU General Public
   License version 3 (GPL-3). See the DESCRIPTION file and
   <https://www.gnu.org/licenses/gpl-3.0.html>. */

/* Offline paired evidence. No network, GSEA, differential tests or shared FDR. */
(() => {
  'use strict';
  const D = JSON.parse(document.getElementById('contrast-evidence-data').textContent), M = D.metadata;
  const el = id => document.getElementById(id), arr = x => Array.isArray(x) ? x : x == null ? [] : [x];
  const number = x => x === null || x === undefined || x === '' || x === 'NA' ? NaN : Number(x);
  const finite = x => Number.isFinite(number(x));
  function fmt(x) {
    if (!finite(x)) return 'not available';
    const n = number(x), magnitude = Math.abs(n);
    if (n === 0) return '0';
    if (magnitude < .01 || magnitude >= 1000) return n.toExponential(2).replace(/\.?0+e/, 'e').replace('e+', 'e');
    return Number(n.toPrecision(3)).toString();
  }
  const pct = x => finite(x) ? `${number(x).toFixed(1)}%` : 'not available';
  const stateText = x => ({non_evaluable:'Not evaluable',non_significant:'Does not pass FDR cutoff',significant:'Passes FDR cutoff',not_recorded:'Not recorded',recorded:'Recorded',empty:'Recorded empty',not_available:'Not available',empty_union:'Undefined: both leading edges empty',defined:'Available',not_in_standardized_DE:'Not in this DE table',both:'Significant in both',A_only:'Significant only in A',B_only:'Significant only in B',neither:'Significant in neither'}[x] || String(x == null ? 'Not available' : x).replace(/_/g, ' '));
  const color = x => !finite(x) ? '#9ca3af' : number(x) > 0 ? '#b2182b' : number(x) < 0 ? '#2166ac' : '#64748b';
  const bool = x => x === true || x === 'TRUE' || x === 'true';
  const sideLabel = side => M['analysis_label_'+side] || M['analysis_'+side];
  const byID = (a, b) => String(a) < String(b) ? -1 : String(a) > String(b) ? 1 : 0;
  const short = (s, n = 47) => String(s).length > n ? String(s).slice(0, n - 3) + '...' : String(s);
  const escTSV = x => x == null ? 'NA' : String(x).replace(/[\t\r\n]/g, ' ');
  const categories = new Map(D.categories.map(x => [String(x.category_id), x]));
  const genes = new Map(D.genes.map(x => [String(x.symbol), x]));
  const conservation = new Map(D.conservation.map(x => [String(x.pathway), x]));
  const setsByCategory = new Map(), le = {a: new Map(), b: new Map()};
  for (const s of D.sets) {const id = String(s.category_id); if (!setsByCategory.has(id)) setsByCategory.set(id, []); setsByCategory.get(id).push(s);}
  for (const e of D.leading_edges) {const p = String(e.pathway); if (!le[e.side].has(p)) le[e.side].set(p, new Set()); le[e.side].get(p).add(String(e.symbol));}
  const edges = (side, p) => le[side].get(String(p)) || new Set();
  let active = null, selectedGene = '', current = null;
  const svgNS = 'http://www.w3.org/2000/svg';
  function svg(tag, attrs = {}, text) {const n = document.createElementNS(svgNS, tag); for (const [k, v] of Object.entries(attrs)) n.setAttribute(k, String(v)); if (text !== undefined) n.textContent = text; return n;}
  function txt(root, x, y, text, attrs = {}) {const t = svg('text', {x, y, 'font-size':12, fill:'#183244', ...attrs}, text); root.appendChild(t); return t;}
  function rootSVG(width, height, label) {const n = svg('svg', {xmlns:svgNS, width, height, viewBox:`0 0 ${width} ${height}`, 'font-family':'system-ui, sans-serif', role:'img', 'aria-label':label}); n.appendChild(svg('rect', {width, height, fill:'#fff'})); return n;}
  function tooltip(n, text) {n.appendChild(svg('title', {}, text)); return n;}
  function point(root, side, x, y, value, size = 5, hollow = false) {const n = side === 'a' ? svg('circle', {cx:x,cy:y,r:size,fill:color(value)}) : svg('rect', {x:x-size, y:y-size, width:size*2, height:size*2, fill:color(value)}); n.setAttribute('data-side',side); if(hollow){n.setAttribute('fill','white');n.setAttribute('stroke',color(value));n.setAttribute('stroke-width','2');} root.appendChild(n); return n;}
  function empty(id, text) {el(id).replaceChildren(); const p = document.createElement('p'); p.className = 'empty'; p.textContent = text; el(id).appendChild(p);}
  function table(id, heads, rows) {const t = document.createElement('table'), h = document.createElement('thead'), tr = document.createElement('tr'); heads.forEach(x => {const th = document.createElement('th'); th.scope = 'col'; th.textContent = x; tr.appendChild(th);}); h.appendChild(tr); t.appendChild(h); const b = document.createElement('tbody'); rows.forEach(row => {const r = document.createElement('tr'); row.forEach(v => {const c = document.createElement('td'); if (v instanceof Node) c.appendChild(v); else c.textContent = v == null ? 'not available' : String(v); r.appendChild(c);}); b.appendChild(r);}); t.appendChild(b); el(id).replaceChildren(t);}
  function download(name, text, type = 'text/plain;charset=utf-8') {const url = URL.createObjectURL(new Blob([text], {type})); const a = document.createElement('a'); a.href = url; a.download = name; a.click(); setTimeout(() => URL.revokeObjectURL(url), 500);}
  function downloadSVG(id, name, includeCurrentState = false) {
    const original = el(id).querySelector('svg'); if (!original) return;
    const n = original.cloneNode(true);
    // Keep a current-view SVG auditable without changing the visible plot.
    if (includeCurrentState && active && current) n.appendChild(svg('metadata', {}, JSON.stringify({settings:settings(), source_tsv:source()})));
    download(name, new XMLSerializer().serializeToString(n), 'image/svg+xml');
  }
  function syncShell() {
    if (!window.lisaReportShell) return;
    const ids={contrast_id:M.contrast_id,analysis_a:M.analysis_a,analysis_b:M.analysis_b,collection:M.collection,tier:M.tier};
    if(active)ids.category_id=active;if(selectedGene)ids.gene=selectedGene;
    window.lisaReportShell.updateContext({analysis:M.contrast_label,collection:M.collection,
      selection:active?categories.get(active).category_display_name+(selectedGene?' / '+selectedGene:''):'All categories',
      direction:'A minus B (descriptive)',cutoff:`GSEA FDR <= ${M.gsea_padj_cutoff}`,exact_ids:ids});
  }
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
    const extensionArtifact=parts.length>4 && parts.slice(0,3).every(part=>part==='..') &&
      parts[3]==='artifacts' && parts.slice(4).every(safePart);
    const runArtifact=parts.length>5 && parts.slice(0,4).every(part=>part==='..') &&
      (parts[4]==='outputs'||parts[4]==='artifacts') && parts.slice(5).every(safePart);
    return extensionArtifact||runArtifact;
  }
  function localLink(label, href) {const a=document.createElement('a');a.href=href;a.download='';a.textContent=label;return a;}
  // Ported from inst/report_assets/lisa_shell.js (startPending/rescroll): the
  // saved-product images are lazy and have zero height until they load, so a
  // one-shot scrollIntoView lands in the wrong place. Re-apply the scroll on
  // the next frame, on window "load" and on every document.body resize, until
  // the first user-initiated input or a bounded deadline.
  let pendingTarget=null, pendingDeadline=0, pendingObserver=null, pendingPoller=null, interactionBound=false;
  function stopPending() {
    pendingTarget=null;pendingDeadline=0;
    if(pendingObserver){pendingObserver.disconnect();pendingObserver=null;}
    if(pendingPoller){clearInterval(pendingPoller);pendingPoller=null;}
  }
  function rescroll() {
    if(!pendingTarget||!pendingTarget.isConnected){stopPending();return;}
    pendingTarget.scrollIntoView({block:'start',behavior:'instant'});
    if(Date.now()>=pendingDeadline)stopPending();
  }
  function startPending(target) {
    stopPending();
    if(!target)return;
    pendingTarget=target;pendingDeadline=Date.now()+2000;
    if(!interactionBound){
      interactionBound=true;
      ['wheel','keydown','touchstart'].forEach(type=>window.addEventListener(type,()=>{if(pendingTarget)stopPending();},{capture:true,passive:true}));
    }
    requestAnimationFrame(rescroll);
    if(document.readyState!=='complete')
      window.addEventListener('load',function onLoad(){window.removeEventListener('load',onLoad);if(pendingTarget)requestAnimationFrame(rescroll);});
    if(window.ResizeObserver){pendingObserver=new ResizeObserver(()=>{if(pendingTarget)rescroll();});pendingObserver.observe(document.body);}
    else pendingPoller=setInterval(()=>{if(pendingTarget)rescroll();else{clearInterval(pendingPoller);pendingPoller=null;}},200);
    setTimeout(()=>{if(pendingTarget)stopPending();},2000);
  }
  // The report generator injects a static announcement block near the top of
  // the page; only its per-category count is client-side state.
  function updateProductsAnchor(count) {
    const anchor=el('evidence-products-anchor');
    if(!anchor)return;
    const current=anchor.querySelector('[data-lisa-products-current]');
    if(current)current.textContent=count?`This category: ${count} saved product${count===1?'':'s'}.`:'This category has no saved products.';
    const jump=anchor.querySelector('.evidence-products-jump');
    if(jump)jump.hidden=!count;
  }
  function bindProductsJump() {
    const jump=document.querySelector('.evidence-products-jump');
    if(!jump)return;
    // Deliberately not native fragment navigation: this page's hashchange
    // handler reads the hash as its scope/category selector, so letting the
    // browser set #category-products would reset the displayed category.
    jump.addEventListener('click',event=>{
      const target=el('category-products');
      if(!target||target.hidden)return;
      event.preventDefault();
      startPending(target);
    });
  }
  function svgStatus(target,product,assets){const requested=arr(D.formats).map(f=>String(f).toLowerCase()).includes('svg');const fmts=assets.map(a=>String((a&&a.format)||'').toLowerCase());if(!requested||!fmts.includes('png')||fmts.includes('svg'))return;const pill=document.createElement('span');pill.className='svg-status';pill.textContent=/kegg_pathway_map|contrast_kegg_map|painted/i.test(String(product.product||'')+' '+assets.map(a=>a.href).join(' '))?'SVG not available for this product type':'SVG missing';target.appendChild(pill);}
  function renderCategoryProducts() {
    const root=el('category-products'), cards=el('category-product-cards'); cards.replaceChildren();
    // Attachments are accepted only for the selected exact category ID.  The
    // report assembler owns global products; this page intentionally does not
    // infer category membership from a filename or render global galleries.
    const products=arr(D.category_products).filter(product=>product && product.product !== 'kegg' && String(product.category_id)===String(active) &&
      Array.isArray(product.assets) && product.assets.some(asset=>asset && safeCategoryAsset(asset.href)));
    root.hidden=!products.length;
    updateProductsAnchor(products.length);
    products.forEach(product=>{
      const card=document.createElement('article');card.className='category-product-card';
      const title=document.createElement('h3');title.textContent=product.label||product.product||'Saved product';card.appendChild(title);
      const downloads=document.createElement('div');downloads.className='downloads';
      const assets=product.assets.filter(asset=>asset && safeCategoryAsset(asset.href));
      assets.forEach(asset=>{const a=localLink((/^r$/i.test(asset.format) ? 'R script' : String(asset.format||'file').toUpperCase()),asset.href);
        if(typeof asset.source_name==='string' && !/[\\/\x00-\x1f]/.test(asset.source_name)){a.download=asset.source_name;a.title=asset.source_name;}
        downloads.appendChild(a);});card.appendChild(downloads);svgStatus(downloads,product,assets);
      const image=assets.find(asset=>/^(png|svg)$/i.test(String(asset.format||'')));
      // The thumb container reserves its box through CSS aspect-ratio before
      // the lazy image loads, so the panel cannot grow underneath an in-flight
      // jump to #category-products.
      if(image){const box=document.createElement('a');box.href=image.href;box.target='_blank';box.rel='noopener';box.className='category-product-thumb';const img=document.createElement('img');img.src=image.href;img.alt=product.label||product.product||'Saved category product';img.loading='lazy';img.decoding='async';box.appendChild(img);card.appendChild(box);}
      if(!image){const pdf=assets.find(asset=>/^pdf$/i.test(String(asset.format||'')));if(pdf){const preview=document.createElement('object');preview.data=pdf.href;preview.type='application/pdf';preview.style.width='100%';preview.style.height='360px';preview.appendChild(localLink('Open PDF',pdf.href));card.appendChild(preview);}}
      cards.appendChild(card);
    });
  }
  function hashFor(category) {const p = new URLSearchParams({contrast_id:String(M.contrast_id), analysis_a:String(M.analysis_a), analysis_b:String(M.analysis_b), collection:String(M.collection), tier:String(M.tier), category_id:String(category)}); if (selectedGene) p.set('gene', selectedGene); return '#' + p.toString();}
  function minFDR(s) {return Math.min(...['a','b'].map(side => bool(s['evaluable_'+side]) && finite(s['gsea_fdr_'+side]) ? number(s['gsea_fdr_'+side]) : Infinity));}
  function largest(s,key) {return Math.max(...['a','b'].map(side => finite(s[key+'_'+side]) ? Math.abs(number(s[key+'_'+side])) : -Infinity));}
  function makeView() {
    const full = setsByCategory.get(active)||[], mode = el('selection').value;
    let ss = full.filter(s => mode==='intersection' ? bool(s.significant_a)&&bool(s.significant_b) : mode==='evaluable' ? bool(s.evaluable_a)||bool(s.evaluable_b) : bool(s.significant_a)||bool(s.significant_b));
    const sorter = el('set-order').value;
    ss.sort((a,b) => (sorter==='fdr'?minFDR(a)-minFDR(b):sorter==='effect'?largest(b,'NES')-largest(a,'NES'):0)||byID(a.pathway,b.pathway));
    const n = ss.length; ss = ss.slice(0,number(M.max_sets));
    const recurrence = new Map();
    full.forEach(s => ['a','b'].forEach(side => {if(bool(s['significant_'+side])) edges(side,s.pathway).forEach(g => recurrence.set(g,(recurrence.get(g)||0)+1));}));
    const all = new Set(); ss.forEach(s => ['a','b'].forEach(side => edges(side,s.pathway).forEach(g => all.add(g))));
    let gs = [...all].map(symbol => genes.get(symbol)||{symbol});
    const order = el('gene-order').value, geneFDR = g => Math.min(...['a','b'].map(side => finite(g['de_fdr_'+side]) ? number(g['de_fdr_'+side]) : Infinity));
    gs.sort((a,b)=>(order==='recurrence'?(recurrence.get(b.symbol)||0)-(recurrence.get(a.symbol)||0):order==='fdr'?geneFDR(a)-geneFDR(b):order==='effect'?largest(b,'log2FC')-largest(a,'log2FC'):0)||byID(a.symbol,b.symbol));
    const nGenes = gs.length; gs = gs.slice(0,number(M.max_genes));
    return {sets:ss,genes:gs,n_matching_sets:n,n_available_genes:nGenes,recurrence,selection:mode};
  }
  function matrix() {
    const {sets:ss,genes:gg} = current, left=330, cell=23, annotation=175, gap=55, panel=Math.max(1,gg.length)*cell+annotation;
    const width=left+2*panel+gap+30, top=294, row=26, height=top+Math.max(1,ss.length)*row+145;
    const r=rootSVG(width,height,`${categories.get(active).category_display_name}: A/B aligned leading edges`);
    txt(r,16,25,categories.get(active).category_display_name,{'font-size':17,'font-weight':700});
    txt(r,16,48,`${M.contrast_label} | ${M.collection} | ${M.tier} | ${current.selection}`);
    txt(r,16,70,`${ss.length}/${current.n_matching_sets} sets; ${gg.length}/${current.n_available_genes} genes. Same order in A and B.`);
    for(const side of ['a','b']) {
      const offset=left+(side==='a'?0:panel+gap);
      txt(r,offset,99,`${side.toUpperCase()}: ${sideLabel(side)}`,{'font-size':14,'font-weight':700});
      txt(r,offset,122,'Gene log2FC (top) and DE FDR (below)',{'font-size':11});
      gg.forEach((g,j)=>{
        const x=offset+j*cell+cell/2;
        if(selectedGene===g.symbol) r.appendChild(svg('rect',{x:x-cell/2,y:top-17,width:cell,height:Math.max(1,ss.length)*row+20,fill:'#fff1c8'}));
        // Vertical numeric labels fit the same gene columns without hiding data.
        // Tooltips retain the full stored number; formatting affects display only.
        for (const [field,y,size] of [['log2FC',190,9],['de_fdr',254,8]]) {
          const value=g[field+'_'+side], available=finite(value);
          const label=txt(r,x,y,available?fmt(value):'NA',{'font-size':size,'text-anchor':'start',
            transform:`rotate(-90 ${x} ${y})`,fill:field==='log2FC'?color(value):'#183244',
            class:'gene-numeric-annotation','data-gene':g.symbol,'data-side':side,'data-field':field});
          tooltip(label,`${g.symbol}: ${field==='log2FC'?'log2FC':'FDR DE'} ${available?String(value):'not available'}`);
        }
        const t=txt(r,x,top+ss.length*row+20,g.symbol,{'font-size':10,transform:`rotate(55 ${x} ${top+ss.length*row+20})`});
        t.setAttribute('tabindex','0');t.style.cursor='pointer';t.addEventListener('click',()=>selectGene(g.symbol));t.addEventListener('keydown',ev=>{if(ev.key==='Enter')selectGene(g.symbol);});
      });
      const ann=offset+Math.max(1,gg.length)*cell+12;
      txt(r,ann,top-16,'NES',{'font-size':11});txt(r,ann+48,top-16,'FDR GSEA',{'font-size':11});txt(r,ann+115,top-16,'LE shown/total',{'font-size':10});
      ss.forEach((s,i)=>{
        const y=top+i*row;
        r.appendChild(svg('rect',{x:offset,y:y-14,width:Math.max(1,gg.length)*cell,height:row-1,fill:!bool(s['evaluable_'+side])?'#e2e5e8':i%2?'#f1f5f9':'#fff',opacity:.7}));
        if(side==='a') tooltip(txt(r,12,y,short(s.pathway,45),{'font-size':10}),s.pathway);
        let shown=0; const member=edges(side,s.pathway);
        gg.forEach((g,j)=>{if(member.has(g.symbol)){shown++;tooltip(point(r,side,offset+j*cell+cell/2,y-3,s['NES_'+side],4),`${g.symbol} in ${side.toUpperCase()} leading edge of ${s.pathway}`);}});
        txt(r,ann,y,fmt(s['NES_'+side]),{'font-size':10,fill:color(s['NES_'+side])});txt(r,ann+48,y,fmt(s['gsea_fdr_'+side]),{'font-size':10});
        const recorded=bool(s['evaluable_'+side])&&['recorded','empty'].includes(s['leading_edge_state_'+side]);
        txt(r,ann+115,y,recorded?`${shown}/${member.size}`:'not available',{'font-size':9});
      });
    }
    if(!ss.length)txt(r,20,top,'No gene sets match this selection; complete states remain in the tables.');
    else if(!gg.length)txt(r,left,top+ss.length*row+30,'No recorded leading-edge genes in the selected sets.');
    el('matrix').replaceChildren(r);
  }
  function support() {
    const c=categories.get(active), rows=[];
    const fields=[['Mapped sets','n_mapped',String],['Evaluable sets','n_evaluable',String],['Significant sets (denominator n)','n_genesets_significant',String],['Positive NES','positive_pct',pct],['Negative NES','negative_pct',pct],['Exactly zero NES','zero_pct',pct],['Mean NES (significant sets)','mean_NES',fmt],['Displayed contrast endpoint','display_mean_NES',fmt],['Median NES','median_NES',fmt],['P25 - P75 (dispersion, not CI)',null,null],['Direction concordance with own mean','same_direction_pct',pct]];
    fields.forEach(([label,key,format])=>rows.push([label,...['a','b'].map(side=>{if(['positive_pct','negative_pct','zero_pct'].includes(key)){const count={positive_pct:'n_pos_genesets',negative_pct:'n_neg_genesets',zero_pct:'n_zero_genesets'}[key];return `${c[count+'_'+side]} / ${c['n_genesets_significant_'+side]} (${pct(c[key+'_'+side])})`;}return key?format(c[key+'_'+side]):`${fmt(c['p25_NES_'+side])} - ${fmt(c['p75_NES_'+side])}`;})]));
    table('support',['Support / summary',`A: ${sideLabel('a')}`,`B: ${sideLabel('b')}`],rows);
    el('composition').textContent=`Significant in both: ${c.n_significant_both}; only A: ${c.n_significant_A_only}; only B: ${c.n_significant_B_only}. Mean difference A-B: ${fmt(c.delta_mean_NES)}. ${c.delta_state==='not_in_canonical_contrast'?'This category has no row in the original contrast table; no delta has been added. ':''}Significant means use each side's own support. A hollow endpoint uses original non-significant contextual support, not zero. Changes can reflect different contributing sets.`;
  }
  function jaccard(a,b,available){if(!available)return {intersection_n:null,union_n:null,jaccard:null,state:'not available'};let n=0;a.forEach(g=>{if(b.has(g))n++;});const u=a.size+b.size-n;return {intersection_n:n,union_n:u,jaccard:u?n/u:null,state:u?'defined':'empty union'};}
  function overlapTables() {
    table('conservation',['Gene set','LE genes A','LE genes B','Shared','Union','Jaccard A/B','State'],current.sets.map(s=>{const c=conservation.get(s.pathway)||{};return [s.pathway,c.n_a,c.n_b,c.intersection_n,c.union_n,fmt(c.jaccard),stateText(c.overlap_state)];}));
    const prior=el('overlap-set').value;el('overlap-set').replaceChildren();current.sets.forEach(s=>{const o=document.createElement('option');o.value=s.pathway;o.textContent=s.pathway;el('overlap-set').appendChild(o);});if(current.sets.some(s=>s.pathway===prior))el('overlap-set').value=prior;
    withinOverlap();
  }
  function withinOverlap(){const reference=current.sets.find(s=>s.pathway===el('overlap-set').value), rows=[];if(!reference){empty('overlaps','No gene sets selected.');return;}
    current.sets.filter(s=>s.pathway!==reference.pathway).forEach(s=>['a','b'].forEach(side=>{
      const available=[s,reference].every(v=>bool(v['evaluable_'+side])&&['recorded','empty'].includes(v['leading_edge_state_'+side]));
      const o=jaccard(edges(side,reference.pathway),edges(side,s.pathway),available);rows.push([side.toUpperCase(),s.pathway,o.intersection_n,o.union_n,fmt(o.jaccard),o.state]);}));
    table('overlaps',['Side','Compared gene set','Shared LE genes','Union','Jaccard','State'],rows);}
  function geneExplorerHref(side) {
    const p = new URLSearchParams({gene:selectedGene, analysis_id:M['analysis_'+side], collection:M.collection,
      tier:M.tier, return_contrast_id:M.contrast_id, return_category_id:active});
    // Only report-generated paths establish a portable return route. The
    // directory is the actual contrast/output scope, not merely contrast_id.
    let parts = [];
    try { parts = location.pathname.split('/').slice(-5).map(decodeURIComponent); } catch (_) { /* no local return route */ }
    if (parts.length === 5 && parts[0] === 'report_pages' && parts[1] === 'contrast_evidence' &&
        parts[3] === String(M.collection) && parts[4] === 'index.html' && parts[2] &&
        !['.', '..'].includes(parts[2]) && !/[\/\\\u0000-\u001f\u007f]/.test(parts[2])) {
      p.set('return_contrast_scope', parts[2]);
      p.set('return_analysis_a', M.analysis_a); p.set('return_analysis_b', M.analysis_b);
      p.set('return_side', side);
    }
    // Gene search reads location.search; contrast-local category navigation uses a hash.
    return M.gene_explorer_href + '?' + p.toString();
  }
  function geneDetails(){const rows=[], selected=selectedGene?[genes.get(selectedGene)].filter(Boolean):current.genes;
    selected.forEach(g=>['a','b'].forEach(side=>rows.push([g.symbol,side.toUpperCase(),sideLabel(side),fmt(g['log2FC_'+side]),fmt(g['de_fdr_'+side]),stateText(g['de_state_'+side]||'not_recorded'),current.recurrence.get(g.symbol)||0])));
    table('gene-table',['Gene','Side','Analysis','log2FC','FDR DE','DE state','LE recurrence (A+B connections)'],rows);
    el('gene-notice').textContent=selectedGene&&!current.genes.some(g=>g.symbol===selectedGene)?(genes.has(selectedGene)?`${selectedGene} exists in the evidence but is outside this displayed matrix; its original side-specific DE is shown below.`:`No exact symbol ${selectedGene} occurs in this evidence.`):selectedGene?`Selected gene: ${selectedGene}`:'';
    el('gene-links').replaceChildren();
    if(selectedGene&&M.gene_explorer_href)['a','b'].forEach(side=>{const a=document.createElement('a');a.href=geneExplorerHref(side);a.textContent=`Explore this gene in ${side.toUpperCase()}: ${sideLabel(side)}`;el('gene-links').append(a,document.createElement('br'));});
  }
  function setChart(){const ss=current.sets, width=1150,left=450,span=360,row=47,top=75,r=rootSVG(width,top+Math.max(1,ss.length)*row+40,'Aligned original gene-set NES for A and B');
    const all=setsByCategory.get(active)||[], lim=Math.max(1,...all.flatMap(s=>[s.NES_a,s.NES_b]).filter(finite).map(x=>Math.abs(number(x))))*1.1, x=v=>left+span/2+number(v)/lim*span/2;
    txt(r,15,24,`${categories.get(active).category_display_name}: gene-set NES`,{'font-size':16,'font-weight':700});
    txt(r,845,52,'A: NES / FDR');txt(r,990,52,'B: NES / FDR');
    r.appendChild(svg('line',{x1:x(0),x2:x(0),y1:55,y2:top+ss.length*row,stroke:'#b4c2cc','stroke-dasharray':'3 3'}));
    ss.forEach((s,i)=>{const y=top+i*row;tooltip(txt(r,15,y,short(s.pathway,60),{'font-size':11}),s.pathway);
      if(finite(s.NES_a)&&finite(s.NES_b))r.appendChild(svg('line',{x1:x(s.NES_a),x2:x(s.NES_b),y1:y-5,y2:y+5,stroke:'#8394a1'}));
      ['a','b'].forEach(side=>{if(finite(s['NES_'+side]))point(r,side,x(s['NES_'+side]),y+(side==='a'?-5:5),s['NES_'+side]);});
      txt(r,845,y,`${fmt(s.NES_a)} / ${fmt(s.gsea_fdr_a)}`,{'font-size':11});txt(r,990,y,`${fmt(s.NES_b)} / ${fmt(s.gsea_fdr_b)}`,{'font-size':11});});el('set-chart').replaceChildren(r);
    table('set-table',['Gene set','A NES','A GSEA FDR','A result','A leading edge','B NES','B GSEA FDR','B result','B leading edge','Significance membership'],ss.map(s=>[s.pathway,fmt(s.NES_a),fmt(s.gsea_fdr_a),stateText(s.state_a||(!bool(s.present_a)?'not mapped':'non_evaluable')),stateText(s.leading_edge_state_a||'not_recorded'),fmt(s.NES_b),fmt(s.gsea_fdr_b),stateText(s.state_b||(!bool(s.present_b)?'not mapped':'non_evaluable')),stateText(s.leading_edge_state_b||'not_recorded'),stateText(s.significance_membership)]));
  }
  function settings(){const c=categories.get(active);return {schema_version:'1.0',contrast_id:M.contrast_id,analysis_a:M.analysis_a,analysis_b:M.analysis_b,collection:M.collection,tier:M.tier,category_id:active,selected_gene:selectedGene,selection:el('selection').value,set_order:el('set-order').value,gene_order:el('gene-order').value,max_sets:M.max_sets,max_genes:M.max_genes,gsea_padj_cutoff:M.gsea_padj_cutoff,de_padj_cutoff:M.de_padj_cutoff,delta_mean_NES:c.delta_mean_NES,delta_state:c.delta_state,displayed_set_ids:current.sets.map(s=>s.pathway),displayed_gene_symbols:current.genes.map(g=>g.symbol),n_matching_sets:current.n_matching_sets,n_available_genes:current.n_available_genes};}
  function source(){const columns=['row_type','contrast_id','contrast_label','analysis_a','analysis_b','collection','tier','category_id','category_display_name','selection','n_matching_sets','n_available_genes','max_sets','max_genes','side','pathway','symbol','set_order','gene_order','NES','gsea_fdr','evaluable','significant','leading_edge_state','log2FC','de_fdr','de_state'];
    const base={...M,...settings(),category_display_name:categories.get(active).category_display_name}, rows=[{...base,row_type:'metadata'}];
    ['a','b'].forEach(side=>{current.sets.forEach((s,i)=>{const row={...base,row_type:'set',side,pathway:s.pathway,set_order:i+1};['NES','gsea_fdr','evaluable','significant','leading_edge_state'].forEach(k=>row[k]=s[k+'_'+side]);rows.push(row);});current.genes.forEach((g,j)=>{const row={...base,row_type:'gene',side,symbol:g.symbol,gene_order:j+1};['log2FC','de_fdr','de_state'].forEach(k=>row[k]=g[k+'_'+side]);rows.push(row);});current.sets.forEach((s,i)=>current.genes.forEach((g,j)=>{if(edges(side,s.pathway).has(g.symbol))rows.push({...base,row_type:'membership',side,pathway:s.pathway,symbol:g.symbol,set_order:i+1,gene_order:j+1});}));});
    return [columns.join('\t'),...rows.map(row=>columns.map(k=>escTSV(row[k]===undefined?'':row[k])).join('\t'))].join('\n')+'\n';}
  function updateLocation(){if(!active)return;try{history.replaceState(null,'',hashFor(active));}catch(_){/* Offline viewers may disable history. */}syncShell();}
  function selectGene(g){selectedGene=String(g||'');el('gene-search').value=selectedGene;redraw();updateLocation();}
  function renderStaticFigure(){
    const ix=arr(D.figure_index).find(x=>String(x.category_id)===String(active)), downloads=el('static-downloads'), figure=el('static-figure');
    downloads.replaceChildren();figure.replaceChildren();if(!ix)return;
    arr(D.formats).forEach(format=>downloads.appendChild(localLink(`Default paired matrix ${String(format).toUpperCase()}`,`figures/${ix.stem}.${format}`)));
    if(D.source_data !== false) downloads.appendChild(localLink('Default figure source TSV',`figures/${ix.stem}_source.tsv`));
    const imageFormat=arr(D.formats).find(format=>['png','svg'].includes(String(format).toLowerCase()));
    if(imageFormat){const img=document.createElement('img');img.src=`figures/${ix.stem}.${imageFormat}`;img.alt=`Default paired matrix for ${categories.get(active).category_display_name}`;img.loading='lazy';figure.appendChild(img);}
  }
  function redraw(){if(!active)return;current=makeView();el('view-counts').textContent=`Showing ${current.sets.length} of ${current.n_matching_sets} matching sets and ${current.genes.length} of ${current.n_available_genes} available leading-edge genes. All data remain in the complete tables.`;matrix();overlapTables();geneDetails();setChart();}
  function openCategory(id, writeLocation=true){
    if(!categories.has(id))return;active=id;selectedGene=selectedGene||'';el('category').value=id;
    const c=categories.get(id), info=el('category-info');info.replaceChildren();const title=document.createElement('h2');title.id='category-title';title.textContent=c.category_display_name;info.appendChild(title);
    if(c.category_description){const p=document.createElement('p');p.textContent=c.category_description;info.appendChild(p);}
    const group=c.macrogroup_name||c.macrogroup_id||'Macrogroup not recorded', summary=document.createElement('p');summary.className='note';summary.textContent=`${group}. A minus B: ${fmt(c.delta_mean_NES)} (${stateText(c.delta_state)}).`;info.appendChild(summary);
    const cards=[['Mapped sets',`${c.n_mapped_a} / ${c.n_mapped_b}`],['Evaluable sets',`${c.n_evaluable_a} / ${c.n_evaluable_b}`],['Significant sets',`${c.n_genesets_significant_a} / ${c.n_genesets_significant_b}`],['Significant in both / A only / B only',`${c.n_significant_both} / ${c.n_significant_A_only} / ${c.n_significant_B_only}`]];
    el('counts').replaceChildren(...cards.map(([label,value])=>{const card=document.createElement('div');card.className='count';const strong=document.createElement('strong');strong.textContent=value;const span=document.createElement('span');span.textContent=label;card.append(strong,span);return card;}));
    support();renderCategoryProducts();renderStaticFigure();redraw();if(writeLocation)updateLocation();else syncShell();
  }
  function clearSelection(){selectedGene='';el('gene-search').value='';redraw();updateLocation();}
  function clearDisplayedCategory(){
    active=null;current=null;selectedGene='';el('category').value='';el('gene-search').value='';
    ['category-info','counts','support','category-product-cards','matrix','conservation','overlaps','gene-table','gene-links','set-chart','set-table','static-downloads','static-figure'].forEach(id=>el(id).replaceChildren());
    el('category-products').hidden=true;el('overlap-set').replaceChildren();el('view-counts').textContent='';el('composition').textContent='';el('gene-notice').textContent='';
  }
  function fromHash(){
    const hash=new URLSearchParams(location.hash.slice(1)), query=new URLSearchParams(location.search), get=key=>hash.has(key)?hash.get(key):query.get(key);
    for(const k of ['contrast_id','analysis_a','analysis_b','collection','tier'])if(get(k)!=null&&get(k)!==String(M[k])){clearDisplayedCategory();el('scope-error').hidden=false;el('scope-error').textContent='This link belongs to a different contrast, collection or tier. No category is displayed under the wrong scope.';return;}
    const id=get('category_id')||query.get('category')||D.categories[0]?.category_id||'';
    if(id&&!categories.has(String(id))){clearDisplayedCategory();el('scope-error').hidden=false;el('scope-error').textContent=`No exact category ${id} occurs in this contrast. No fallback category was substituted.`;return;}
    el('scope-error').hidden=true;selectedGene=get('gene')||'';el('gene-search').value=selectedGene;openCategory(String(id),false);
  }
  el('side-context').textContent=`A: ${sideLabel('a')} | B: ${sideLabel('b')}. Positive direction in A: ${M.positive_contrast_a||'not recorded'}; in B: ${M.positive_contrast_b||'not recorded'}.`;
  D.categories.forEach(c=>{const o=document.createElement('option');o.value=c.category_id;o.textContent=`${c.category_display_name} [${c.category_id}]`;el('category').appendChild(o);});
  ['categories','sets','genes','leading_edges','conservation','metadata','figure_index'].forEach(k=>{const li=document.createElement('li'),a=document.createElement('a');a.href=`tables/${k}.tsv`;a.download='';a.textContent=k.replace(/_/g,' ');li.appendChild(a);el('downloads').appendChild(li);});
  ['selection','set-order','gene-order'].forEach(id=>el(id).addEventListener('change',redraw));el('overlap-set').addEventListener('change',withinOverlap);
  el('category').addEventListener('change',()=>{selectedGene='';el('gene-search').value='';openCategory(el('category').value);});
  el('gene-search').addEventListener('change',()=>selectGene(el('gene-search').value.trim()));el('clear-selection').addEventListener('click',clearSelection);
  el('reset').addEventListener('click',()=>{el('selection').value='union';el('set-order').value='fdr';el('gene-order').value='recurrence';selectedGene='';el('gene-search').value='';redraw();updateLocation();});
  for (const [id, enabled] of [['export-matrix', arr(D.formats).includes('svg')], ['export-set-chart', arr(D.formats).includes('svg')], ['export-source', D.source_data !== false], ['export-settings', D.recipes !== false]]) { el(id).hidden = !enabled; el(id).disabled = !enabled; }
  el('export-matrix').addEventListener('click',()=>downloadSVG('matrix','contrast_evidence_current.svg',true));
  el('export-set-chart').addEventListener('click',()=>downloadSVG('set-chart','contrast_gene_sets_current.svg',true));
  el('export-source').addEventListener('click',()=>download('contrast_evidence_current_source.tsv',source(),'text/tab-separated-values;charset=utf-8'));
  el('export-settings').addEventListener('click',()=>download('contrast_evidence_current_settings.json',JSON.stringify(settings(),null,2),'application/json'));
  bindProductsJump();
  document.addEventListener('lisa:shell-ready',syncShell);addEventListener('hashchange',fromHash);fromHash();
})();
