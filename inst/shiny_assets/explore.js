/* Injects the per-category exploration control into the native report shown in
   the iframe. The report remains the presentation; this is not a dashboard and
   it never fetches or generates anything.

   Three properties this file is responsible for:

   G1-4 placement. A block belongs to the category the reader has actually
   opened, not to the foot of a collection page. The route now resolves to the
   evidence sheet (`report_pages/evidence/...`, `report_pages/contrast_evidence/...`)
   which renders one category at a time and mirrors the choice into `?category=`
   and `#category-<id>`; rows are filtered by that active category, so category A
   can never show B's figure. Insertion uses the same rule as the static
   exporter in R/explore_presentation.R -- before the first `section.panel`
   inside `main`, else at the end of `main` -- so the running app and the
   exported bundle put the block in the same place.

   G1-4 focus. The server sends the payload once a second. Rebuilding the
   injected DOM on every message discarded keyboard focus and reset a pressed
   button once a second. Everything rendered here is a pure function of
   (route, active category, rows, feedback, pending), so that tuple is
   fingerprinted and the DOM is touched only when it changes.

   G1-3 feedback. A submission can be refused -- one worker at a time, a lock
   held elsewhere, a request the engine will not accept. The refusal is shown in
   the category the reader is looking at and the button becomes usable again. A
   refused submission never reports success. */
(function () {
  var rows = [], feedback = {}, pending = {}, selections = {}, applied = null,
      observer = null;

  function frameDoc() {
    var frame = document.getElementById('lisa-explore-report');
    if (!frame) { return null; }
    try { return frame.contentDocument || null; } catch (error) { return null; }
  }
  // An iframe is a separate document and does not inherit the shell's CSS.
  // Reuse the already-served, same-origin asset URL rather than copying styles
  // into the native report or introducing a second stylesheet. The guard lives
  // in the child document, so each report navigation gets one link of its own.
  function installStyles(doc) {
    var parentStyle = document.getElementById('lisa-explore-stylesheet');
    if (!doc.head || !parentStyle || !parentStyle.href ||
        doc.getElementById('lisa-explore-iframe-stylesheet')) { return; }
    var style = doc.createElement('link');
    style.id = 'lisa-explore-iframe-stylesheet';
    style.rel = 'stylesheet';
    style.href = parentStyle.href;
    doc.head.appendChild(style);
  }
  function routeFor(doc) {
    // `/lisa-explore-report-<nonce>/report_pages/evidence/<o>/<c>/index.html`
    // -> `report_pages/evidence/<o>/<c>/index.html`
    var path = (doc.location && doc.location.pathname) || '';
    try { path = decodeURIComponent(path); } catch (error) { /* keep raw */ }
    return path.replace(/^.*?(report_pages\/)/, '$1');
  }
  function scopedBy(doc) { return !!doc.getElementById('category-info'); }
  function activeCategory(doc) {
    var select = doc.getElementById('category');
    if (select && select.value) { return select.value; }
    var view = doc.defaultView;
    try {
      var query = new URLSearchParams(view.location.search).get('category');
      if (query) { return query; }
    } catch (error) { /* no search string available */ }
    var hash = (view && view.location.hash) || '';
    if (hash.indexOf('#category-') === 0) {
      try { return decodeURIComponent(hash.slice(10)); } catch (error) { return hash.slice(10); }
    }
    return '';
  }
  function visibleRows(doc) {
    var route = routeFor(doc), scoped = scopedBy(doc), category = activeCategory(doc);
    return rows.filter(function (row) {
      if (row.route !== route) { return false; }
      // H3. A collection-wide row has no category at all, so filtering it by the
      // open category would hide it everywhere. It belongs to the owner's own
      // collection block on the native page, which is a different route anyway.
      if (row.scope === 'collection') { return true; }
      // On a sheet that routes by category, only the open category's row is its
      // own. On the navigator fallback every category is on the page already.
      return !scoped || row.category_id === category;
    });
  }
  function isMap(row) {
    return row.product === 'kegg_pathway_map' || row.product === 'contrast_kegg_map';
  }
  function groupKey(row) {
    return [row.route, row.scope, row.context, row.nav_collection,
            row.category_id, row.product].join('\u001f');
  }
  function groupedRows(shown) {
    var groups = [], byKey = {};
    shown.forEach(function (row) {
      if (!isMap(row)) {
        groups.push({ key: row.request_id, rows: [row] });
        return;
      }
      var key = groupKey(row);
      if (!byKey[key]) {
        byKey[key] = { key: key, rows: [] };
        groups.push(byKey[key]);
      }
      byKey[key].rows.push(row);
    });
    groups.forEach(function (group) {
      if (!isMap(group.rows[0])) { return; }
      group.rows.sort(function (left, right) {
        var lr = Number(left.entity_rank || Number.MAX_SAFE_INTEGER);
        var rr = Number(right.entity_rank || Number.MAX_SAFE_INTEGER);
        if (lr !== rr) { return lr - rr; }
        return String(left.entity).localeCompare(String(right.entity));
      });
    });
    return groups;
  }
  function selectedRow(group) {
    var selected = selections[group.key];
    var row = group.rows.find(function (candidate) {
      return candidate.request_id === selected;
    }) || group.rows[0];
    selections[group.key] = row.request_id;
    return row;
  }
  function fingerprint(doc, groups) {
    return JSON.stringify([routeFor(doc), activeCategory(doc), groups.map(function (group) {
      var selected = selectedRow(group);
      return [group.key, selected.request_id, group.rows.map(function (row) {
        var note = feedback[row.request_id];
        return [row.request_id, row.state, row.entity_title || '', row.entity_rank || '',
                row.block ? 1 : 0, pending[row.request_id] ? 1 : 0,
                note ? note.message : ''];
      })];
    })]);
  }
  // H2 shows four products, each with its own control, so a heading that says
  // only "Volcano" or "Extended figure" would mislabel three of them. The exact
  // selector is part of the name: two heatmaps of one category differ only by
  // their scale, and two pathway maps only by their pathway.
  var PRODUCT_TITLES = {
    volcano: 'Volcano',
    gene_cards: 'Prioritized gene card',
    heatmap: 'Gene heatmap',
    kegg_pathway_map: 'KEGG pathway map',
    // H3. The five native products, with the section names H3_SCOPE fixes.
    de_recurrent_genes: 'Recurrent genes',
    contrast_gene_card: 'Contrast gene card',
    contrast_paired_heatmap: 'Paired gene heatmaps',
    contrast_gene_category_network: 'Gene-category network',
    contrast_kegg_map: 'KEGG maps / painted pathways'
  };
  function selectorText(row) {
    var parts = [];
    // The network's variant is a fixed exploration identity, not a scale, so it
    // is not suffixed with the word "scale" the heatmap uses.
    if (row.entity) {
      parts.push(row.entity_title ? row.entity_title + ' (' + row.entity + ')' : row.entity);
    }
    if (row.variant) {
      parts.push(row.product === 'heatmap' ? row.variant + ' scale' : row.variant);
    }
    return parts.join(', ');
  }
  function subjectText(row) {
    if (row.scope === 'collection') {
      return [row.context || '', row.nav_collection || row.collection || '']
        .filter(function (part) { return !!part; }).join(' ');
    }
    return row.category_id;
  }
  function heading(row) {
    var title = PRODUCT_TITLES[row.product] || 'Extended figure';
    var selector = selectorText(row);
    if (selector) { title += ' (' + selector + ')'; }
    var subject = subjectText(row);
    return subject ? title + ' · ' + subject : title;
  }
  function mapOptionLabel(row) {
    var title = row.entity_title || 'Untitled pathway';
    var rank = row.entity_rank ? 'rank ' + row.entity_rank + ' · ' : '';
    return rank + title + ' (' + row.entity + ') — ' + row.state.replace('_', ' ');
  }
  function control(doc, row, alternatives) {
    var section = doc.createElement('section');
    section.className = 'lisa-explore-control';
    section.dataset.lisaExploreUi = 'true';
    section.dataset.lisaExploreCategory = row.category_id;
    section.dataset.lisaExploreProduct = row.product;
    if (row.entity) { section.dataset.lisaExploreEntity = row.entity; }
    if (row.variant) { section.dataset.lisaExploreVariant = row.variant; }
    section.dataset.lisaExploreScope = row.scope || 'category';
    // A collection-wide control has no category anchor to carry; writing
    // "category-" with nothing after it would name an anchor that exists nowhere.
    section.dataset.lisaAttachAnchor = row.scope === 'collection' ?
      (row.anchor || '') : 'category-' + row.category_id;
    section.innerHTML = '<h2></h2><p>State: <span data-lisa-explore-state></span></p>';
    section.querySelector('h2').textContent = heading(row);
    section.querySelector('[data-lisa-explore-state]').textContent = row.state.replace('_', ' ');
    if (isMap(row)) {
      var label = doc.createElement('label');
      label.className = 'lisa-explore-map-selector';
      label.textContent = 'Native KEGG pathway (title and ID)';
      var select = doc.createElement('select');
      select.setAttribute('aria-label', 'Native KEGG pathway title and ID');
      alternatives.forEach(function (candidate) {
        var option = doc.createElement('option');
        option.value = candidate.request_id;
        option.textContent = mapOptionLabel(candidate);
        option.selected = candidate.request_id === row.request_id;
        select.appendChild(option);
      });
      select.addEventListener('change', function () {
        selections[groupKey(row)] = select.value;
        // Selection changes presentation only. It deliberately sends no Shiny
        // input; a map is generated only by the explicit button below.
        applied = null;
        apply();
      });
      label.appendChild(select);
      section.insertBefore(label, section.querySelector('p'));
    }
    // A product that cannot be drawn says why, in the reader's own words. A
    // disabled control with no explanation is indistinguishable from a bug.
    if (row.reason) {
      var why = doc.createElement('p');
      why.className = 'lisa-explore-reason';
      why.textContent = row.reason;
      section.appendChild(why);
    }
    if (row.state === 'available' && row.block) {
      var available = doc.createElement('div');
      available.innerHTML = row.block;
      available.querySelectorAll('[data-lisa-explore-attachment]').forEach(function (node) {
        node.hidden = false;
      });
      section.appendChild(available);
    } else if (row.state === 'ungenerated' || row.state === 'failed') {
      var button = doc.createElement('button');
      button.type = 'button';
      var waiting = !!pending[row.request_id];
      button.disabled = waiting;
      // "Requesting" is a request, not a result. The success wording only ever
      // comes from a state the server actually reported.
      button.textContent = waiting ? 'Requesting…' : (row.label || 'Generate figure');
      button.addEventListener('click', function () {
        var api = window.parent && window.parent.Shiny;
        if (!api) { return; }
        pending[row.request_id] = true;
        delete feedback[row.request_id];
        button.disabled = true;
        button.textContent = 'Requesting…';
        // entity/variant are sent with the click. Without them the server would
        // rebuild a request for the product's default and draw the wrong figure.
        api.setInputValue('explore_generate', {
          unit_type: row.unit_type, analysis_id: row.analysis_id,
          contrast_id: row.contrast_id, collection: row.collection,
          category_id: row.category_id, product: row.product,
          entity: row.entity || '', variant: row.variant || '',
          nonce: Date.now()
        }, { priority: 'event' });
      });
      section.appendChild(button);
    }
    var note = feedback[row.request_id];
    if (note && note.message) {
      var error = doc.createElement('p');
      error.className = 'lisa-explore-error';
      error.setAttribute('role', 'alert');
      error.textContent = note.message;
      section.appendChild(error);
    }
    return section;
  }
  // The one placement rule, matching R/explore_presentation.R exactly: insert
  // before the first `section.panel` of the target context, else append to it.
  function insert(doc, container, node) {
    var panel = container.querySelector('section.panel');
    if (panel) { container.insertBefore(node, panel); } else { container.appendChild(node); }
  }
  // H3. The target context of a collection-wide row is the owner's own
  // collection block, addressed by `data-lisa-nav-context` and
  // `data-lisa-nav-collection` -- never the first panel on the page. A row whose
  // context is genuinely absent is skipped rather than dropped somewhere else.
  function containerFor(doc, main, row) {
    if (row.scope !== 'collection') { return main; }
    var context = doc.querySelector('[data-lisa-nav-context="' +
      (row.context || '').replace(/"/g, '\\"') + '"]');
    if (!context) { return null; }
    var block = context.querySelector('[data-lisa-nav-collection="' +
      (row.nav_collection || row.collection || '').replace(/"/g, '\\"') + '"]');
    if (!block) { return null; }
    return block.querySelector('.collection-inner') || block;
  }
  function apply() {
    var doc = frameDoc();
    if (!doc) { return; }
    installStyles(doc);
    var main = doc.querySelector('main');
    if (!main) { return; }
    var shown = visibleRows(doc);
    var groups = groupedRows(shown);
    var signature = fingerprint(doc, groups);
    // Nothing the reader can see has changed, so touch nothing: a focused
    // button keeps focus through every idle second.
    if (signature === applied) { return; }
    applied = signature;
    doc.querySelectorAll('[data-lisa-explore-ui="true"]').forEach(function (node) { node.remove(); });
    groups.forEach(function (group) {
      var row = selectedRow(group);
      var container = containerFor(doc, main, row);
      if (!container) { return; }
      if (!isMap(row) && row.state === 'available' && row.block) {
        var wrap = doc.createElement('div');
        wrap.dataset.lisaExploreUi = 'true';
        wrap.dataset.lisaExploreCategory = row.category_id;
        wrap.innerHTML = row.block;
        // The block is emitted `hidden` for the static bundle, where a script
        // reveals the active category. Here the row list is already filtered to
        // the open category, so it is shown directly.
        wrap.querySelectorAll('[data-lisa-explore-attachment]').forEach(function (node) {
          node.hidden = false;
        });
        insert(doc, container, wrap);
      } else {
        insert(doc, container, control(doc, row, group.rows));
      }
    });
  }
  function watchFrame() {
    var doc = frameDoc();
    if (!doc) { return; }
    applied = null;
    if (observer) { observer.disconnect(); observer = null; }
    var info = doc.getElementById('category-info');
    if (info && window.MutationObserver) {
      // The sheet re-renders #category-info whenever the reader changes
      // category, and records the change with history.replaceState -- which
      // fires neither popstate nor hashchange. This is the only reliable signal
      // that the open category changed, and it costs nothing while it does not.
      observer = new MutationObserver(function () { apply(); });
      observer.observe(info, { childList: true });
    }
    var select = doc.getElementById('category');
    if (select) { select.addEventListener('change', function () { apply(); }); }
    apply();
  }
  window.addEventListener('load', function () {
    var frame = document.getElementById('lisa-explore-report');
    if (frame) { frame.addEventListener('load', watchFrame); }
    watchFrame();
  });
  if (window.Shiny) {
    window.Shiny.addCustomMessageHandler('lisa-explore-controls', function (payload) {
      var message = payload || {};
      rows = message.rows || [];
      // An empty R list arrives as `[]`, not `{}`; either way there are no notes.
      feedback = (message.feedback && !Array.isArray(message.feedback)) ? message.feedback : {};
      // A request the server has answered -- with a state or with a refusal --
      // is no longer pending, so its button stops saying "Requesting".
      rows.forEach(function (row) {
        if (feedback[row.request_id] || row.state !== 'ungenerated') {
          delete pending[row.request_id];
        }
      });
      apply();
    });
  }
})();
