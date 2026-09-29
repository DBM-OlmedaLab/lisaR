/* Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
   Científicas (CSIC). Author: David Olmeda Casadomé.
   This file is part of lisaR, free software under the GNU General Public
   License version 3 (GPL-3). See the DESCRIPTION file and
   <https://www.gnu.org/licenses/gpl-3.0.html>. */

/* Offline enhancement. Navigation follows the emitted product inventory. */
(function () {
  "use strict";
  var fields = ["study", "analysis", "collection", "selection", "direction", "cutoff"];
  var reportNavigation = null;
  function shellNode() { return document.querySelector("[data-lisa-shell]"); }
  function updateHeight() {
    var shell = shellNode();
    if (!shell) return;
    var height = window.getComputedStyle(shell).position === "sticky" ? shell.getBoundingClientRect().height : 0;
    document.documentElement.style.setProperty("--lisa-shell-height", Math.ceil(height) + "px");
  }
  function updateContext(update) {
    var shell = shellNode();
    if (!shell || !update || typeof update !== "object") return false;
    fields.forEach(function (key) {
      if (!Object.prototype.hasOwnProperty.call(update, key)) return;
      var node = shell.querySelector('[data-lisa-context="' + key + '"]');
      var row = shell.querySelector('[data-lisa-context-field="' + key + '"]');
      var text = update[key] == null ? "" : String(update[key]);
      if (node) node.textContent = text;
      if (row) row.hidden = !text;
    });
    if (Object.prototype.hasOwnProperty.call(update, "exact_ids")) {
      var ids = update.exact_ids || {};
      var list = shell.querySelector("[data-lisa-exact-ids]");
      var details = shell.querySelector(".lisa-shell-identifiers");
      if (list) {
        while (list.firstChild) list.removeChild(list.firstChild);
        Object.keys(ids).forEach(function (key) {
          var row = document.createElement("div");
          var term = document.createElement("dt");
          var description = document.createElement("dd");
          var code = document.createElement("code");
          term.textContent = key;
          code.textContent = ids[key] == null ? "" : String(ids[key]);
          description.appendChild(code);
          row.appendChild(term);
          row.appendChild(description);
          list.appendChild(row);
        });
        if (details) details.hidden = !list.children.length;
      }
    }
    var context = shell.querySelector(".lisa-shell-context");
    if (context) {
      context.hidden = !Array.prototype.some.call(shell.querySelectorAll("[data-lisa-context-field]"),
        function (row) { return !row.hidden; }) &&
        !shell.querySelector(".lisa-shell-identifiers:not([hidden])");
    }
    updateHeight();
    if (reportNavigation) reportNavigation.syncContext(update);
    return true;
  }
  function setActive(key) {
    var shell = shellNode();
    if (!shell) return false;
    var links = shell.querySelectorAll("[data-lisa-route]");
    var match = Array.prototype.filter.call(links, function (link) {
      return link.getAttribute("data-lisa-route") === key;
    });
    if (key != null && match.length !== 1) return false;
    Array.prototype.forEach.call(links, function (link) { link.removeAttribute("aria-current"); });
    if (match.length) match[0].setAttribute("aria-current", "page");
    return true;
  }
  function initReportNavigation(shell) {
    var script = document.getElementById("lisa-report-navigation");
    var bar = shell.querySelector("[data-lisa-explorer]");
    if (!script || !bar) return null;
    var inventory;
    try { inventory = JSON.parse(script.textContent); }
    catch (error) { console.error("Cannot read report navigation inventory", error); return null; }
    if (!inventory || inventory.version !== 1 || !Array.isArray(inventory.contexts)) return null;
    var external = inventory.mode === "external";
    var contexts = inventory.contexts;
    var controls = {};
    ["context", "collection", "section"].forEach(function (key) {
      controls[key] = bar.querySelector('[data-lisa-nav-select="' + key + '"]');
    });
    var state = { context: "", collection: "", section: "" };
    var changing = false;
    var scrolling = false;
    var scrollTimer = null;
    // A hash-restored or programmatic navigation may target a section whose
    // final position has not settled yet (images/figures above it still
    // loading). While a restoration is pending we keep re-applying
    // scrollIntoView (on the next frame, on window "load", and on any
    // document.body size change) until either the first user-initiated
    // scroll/wheel/key/touch input or a short deadline, and we hold the
    // `scrolling` guard the whole time so the debounced scroll listener
    // below cannot observe a not-yet-settled position and overwrite state.
    var pendingTarget = null;
    var pendingMoveFocus = false;
    var pendingDeadline = 0;
    var pendingGuardTimer = null;
    var pendingObserver = null;
    var pendingPoller = null;
    var interactionBound = false;
    function releaseScrollGuard() {
      if (pendingGuardTimer) window.clearTimeout(pendingGuardTimer);
      // Longer than the 100 ms scroll-listener debounce, so its handler
      // cannot fire in the gap between the guard clearing and the last
      // programmatic scroll settling.
      pendingGuardTimer = window.setTimeout(function () { pendingGuardTimer = null; scrolling = false; }, 150);
    }
    function stopPending() {
      var hadPending = !!pendingTarget;
      pendingTarget = null; pendingMoveFocus = false; pendingDeadline = 0;
      if (pendingObserver) { pendingObserver.disconnect(); pendingObserver = null; }
      if (pendingPoller) { window.clearInterval(pendingPoller); pendingPoller = null; }
      if (hadPending) releaseScrollGuard();
    }
    function onUserInteraction() { if (pendingTarget) stopPending(); }
    function bindInteractionListeners() {
      if (interactionBound) return;
      interactionBound = true;
      ["wheel", "keydown", "touchstart"].forEach(function (type) {
        window.addEventListener(type, onUserInteraction, { capture: true, passive: true });
      });
    }
    function rescroll() {
      if (!pendingTarget || !pendingTarget.isConnected) { stopPending(); return; }
      scrolling = true;
      if (pendingGuardTimer) { window.clearTimeout(pendingGuardTimer); pendingGuardTimer = null; }
      updateHeight();
      pendingTarget.scrollIntoView({ block: "start", behavior: "instant" });
      if (pendingMoveFocus) {
        if (!pendingTarget.hasAttribute("tabindex")) pendingTarget.setAttribute("tabindex", "-1");
        pendingTarget.focus({ preventScroll: true });
      }
      if (Date.now() >= pendingDeadline) stopPending();
    }
    function startPending(target, moveFocus) {
      stopPending();
      pendingTarget = target; pendingMoveFocus = !!moveFocus; pendingDeadline = Date.now() + 2000;
      bindInteractionListeners();
      window.requestAnimationFrame(rescroll);
      if (document.readyState !== "complete") {
        window.addEventListener("load", function onLoad() {
          window.removeEventListener("load", onLoad);
          if (pendingTarget) window.requestAnimationFrame(rescroll);
        });
      }
      if (window.ResizeObserver) {
        pendingObserver = new ResizeObserver(function () { if (pendingTarget) rescroll(); });
        pendingObserver.observe(document.body);
      } else {
        pendingPoller = window.setInterval(function () {
          if (pendingTarget) rescroll(); else { window.clearInterval(pendingPoller); pendingPoller = null; }
        }, 200);
      }
      window.setTimeout(function () { if (pendingTarget) stopPending(); }, 2000);
    }
    var has = function (x) { return typeof x === "string" && x.length > 0; };
    var unique = function (rows) {
      if (!Array.isArray(rows) || !rows.every(function (x) { return x && typeof x === "object"; })) return false;
      var ids = rows.map(function (x) { return x.id; });
      return ids.every(has) && new Set(ids).size === ids.length;
    };
    var safeHref = function (href) {
      return has(href) && !/^[\s/]|[\s]$|[\\\x00-\x1f]|^[A-Za-z][A-Za-z0-9+.-]*:/.test(href);
    };
    var sectionHref = function (section) { return section.href || "#" + encodeURIComponent(section.id); };
    var valid = unique(contexts) && contexts.every(function (context) {
      return has(context.label) && Array.isArray(context.collections) && unique(context.collections) &&
        context.collections.every(function (collection) {
          return has(collection.label) && Array.isArray(collection.sections) && unique(collection.sections) &&
            collection.sections.every(function (section) {
              return has(section.label) && safeHref(sectionHref(section)) &&
                (external || !!document.getElementById(section.id));
            });
        });
    });
    if (!valid) {
      console.error("Report navigation inventory does not match its available sections");
      return null;
    }
    // An overview or standalone report may have no generated result sections.
    // Leave the contextual controls hidden without fabricating a destination.
    if (!contexts.length) return null;
    var find = function (rows, id) { return rows.find(function (x) { return x.id === id; }); };
    function selected() {
      var context = find(contexts, state.context);
      var collection = context && find(context.collections, state.collection);
      var section = collection && find(collection.sections, state.section);
      return { context: context, collection: collection, section: section };
    }
    function options(control, rows, value, placeholder) {
      while (control.firstChild) control.removeChild(control.firstChild);
      if (!value) {
        var empty = document.createElement("option"); empty.value = "";
        empty.textContent = placeholder; control.appendChild(empty);
      }
      rows.forEach(function (row) {
        var option = document.createElement("option"); option.value = row.id;
        option.textContent = row.label; control.appendChild(option);
      });
      control.value = value || "";
      control.disabled = !rows.length;
      var selectedOption = rows.find(function (row) { return row.id === value; });
      control.title = selectedOption ? selectedOption.label : placeholder;
    }
    function refresh() {
      var choice = selected();
      if (external) shell.setAttribute("data-navigation-aggregate", String(!state.context || !state.collection));
      options(controls.context, contexts, state.context, "Choose an analysis or contrast");
      options(controls.collection, choice.context ? choice.context.collections : [], state.collection, "Choose a collection");
      options(controls.section, choice.collection ? choice.collection.sections : [], state.section, "Choose a section");
      bar.querySelector("[data-lisa-nav-location]").textContent = [choice.context, choice.collection, choice.section]
        .filter(Boolean).map(function (x) { return x.label; }).join(" / ");
      var kinds = new Set(contexts.map(function (context) { return context.kind; }));
      bar.querySelector("[data-lisa-nav-context-label]").textContent = kinds.size === 1 && kinds.has("contrast") ?
        "Contrast" : kinds.size === 1 && kinds.has("analysis") ? "Analysis" : "Analysis / contrast";
      updateHeight();
    }
    function choose(contextId, collectionId, sectionId, sectionKind) {
      var context = find(contexts, contextId);
      if (!context) return null;
      var collection = find(context.collections, collectionId) || context.collections[0];
      if (!collection) return null;
      var section = find(collection.sections, sectionId) || collection.sections.find(function (x) {
        return has(sectionKind) && x.kind === sectionKind;
      }) || collection.sections[0];
      if (!section) return null;
      return { context: context.id, collection: collection.id, section: section.id };
    }
    function targetForChoice(choice) { return document.getElementById(choice.section); }
    function choiceForTarget(target) {
      if (!target) return null;
      for (var i = 0; i < contexts.length; i++) {
        var context = contexts[i];
        for (var j = 0; j < context.collections.length; j++) {
          var collection = context.collections[j];
          for (var k = 0; k < collection.sections.length; k++) {
            var section = collection.sections[k]; var node = document.getElementById(section.id);
            if (node && (node === target || node.contains(target)))
              return { context: context.id, collection: collection.id, section: section.id };
          }
        }
      }
      // Existing links may point to an analysis wrapper rather than a section.
      var contextWrapper = target.closest("[data-lisa-nav-context]");
      var collectionWrapper = target.closest("[data-lisa-nav-collection]");
      if (contextWrapper) return choose(contextWrapper.dataset.lisaNavContext,
        collectionWrapper ? collectionWrapper.dataset.lisaNavCollection : "", "", "");
      return null;
    }
    function hashChoice() {
      var id;
      try { id = decodeURIComponent(window.location.hash.slice(1)); } catch (_) { return null; }
      return id ? choiceForTarget(document.getElementById(id)) : null;
    }
    function writeHistory(section, replace) {
      var url = new URL(sectionHref(section), window.location.href);
      if (url.href === window.location.href) return;
      if (url.pathname !== window.location.pathname || url.search !== window.location.search) return;
      try { window.history[replace ? "replaceState" : "pushState"](null, "", url.href); }
      catch (_) { window.location.hash = url.hash; }
    }
    function apply(choice, navigate, historyMode, moveFocus) {
      if (!choice) return false;
      changing = true; state = choice; refresh();
      var result = selected();
      if (external) {
        changing = false;
        if (navigate) window.location.assign(sectionHref(result.section));
        return true;
      }
      var target = targetForChoice(choice);
      if (!target) { changing = false; return false; }
      var activeContext = target.closest("[data-lisa-nav-context]");
      var activeCollection = target.closest("[data-lisa-nav-collection]");
      document.querySelectorAll("details[data-lisa-nav-context]").forEach(function (node) { node.open = node === activeContext; });
      document.querySelectorAll("details[data-lisa-nav-collection]").forEach(function (node) { node.open = node === activeCollection; });
      // Full-report products may themselves be collapsible section targets.
      if (target.tagName === "DETAILS") target.open = true;
      for (var ancestor = target.parentElement; ancestor; ancestor = ancestor.parentElement)
        if (ancestor.tagName === "DETAILS") ancestor.open = true;
      var contextNode = target.closest("[data-lisa-context-json]");
      if (contextNode) {
        try { updateContext(JSON.parse(contextNode.dataset.lisaContextJson)); }
        catch (error) { console.error("Cannot read selected report context", error); }
      }
      if (historyMode) writeHistory(result.section, historyMode === "replace");
      if (navigate) startPending(target, moveFocus); else stopPending();
      changing = false;
      return true;
    }
    Object.keys(controls).forEach(function (key) {
      controls[key].addEventListener("change", function () {
        var before = selected(); var kind = before.section && before.section.kind;
        var choice = choose(controls.context.value, controls.collection.value,
          key === "section" ? controls.section.value : "", kind);
        // Keep keyboard focus in the picker, so repeated arrow-key selection
        // remains usable. Target content is already visible after navigation.
        apply(choice, true, "push", false);
      });
    });
    if (!external) {
      document.addEventListener("click", function (event) {
        if (event.defaultPrevented || event.ctrlKey || event.metaKey || event.shiftKey || event.altKey || event.button > 0) return;
        var link = event.target.closest("a[href]");
        if (!link || link.hasAttribute("download") || link.target === "_blank") return;
        var url = new URL(link.href, window.location.href);
        if (url.pathname !== window.location.pathname || url.search !== window.location.search || !url.hash) return;
        var id; try { id = decodeURIComponent(url.hash.slice(1)); } catch (_) { return; }
        var target = document.getElementById(id); var choice = choiceForTarget(target);
        if (!choice) return;
        event.preventDefault(); apply(choice, true, "push", true);
        // Preserve a deeper pre-existing anchor rather than replacing it with
        // the parent section's id (for example a specific figure download).
        if (target && target.id !== choice.section) {
          try { window.history.replaceState(null, "", url.href); } catch (_) { /* local file history may be unavailable */ }
          window.requestAnimationFrame(function () { target.scrollIntoView({ block: "start" }); });
        }
      });
      var restore = function () { var choice = hashChoice(); if (choice) apply(choice, true, null, false); };
      window.addEventListener("hashchange", restore); window.addEventListener("popstate", restore);
      document.addEventListener("toggle", function (event) {
        var node = event.target;
        if (changing || !node.open || !node.matches("details[data-lisa-nav-context],details[data-lisa-nav-collection]")) return;
        var contextNode = node.closest("[data-lisa-nav-context]");
        var collectionNode = node.closest("[data-lisa-nav-collection]");
        if (!contextNode) return;
        if (contextNode.dataset.lisaNavContext === state.context &&
            (!collectionNode || collectionNode.dataset.lisaNavCollection === state.collection)) return;
        apply(choose(contextNode.dataset.lisaNavContext,
          collectionNode ? collectionNode.dataset.lisaNavCollection : "", "", ""), false, "replace", false);
      }, true);
      window.addEventListener("scroll", function () {
        if (scrollTimer) return;
        scrollTimer = window.setTimeout(function () {
          scrollTimer = null;
          if (changing || scrolling) return;
          var choice = selected(); if (!choice.collection) return;
          var top = shell.getBoundingClientRect().bottom + 24;
          var sections = choice.collection.sections.map(function (section) {
            return { section: section, node: document.getElementById(section.id) };
          }).filter(function (row) { return row.node && row.node.getClientRects().length; });
          var above = sections.filter(function (row) { return row.node.getBoundingClientRect().top <= top; });
          var current = above.length ? above[above.length - 1] : sections[0];
          if (current && current.section.id !== state.section) {
            state.section = current.section.id; refresh(); writeHistory(current.section, true);
          }
        }, 100);
      }, { passive: true });
    }
    var initial = inventory.current;
    if (external) {
      if (initial && initial.context) apply(choose(initial.context, initial.collection,
        initial.section, initial.section), false, null, false);
      else refresh();
    } else {
      var fromHash = hashChoice();
      apply(fromHash || (initial && choose(initial.context, initial.collection, initial.section, initial.section)) ||
        choose(contexts[0].id, "", "", ""), !!fromHash, null, false);
    }
    bar.hidden = false; shell.setAttribute("data-navigation-ready", "true"); updateHeight();
    return {
      getState: function () { return { context: state.context, collection: state.collection, section: state.section }; },
      select: function (context, collection, section) { return apply(choose(context, collection, section, section), true, "push", false); },
      syncContext: function (update) {
        if (changing || !external || !update.exact_ids) return;
        var ids = update.exact_ids;
        var kind = ids.contrast_id ? "contrast" : "analysis";
        var id = ids.contrast_id || ids.analysis_id;
        // Gene evidence can explicitly return to all analyses or collections.
        // An empty filter is not the previous scope or the first scope.
        if (!id) {
          if (Object.prototype.hasOwnProperty.call(ids, "analysis_id") ||
              Object.prototype.hasOwnProperty.call(ids, "contrast_id")) {
            state = { context: "", collection: "", section: "" }; refresh();
          }
          return;
        }
        var context = contexts.find(function (row) {
          return row.kind === kind && (row.scientific_id || row.id) === String(id);
        });
        if (!context) return;
        var hasCollection = Object.prototype.hasOwnProperty.call(ids, "collection") ||
          Object.prototype.hasOwnProperty.call(update, "collection");
        var collection = Object.prototype.hasOwnProperty.call(ids, "collection") ? ids.collection : update.collection;
        if (hasCollection && (collection == null || collection === "")) {
          state = { context: context.id, collection: "", section: "" }; refresh(); return;
        }
        var before = selected(); var sectionKind = before.section ? before.section.kind : "";
        var choice = choose(context.id, String(collection || ""), "", sectionKind);
        if (choice) { state = choice; refresh(); }
      }
    };
  }
  function init() {
    var shell = shellNode();
    if (!shell || shell.getAttribute("data-shell-initialized") === "true") return;
    shell.setAttribute("data-shell-initialized", "true");
    var button = shell.querySelector(".lisa-shell-menu");
    var menu = shell.querySelector(".lisa-shell-nav");
    if (button && menu) {
      var close = function () {
        shell.removeAttribute("data-menu-open");
        button.setAttribute("aria-expanded", "false");
        updateHeight();
      };
      shell.setAttribute("data-menu-ready", "true");
      button.addEventListener("click", function () {
        var open = button.getAttribute("aria-expanded") !== "true";
        button.setAttribute("aria-expanded", String(open));
        shell.setAttribute("data-menu-open", String(open));
        updateHeight();
      });
      menu.addEventListener("click", function (event) {
        if (event.target.closest("a")) close();
      });
      shell.addEventListener("keydown", function (event) {
        if (event.key === "Escape" && button.getAttribute("aria-expanded") === "true") {
          close();
          button.focus();
        }
      });
    }
    var skip = document.querySelector(".lisa-shell-skip");
    if (skip) skip.addEventListener("click", function () {
      var target = document.getElementById(skip.getAttribute("href").slice(1));
      if (target) {
        if (!target.hasAttribute("tabindex")) target.setAttribute("tabindex", "-1");
        target.focus();
      }
    });
    if (window.ResizeObserver) new ResizeObserver(updateHeight).observe(shell);
    window.addEventListener("resize", updateHeight);
    reportNavigation = initReportNavigation(shell);
    window.lisaReportShell.navigation = reportNavigation;
    updateHeight();
    document.dispatchEvent(new CustomEvent("lisa:shell-ready"));
  }
  window.lisaReportShell = { updateContext: updateContext, setActive: setActive, init: init };
  document.addEventListener("lisa:context", function (event) { updateContext(event.detail); });
  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", init);
  else init();
}());
// Shared offline help, including dynamically rendered category sheets.
(function () {
  'use strict';
  function initHommelHelp() {
    const dialog = document.getElementById('hommel-support-help');
    if (!dialog || dialog.dataset.bound) return;
    dialog.dataset.bound = 'true';
    let opener = null;
    document.addEventListener('click', event => {
      const button = event.target.closest('[data-hommel-help]');
      if (button) { opener = button; dialog.showModal(); dialog.querySelector('[data-hommel-close]').focus(); }
      if (event.target.closest('[data-hommel-close]')) dialog.close();
    });
    dialog.addEventListener('close', () => { if (opener && opener.isConnected) opener.focus(); });
    dialog.addEventListener('click', event => {
      if (event.target !== dialog) return;
      const rect = dialog.getBoundingClientRect();
      if (event.clientX < rect.left || event.clientX > rect.right || event.clientY < rect.top || event.clientY > rect.bottom) dialog.close();
    });
  }
  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', initHommelHelp);
  else initHommelHelp();
}());
