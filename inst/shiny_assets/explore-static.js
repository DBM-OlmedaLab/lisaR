/* Reveals an attached exploration figure inside the category the reader has
   actually opened, in an exported bundle with R and Shiny closed.

   The evidence sheet renders one category at a time and mirrors the reader's
   choice into `?category=<id>` and `#category-<id>`. The exporter therefore
   emits one block per category, each carrying `data-lisa-explore-category` and
   the `hidden` attribute, and this script shows exactly the active one. Without
   it every category's figure would be visible on every category, which is the
   defect G1-4 describes.

   It reads the page; it never fetches, generates, or rewrites content. */
(function () {
  function activeCategory(doc) {
    var select = doc.getElementById('category');
    if (select && select.value) { return select.value; }
    var view = doc.defaultView || window;
    try {
      var query = new URLSearchParams(view.location.search).get('category');
      if (query) { return query; }
    } catch (error) { /* Some file viewers expose no search string. */ }
    var hash = (view.location.hash || '');
    if (hash.indexOf('#category-') === 0) {
      try { return decodeURIComponent(hash.slice(10)); } catch (error) { return hash.slice(10); }
    }
    return '';
  }
  function apply(doc) {
    // H3. A collection-wide block belongs to a whole analysis/collection or
    // contrast/collection and carries no category, so it is never revealed or
    // hidden by the active category. In practice it lands only on the native
    // collection pages, which are not category-scoped and do not load this
    // script at all; excluding it here makes that independent of where it lands.
    var blocks = doc.querySelectorAll(
      '[data-lisa-explore-category]:not([data-lisa-explore-scope="collection"])');
    if (!blocks.length) { return; }
    var category = activeCategory(doc);
    Array.prototype.forEach.call(blocks, function (node) {
      var mine = node.getAttribute('data-lisa-explore-category') === category;
      // Only write when the state actually changes, so switching category does
      // not touch the nodes that were already correct.
      if (node.hidden !== !mine) { node.hidden = !mine; }
    });
  }
  function start() {
    var doc = document;
    apply(doc);
    var info = doc.getElementById('category-info');
    if (info && window.MutationObserver) {
      // The sheet re-renders #category-info on every category change, including
      // the first one, and it records the change with history.replaceState --
      // which fires neither popstate nor hashchange. Observing that one element
      // is the only reliable signal, and it costs nothing while the reader stays
      // on one category.
      new MutationObserver(function () { apply(doc); }).observe(info, { childList: true });
    }
    var select = doc.getElementById('category');
    if (select) { select.addEventListener('change', function () { apply(doc); }); }
    window.addEventListener('hashchange', function () { apply(doc); });
    window.addEventListener('popstate', function () { apply(doc); });
  }
  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', start);
  } else {
    start();
  }
})();
