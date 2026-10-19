// Open every content link in a new tab (in-page anchors and heading permalinks excluded).
(function () {
  function apply() {
    document.querySelectorAll('.md-content a[href]').forEach(function (a) {
      var href = a.getAttribute('href');
      if (!href || href.charAt(0) === '#' || a.classList.contains('headerlink')) return;
      a.setAttribute('target', '_blank');
      a.setAttribute('rel', 'noopener');
    });
  }
  if (typeof document$ !== 'undefined' && document$.subscribe) {
    document$.subscribe(apply);
  } else {
    document.addEventListener('DOMContentLoaded', apply);
  }
})();
