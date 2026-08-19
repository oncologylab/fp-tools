(() => {
  const setPageClasses = () => {
    const path = window.location.pathname.replace(/\/+$/, "/");
    const classes = {
      "fp-page-api-reference": /\/api\/$/.test(path),
      "fp-page-get-started": /\/get-started\//.test(path),
      "fp-page-tool-overview": /\/get-started\/tool-overview\/$/.test(path),
      "fp-page-command-guide": /\/get-started\/commands\/[^/]+\/$/.test(path),
    };

    Object.entries(classes).forEach(([name, enabled]) => {
      document.body.classList.toggle(name, enabled);
    });
  };

  setPageClasses();
  if (typeof document$ !== "undefined") {
    document$.subscribe(setPageClasses);
  }
})();
