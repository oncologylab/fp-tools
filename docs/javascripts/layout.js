(() => {
  const normalizePath = (path) => {
    const withoutIndex = path.replace(/\/index\.html$/, "/");
    return withoutIndex.replace(/\/+$/, "/") || "/";
  };

  const setPageClasses = () => {
    const path = normalizePath(window.location.pathname);
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

  const setHeaderTabState = () => {
    const items = [...document.querySelectorAll(".fp-header-tabs .md-tabs__item")];
    if (!items.length) return;

    const currentPath = normalizePath(window.location.pathname);
    const links = items.map((item) => item.querySelector(".md-tabs__link"));
    const destinations = links.map((link) =>
      link ? normalizePath(new URL(link.href, window.location.href).pathname) : ""
    );
    let activeIndex = destinations.findIndex(
      (destination, index) => index > 0 && currentPath.startsWith(destination)
    );
    if (activeIndex < 0) activeIndex = 0;

    items.forEach((item, index) => {
      const active = index === activeIndex;
      item.classList.toggle("md-tabs__item--active", active);
      const link = links[index];
      if (!link) return;
      if (active) link.setAttribute("aria-current", "page");
      else link.removeAttribute("aria-current");
    });
  };

  const updateLayout = () => {
    setPageClasses();
    setHeaderTabState();
  };

  updateLayout();
  if (typeof document$ !== "undefined") {
    document$.subscribe(updateLayout);
  }
})();
