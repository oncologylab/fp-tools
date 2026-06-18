# GUI Demo

`fp-tools-gui` is a Streamlit application. It is useful for interactive local exploration, but it requires a Python runtime and cannot run directly on GitHub Pages.

## Local GUI

```bash
pip install "fp-tools-bio[gui]"
fp-tools-gui
```

## Recommended Live Demo Hosting

Use one of these external app hosts and link to it from this page:

- **Streamlit Community Cloud** for the simplest Streamlit deployment.
- **Hugging Face Spaces** for a public app with pinned dependencies and example data.
- **Binder** for a notebook-like temporary environment.

For the live demo, keep the dataset small and read-only. The demo should show how to configure a run, inspect existing example outputs, and open static HTML reports.

## Static GitHub Pages Content

GitHub Pages should host:

- screenshots and workflow diagrams;
- standalone report HTML files;
- command/API documentation;
- links to a hosted GUI app.

GitHub Pages should not run full `fp-tools-gui` or compute-heavy workflows.

![fp-tools pseudobulk workflow](assets/fp-tools-pseudo-bulk.png)
