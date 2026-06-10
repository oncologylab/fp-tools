# Cover Letter

Dear Editor,

On behalf of my co-author, I am pleased to submit our manuscript, **"fp-tools: A
Reproducible Command-Line Platform for Classical and Multiscale ATAC-seq
Footprinting, Supervised TF-Binding Prediction, Motif Discovery, and Variant
Scoring,"** for consideration as an Article in *BioMedInformatics*.

Transcription-factor footprinting from ATAC-seq is a central technique in
regulatory genomics, but its surrounding software ecosystem is fragmented across
single-purpose tools with heterogeneous interfaces, inconsistent reproducibility,
and uneven test coverage. This raises the practical cost of building and
benchmarking end-to-end analyses. Our manuscript presents `fp-tools`, a standalone,
pip-installable package that unifies the core footprinting workflows behind stable
command-line entry points and adds a series of explicitly opt-in scientific
modules—multiscale and nucleosome-aware scoring, supervised TF-binding prediction,
motif-relaxed/motif-free recovery, de novo motif discovery, variant scoring,
single-cell pseudobulk aggregation, replicate-aware differential-binding
uncertainty reporting, and a competition-aware decomposition of overlapping
TF-scale and nucleosome-scale footprint signals—while leaving the classical
commands unchanged when they are not invoked.

Using public ENCODE data across multiple transcription factors and cell types, we
demonstrate that the package reproducibly recovers TF occupancy and that integrating
sequence, accessibility, and cut-site footprint evidence improves discrimination
over single-signal baselines. We believe this work fits the scope of
*BioMedInformatics* because it emphasizes reproducible, well-engineered informatics
infrastructure for biomedical data analysis: deterministic regression testing, a
tiered public-data benchmark design with standard classification and calibration
metrics, and publication-ready figure generation, all built on open data and
open-source software.

We confirm that this manuscript is original, has not been published previously, and
is not under consideration for publication elsewhere. All authors have read and
approved the manuscript and agree to its submission. The authors declare no
conflicts of interest.

Thank you for your consideration.

Sincerely,

Chunling Yi, Ph.D. (Corresponding Author)
Professor, Department of Oncology
Lombardi Comprehensive Cancer Center
Georgetown University Medical Center, Washington, DC, USA
chunling.yi@georgetown.edu
