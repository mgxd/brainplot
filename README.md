##README

1. Clone repo (including 32k Atlas)
2. Change paths at top of `om_py_plot_surfs.py` (add argparser to fix this)
3. Install [CIFTI](https://github.com/satra/nibabel/tree/enh/cifti2) version of nibabel (DISCLAIMER: best to make dedicated environment for this script)
4. Before running: `export QT_API=pyqt` in terminal
4. Run with commented command at top of script


###TO-DO
1. Add `mask` input to show brain coverage
2. Mutable image size
3. Handle flexible directory structure - Subject Level QA, Group Analysis, Etc
