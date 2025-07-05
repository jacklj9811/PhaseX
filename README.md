# PhaseX

This repository contains MATLAB implementations of the 3D-GESPAR algorithm in the `3DGSPR` directory.
A lightweight Python port is available in the `py3dgespar` package. You can run a small
demonstration using:

```bash
python -m py3dgespar
```

Use the `--plot` flag to display a 1D comparison plot of the flattened
ground-truth signal and the recovered result (requires `matplotlib`):

```bash
python -m py3dgespar --plot
```

The Python code provides basic functions for experimenting with the 3D-GESPAR
approach. The 1D variant assumes an even signal length; providing an odd length
will raise a ``ValueError``.
