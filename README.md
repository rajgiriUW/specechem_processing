# Spectroelectrochemistry processing

```
Contact info:

Rajiv Giridharagopal, Ph.D.
University of Washington
Department of Chemistry
rgiri@uw.edu
```

Python tools for processing spectroelectrochemistry data, used mainly by the Ginger Lab at the University of Washington. This code was split out of [OECT_processing](https://github.com/rajgiriUW/OECT_processing), which now handles OECT device processing only.

## Modules
* `uvvis`: the `UVVis` class for time-dependent UV-Vis spectra recorded at stepped potentials (absorbance vs. voltage, single-wavelength kinetics, charge and current vs. time, banded exponential fits)
* `uvvis_h5`: load spectroelectrochemistry `.h5` files into `UVVis` objects
* `uvvis_plot`: plotting helpers for `UVVis` data
* `read_files`: find and sort the `steps(N).txt` / `spectra(N).txt` files in a data folder
* `cv`: cyclic voltammetry loading and analysis
* `eis`: electrochemical impedance spectroscopy loading and plotting

## Requirements
Python 3.11+ and pandas 3 (see `requirements.txt`).

## Installation
The test data are stored with [Git LFS](https://git-lfs.com), so install it once (`git lfs install`) before cloning:
```
git clone https://github.com/rajgiriUW/specechem_processing.git
cd specechem_processing
pip install -e .
```

## Quick start
Load the `.h5` files saved by the Python spectroelectrochemistry software:
```python
from specechem_processing import uvvis_h5

data = uvvis_h5.convert_h5_python('sample1_run3_doping.h5')
dedata = uvvis_h5.convert_h5_python('sample1_run3_dedoping.h5')  # pairs with sample1_run3_doping.h5

data.single_wl_time(potential=0.6, wavelength=800)
data.spec_echem_voltage(time=0)
data.banded_fits(wl_start=700, wl_stop=900, voltage=0.6)
```
`spectra_vs_time` is keyed by potential (V). Dedoping steps are all held at the same potential, so each is keyed by the doping step it follows, read from the paired doping file (pass `doping_path=` if it isn't next to the dedoping file with the same name). `current` is in A and `charge` is in C.

A worked example is in `notebooks/UV-Vis Notebook_updated.ipynb`.

## Tests
Run from the repository root:
```
pytest tests/
```
The tests use `tests/test_h5/sample1_run3_doping.h5` and `sample1_run3_dedoping.h5` and check results against the raw arrays in those files.
