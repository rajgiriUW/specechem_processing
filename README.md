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

## Installation
```
git clone https://github.com/rajgiriUW/specechem_processing.git
cd specechem_processing
pip install -e .
```

## Quick start
```python
from specechem_processing import uvvis, read_files

steps, specs, volts, dedope_steps, dedope_specs = read_files.read_files('path/to/data')
data = uvvis.UVVis(steps, specs, volts)
data.time_dep_spectra(specfiles=specs)
data.spec_echem_voltage(time=0)
```

A worked example is in `notebooks/UV-Vis Notebook_updated.ipynb`.

## Tests
Run from the repository root:
```
pytest tests/
```
See `tests/uvvis_test/README.md` for how the expected values are generated.
