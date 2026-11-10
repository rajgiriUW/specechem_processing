'''
Tests for specechem_processing using h5 files saved by the Python
spectroelectrochemistry software (sample1_run3, PgbTTT in KTFSI/AcN).

The data are converted with uvvis_h5.convert_h5_python, and results are checked
against the raw arrays read directly from the h5 files, so no stored
expected-value files are needed.

Run from the repository root:  pytest tests/
'''

import hashlib
import shutil
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import pytest
from scipy import integrate as spint

from specechem_processing import uvvis_h5

DATA = Path(__file__).parent / 'test_h5'
FILES = {'doping': DATA / 'sample1_run3_doping.h5',
         'dedoping': DATA / 'sample1_run3_dedoping.h5'}
DOPING_POTENTIALS = [-0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]


def _sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def _raw(path):
    '''Raw arrays per step, in numeric group order, read straight from the h5.'''
    with h5py.File(path, 'r') as f:
        keys = sorted((k for k in f if k.isdigit()), key=int)
        return {'wavelength': f['wavelength'][()],
                'absorbance': [f[k]['absorbance_vs_time'][()] for k in keys],
                'time': [f[k]['time'][()] for k in keys],
                'echem_time': [f[k]['echem']['time'][()] for k in keys],
                'current': [f[k]['echem']['current'][()] for k in keys]}


@pytest.fixture(scope='module')
def hashes_before():
    return {name: _sha256(p) for name, p in FILES.items()}


@pytest.fixture(scope='module')
def converted(hashes_before):
    return {name: uvvis_h5.convert_h5_python(p) for name, p in FILES.items()}


@pytest.fixture(scope='module')
def raw():
    return {name: _raw(p) for name, p in FILES.items()}


@pytest.fixture(params=['doping', 'dedoping'])
def name(request):
    return request.param


# --- conversion ----------------------------------------------------------------

def test_conversion_does_not_modify_h5(converted, hashes_before):
    for name, p in FILES.items():
        assert _sha256(p) == hashes_before[name]


def test_potential_keys(converted, name):
    # Doping steps are keyed by set potential. Dedoping steps are all at -0.5 V,
    # so each is keyed by the doping potential it follows.
    data = converted[name]
    assert list(data.spectra_vs_time) == pytest.approx(DOPING_POTENTIALS)
    assert list(data.potentials) == pytest.approx(DOPING_POTENTIALS)
    assert list(data.current.columns) == pytest.approx(DOPING_POTENTIALS)


def test_dedoping_potentials_from_doping_file(converted):
    # Explicit doping_path gives the same keys as the automatic sibling lookup
    explicit = uvvis_h5.convert_h5_python(FILES['dedoping'], doping_path=FILES['doping'])
    assert list(explicit.spectra_vs_time) == list(converted['dedoping'].spectra_vs_time)
    assert list(explicit.spectra_vs_time) == list(converted['doping'].spectra_vs_time)


def test_dedoping_without_doping_file_uses_metadata(tmp_path):
    # Copy of the dedoping file alone in a temp folder: no doping file to pair with
    lone = tmp_path / FILES['dedoping'].name
    shutil.copyfile(FILES['dedoping'], lone)
    data = uvvis_h5.convert_h5_python(lone)
    assert list(data.spectra_vs_time) == pytest.approx(DOPING_POTENTIALS)


def test_spectra_match_raw_h5(converted, raw, name):
    data, r = converted[name], raw[name]
    for df, absorbance, t in zip(data.spectra_vs_time.values(), r['absorbance'], r['time']):
        assert df.shape == absorbance.shape
        np.testing.assert_array_equal(df.values, absorbance)
        np.testing.assert_array_equal(df.index.values, r['wavelength'])
        np.testing.assert_array_equal(df.columns.values, np.round(t, 2))


def test_current_and_charge(converted, raw, name):
    data, r = converted[name], raw[name]
    tx = np.round(r['echem_time'][0], 2)
    np.testing.assert_array_equal(data.current.index.values, tx)
    for p, current in zip(data.current.columns, r['current']):
        np.testing.assert_array_equal(data.current[p].values, current)
        # charge in coulombs
        assert data.charge[p].iloc[0] == pytest.approx(spint.trapezoid(current, x=tx))


# --- analysis ------------------------------------------------------------------

def test_time_index(converted, raw, name):
    data = converted[name]
    data.time_index()
    np.testing.assert_array_equal(data.tx, np.round(raw[name]['time'][0], 2))


def test_single_wl_time(converted, name):
    data = converted[name]
    df = data.spectra_vs_time[0.4]
    data.single_wl_time(potential=0.4, wavelength=800, smooth=3)

    row = df.iloc[df.index.searchsorted(800)]
    pd.testing.assert_series_equal(data.time_spectra, row)
    assert data.time_spectra_norm.min() == pytest.approx(0)
    assert data.time_spectra_norm.max() == pytest.approx(1)
    assert len(data.time_spectra_sm) == len(row)


def test_spec_echem_voltage(converted, name):
    data = converted[name]
    data.spec_echem_voltage(time=0, smooth=3)
    assert list(data.spectra.columns) == pytest.approx(DOPING_POTENTIALS)
    for p, df in data.spectra_vs_time.items():
        np.testing.assert_array_equal(data.spectra[p].values, df.iloc[:, 0].values)
    assert data.spectra_sm.shape == data.spectra.shape


def test_abs_vs_voltage(converted, name):
    data = converted[name]
    data.time_index()
    data.abs_vs_voltage(wavelength=800, time=0)
    assert list(data.vt.index) == pytest.approx(DOPING_POTENTIALS)
    for p, df in data.spectra_vs_time.items():
        assert data.vt.loc[p, 'Abs'] == pytest.approx(df.iloc[df.index.searchsorted(800), 0])


@pytest.mark.parametrize('fittype', ['exp', 'stretched'])
def test_banded_fits(converted, fittype):
    data = converted['doping']
    data.single_wl_time(potential=0.6, wavelength=800)
    data.banded_fits(wl_start=700, wl_stop=900, voltage=0.6, fittype=fittype)

    band = data.spectra_vs_time[0.6].loc[700:900]
    assert band.index.min() >= 700 and band.index.max() <= 900
    # banded_fits skips the first wavelength in the band
    assert len(data.fits) == len(band) - 1
    assert np.all(np.isfinite(np.asarray(data.fits, dtype=float)))


def test_banded_fits_bad_fittype(converted):
    data = converted['doping']
    data.single_wl_time(potential=0.6, wavelength=800)
    with pytest.raises(ValueError):
        data.banded_fits(voltage=0.6, fittype='linear')
