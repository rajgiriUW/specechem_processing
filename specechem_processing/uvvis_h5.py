import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import integrate as spint
from .uvvis import UVVis


def save_h5(data, filename):
    '''
    Saves UVVis data to an HDF5 file (.h5).

    Parameters
    ----------
    data : UVVis
        UVVis object containing spectra_vs_time, current, charge, and potentials.
    filename : str or Path
        File path to save to.
    '''

    if isinstance(filename, str):
        filename = Path(filename)

    with h5py.File(filename, 'w') as f:
        f.create_dataset('potentials', data=data.potentials)

        try:
            f.create_dataset('charge', data=data.charge.values[0])
        except:
            pass

        for p in data.spectra_vs_time:
            df = data.spectra_vs_time[p]
            grp = f.create_group(str(p))
            grp.create_dataset('data', data=df.values)
            grp.create_dataset('index', data=df.index.values)
            grp.create_dataset('columns', data=df.columns.values.astype(float))

        try:
            grp = f.create_group('current')
            grp.create_dataset('data', data=data.current.values)
            grp.create_dataset('index', data=data.current.index.values)
            grp.create_dataset('columns', data=data.current.columns.values.astype(float))
        except:
            pass

    return


def convert_h5(h5file):
    '''
    Loads a saved HDF5 file and returns a UVVis object.

    Parameters
    ----------
    h5file : str or Path
        Path to the HDF5 file to load.

    Returns
    -------
    UVVis
        UVVis object populated with spectra_vs_time, current, charge, and potentials.
    '''
    data = UVVis(None, None, None)

    with h5py.File(h5file, 'r') as file:
        data.potentials = file['potentials'][()]

        non_data = {'current', 'charge', 'potentials'}
        folders = [k for k in file.keys() if k not in non_data]

        try:
            folders_num = [float(p) for p in folders]
        except:  # for old 'x-1.0V' style, crops 'x'
            folders_num = [float(p[1:]) for p in folders]

        df_dict = {}
        for v, n in zip(folders, folders_num):
            grp = file[v]
            df = pd.DataFrame(data=grp['data'][()],
                               index=grp['index'][()],
                               columns=grp['columns'][()])
            df.index.name = 'Wavelength (nm)'
            df.columns.name = 'Time (s)'
            df_dict[n] = df
            data.tx = np.round(grp['columns'][()], 2)

        data.spectra_vs_time = df_dict

        grp = file['current']
        data.current = pd.DataFrame(data=grp['data'][()],
                                    index=grp['index'][()],
                                    columns=grp['columns'][()])
        data.current.index.name = 'Time (s)'
        data.current.columns.name = 'Potential (V)'

        data.charge = pd.DataFrame(data=file['charge'][()],
                                   index=data.potentials.T)

    return data

def _step_keys(h5_file):
    '''Numbered step groups, sorted numerically so 10 comes after 9, not after 1.'''
    return sorted((k for k in h5_file.keys() if k.isdigit()), key=int)


def _set_potentials(h5_file):
    '''Set potential of each step, from the potential_set attribute if present.'''
    potentials = []
    for k in _step_keys(h5_file):
        attrs = h5_file[k].attrs
        if 'potential_set' in attrs:
            v = float(attrs['potential_set'])
        else:
            v = float(h5_file[k]['echem']['potential'][0])
        potentials.append(np.round(v, 2) + 0.0)  # + 0.0 turns -0.0 into 0.0
    return potentials


def _doping_schedule(h5_file, n):
    '''The n doping potentials from the run's metadata_json, or None if absent.'''
    try:
        settings = json.loads(h5_file.attrs['metadata_json'])['settings']
        start = settings['doping_potential_start']
        step = settings['doping_potential_step']
    except (KeyError, TypeError, ValueError):
        return None

    return [np.round(start + i * step, 2) + 0.0 for i in range(n)]


def _paired_doping_path(h5_path):
    '''sample_dedoping.h5 -> sample_doping.h5 in the same folder, if it exists.'''
    h5_path = Path(h5_path)
    if 'dedoping' not in h5_path.name:
        return None
    candidate = h5_path.with_name(h5_path.name.replace('dedoping', 'doping'))
    return candidate if candidate.exists() else None


def convert_h5_python(h5_path, doping_path=None):
    '''
    Converts an h5 file saved by the Python version of the spectroelectrochemistry
    software into a UVVis object. Files are only ever opened read-only.

    The file format looks like this (one numbered group per potential step):
    /
    ├ 0                       attrs: potential_set, potential_measured, label, ...
      ├ absorbance_vs_time    (wavelength x time)
      ├ counts_vs_time
      ├ dark
      ├ echem
        ├ current             (A)
        ├ potential           (V)
        ├ time                (potentiostat clock, s)
      ├ reference
      ├ time                  (spectrometer time for each spectrum, s)
      ├ time_spectrometer
    ├ 1 ...
    ├ wavelength              (nm)

    Parameters
    ----------
    h5_path : str or Path
        Path to the .h5 file on disk.
    doping_path : str or Path, optional
        For a dedoping file: the paired doping .h5. Every dedoping step is held at
        the same potential, so each step is keyed by the potential of the doping
        step it follows, read from this file. If None, a sibling file with
        'dedoping' replaced by 'doping' in its name is used if it exists, then
        the doping schedule in the run metadata, then the step number.

    Returns
    -------
    UVVis
        * spectra_vs_time[potential]: index = wavelength (nm), columns = spectrum
          time (s, rounded to 0.01 s), values = absorbance.
        * potentials: the keys of spectra_vs_time, in step order.
        * current: index = potentiostat time (s), one column per potential (A).
        * charge: integrated current per potential (C).
    '''

    data = UVVis()
    spectra = {}
    current = {}

    with h5py.File(h5_path, 'r') as h5_file:
        keys = _step_keys(h5_file)
        wavelength = h5_file['wavelength'][()]

        potentials = _set_potentials(h5_file)
        if len(set(potentials)) < len(potentials):
            # e.g. dedoping: every step is at the same potential
            potentials = None
            doping_path = doping_path or _paired_doping_path(h5_path)
            if doping_path is not None:
                with h5py.File(doping_path, 'r') as doping_file:
                    potentials = _set_potentials(doping_file)
                if len(potentials) != len(keys):
                    raise ValueError(f'{doping_path} has {len(potentials)} steps but '
                                     f'{h5_path} has {len(keys)}')
            if potentials is None:
                potentials = _doping_schedule(h5_file, len(keys)) or list(range(len(keys)))

        tx = None
        for k, potential in zip(keys, potentials):
            times = np.round(h5_file[k]['time'][()], 2)
            spectra[potential] = pd.DataFrame(h5_file[k]['absorbance_vs_time'][()],
                                              index=wavelength, columns=times)

            if tx is None:
                tx = np.round(h5_file[k]['echem']['time'][()], 2)
            current[potential] = h5_file[k]['echem']['current'][()]

    data.spectra_vs_time = spectra
    data.potentials = np.array(potentials)
    data.current = pd.DataFrame(current, index=tx)

    data.charge = pd.DataFrame({p: [spint.trapezoid(data.current[p].values, x=tx)]
                                for p in data.current}, index=[0])
    data.charge.columns.name = 'Potential (V)'

    return data
