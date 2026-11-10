from setuptools import setup, find_packages

from os import path
this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='specechem_processing',
    version='0.1.0',
    description='Spectroelectrochemistry (UV-Vis), cyclic voltammetry, and EIS processing',
    long_description=long_description,
    long_description_content_type='text/markdown',

    author='Rajiv Giridharagopal',
    author_email='rgiri@uw.edu',
    license='MIT',
    url='https://github.com/rajgiriUW/specechem_processing',

    packages=find_packages(exclude=["*.tests", "*.tests.*", "tests.*", "tests"]),
    install_requires=['numpy >=1.18',
                      'pandas >=2.0.1',
                      'scipy >=1.14',
                      'matplotlib >=3.6.0',
                      'seaborn',
                      'h5py >= 3.2.0',
                      ],
)
