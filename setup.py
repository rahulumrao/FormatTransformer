# setup.py
from setuptools import setup, find_packages
'''
setup file
'''

setup(
    name='FormatTransformer',                 # package name
    version='0.1',                             # version
    description='Convert XYZ trajectories to VASP POSCAR files and center structures',
    long_description=open('README.md').read(),
    author='Rahul Verma',
    author_email='rverma7@ncsu.edu',
    packages=find_packages(include=["src"]),
    install_requires=['numpy', 'ase', 'MDAnalysis', 'tqdm'],
    classifiers = [
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',
        ],
    python_requires = '>=3.0',
    keywords = 'FormatTransformer',
    url = '',
    docs_url = '',
    entry_points={                             # Command-line scripts
        'console_scripts': [
            'xyzcenter = src.center_molecule:main',
            'xyzconverter = src.writer:main',
        ],
    },
)
