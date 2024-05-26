from pathlib import Path
from setuptools import find_packages, setup

# Load version number
__version__ = ''
version_file = Path(__file__).parent.absolute() / 'chemfunc' / '_version.py'

with open(version_file) as fd:
    exec(fd.read())

# Load README
with open('README.md', encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='chemfunc',
    version=__version__,
    author='Kyle Swanson',
    author_email='swansonk.14@gmail.com',
    description='Chem Func',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/swansonk14/chemfunc',
    download_url=f'https://github.com/swansonk14/chemfunc/v_{__version__}.tar.gz',
    project_urls={
        'Source': 'https://github.com/swansonk14/chemfunc',
        'PyPi': 'https://pypi.org/project/chemfunc/'
    },
    license='MIT',
    packages=find_packages(),
    package_data={'chemfunc': ['py.typed']},
    entry_points={
        'console_scripts': [
            'chemfunc=chemfunc.main:main',
        ]
    },
    install_requires=[
        'descriptastorus',
        'matplotlib',
        'numpy',
        'pandas',
        'rdkit',
        'scikit-learn',
        'tqdm>=4.66.3',
        'typed-argument-parser>=1.9.0'
    ],
    python_requires='>=3.10',
    classifiers=[
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent'
    ],
    keywords=[
        'computational chemistry'
    ]
)
