from setuptools import find_packages, setup

__version__ = '0.0.1'

# Load README
with open('README.md', encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='chem_utils',
    author='Kyle Swanson',
    author_email='swansonk.14@gmail.com',
    description='Chem Utils',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/swansonk14/chem_utils',
    download_url=f'https://github.com/swansonk14/chem_utils/v_{__version__}.tar.gz',
    project_urls={
        'Source': 'https://github.com/swansonk14/chem_utils',
        'PyPi': 'https://pypi.org/project/chem_utils/'
    },
    license='MIT',
    packages=find_packages(),
    package_data={'chem_utils': ['py.typed']},
    install_requires=[
        'matplotlib',
        'numpy',
        'pandas',
        'scikit-learn',
        'scikit-learn-intelex',
        'tqdm',
        'typed-argument-parser',
        'umap-learn'
    ],
    python_requires='>=3.10',
    classifiers=[
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.10',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent'
    ],
    keywords=[
        'computational chemistry'
    ]
)
