import os
import codecs
import setuptools


# https://packaging.python.org/guides/single-sourcing-package-version/
def read(rel_path):
    here = os.path.abspath(os.path.dirname(__file__))
    with codecs.open(os.path.join(here, rel_path), 'r') as fp:
        return fp.read()


def get_version(rel_path):
    for line in read(rel_path).splitlines():
        if line.startswith('__version__'):
            delim = '"' if '"' in line else "'"
            return line.split(delim)[1]
    else:
        raise RuntimeError("Unable to find version string.")


with open("README.md", "r") as fh:
    long_description = fh.read()

name = 'planetMagFields'
version = get_version('planetmagfields/__init__.py')
description = '''Routines to plot magnetic fields
 of planets in our solar system '''
copyright = '2025 Ankit Barik'


setuptools.setup(
    name=name,
    version=version,
    author='Barik, Ankit',
    author_email='abarik@jhu.edu ',
    packages=['planetmagfields'],
    description=description,
    long_description=long_description,
    install_requires=[
        'numpy>=1.18',
        'scipy>=1.5.4',
        'matplotlib>=3',
        ],
    package_data={'planetmagfields': ['data/*']},
    include_package_data=True,
)
