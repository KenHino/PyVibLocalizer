"""Minimal setup file for tasks project."""

from setuptools import setup, find_packages
def _requires_from_file(filename):
    return open(filename).read().splitlines()

setup(
    name='pyvib',
    version='0.0.0',

    author='Kentaro Hino',
    url='https://kenhino.github.io/PyVibLocalizer/index.html',

    packages=find_packages(where='.'),
    package_dir={'': '.'},
    install_requires=_requires_from_file('requirements.txt')
)
