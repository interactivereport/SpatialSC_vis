from setuptools import setup
from setuptools import find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name = 'SpatialSCvis',
    version = '0.0.3',
    description = 'A package to enable spatial visualization for single-cell data',
    long_description = long_description,
    long_description_content_type = "text/markdown",
    classifiers = [
    "Programming Language :: Python :: 3",
    "License :: OSI Approved :: MIT License",
    "Operating System :: OS Independent",
],
    url = 'https://github.com/interactivereport/SpatialSC_vis',
    author = 'Wanli Wang',
    author_email = 'wendywangjk@gmail.com',
    license = 'MIT',
    install_requires = [
    "setuptools>=61.0",
    "numpy==1.20.0",
    "scanpy==1.8.2",
    "pandas==1.4.2",],
    packages=find_packages(),
    py_modules= ['SpatialSCvis'],
    package_dir = {'': 'src'},

)

