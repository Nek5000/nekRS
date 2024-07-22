# nekRS Documentation

This subdirectory contains documentation of the [nekRS](https://github.com/Nek5000/nekRS) project
using the [Sphinx](http://www.sphinx-doc.org/) documentation framework. A read the docs website
is hosted [here](https://nekrs.readthedocs.io/en/latest/).

## How to build locally

If you are developing the documentation you can preview the website locally before a pull request is
merged.

### Dependencies

The documentation requires a python 3 installation and a few pip packages including 
[Sphinx](https://pypi.python.org/pypi/Sphinx) and [sphinx_rtd_theme](https://pypi.python.org/pypi/sphinx_rtd_theme) Python packages.

The recommended way to setup the build environment is to create a local venv and 
then install the packages via the requirement.txt file.

```sh
cd doc
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

This environment can be disabled with `deactivate`. It can re-enabled from any 
terminal with `source $NEKRS_HOME/doc/venv/bin/activate`

### Building

Once the dependencies have been acquired `make html` (from the doc directory) builds 
the user documentation as a set of interlinked HTML and image files. The top-level 
webpage is `build/html/index.html`. To view this documentation as a navigable HTML
web page, simply navigate to the `build/html/index.html` file in your file system
and open with a web browser.

## How to contribute

Please create a fork of the repository and make pull/merge requests. Keep in 
mind that the number of binary files should be kept minimal. The Makefile should be 
adapted to any special build requirements.

New issues or requests are welcome to be reported.
