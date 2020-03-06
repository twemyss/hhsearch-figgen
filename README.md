# hhsearch-figgen

This package generates PDF figures from HHSuite data. It requires a HHSearch log file, and HHSuite HMMs for each group of proteins mentioned in the log file. Configuration is via a JSON file.

## Getting Started

These instructions discuss the basic process of installation

### Installing

hhsearch-figgen requires Python 3, PyCairo, and NumPy.

The Python package manager `pip` is the easiest way to install the required modules.

If you don't have `pip` installed, you can install it by following the instructions [here](https://pip.pypa.io/en/stable/installing/).

Installing PyCairo and Numpy is then as simple as running the following commands.

```
pip3 install pycairo
pip3 install numpy
```

This may fail on some systems. If so, you should install Cairo, and try again. This can be done in the following ways:

MacOS: [Homebrew](https://brew.sh/) is the easiest way to install Cairo. Once Homebrew [is installed](https://brew.sh/), Cairo can be installed as follows:

```
brew install cairo
brew install py3cairo
```
Once this is done, follow the `pip3` instructions above and the program should install successfully.

### Running HHSearch and making HMMs


### 