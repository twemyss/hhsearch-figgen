# hhsearch-figgen

This package generates PDF figures from HHSuite data. It requires a HHSearch log file, and HHSuite HMMs for each group of proteins mentioned in the log file. Configuration is via a JSON file.

## Getting Started

These instructions discuss the basic process of installation

### Installing

`hhsearch-figgen` requires Python 3, PyCairo, and NumPy.

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

You must run two HHsearches for this program, one using the PfamA 32.0 database, and the second using your database (in all our examples, called `eukarya`).

Example commands are:

```
hhsearch -i hmms/hmm_og.hhm -d database -cpu 8 -o results/search_og_database.txt -oa3m results_dir/search_og_database.a3m -all
```

You must run hhsearch for your protein group of interest, and have access to HMMs for all the hits that you want to include in your results.

### Generating a configuration file

`hhsearch-figgen` requires a configuration file to tell it about the figure you wish to generate. An example configuration file is located in `test-data/kkt17.json`.

The file is split into 6 sections:

##### 1: Master

This section governs the features of the main "master" HMM used in your figure.

`name` is a plain-text display name for the master HMM

`hmm_file` is the name of the master HMM file. It should be an absolute path (e.g. "/home/thomas/KKT17.hmm")

##### 2: Searches

These files are the HHSearch results files for your master.

`pfam` is the absolute path of the HHSearch against the PFAM database

`hits` is the absolute path of the HHSearch in your custom database

##### 3: Domains

It is rare to want to include all of the PFAM domains - many overlap and/or might not be relevant.

You should examine the hits file, and then choose the hit numbers of the domains that you wish to show on the figure.

This configuration section holds a JSON array of each domain you wish to include

```
"pfam_hit_number": [4, 5],
"name": "Coiled-Coil",
"colour": [0.9, 0.1, 0.1, 0.5]
```

In this example, hits number 4 and 5 would be labelled as "Coiled-Coil" and shown in a semi-transparent red colour.

##### 4: Page

The program is not currently smart when it comes to positioning the figure on the page or sizing the page. Tuning these parameters should help to get your figure shown fully.

The first two parameters define the page size: `horizontal_padding` sets how much space (left/right) of your master HMM to have on the page, whereas 
`height` sets the page height.

The second two parameters, `padding_left` and `padding_top`, set where the master HMM is drawn on the page. **You will likely need to adjust `padding_top` when you change the location of the master HMM between normal and split figure views (see end of section 5).**

##### 5: Output

`file_name` is the name of the PDF figure to save. It should typically end in '.pdf'.

`max_hits` is the maximum number of hits to show on each page./

`cutoff` is the E-value cutoff for hits to be shown in green (rather than grey)

The program can generate two types of figures. One has all the hits shown below the master ("normal view"), and the other has hits that hit before a certain amino acid position on the master above the master, and those that hit after that amino acid position below the master. This is called "split view", and can be enabled by setting the parameters `split` and `split_at`. The first is a boolean, false for normal and true for split. The second is the position threshold for the start of the hit. 

##### 6: Colours

Colours defines the colours to use in the profile HMM bitscore plots. It should typically not be changed - the current theme is based on the ClustalX colour scheme.