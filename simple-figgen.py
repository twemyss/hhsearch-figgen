#!/usr/bin/env python

import sys
import json
import ParsedHMM
import OutputFigure

# Check the arguments specified to the program
if len(sys.argv) is not 2:
    print("Invalid arguments. You must specify a configuration file as an argument to the program.")
    print("Usage: python hhsearch-figgen.py [config.json file]")
    sys.exit(1)

# Load the configuration
with open(sys.argv[1]) as json_file:
    config = json.load(json_file)

# Load the HMM of the target sequence, and extract the length


# Generate the output figure
output_figure = OutputFigure.OutputFigure(config)

