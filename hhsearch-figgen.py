#!/usr/bin/env python
import sys
import json
import ParsedHMM
from OutputFigure import OutputFigure
from argparse import ArgumentParser

if __name__ == "__main__":
    """
    Main command line interface for the software. This uses argparse to validate the arguments
    and pass them to the figure generation script. It requires a configuration JSON file that
    contains details of the hits, regions to include, and the secondary structure data. An example
    of this can be seen in test-data/kkt17.json.
    """

    # Use argparse to parse the arguments
    parser = ArgumentParser(description="Generate plots of hhsearch results and secondary structure data")
    parser.add_argument('config', help="The path to the JSON configuration file. See test-data/kkt17.json for an example.")
    arguments= parser.parse_args()

    # Load the configuration
    with open(arguments.config) as json_file:
        config = json.load(json_file)

    # Generate the output figure
    output_figure = OutputFigure(config, True)

