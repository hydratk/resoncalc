# -*- coding: utf-8 -*
"""Input

.. module:: input
   :synopsis: Module for processing input

"""

import resoncalc.output as output
import resoncalc.detection as detection
from argparse import ArgumentParser
from os import path
from json import load, JSONDecodeError

def process_command():
    """Process command

    usage: resoncalc [-h] [-o OUTPUT] [-v] [-s] [-g] [-t TITLE] input

    Calculate bound states and resonances for potential

    positional arguments:
      input                 input file with computation settings

    options:
      -h, --help            show this help message and exit
      -o OUTPUT, --output OUTPUT
                            output directory
      -v, --verbose         verbose mode
      -s, --silent          silent mode
      -g, --generate        generate graphs from data
      -t, --title           output title, used for generate

    """     

    # command options
    parser = ArgumentParser(description='Calculate bound states and resonances for potential')
    parser.add_argument('string', metavar='input', help='input file with computation settings')
    parser.add_argument('-o', '--output', help='output directory')
    parser.add_argument('-v', '--verbose', action='store_true', help='verbose mode')
    parser.add_argument('-s', '--silent', action='store_true', help='silent mode')
    parser.add_argument('-g', '--generate', action='store_true', help='generate graphs from data')
    parser.add_argument('-t', '--title', help='output title, used for generate')
    args = parser.parse_args()

    # verbose
    verbose = args.verbose
    if (verbose):
        output.log_level = 2

    # silent
    silent = args.silent
    if (silent):
        output.stdout = False

    # generate
    generate = args.generate

    # title
    title = args.title

    # output
    outdir = args.output 
    if (outdir is not None):
        if (path.exists(outdir)):
            output.outdir = outdir
        else:
            print('Output directory {0} not found'.format(outdir))
            return -1

    # input
    infile = args.string
    computations = None
    if (path.exists(infile)):
        if (generate):
            output.generate_graphs(infile, title)
            return 0
        else:
            computations = load_input(infile)
    else:
        print('Input file {0} not found'.format(infile))
        return -1

    # computation processing
    if (computations is not None):
        try:
            process_computations(computations)
            return 0
        except KeyboardInterrupt as ex:
            print('Program terminated by user')
            return -1

def load_input(fname):
    """Load input computations file

    Args:
        fname (str): input filename

    Returns:
        list: computations

    """     

    try:
        with open(fname, 'r') as f:
            computations = load(f)
            if (type(computations) is not list):
                computations = [computations]

            return computations
        
    except JSONDecodeError as ex:
        print('Failed to parse input file {0}: {1}'.format(fname, ex))
        return None

def process_computations(computations):
    """Process computations

    Args:
        computations (list): computations settings

    """        

    for computation in computations:
        detection.perform_detection_loop(computation)
