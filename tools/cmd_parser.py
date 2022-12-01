# -*- coding: utf-8 -*-
'''
Python 3.10.6
12 September 2022
'''
import argparse

class CommandLine:
    def __init__(self):
        self.args = {'fasta':None, 'outdir':None}
        parser = argparse.ArgumentParser(description = "PfaSTer - pneumococcal fasta serotyper")
        parser.add_argument("-f", "--fasta", help = "path to genome (.fasta)", required = True, default = "")
        parser.add_argument("-o", "--output", help = "output directory (required to save results)", required = False, default = "")

        argument = parser.parse_args()
        if argument.fasta: self.args['fasta'] = argument.fasta
        if argument.output: self.args['outdir'] = argument.output
