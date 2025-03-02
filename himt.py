#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse
import os
import datetime
from version import get_versions
from filter import filter
from assemble import assemble
from assess import assess
from compare import compare

class Logger():
    def __init__(self, filename="log.txt"):
        self.terminal = sys.stdout
        self.log = open(filename, "w", encoding='utf-8')  # 指定编码方式为 UTF-8
    def write(self, message):
        self.log.write(message)
        self.terminal.write(message)
        self.log.flush()
    def flush(self):
        pass

def main():
    version = get_versions()
    parser = argparse.ArgumentParser(description="a toolkit for assembling mitochondrial genome" + \
    "\n" + "version " + version + "\n" + \
    "for any questions, please put them forward on https://github.com/tang-shuyuan/himt or https://bioanno.com.",\
    usage="himt [function] [argument]",formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-v', '--version', action='version', version='HiMT ' + version)
    function_description="""
    assemble      Assemble mitochondrial genome with HiFi data 
    assess        Assess the assembly quality of the mitochondrial genome
    filter        Filter low-depth nuclear genome sequencing reads
    compare       Visualize the collinearity between two genomes
     """

    subparsers = parser.add_subparsers(title='function', description=function_description, metavar='')
    parser_assemble = subparsers.add_parser('assemble', description="Assembling mitochondrail genome with HiFi\
     sequcncing data ",usage="himt assemble [argument]"+"\n" + \
    "please use 'himt assemble -h or --help' to show help information",add_help=False)

    required_group = parser_assemble.add_argument_group('Required arguments')

    required_group.add_argument('-i', '--input_file', required=True, \
                                help='input a fasta or fastq file,gz compressed files are supported.')
    required_group.add_argument('-o', '--output_dir', required=True, \
                                help='output directory.')
    optional_group=parser_assemble.add_argument_group('Optional arguments')

    optional_group.add_argument('-h', '--help', action='help',default=argparse.SUPPRESS,
                                help='Show this help message and exit')
    optional_group.add_argument('-n', '--head_number', type=int, default=4, \
                                help='default=4,The number of kmer species randomly selected.')
    optional_group.add_argument('-t', '--thread', type=int, default=2, \
                                help='default=2 The number of thread used during code execution.')
    optional_group.add_argument('-b', '--base_number', type=int, default=3, choices=[3, 4], \
                                help='default=3,only can be 3 and 4,The number of bases at the beginning of kmer.')
    optional_group.add_argument('-fd', '--filter_depth', type=int, default=0, \
                                help='read depths below this value will be filtered.\
                                You must input the -p parameter to enable the use of the -fd parameter')
    optional_group.add_argument('-fp', '--filter_percentage', type=float, default=0.3, \
                                help='default=0.3,The depth of the mitochondrial genome obtained by blast, \
                                    the proportion adjusted downwards on this value.')
    optional_group.add_argument('-p', '--proportion', type=float, default=0, \
                                help='default=1,The percentage of the selected dataset from the entire file,\
                                    choose a value from 0-1.')
    optional_group.add_argument('-c', '--accuracy', type=float, default=0.8,
                                help='default=0.8,If one read has a high-frequency kmer ratio exceeding this value, \
                                    it will be considered as a high-frequency read,choose a value from 0-1.')
    optional_group.add_argument("-s", "--species", default="plant", choices=["plant", "animal"], \
                                help="default=plant,Species category,only can be plant or animal.")
    optional_group.add_argument('--no_flye_meta', action='store_true', help="By default, \
                                    we use flye Meta pattern to assemble the mitochondrial genome. \
                                    If you don't want to use meta pattern, add this parameter.")
    optional_group.add_argument('-x', "--normalize_depth", type=int, default=0, help="Normalize the \
                mitochondrial genome depth to a value.If the input value exceeds the mitochondrial genome depth, \
                retain the maximum mitochondrial genome depth,the default mitochondrial genome depth ranges between 15 and 50.\
                input a value less than 0 (such as:-1) to retain the maximum mitogenome depth")

    parser_assemble.set_defaults(func=assemble)

    parser_filter = subparsers.add_parser('filter',description="Filtering the low-depth nuclear genome reads",\
    usage="himt filter [argument]"+"\n"+"please use 'himt filter -h or --help' to show help information",add_help=False)
    required_group = parser_filter.add_argument_group('Required arguments')
    required_group.add_argument('-i', '--input_file', required=True, \
                                help='input a fasta or fastq file,gz compressed files are supported.')
    required_group.add_argument('-o', '--output_dir', required=True, \
                                help='output directory.')

    optional_group = parser_filter.add_argument_group('Optional arguments')
    optional_group.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                                help='Show this help message and exit')
    optional_group.add_argument('-n', '--head_number', type=int, default=4, \
                                help='default=4,The number of kmer species randomly selected.')
    optional_group.add_argument('-t', '--thread', type=int, default=2, \
                                help='default=2 The number of thread used during code execution.')
    optional_group.add_argument('-b', '--base_number', type=int, default=3, choices=[3, 4], \
                                help='default=3,only can be 3 and 4,The number of bases at the beginning of kmer.')
    optional_group.add_argument('-fd', '--filter_depth', type=int, default=0, \
                                help='read depths below this value will be filtered.\
                                    You must input the -p parameter to enable the use of the -fd parameter')
    optional_group.add_argument('-fp', '--filter_percentage', type=float, default=0.3, \
                                help='default=0.3,The depth of the mitochondrial genome obtained by blast,\
                                    the proportion adjusted downwards on this value.')
    optional_group.add_argument('-p', '--proportion', type=float, default=0,
                                help='The percentage of the selected dataset from the entire file,\
                                    choose a value from 0-1.')
    optional_group.add_argument('-c', '--accuracy', type=float, default=0.8,
                                help='default=0.8,If one read has a high-frequency kmer ratio exceeding this value, \
                                    it will be considered as a high-frequency read,choose a value from 0-1.')
    optional_group.add_argument("-s", "--species", default="plant", choices=["plant", "animal"], \
                                help="default=plant,Species category,only can be plant or animal.")
    optional_group.add_argument('-x', "--normalize_depth", type=int, default=0, help="Normalize the \
            mitochondrial genome depth to a value.If the input value exceeds the mitochondrial genome depth, \
            retain the maximum mitochondrial genome depth,the default mitochondrial genome depth ranges between 15 and 50.\
            input a value less than 0 (such as:-1) to retain the maximum mitogenome depth")
    parser_filter.set_defaults(func=filter)



    parser_compare = subparsers.add_parser('compare', description= \
        "compare the collinearity between two genomes based on alignment results from minimap2", \
        usage="himt compare [argument]" + "\n" + "please use 'himt compare -h or --help' to show help information",
        add_help=False)
    required_group = parser_compare.add_argument_group('Required arguments')
    required_group.add_argument('-r', '--reference', required=True, help='input reference genome')
    required_group.add_argument('-q', '--query', required=True, help='input query genome')
    required_group.add_argument('-o', '--output_dir', required=True, help='output directory')
    optional_group = parser_compare.add_argument_group('Optional arguments')
    optional_group.add_argument('-c', '--category', default='mitochondrial', \
    choices=['mitochondrial', 'chloroplast', 'other'], help='default=mitochondrial,choose the category of \
    genome you want to compare')
    # optional_group.add_argument('-f', '--filter',type=int, help='filter synteny block length less  \
    # than this value')
    optional_group.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                                help='Show this help message and exit')
    parser_compare.set_defaults(func=compare)

    parser_assess = subparsers.add_parser("assess",description="Assessing the assembly quality\
     of the mitochondrial genome",usage="himt assess [argument]"+
    "\n"+"please use 'himt assess -h or --help' to show help information")
    parser_assess.add_argument('-c','--category',default='mitochondrial',\
    choices=['mitochondrial','chloroplast'],help='default=mitochondrial,Choose the category of organelles you want to assess.')

    required_group = parser_assess.add_argument_group('Required arguments')

    required_group.add_argument('-i', '--input_file', required=True, \
                                help='input a fasta or gfa file.')
    required_group.add_argument('-o', '--output_dir', required=True, \
                                help='output directory.')
    parser_assess.set_defaults(func=assess)

    args = parser.parse_args()
    if len(vars(args)) == 0:
        parser.print_help()
        sys.exit(1)
    os.makedirs(args.output_dir, exist_ok=True)
    logger = Logger(os.path.join(args.output_dir, "HiMT.log"))
    sys.stdout = logger
    sys.stderr = logger
    print(f"({datetime.datetime.now()}) your running command: ", ' '.join(sys.argv))
    args.func(args)

if __name__ == "__main__":
    main()