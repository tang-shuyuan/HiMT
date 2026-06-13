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
from refassemble import refassemble
from view import view

class Logger():
    def __init__(self, filename="log.txt"):
        self.terminal = sys.stdout
        self.log = open(filename, "w", encoding='utf-8')
    def write(self, message):
        self.log.write(message)
        self.terminal.write(message)
        self.log.flush()
    def flush(self):
        pass


def main():
    version = get_versions()
    parser = argparse.ArgumentParser(description="An Integrative Toolkit for Assembling Organelle Genomes " + \
    "\n" + "version " + version + "\n" + \
    "For any questions, please submit them via https://github.com/tang-shuyuan/HiMT or https://bioanno.com.",\
    usage="himt [function] [argument]",formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-v', '--version', action='version', version='HiMT ' + version)
    function_description="""
    assemble      Assemble mitochondrial genome with HiFi data 
    assess        Assess the assembly quality of the mitochondrial genome
    filter        Filter low-depth nuclear genome sequencing reads
    compare       Visualize the collinearity between two genomes
    refassemble   Reference-based extraction and assembly of reads
    view          Visualize the graph of a GFA file
     """

    subparsers = parser.add_subparsers(title='function', description=function_description, metavar='')
    parser_assemble = subparsers.add_parser('assemble', description="Assemblie mitochondrail genome with HiFi\
     sequcncing data ",usage="himt assemble [argument]"+"\n" + \
    "please use 'himt assemble -h or --help' to show help information",add_help=False)

    required_group = parser_assemble.add_argument_group('Required arguments')

    required_group.add_argument('-i', '--input_file', required=True, \
                                help='input a fasta or fastq file,gz compressed files are supported.')
    required_group.add_argument('-o', '--output_dir', required=True, \
                                help='output directory.')

    optional_group=parser_assemble.add_argument_group('Optional arguments')

    optional_group.add_argument('-h', '--help', action='help',
                                help='Show this help message and exit')
    optional_group.add_argument("-s", "--species", default="plant", choices=["plant", "animal"], \
                                help="default=plant,Species category,only can be plant or animal.")
    optional_group.add_argument('-d', '--data_type', default='HiFi', choices=["HiFi", "CLR", "ONT"], \
                                help=" default=HiFi,Choose your sequencing data type (HiFi,CLR or ONT)")
    optional_group.add_argument('-k', '--kmer_length', type=int,default=21,
                                help="default=21")
    optional_group.add_argument('-n', '--head_number', type=int, default=4, \
                                help='default=4,The number of kmer species randomly selected.')
    optional_group.add_argument('-t', '--thread', type=int, default=2, \
                                help='default=2 The number of threads used for flye assembly.')
    optional_group.add_argument('-e', '--extract_parallel', type=int, default=2,
                                help='default=2,The number of k-mers processed in parallel.')

    optional_group.add_argument('-b', '--base_number', type=int, default=3, choices=[3, 4], \
                                help='default=3,only can be 3 and 4,The number of bases at the beginning of kmer.')
    optional_group.add_argument('-fd', '--filter_depth', type=int, default=0, \
                                help='Read depths below this value will be filtered.\
                                You must input the -p parameter to enable the use of the -fd parameter')
    optional_group.add_argument('-fp', '--filter_percentage', type=float, default=0.3, \
                                help='default=0.3,The depth of the mitochondrial genome obtained by blast, \
                                    the proportion adjusted downwards on this value.')
    optional_group.add_argument('-p', '--proportion', type=float, default=0, \
                                help='default=0,The percentage of the selected dataset from the entire file,\
                                    choose a value from 0-1.')
    optional_group.add_argument('-c', '--accuracy', type=float, default=0.8,
                                help="default=0.8,choose a value between 0 and 1 .If a read’s ratio of high-frequency \
                                k-mers exceeds the threshold , it will be classified as a high-frequency read.")
    optional_group.add_argument('--no_flye_meta', action='store_true', help="By default, \
                                    we use flye Meta pattern to assemble the mitochondrial genome. \
                                    If you don't want to use meta pattern, add this parameter.")
    optional_group.add_argument('-x', "--normalize_depth", type=int, default=0, help="Normalize the \
                mitochondrial genome depth to a value.If the input value exceeds the mitochondrial genome depth, \
                retain the maximum mitochondrial genome depth,the default mitochondrial genome depth ranges between 15 and 50.\
                input a value less than 0 (such as:-1) to retain the maximum mitogenome depth")
    optional_group.add_argument('--__internal_seed', type=int, default=0, help=argparse.SUPPRESS)

    parser_assemble.set_defaults(func=assemble)

    parser_filter = subparsers.add_parser('filter',description="Filter the low-depth nuclear genome reads",\
    usage="himt filter [argument]"+"\n"+"please use 'himt filter -h or --help' to show help information",add_help=False)
    required_group = parser_filter.add_argument_group('Required arguments')
    required_group.add_argument('-i', '--input_file', required=True, \
                                help='input a fasta or fastq file,gz compressed files are supported.')
    required_group.add_argument('-o', '--output_dir', required=True, \
                                help='output directory.')


    optional_group = parser_filter.add_argument_group('Optional arguments')
    optional_group.add_argument('-h', '--help', action='help',
                                help='Show this help message and exit')
    optional_group.add_argument('-d', '--data_type', default='HiFi', choices=["HiFi", "CLR", "ONT"],\
                                help=" default=HiFi,Choose your sequencing data type (HiFi,CLR or ONT)")
    optional_group.add_argument('-k', '--kmer_length', type=int, default=21,
                                help=" default=21")
    optional_group.add_argument('-n', '--head_number', type=int, default=4, \
                                help='default=4,The number of kmer species randomly selected.')
    optional_group.add_argument('-t', '--thread', type=int, default=2, \
                                help='default=2 The number of threads used for flye assembly.')
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
                                help="default=0.8,choose a value between 0 and 1 .If a read’s ratio of high-frequency \
                                    k-mers exceeds the threshold , it will be classified as a high-frequency read.")
    optional_group.add_argument("-s", "--species", default="plant", choices=["plant", "animal"], \
                                help="default=plant,Species category,only can be plant or animal.")
    optional_group.add_argument('-x', "--normalize_depth", type=int, default=0, help="Normalize the \
            mitochondrial genome depth to a value.If the input value exceeds the mitochondrial genome depth, \
            retain the maximum mitochondrial genome depth,the default mitochondrial genome depth ranges between 15 and 50.\
            input a value less than 0 (such as:-1) to retain the maximum mitogenome depth")
    optional_group.add_argument('-e', '--extract_parallel', type=int, default=2,
                                help='default=2,The number of k-mers processed in parallel.')
    optional_group.add_argument('--__internal_seed', type=int, default=0,help=argparse.SUPPRESS)


    parser_filter.set_defaults(func=filter)



    parser_compare = subparsers.add_parser('compare', description= \
        "Visualize the collinearity between two genomes", \
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
    optional_group.add_argument('-h', '--help', action='help',
                                help='Show this help message and exit')

    parser_compare.set_defaults(func=compare)

    parser_assess = subparsers.add_parser("assess",description="Assess the assembly quality\
     of the mitochondrial genome",usage="himt assess [argument]"+
    "\n"+"please use 'himt assess -h or --help' to show help information",add_help=False)

    required_group = parser_assess.add_argument_group('Required arguments')

    required_group.add_argument('-i', '--input_file', required=True, \
                                help='input a fasta or gfa file.')
    required_group.add_argument('-o', '--output_dir', required=True, \
                                help='output directory.')

    optional_group = parser_assess.add_argument_group('Optional arguments')
    optional_group.add_argument('-h', '--help', action='help',
                                help='Show this help message and exit')
    optional_group.add_argument("-s", "--species", default="plant", choices=["plant", "animal"],
                                help="default=plant,Species category,only can be plant or animal.")
    optional_group.add_argument('-c', '--category', choices=['mitochondrial', 'chloroplast'],
                                help='Choose the organelle type for analysis. \
                                HiMT will automatically identify it if left unspecified.')

    parser_assess.set_defaults(func=assess)

    parser_refassemble = subparsers.add_parser('refassemble',
    description= "Obtain data via reference genome and assemble",
    usage="himt compare [argument]" + "\n" + "please use 'himt refassemble -h or --help' to show help information",
                                           add_help=False)
    required_group = parser_refassemble.add_argument_group('Required arguments')
    required_group.add_argument('-i','--input_file',required=True,help='Input file')
    required_group.add_argument('-r', '--ref_genome', required=True, help='Input reference genome')
    required_group.add_argument('-o','--output_dir',required=True,help='Output directory')

    optional_group = parser_refassemble.add_argument_group('Optional arguments')
    optional_group.add_argument('-h', '--help', action='help',
                                help='Show this help message and exit')
    optional_group.add_argument('-d', '--data_type', default='HiFi', choices=["HiFi", "CLR", "ONT"],\
                                help=" default=HiFi,Choose your sequencing data type (HiFi,CLR or ONT)")

    optional_group.add_argument('-c', '--category', choices=['mitochondrial', 'chloroplast'],
                                help='select the category of organelles')

    optional_group.add_argument("-s", "--species", default="plant", choices=["plant", "animal"],\
                                help="Species category,only can be plant or animal (default=plant)")
    optional_group.add_argument('-l','--length',type=int,default=0,
                                help='filter reads with mapping length below this value (default: 0)')
    optional_group.add_argument('-p', '--percent',type=float,default=0.6,
                                help='extract reads with mapping length pencentage >= this value (default: 0.6)')
    optional_group.add_argument('-t', '--thread', type=int, default=2, \
                                help='default=2 ,The number of threads used for flye assembly')
    optional_group.add_argument('--no_flye_meta', action='store_true', help="By default, \
                                        we use flye Meta pattern to assemble the mitochondrial genome. \
                                        If you don't want to use meta pattern, add this parameter.")
    parser_refassemble.set_defaults(func=refassemble)

    parser_view = subparsers.add_parser(
        'view',
        description="Visualize the graph of a GFA file",
        usage="himt view [argument]" + "\n" + "please use 'himt view -h or --help' to show help information",
        add_help=False
    )
    required_group = parser_view.add_argument_group('Required arguments')
    required_group.add_argument('-i', '--input_file', required=True, help='input a gfa file.')
    required_group.add_argument('-o', '--output_dir', required=True, help='output directory.')

    optional_group = parser_view.add_argument_group('Optional arguments')
    optional_group.add_argument('-h', '--help', action='help',
                                help='Show this help message and exit')
    optional_group.add_argument('-f', '--image_format', default='svg', choices=['svg', 'png'],
                                help='default=svg,output image format.')
    optional_group.add_argument('-n', '--output_name', default='',
                                help='output image name. By default, the input file name is used.')
    optional_group.add_argument('-W', '--width', type=int, default=2200,
                                help='default=2200,output image width.')
    optional_group.add_argument('-H', '--height', type=int, default=1600,
                                help='default=1600,output image height.')
    optional_group.add_argument('--node_label', default='name-depth',
                                choices=['none', 'name', 'depth', 'length', 'name-depth', 'name-length', 'name-depth-length'],
                                help='default=name-depth,node label display mode.')
    optional_group.add_argument('--edge_label', default='name', choices=['none', 'name', 'derived'],
                                help='default=name,edge label display mode.')
    optional_group.add_argument('--font_size', type=int,
                                help='set node and edge label font size.')
    optional_group.add_argument('--node_font_size', type=int,
                                help='set node label font size.')
    optional_group.add_argument('--edge_font_size', type=int,
                                help='set edge label font size.')
    optional_group.add_argument('--font_color',
                                help='set label color, such as black or #000000.')
    optional_group.add_argument('--seed', type=int,
                                help='random seed for stable layout and colors.')
    parser_view.set_defaults(func=view)


    args = parser.parse_args()
    if len(vars(args)) == 0:
        parser.print_help()
        sys.exit(1)
    os.makedirs(args.output_dir, exist_ok=True)
    logger = Logger(os.path.join(args.output_dir, "HiMT.log"))
    sys.stdout = logger
    sys.stderr = logger
    print("conda",version)
    print(f"({datetime.datetime.now()}) your running command: ", ' '.join(sys.argv))
    args.func(args)

if __name__ == "__main__":
    main()
