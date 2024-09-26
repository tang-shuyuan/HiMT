import sys
import argparse
import os
import datetime
from version import get_versions
from fliter import fliter
from assembly import assembly
from assess import assese

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
    "for any questions, please put them forward on https://github.com/tang-shuyuan/EeayMT or https://bioanno.com.",\
    usage="himt [function] [argument]",formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-v', '--version', action='version', version='HiMT ' + version)
    function_description="""
    filter        Flitering the low-depth nuclear genome reads
    assembly      Assembling mitochondrail genome with HiFi sequcncing data 
    assese        Assessing the assembly quality of the mitochondrial genome """

    subparsers = parser.add_subparsers(title='function', description=function_description, metavar='')

    parser_fliter = subparsers.add_parser('fliter',description="Flitering the low-depth nuclear genome reads",\
                                          usage="himt fliter [argument]")
    optional_group = parser_fliter._action_groups.pop()
    required_group = parser_fliter.add_argument_group('Required arguments')

    # optional_group = parser_fliter.add_argument_group('Optional arguments')

    required_group.add_argument('-i', '--input_file', required=True, \
                                help='input a fasta or fastq file,gz compressed files are supported')
    required_group.add_argument('-o', '--output_dir', required=True, \
                                help='output directory')
    optional_group.add_argument('-n', '--head_number', type=int, default=4, \
                                help='default=4,The number of kmer species randomly selected')
    optional_group.add_argument('-t', '--thread', type=int, default=2, \
                                help='default=2 The number of thread used during code execution.')
    optional_group.add_argument('-b', '--base_number', type=int, default=3, choices=[3, 4], \
                                help='default=3,only can be 3 and 4,The number of bases at the beginning of kmer')
    optional_group.add_argument('-fd', '--fliter_depth', type=int, default=0, \
                                help='read depths below this value will be filtered')
    optional_group.add_argument('-fp', '--fliter_percentage', type=float, default=0.3, \
                                help='default=0.3,The depth of the mitochondrial genome obtained by blast,\
                                the proportion adjusted downwards on this value.')
    optional_group.add_argument('-p', '--proportion', type=float, default=1,
                                help='default=1,The percentage of the selected dataset from the entire file,\
                                choose a value from 0-1')
    optional_group.add_argument('-c', '--accuracy', type=float, default=0.8,
                                help='default=0.8,If one read has a high-frequency kmer ratio exceeding this value, \
                                it will be considered as a high-frequency read,choose a value from 0-1')
    optional_group.add_argument("-s", "--species", default="plant", choices=["plant", "animal"], \
                                help="default=plant,Species category,only can be plant or animal")
    optional_group.add_argument('-x','--normalize_depth',type=int,default=0,\
                                help="Normalize the mitochondrial genome depth to a value.\
                                Default returns the maximum mitochondrial genome depth,If the input value \
                                exceeds the mitochondrial genome depth, retain the maximum mitochondrial genome depth")
    parser_fliter._action_groups.append(optional_group)

    parser_fliter.set_defaults(func=fliter)



    parser_assembly = subparsers.add_parser('assembly', description="Assembling mitochondrail genome with HiFi\
     sequcncing data ",usage="himt assembly [argument]")
    optional_group = parser_assembly._action_groups.pop()
    required_group = parser_assembly.add_argument_group('Required arguments')

    required_group.add_argument('-i', '--input_file', required=True, \
                                help='input a fasta or fastq file,gz compressed files are supported')
    required_group.add_argument('-o', '--output_dir', required=True, \
                                help='output directory')

    optional_group.add_argument('-n', '--head_number', type=int, default=4, \
                                help='default=4,The number of kmer species randomly selected')
    optional_group.add_argument('-t', '--thread', type=int, default=2, \
                                help='default=2 The number of thread used during code execution.')
    optional_group.add_argument('-b', '--base_number', type=int, default=3, choices=[3, 4], \
                                help='default=3,only can be 3 and 4,The number of bases at the beginning of kmer')
    optional_group.add_argument('-fd', '--fliter_depth', type=int, default=0, \
                                help='read depths below this value will be filtered')
    optional_group.add_argument('-fp', '--fliter_percentage', type=float, default=0.3, \
                                help='default=0.3,The depth of the mitochondrial genome obtained by blast, \
                                the proportion adjusted downwards on this value.')
    optional_group.add_argument('-p', '--proportion', type=float, default=1, \
                                help='default=1,The percentage of the selected dataset from the entire file,\
                                choose a value from 0-1')
    optional_group.add_argument('-c', '--accuracy', type=float, default=0.8,
                                help='default=0.8,If one read has a high-frequency kmer ratio exceeding this value, \
                                it will be considered as a high-frequency read,choose a value from 0-1')
    optional_group.add_argument("-s", "--species", default="plant", choices=["plant", "animal"], \
                                help="default=plant,Species category,only can be plant or animal")
    optional_group.add_argument('--no_flye_meta', action='store_true', help="By default, \
                                we use flye Meta pattern to assemble the mitochondrial genome. \
                                If you don't want to use meta pattern, add this parameter")
    optional_group.add_argument('-x',"--normalize_depth" ,type=int,default=0, help="Normalize the \
    mitochondrial genome depth to a value.If the input value exceeds the mitochondrial genome depth, \
    retain the maximum mitochondrial genome depth,the default mitochondrial genome depth ranges between 15 and 50,\
    input a value less than 0 (such as:-1) to retain the maximum mitogenome depth.")
    parser_assembly._action_groups.append(optional_group)
    parser_assembly.set_defaults(func=assembly)

    parser_assese = subparsers.add_parser("assese",description="Assessing the assembly quality\
     of the mitochondrial genome",usage="himt assese [argument]")
    required_group = parser_assese.add_argument_group('Required arguments')

    required_group.add_argument('-i', '--input_file', required=True, \
                                help='input a fasta or gfa file')
    required_group.add_argument('-o', '--output_dir', required=True, \
                                help='output directory')
    parser_assese.set_defaults(func=assese)

    args = parser.parse_args()
    if len(vars(args)) == 0:
        parser.print_help()
        sys.exit(1)
    os.makedirs(args.output_dir, exist_ok=True)
    logger = Logger(os.path.join(args.output_dir, "easymt.log"))
    sys.stdout = logger
    sys.stderr = logger
    print(f"({datetime.datetime.now()}) your running command: ", ' '.join(sys.argv))
    args.func(args)

if __name__ == "__main__":
    main()