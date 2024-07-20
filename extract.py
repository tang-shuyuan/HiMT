#!/usr/bin/env python
# -*- coding: utf-8 -*-

import gzip
import os
import datetime
import random
from multiprocessing import Pool
import sys
import argparse
from version import get_versions


def process_input_file(infile):
    try:
        with gzip.open(infile, 'rb') as f_in:
            first_item = f_in.read(1).decode("utf-8")
            print(f"({datetime.datetime.now()}) uncompressing file")
            if first_item == ">":
                with open(os.path.join(output_dir, "process.fa"), "wb") as f_out:
                    f_out.write(first_item.encode("utf-8"))
                    f_out.writelines(line for line in f_in)
                return os.path.join(output_dir, "process.fa")
            elif first_item == '@':
                print(f"({datetime.datetime.now()}) transfoming fq to fa")
                with open(os.path.join(output_dir, "process.fa"), 'wb') as f_out:
                    f_in.seek(0)
                    for i, line in enumerate(f_in):
                        if i % 4 == 0:
                            f_out.write(">".encode("utf-8") + line[1:])
                        elif i % 4 == 1:
                            f_out.write(line)
                return os.path.join(output_dir, "process.fa")
            else:
                print("file_type is unknow")
                sys.exit()
    except OSError:
        with open(infile, "r") as f:
            first_item = f.read(1)
            if first_item == ">":
                return infile
            elif first_item == '@':
                print(f"({datetime.datetime.now()}) transfoming fq to fa")
                with open(os.path.join(output_dir, "process.fa"), "w") as f_out:
                    f.seek(0)
                    for i, line in enumerate(f):
                        if i % 4 == 0:
                            f_out.write(">"+line[1:])
                        elif i % 4 == 1:
                            f_out.write(line)
                return os.path.join(output_dir, "process.fa")
            else:
                print("file_type unknown")
                sys.exit()

##数据是一行id多行序列的情况要矫正一下
def stardan_format(infile):
    with open(infile, "r") as f:
        for i, line in enumerate(f):
            if i in [0, 2, 4, 6, 8, 10]:
                first_item = line[0]
                if first_item != ">":
                    unstardan_format = True
                    break
                    print("modifying format")
                else:
                    unstardan_format = False
            elif i == 11:
                break
    if unstardan_format:
        with open(infile, "r") as f, open(os.path.join(output_dir, "stardand_process.fa"), "w") as f_w:
            for i, line in enumerate(f):
                if line.startswith(">"):
                    if i == 0:
                        f_w.write(line)
                    else:
                        f_w.write("\n" + line)
                else:
                    f_w.write(line.strip())
        return os.path.join(output_dir, "stardand_process.fa")
    else:
        return infile


def get_proportion_file(infile,outfile,proportion):
    with open(infile, "r") as f:
        line_number = sum(1 for line in f)
    base_line = list(range(1, line_number, 2))
    select_line = set(random.sample(base_line, int(len(base_line) * proportion)))
    with open(infile,"r") as f,open(outfile,"w") as f_w:
        for i,line in enumerate(f):
            if i + 1 in select_line or i in select_line:
                f_w.write(line)
    return os.path.join(output_dir,"proportion.fa")


def generate_header(n):
    Base = ["A", "T", "G", "C"]
    header=[]
    if n==1:
        return Base
    else:
        for i in Base:
            for j in generate_header(n-1):
                header.append(i + j)
        return header

###调用blast确定线粒体覆盖深度，然后确定提取深度，怎么定这个参数？
def determining_motigenome_depth(assemble_species,process_file):
    if assemble_species=="animal":
        prot_sequence=os.path.join(config_path,"Drosophila_gunungcola_CM045947.fasta")
    else:
        prot_sequence = os.path.join(config_path,"Arabidopsis_protein.fasta")
    # prot_sequence = os.path.join(config_path,"sequence.txt")
    # prot_sequence = "/mnt/g/test2/human.txt"
    out_blast_db = os.path.join(output_dir, "database")
    blast_result = os.path.join(output_dir, "blast_result")
    command1 = f"makeblastdb -in {process_file} -dbtype nucl -out {out_blast_db}"
    os.system(command1)
    command2 = f"tblastn -num_descriptions 100000 -num_alignments 100 -num_threads {manual_input_process_number}\
    -db {out_blast_db} -query {prot_sequence} -evalue 1e-10 -out {blast_result}"
    os.system(command2)
    mitochondrial_depth = {}
    with open(blast_result, "r") as f:
        for i, line in enumerate(f):
            if line.strip().startswith("Query="):
                line_number = i
                current_gene = line.strip()[7:]
                mitochondrial_depth[current_gene] = 0
                first = True
            elif line.strip().startswith(">") and first == True:
                blast_reads = i - line_number - 8
                mitochondrial_depth[current_gene] += blast_reads
                first = False
            else:
                continue
    print(f"({datetime.datetime.now()}) The minimum depth is{min(mitochondrial_depth.values())}")
    reads_depth = sorted(list(mitochondrial_depth.values()))
    print(reads_depth)
    for i, depth in enumerate(reads_depth):
        if depth != 0:
            min_reads_depth = depth
            break
    try:
        next_depth = reads_depth[i + 1]
    except IndexError:
        print("Check that you have entered the correct species category")
        sys.exit()
    if next_depth > 2 * min_reads_depth :
        min_reads_depth = next_depth
    print(f"({datetime.datetime.now()}) The depth of actual application is {min_reads_depth}")
    return min_reads_depth

##找某种固定header的kmer，不用将数据写入外存，直接得到高频kmer。
def process_header(head,threthold,file):
    dic = {}
    with open(file, 'r') as f_in:
        for i, line in enumerate(f_in):
            if i % 2 == 1:
                start = 0
                while True:
                    start = line.strip().find(head, start)
                    if start == -1 or start + kmer_length > len(line.strip()):
                        break
                    kmer = line.strip()[start:start + 21]
                    start += 1
                    if kmer not in dic:
                        dic[kmer] = 1
                    else:
                        dic[kmer] += 1
    high_frequency_kmer=set()
    for key,value in dic.items():
        if value>= threthold:
            high_frequency_kmer.add(key)
    print(f"{head} consume memory {sys.getsizeof(dic) / (1024.0 ** 3)} GB ")
    del dic
    return high_frequency_kmer

###提取reads
def extract_reads(process_file,high_frequency_kmer):
    if min_depth <=50:
        n=1
    elif min_depth<100:
        n = 15/min_depth
    elif min_depth <=500:
        n=0.1
    else:
        n=30/min_depth
    with open(process_file, "r") as f:
        line_number = sum(1 for line in f)
    base_line = list(range(1, line_number, 2))
    select_line = set(random.sample(base_line, int(len(base_line) * n)))
    with open(process_file, 'r') as f_in,open(os.path.join(output_dir,"extract.fa"),"w") as w:
        for i, line in enumerate(f_in):
            if i + 1 in select_line or i in select_line:
                if i % 2 ==0:
                    id_line=line
                elif i % 2 == 1:
                    if len(line.strip()) > kmer_length:
                        kmers = set()
                        # for j in range(len(line.strip()) - kmer_length + 1):
                        #     kmer = line.strip()[j:j + kmer_length]
                        #     if kmer[0:3] in header:
                        #         kmers.add(kmer)
                        for head in header:
                            start = 0
                            while True:
                                start = line.strip().find(head, start)
                                if start == -1 or start + kmer_length > len(line.strip()):
                                    break
                                kmer = line.strip()[start:start + 21]
                                start += 1
                                kmers.add(kmer)
                        if len(kmers) > 0 and len(kmers & high_frequency_kmer) / len(kmers) >= 0.8:
                            w.write(id_line)
                            w.write(line)
def simple_gfa(assemble_species,input_file):
    with open(input_file, "r") as f_in, \
            open(os.path.join(output_dir, "out.fa"), "w") as f_w:
        for line in f_in:
            if line.startswith("S"):
                line_content = line.strip().split("\t")
                f_w.write(">" + line_content[1] + "\n")
                f_w.write(line_content[2] + "\n")
    fa_file = os.path.join(output_dir, "out.fa")
    if assemble_species == "animal":
        prot_sequence = os.path.join(config_path, "Drosophila_gunungcola_CM045947.fasta")
    else:
        prot_sequence = os.path.join(config_path, "Arabidopsis_protein.fasta")
    out_blast_db = os.path.join(output_dir, "temp_db")
    blast_result = os.path.join(output_dir, "temp_blast_result")

    command4 = f"makeblastdb -in {fa_file} -dbtype nucl -out {out_blast_db}"
    os.system(command4)
    command5 = f"tblastn -num_descriptions 100000 -num_alignments 100 -num_threads {manual_input_process_number} \
    -db {out_blast_db} -query {prot_sequence} -evalue 1e-10 -out {blast_result}"
    os.system(command5)
    with open(blast_result, "r") as f:
        id = set()
        need = False
        for i, line in enumerate(f):
            if line.startswith("Sequences producing"):
                need = True
            elif line.startswith(">"):
                need = False
            else:
                if need == True and line.strip():
                    line_content = line.strip().split("\x20")
                    id.add(line_content[0])
    with open(input_file, "r") as f, open(os.path.join(output_dir, "extract.gfa"), "w") as f_w:
        id_new = id
        while True:
            for line in f:
                if line.startswith("L"):
                    line_content = line.strip().split("\t")
                    if line_content[1] in id_new or line_content[3] in id_new:
                        id_new.add(line_content[1])
                        id_new.add(line_content[3])
            if id_new == id:
                break
            else:
                id = id_new
                f.seek(0)
        f.seek(0)
        for line in f:
            if line.startswith("H"):
                f_w.write(line)
            elif line.startswith("S") or line.startswith("L") or line.startswith("P"):
                if line.split("\t")[1] in id:
                    f_w.write(line)

def get_hard_disk(path):
    total_size = 0
    for dirpath, dirnames, filenames in os.walk(path):
        for f in filenames:
            fp = os.path.join(dirpath, f)
            total_size += os.path.getsize(fp)
    print(f"The size of the folder is: {total_size/(1024**3)} GB")

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

class MyParser(argparse.ArgumentParser):
    def error(self, message):
        self.print_help()
        self.exit(2, '%s: error: %s\n' % (self.prog, message))

def main():
    version=get_versions()
    github = "If you have any questions, please put them forward on https://github.com/tang-shuyuan/EeayMT or https://bioanno.com."
    parser = MyParser(description="this is a tools for extracting high depth reads and \
calling flye to assemble mitochondrial genome"+"\n"+version+"\n"+github,formatter_class=argparse.RawTextHelpFormatter)
#     parser = argparse.ArgumentParser(description= "this is a tools for extracting high depth reads and \
# calling flye to assemble mitochondrial genome"+"\n"+version+"\n"+github,formatter_class=argparse.RawTextHelpFormatter)

    required_group = parser.add_argument_group('Required arguments')
    required_group.add_argument('-i', '--input_file', required=True, \
                                help='input a fasta or fastq file,gz compressed files are supported')
    required_group.add_argument('-o', '--output_dir', required=True, \
                                help='output directory')
    optional_group = parser.add_argument_group('Options arguments')
    optional_group.add_argument('-n', '--head_number', type=int, default=4, \
                                help='default=4,The number of kmer species randomly selected')
    optional_group.add_argument('-t', '--thread', type=int, default=2, \
                                help='default=2 The number of thread used during code execution.')
    optional_group.add_argument('-b', '--base_number', type=int, default=3,choices=[3,4],\
                                help='default=3,only can be 3 and 4,The number of bases at the beginning of kmer')
    optional_group.add_argument('-fd', '--fliter_depth', type=int, default=0, \
                                help='read depths below this value will be filtered')
    optional_group.add_argument('-fp', '--fliter_percentage', type=float, default=0.3, \
                                help='default=0.3,The depth of the mitochondrial genome obtained by blast, the proportion adjusted downwards on this value.')
    optional_group.add_argument('-p', '--proportion', type=float, default=1, \
                                help='default=1,The percentage of the selected dataset from the entire file')
    optional_group.add_argument("-s","--species",default="plant",choices=["plant","animal"],\
                                help="default=plant,Species category,only can be plant or animal")
    parser.add_argument('-v', '--version', action='version', version='EasyMT ' + version)
    args = parser.parse_args()
    # parser._action_groups.reverse()
    global kmer_length,config_path,output_dir,manual_input_process_number,min_depth,header
    kmer_length = 21
    input_file = args.input_file
    config_path=os.path.dirname(os.path.abspath(sys.argv[0]))
    output_dir = args.output_dir
    title_base_number = args.base_number
    manual_input_head_number = args.head_number
    manual_input_process_number = args.thread
    manaul_input_fliter_threshold = args.fliter_depth
    manaul_input_fliter_percentage= args.fliter_percentage
    manual_input_file_proportion =args.proportion
    assmble_species=args.species
    os.makedirs(output_dir, exist_ok=True)
    logger = Logger(os.path.join(output_dir,"easymt.log"))
    sys.stdout = logger
    sys.stderr = logger
    print(f"({datetime.datetime.now()}) your running command: ",' '.join(sys.argv))
    if title_base_number==3:
        if manual_input_head_number<4 or manual_input_head_number>64:
            print("if base number is 3,head number only can be between 4 to 64 ")
            sys.exit()
    elif title_base_number==4:
        if manual_input_head_number==4:
            manual_input_head_number = 16
        if manual_input_head_number<16 or manual_input_head_number>256:
            print("if base number is 4,head number only can be be tween 16 to 256")
            sys.exit()
    else:
        print("title base number only can be 3 and 4")


    print(f"({datetime.datetime.now()}) processing file")
    if title_base_number==3 and manual_input_head_number==4:
        header={x+y for x,y in zip(generate_header(1),random.sample(generate_header(2),4))}
    else:
        header =set(random.sample(generate_header(title_base_number), manual_input_head_number))
    print(f"({datetime.datetime.now()}) random head {header}")
    #判断文件类型，改解压或者转fa格式
    process_file = process_input_file(input_file)
    process_file = stardan_format(process_file)
    #只取一部分数据进行提取，会加快速度
    if manual_input_file_proportion !=1:
        process_file=get_proportion_file(process_file,os.path.join(output_dir,"proportion.fa"),manual_input_file_proportion)
    if manaul_input_fliter_threshold:
        min_depth=manaul_input_fliter_threshold
    else:
        min_depth = determining_motigenome_depth(assmble_species,process_file)
    high_depth_kmer = set()
    with Pool(processes=manual_input_process_number) as pool:
        results = pool.starmap(process_header, [(h,min_depth*manaul_input_fliter_percentage,process_file) for h in header])
    for result in results:
        high_depth_kmer.update(result)
    print(f"({datetime.datetime.now()}) high depth kmer cosume memory {sys.getsizeof(high_depth_kmer)/(1024.0 ** 3)}")
    extract_reads(process_file,high_depth_kmer)
    print(f"({datetime.datetime.now()}) complete extracting")

    command3=f"flye --pacbio-hifi {os.path.join(output_dir,'extract.fa')} \
        -o {os.path.join(output_dir,'flye_output')} -t {manual_input_process_number} --meta"
    os.system(command3)
    simple_gfa(assmble_species,os.path.join(output_dir,"flye_output/assembly_graph.gfa"))
    get_hard_disk(output_dir)
    print(datetime.datetime.now())
if __name__=="__main__":
    main()