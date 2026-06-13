#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys
import datetime
import subprocess
from filter import filter
from assess import assess

def simple_gfa(args):
    config_path = os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), "config_database")
    input_file=os.path.join(args.output_dir,"flye_output/assembly_graph.gfa")
    os.makedirs(os.path.join(args.output_dir, "blast_output"), exist_ok=True)
    with open(input_file, "r") as f_in, \
            open(os.path.join(args.output_dir,"blast_output","out.fa"), "w") as f_w:
        for line in f_in:
            if line.startswith("S"):
                line_content = line.strip().split("\t")
                f_w.write(">" + line_content[1] + "\n")
                f_w.write(line_content[2] + "\n")
    fa_file = os.path.join(args.output_dir,"blast_output","out.fa")

    if args.species == "animal":
        prot_sequence = os.path.join(config_path, "Drosophila_gunungcola_CM045947.fasta")
    else:
        chr_pro_sequence = os.path.join(config_path,"Arabidopsis_chr_protein.fasta")
        prot_sequence = os.path.join(config_path, "Arabidopsis_protein.fasta")
    out_blast_db = os.path.join(args.output_dir,"blast_output","temp_db")
    blast_result = os.path.join(args.output_dir, "blast_output","temp_blast_result")

    command4 = ["makeblastdb", "-in", fa_file, "-dbtype", "nucl", "-out", out_blast_db]
    result4 = subprocess.run(command4, stdout=subprocess.PIPE,
    stderr=subprocess.STDOUT, text=True, check=True)
    if result4.stdout:
        print(result4.stdout, end="")


    command5 = ["tblastn", "-db", out_blast_db, "-query", prot_sequence, 
                "-evalue", "1e-10", "-out", blast_result, "-outfmt", "6"]
    subprocess.run(command5, stdout=subprocess.PIPE,
    stderr=subprocess.STDOUT, text=True, check=True)


    with open(blast_result, "r") as f:
        mito_id = set()
        mito_edge={}
        for line in f:
            mito_gene=line.strip().split("\t")[0]
            id_name = line.strip().split("\t")[1]
            mito_id.add(id_name)
            if id_name in mito_edge:
                mito_edge[id_name].add(mito_gene)
            else:
                mito_edge[id_name]={mito_gene}
    max_mito_edge=max(mito_edge,key=lambda x: len(mito_edge[x]))

    with open(input_file, "r") as f:
        id_new = mito_id.copy()
        while True:
            for line in f:
                if line.startswith("L"):
                    line_content = line.strip().split("\t")
                    if line_content[1] in id_new or line_content[3] in id_new:
                        id_new.add(line_content[1])
                        id_new.add(line_content[3])
            if id_new == mito_id:
                break
            else:
                mito_id = id_new.copy()
                f.seek(0)
        ###filter short nuclear sequence
        f.seek(0)
        link_seq=set()
        un_link_seq=set()
        for line in f:
            if line.startswith("L"):
                line_content = line.strip().split("\t")
                link_seq.add(line_content[1])
                link_seq.add(line_content[3])
            elif line.startswith("S"):
                if line.strip().split("\t")[1]==max_mito_edge:
                    max_mito_edge_depth=int(line.strip().split("\t")[3].split(":")[2])
        for id in mito_id:
            if id not in link_seq:
                un_link_seq.add(id)
        f.seek(0)
        for line in f:
            if line.startswith("S"):
                line_content = line.strip().split("\t")
                if line_content[1] in un_link_seq:
                    # seq_length=len(line_content[2])
                    seq_depth = int(line_content[3].split(":")[2])
                    if seq_depth<max_mito_edge_depth/2 and len(mito_id)>1:
                        mito_id.remove(line_content[1])


        ###split mitogenome and chloroplast genome
        if args.species=='plant' and args.func.__name__=="assemble":
            chr_blast_result = os.path.join(args.output_dir, "blast_output", "chr_blast_result")
            command6 = ["tblastn", "-db", out_blast_db, "-query", chr_pro_sequence, 
                        "-evalue", "1e-10", "-out", chr_blast_result, "-outfmt", "6"]
            subprocess.run(command6, stderr=subprocess.PIPE, text=True, check=True)

            chr_id = {}
            with open(chr_blast_result, 'r') as f_in:
                for line in f_in:
                    chr_gene = line.strip().split("\t")[0]
                    id_name = line.strip().split("\t")[1]
                    if id_name in chr_id:
                        chr_id[id_name].add(chr_gene)
                    else:
                        chr_id[id_name] = {chr_gene}
            max_edge = max(chr_id,  key=lambda x: len(chr_id[x]))
            chr_id = {max_edge}
            f.seek(0)
            id_new = chr_id.copy()
            while True:
                for line in f:
                    if line.startswith("L"):
                        line_content = line.strip().split("\t")
                        if line_content[1] in id_new or line_content[3] in id_new:
                            id_new.add(line_content[1])
                            id_new.add(line_content[3])
                if id_new == chr_id:
                    break
                else:
                    chr_id = id_new.copy()
                    f.seek(0)
            ##get chloroplast and mitochondrial genome depth
            chr_depth=[]
            mito_depth=[]
            f.seek(0)
            for line in f:
                if line.startswith("S"):
                    line_content = line.strip().split("\t")
                    contig_depth = int(line_content[3].split(":")[2])
                    if line_content[1] in chr_id:
                        chr_depth.append(contig_depth)
                    elif line_content[1] in mito_id-chr_id:
                        mito_depth.append(contig_depth)
            if 0 < len(chr_id) <= 3 and len(chr_depth) > 0 and len(mito_depth)>0:
                if sum(chr_depth)/len(chr_depth) > 3*sum(mito_depth)/len(mito_depth):
                    mito_id=mito_id-chr_id
                    ####write chloroplast genome
                    with open(os.path.join(args.output_dir, "himt_chloroplast.gfa"), "w") as f_w:
                        f.seek(0)
                        for line in f:
                            if line.startswith("H"):
                                f_w.write(line)
                            elif line.startswith("S") or line.startswith("L") or line.startswith("P"):
                                if line.split("\t")[1] in chr_id:
                                    f_w.write(line)

    if args.func.__name__=="refassemble":
        output_gfa = os.path.join(args.output_dir,"himt_refassemble.gfa")
        output_fa = os.path.join(args.output_dir,"himt_refassemble_raw.fa")
    else:
        output_gfa = os.path.join(args.output_dir, "himt_mitochondrial.gfa")
        output_fa = os.path.join(args.output_dir, "himt_mitochondrial_raw.fa")

    with open(input_file, "r") as f,open(output_gfa, "w") as f_w,\
        open(output_fa, "w") as fa_w:
        for line in f:
            if line.startswith("S") :
                if line.split("\t")[1] in mito_id:
                    f_w.write(line)
                    fa_w.write(">"+line.split("\t")[1]+"\n")
                    fa_w.write(line.split("\t")[2]+"\n")
            elif line.startswith("L") or line.startswith("P"):
                if line.split("\t")[1] in mito_id:
                    f_w.write(line)
            else:
                f_w.write(line)

def get_reversed_complement(sequence):
    new_seq=''
    complement_dic={"A":"T","T":"A","G":"C","C":"G"}
    for item in reversed(sequence.upper()):
        try:
            new_seq+=complement_dic[item]
        except KeyError:
            new_seq+=item
    return new_seq
def output_two_chl_hap(chl_gfa_file, output):
    seq_info = {}
    seq_number = {}
    with open(chl_gfa_file, 'r') as f:
        for line in f:
            line_content = line.split("\t")
            if line.startswith("S"):
                seq_info[line_content[1]] = line_content[2]
            elif line.startswith("L"):
                if line_content[1] in seq_number:
                    seq_number[line_content[1]] += 1
                else:
                    seq_number[line_content[1]] = 1
                if line_content[3] in seq_number:
                    seq_number[line_content[3]] += 1
                else:
                    seq_number[line_content[3]] = 1

        for key, value in seq_number.items():
            if value == 4:
                reversed_repeat = key

        if len(seq_info) == 1:
            with open(os.path.join(output, "chloroplast_path1.fa"), 'w') as f_w1:
                f_w1.write(">path1"+"\n")
                for value in seq_info.values():
                    f_w1.write(value)
        elif len(seq_info) == 3:
            f.seek(0)
            for line in f:
                if line.startswith("L"):
                    line_content = line.split("\t")
                    if line_content[2] == "+" and line_content[4] == "+" and reversed_repeat in line_content:
                        if line_content[1] == reversed_repeat:
                            single_copy = line_content[3]
                        else:
                            single_copy = line_content[1]
                        for key in seq_number.keys():
                            if key not in {single_copy, reversed_repeat}:
                                another_single_copy = key
                        with open(os.path.join(output, "chloroplast_path1.fa"), 'w') as f_w1:
                            f_w1.write(">path1"+"\n"+
                                seq_info[single_copy] + seq_info[reversed_repeat] + seq_info[another_single_copy] + \
                                get_reversed_complement(seq_info[reversed_repeat]))
                        with open(os.path.join(output, "chloroplast_path2.fa"), 'w') as f_w2:
                            f_w2.write(">path2"+"\n"+seq_info[single_copy] + seq_info[reversed_repeat] + \
                                       get_reversed_complement(seq_info[another_single_copy]) + get_reversed_complement(
                                seq_info[reversed_repeat]))
                        break

        else:
            print("chloroplast genome seems wrong")

def get_extract_file_info(extract_file):
    with open(extract_file, 'r') as f_in:
        seq_number = 0
        min_seq = float('inf')
        max_seq = 0
        sum_len = 0
        for line in f_in:
            if line.startswith(">"):
                seq_number += 1
            else:
                max_seq = max(max_seq, len(line.strip()))
                min_seq = min(min_seq, len(line.strip()))
                sum_len += len(line.strip())
    return seq_number, sum_len, min_seq, max_seq
def get_gfa_file_info(gfa):
    with open(gfa,'r') as f:
        weighting_total=0
        total_len=0
        for line in f:
            if line.startswith("S"):
                weighting_total += len(line.split("\t")[2])*int(line.strip().split("\t")[3].split(":")[2])
                total_len += len(line.split("\t")[2])
        depth = int(weighting_total/total_len)
    return depth
def file_exist(file_path):
    if os.path.exists(file_path) and os.path.getsize(file_path) > 0:
        return True
    else:
        return False
def get_hard_disk(path):
    total_size = 0
    for dirpath, dirnames, filenames in os.walk(path):
        for f in filenames:
            fp = os.path.join(dirpath, f)
            total_size += os.path.getsize(fp)
    print(f"({datetime.datetime.now()}) The size of the output folder is: {total_size/(1024**3)} GB")

def assemble(args):
    total_seq_number,total_seq_len,total_min,total_max,reduction_radio=filter(args)

    read_type_param={"HiFi":"--pacbio-hifi" , "ONT": "--nano-corr" , "CLR": "--pacbio-corr"}

    extract_fa = os.path.join(args.output_dir, 'extract.fa')
    flye_output = os.path.join(args.output_dir, 'flye_output')
    
    if args.no_flye_meta:
        flye_cmd = ["flye", read_type_param[args.data_type], extract_fa, 
                   "-o", flye_output, "-t", str(args.thread)]
    else:
        flye_cmd = ["flye", read_type_param[args.data_type], extract_fa, 
                   "-o", flye_output, "-t", str(args.thread), "--meta"]

    p = subprocess.Popen(
        flye_cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1  # 行缓冲（对文本模式 + PIPE 生效）
    )

    # 实时读取输出
    for line in p.stdout:
        print(line, end="")

    p.wait()
    if file_exist(os.path.join(args.output_dir,"flye_output/assembly_graph.gfa")):
        simple_gfa(args)
    else:
        print("The flye output file doesn't exist")
        sys.exit(1)

    filter_seq_number,filter_seq_len,filter_min,filter_max,=\
    get_extract_file_info(os.path.join(args.output_dir, 'extract.fa'))

    filter_mito_detph=get_gfa_file_info(os.path.join(args.output_dir, "himt_mitochondrial.gfa"))


    if file_exist(os.path.join(args.output_dir, "himt_chloroplast.gfa")):
        output_two_chl_hap(os.path.join(args.output_dir, "himt_chloroplast.gfa"), args.output_dir)
        filter_chl_depth=get_gfa_file_info(os.path.join(args.output_dir,"himt_chloroplast.gfa"))

    if args.species =="plant":

        args.table3_value = [total_seq_number, total_seq_len, total_min, int(round(total_seq_len// total_seq_number,0)),\
        total_max,int(filter_mito_detph/reduction_radio),filter_seq_number, filter_seq_len,\
        filter_min,int(round(filter_seq_len // filter_seq_number,0)),filter_max,filter_mito_detph]

        args.input_file=os.path.join(args.output_dir, "himt_mitochondrial.gfa")
        args.category='mitochondrial'
        assess(args)

        if file_exist(os.path.join(args.output_dir, "himt_chloroplast.gfa")):
            args.table3_value = args.table3_value = [total_seq_number, total_seq_len, total_min,
            int(round(total_seq_len// total_seq_number,0)),
            total_max,int(filter_chl_depth/reduction_radio),filter_seq_number, filter_seq_len,
            filter_min,int(round(filter_seq_len // filter_seq_number,0)),filter_max,filter_chl_depth]


            args.input_file = os.path.join(args.output_dir, "himt_chloroplast.gfa")
            args.category='chloroplast'
            assess(args)

    elif args.species=="animal":
        args.table3_value = [total_seq_number, total_seq_len, total_min,
                             int(round(total_seq_len // total_seq_number, 0)), \
                             total_max, int(filter_mito_detph / reduction_radio), filter_seq_number, filter_seq_len, \
                             filter_min, int(round(filter_seq_len // filter_seq_number, 0)), filter_max,
                             filter_mito_detph]

        args.input_file = os.path.join(args.output_dir, "himt_mitochondrial.gfa")
        args.category="mitochondrial"
        assess(args)


    get_hard_disk(args.output_dir)