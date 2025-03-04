#!/usr/bin/env python
# -*- coding: utf-8 -*-
import gzip
import os
import datetime
import random
from multiprocessing import Pool
import sys
def process_input_file(input_file,output_dir):
    try:
        with gzip.open(input_file,'rb') as f_in:
            first_item = f_in.read(1).decode('utf-8')
            seq_number = 0
            min_seq = float('inf')
            max_seq = 0
            sum_len = 0
            if first_item == ">":
                for i, line in enumerate(f_in):
                    if i in {2, 4, 6, 8, 10}:
                        item = line.decode("utf-8")[0]
                        if item != ">":
                            unstardan_format = True
                            print(f"({datetime.datetime.now()}) Modifying the file format\
                                     to have one line for each ID followed by one line for its corresponding sequence.")
                            break
                        else:
                            unstardan_format = False
                    elif i == 11:
                        break
                if unstardan_format:
                    f_in.seek(0)
                    with open(os.path.join(output_dir, "process.fa"), "w") as f_w:
                        f_w.write(f_in.readline().decode("utf-8"))
                        seq_number += 1
                        for line in f_in:
                            line = line.decode("utf-8")
                            if line.startswith(">"):
                                f_w.write("\n" + line)
                                seq_number += 1
                            else:
                                f_w.write(line.strip())
                                max_seq = max(max_seq, len(line.strip()))
                                min_seq = min(min_seq, len(line.strip()))
                                sum_len += len(line.strip())
                else:
                    f_in.seek(0)
                    with open(os.path.join(output_dir, "process.fa"), "w") as f_out:
                        for line in f_in:
                            line = line.decode("utf-8")
                            if line.startswith(">"):
                                f_out.write(line)
                                seq_number += 1
                            else:
                                f_out.write(line)
                                max_seq = max(max_seq, len(line.strip()))
                                min_seq = min(min_seq, len(line.strip()))
                                sum_len += len(line.strip())
            elif first_item == '@':
                print(f"({datetime.datetime.now()}) transfoming fq to fa")
                with open(os.path.join(output_dir, "process.fa"), 'w') as f_out:
                    f_in.seek(0)
                    for i, line in enumerate(f_in):
                        line = line.decode("utf-8")
                        if i % 4 == 0:
                            f_out.write(">" + line[1:])
                            seq_number += 1
                        elif i % 4 == 1:
                            f_out.write(line)
                            max_seq = max(max_seq, len(line.strip()))
                            min_seq = min(min_seq, len(line.strip()))
                            sum_len += len(line.strip())
            else:
                print("file type unknow")
                sys.exit()
            return seq_number, sum_len, min_seq, max_seq, os.path.join(output_dir, "process.fa")
    except OSError:
        with open(input_file, 'r') as f_in:
            first_item = f_in.read(1)
            seq_number = 0
            min_seq = float('inf')
            max_seq = 0
            sum_len = 0
            if first_item == ">":
                for i, line in enumerate(f_in):
                    if i in {2, 4, 6, 8, 10}:
                        item = line[0]
                        if item != ">":
                            unstardan_format = True
                            print(f"({datetime.datetime.now()}) Modifying the file format\
                                     to have one line for each ID followed by one line for its corresponding sequence.")
                            break
                        else:
                            unstardan_format = False
                    elif i == 11:
                        break
                if unstardan_format:
                    f_in.seek(0)
                    with open(os.path.join(output_dir, "process.fa"), "w") as f_w:
                        f_w.write(f_in.readline())
                        seq_number += 1
                        for line in f_in:
                            if line.startswith(">"):
                                f_w.write("\n" + line)
                                seq_number += 1
                            else:
                                f_w.write(line.strip())
                                max_seq = max(max_seq, len(line.strip()))
                                min_seq = min(min_seq, len(line.strip()))
                                sum_len += len(line.strip())
                else:
                    f_in.seek(0)
                    for line in f_in:
                        if line.startswith(">"):
                            seq_number += 1
                        else:
                            max_seq = max(max_seq, len(line.strip()))
                            min_seq = min(min_seq, len(line.strip()))
                            sum_len += len(line.strip())
                    return seq_number, sum_len, min_seq, max_seq, input_file
            elif first_item == '@':
                print(f"({datetime.datetime.now()}) transfoming fq to fa")
                with open(os.path.join(output_dir, "process.fa"), 'w') as f_out:
                    f_in.seek(0)
                    for i, line in enumerate(f_in):
                        if i % 4 == 0:
                            f_out.write(">" + line[1:])
                            seq_number += 1
                        elif i % 4 == 1:
                            f_out.write(line)
                            max_seq = max(max_seq, len(line.strip()))
                            min_seq = min(min_seq, len(line.strip()))
                            sum_len += len(line.strip())
            else:
                print("file type unknow")
                sys.exit()
            return seq_number, sum_len, min_seq, max_seq, os.path.join(output_dir, "process.fa")

def generate_header(n):
    Base = ["A", "T", "G", "C"]
    header = []
    if n == 1:
        return Base
    else:
        for i in Base:
            for j in generate_header(n - 1):
                header.append(i + j)
        return header

def determining_mitogenome_depth(assemble_species,config_path,thread,process_file,output_dir):
    if assemble_species == "animal":
        prot_sequence = os.path.join(config_path, "Drosophila_gunungcola_CM045947.fasta")
    else:
        prot_sequence = os.path.join(config_path, "Arabidopsis_protein.fasta")
    out_blast_db = os.path.join(output_dir, "database")
    blast_result = os.path.join(output_dir, "blast_result")
    command1 = f"makeblastdb -in {process_file} -dbtype nucl -out {out_blast_db}"
    os.system(command1)
    command2 = f"tblastn -num_descriptions 100000 -num_alignments 100 -num_threads {thread}\
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
    print(f"({datetime.datetime.now()}) The minimum depth is {min(mitochondrial_depth.values())}")
    reads_depth = sorted(list(mitochondrial_depth.values()))
    print(reads_depth)
    print(mitochondrial_depth)
    for i, depth in enumerate(reads_depth):
        if depth != 0:
            min_reads_depth = depth
            break
    try:
        next_depth = reads_depth[i + 1]
        if next_depth > 2 * min_reads_depth:
            min_reads_depth = next_depth
        print(f"({datetime.datetime.now()}) The depth of actual application is {min_reads_depth}")
    except IndexError:
        print("Check that you have entered the correct species category")
        sys.exit()
    return min_reads_depth
def process_header(head, threthold, file):
    kmer_length = 21
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
    high_frequency_kmer = set()
    for key, value in dic.items():
        if value >= threthold:
            high_frequency_kmer.add(key)
    print(f"{head} consume memory {sys.getsizeof(dic) / (1024.0 ** 3)} GB ")
    del dic
    return high_frequency_kmer

def get_porportion(input,line_number,proportion,output_dir):
    with open(input, "r") as f:
        base_line = list(range(1, line_number, 2))
        select_line = set(random.sample(base_line, int(len(base_line) * proportion)))
        f.seek(0)
        with open(os.path.join(output_dir, "proportion.fa"), "w") as f_w:
            for i, line in enumerate(f):
                if i + 1 in select_line or i in select_line:
                    f_w.write(line)
    return os.path.join(output_dir, "proportion.fa")



def obtain_extract_file(in_file,output_dir,header,high_depth_kmer,kmer_length,accuracy):
    with open(in_file, "r") as f:
        with open(os.path.join(output_dir, "extract.fa"), "w") as w:
            for i, line in enumerate(f):
                if i % 2 == 0:
                    id_line = line
                elif i % 2 == 1:
                    if len(line.strip()) > kmer_length:
                        kmers = set()
                        for head in header:
                            start = 0
                            while True:
                                start = line.strip().find(head, start)
                                if start == -1 or start + kmer_length > len(line.strip()):
                                    break
                                kmer = line.strip()[start:start + 21]
                                start += 1
                                kmers.add(kmer)
                        if len(kmers) > 0 and len(kmers & high_depth_kmer) / len(kmers) >= accuracy:
                            w.write(id_line)
                            w.write(line)

def add_obtain_extract_file(in_file,output_dir,header,high_depth_kmer,kmer_length,accuracy,proportion):
    with open(in_file, "r") as f:
        line_number = sum(1 for line in f)
        base_line = list(range(1, line_number, 2))
        select_line = set(random.sample(base_line, int(len(base_line) * proportion)))
        f.seek(0)
        with open(os.path.join(output_dir, "extract.fa"), "w") as w:
            for i, line in enumerate(f):
                if i + 1 in select_line or i in select_line:
                    if i % 2 == 0:
                        id_line = line
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
                            if len(kmers) > 0 and len(kmers & high_depth_kmer) / len(kmers) >= accuracy:
                                w.write(id_line)
                                w.write(line)
    print(f"({datetime.datetime.now()}) complete extracting")

def filter(args):
    config_path = os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), "config_database")
    input_file=args.input_file
    output_dir=args.output_dir
    os.makedirs(os.path.join(output_dir,"blast_output"), exist_ok=True)
    blast_output=os.path.join(output_dir,"blast_output")
    kmer_length = 21
    print(f"({datetime.datetime.now()}) processing file")

    if args.proportion >1 or args.proportion<0:
        print("The value of --proportion must be 0-1")
        sys.exit()
    if args.accuracy > 1 or args.accuracy < 0:
        print("The value of --accuracy must be 0-1")
        sys.exit()
    seq_number, sum_len, min_seq, max_seq, process_file = process_input_file(input_file, output_dir)
    if args.base_number == 3:
        if args.head_number < 1 or args.head_number > 64:
            print("if base number is 3,head number only can be between 1 to 64 ")
            sys.exit()
    elif args.base_number == 4:
        if args.head_number == 4:
            args.head_number = 16
        if args.head_number < 1 or args.head_number > 256:
            print("if base number is 4,head number only can be be tween 1 to 256")
            sys.exit()
    if args.base_number == 3 and args.head_number == 4:
        header = {x + y for x, y in zip(generate_header(1), random.sample(generate_header(2), 4))}
    else:
        header = set(random.sample(generate_header(args.base_number), args.head_number))
    print(f"({datetime.datetime.now()}) random head {header}")

    high_depth_kmer = set()

    if args.proportion:
        if args.proportion ==1:
            proportion_file=process_file
        else:
            proportion_file = get_porportion(process_file,seq_number*2,args.proportion,output_dir)

        if args.filter_depth:
            estimate_mitogenome_depth="None"

            if args.normalize_depth:
                print("you have manually input an filter threshold ,cann't normalize mitogemome depth")
                sys.exit()

            else:
                min_reads_depth = args.filter_depth
                with Pool(processes=args.thread) as pool:
                    results = pool.starmap(process_header,
                                           [(h, min_reads_depth, proportion_file) for h in header])
                for result in results:
                    high_depth_kmer.update(result)
                obtain_extract_file(proportion_file,output_dir,header,high_depth_kmer,kmer_length,args.accuracy)

        else:
            min_reads_depth=determining_mitogenome_depth(args.species, config_path,args.thread,proportion_file, blast_output)
            estimate_mitogenome_depth=min_reads_depth/args.proportion
            with Pool(processes=args.thread) as pool:
                results = pool.starmap(process_header,
                                       [(h, min_reads_depth * args.filter_percentage, proportion_file) for h in header])
            for result in results:
                high_depth_kmer.update(result)
            ##output designated genome depth
            if args.normalize_depth :
                if args.normalize_depth >min_reads_depth or args.normalize_depth<0:
                    print("retain the maximum mitogenome depth")
                    obtain_extract_file(proportion_file, output_dir, header, high_depth_kmer, kmer_length,args.accuracy)
                else:
                    n=args.normalize_depth/min_reads_depth
                    add_obtain_extract_file(proportion_file, output_dir, header, high_depth_kmer, kmer_length,args.accuracy,n)

            else:
                if min_reads_depth <= 50:
                    n = 1
                elif min_reads_depth <= 150:
                    n = 20 / min_reads_depth
                elif min_reads_depth <= 500:
                    n = 0.1
                else:
                    n = 30 / min_reads_depth

                add_obtain_extract_file(proportion_file,output_dir,header,high_depth_kmer,kmer_length,args.accuracy,n)


    ###auto subsample
    else:
        file_size = os.path.getsize(process_file)/(1024 ** 3)
        if file_size>100:
            auto_proportion= 1/20
        elif file_size > 20:
            auto_proportion= 1/10
        elif file_size >10:
            auto_proportion = 1/5
        elif file_size > 5:
            auto_proportion = 1/2
        else:
            auto_proportion =1

        if args.normalize_depth <0:
            min_reads_depth = determining_mitogenome_depth(args.species, config_path, args.thread, process_file,
                                                           blast_output)
            with Pool(processes=args.thread) as pool:
                results = pool.starmap(process_header,
                                       [(h, min_reads_depth * args.filter_percentage, process_file) for h in header])
            for result in results:
                high_depth_kmer.update(result)
            obtain_extract_file(process_file, output_dir, header, high_depth_kmer, kmer_length, args.accuracy)
            return seq_number, sum_len, min_seq, max_seq, min_reads_depth

        if auto_proportion !=1:
            proportion_file = get_porportion(process_file,seq_number*2,auto_proportion,output_dir)
            print(f"({datetime.datetime.now()}) your file size is {file_size} GB,\
Automatically reduce the data set to {auto_proportion} times of its original size.")
        else:
            proportion_file = process_file
            print(f"({datetime.datetime.now()}) your file size is {file_size} GB")
        min_reads_depth=determining_mitogenome_depth(args.species, config_path,args.thread,proportion_file, blast_output)
        estimate_mitogenome_depth = min_reads_depth / auto_proportion

        ###normalize
        if args.normalize_depth > 0 :
            if args.normalize_depth < min_reads_depth:
                with Pool(processes=args.thread) as pool:
                    results = pool.starmap(process_header,
                    [(h, min_reads_depth * args.filter_percentage, proportion_file) for h in header])
                for result in results:
                    high_depth_kmer.update(result)
                n=args.normalize_depth/min_reads_depth
                add_obtain_extract_file(proportion_file,output_dir,header,high_depth_kmer,kmer_length,args.accuracy,n)
            else:
                if args.normalize_depth < estimate_mitogenome_depth:
                    n = args.normalize_depth * auto_proportion / min_reads_depth
                    proportion_file = get_porportion(process_file, seq_number * 2, n, output_dir)
                else:
                    proportion_file=process_file

                min_reads_depth = determining_mitogenome_depth(args.species, config_path, args.thread, proportion_file,
                                                               blast_output)

                with Pool(processes=args.thread) as pool:
                    results = pool.starmap(process_header,
                                           [(h, min_reads_depth * args.filter_percentage,proportion_file) for h in header])
                for result in results:
                    high_depth_kmer.update(result)
                obtain_extract_file(proportion_file,output_dir,header,high_depth_kmer,kmer_length,args.accuracy)

            return seq_number, sum_len, min_seq, max_seq,estimate_mitogenome_depth

        ####add data
        if min_reads_depth < 40 and auto_proportion != 1:
            ###if depth not exceed 120，tbalstn all
            if estimate_mitogenome_depth< 120 :
                proportion_file=process_file
                min_reads_depth=determining_mitogenome_depth(args.species, config_path,args.thread,proportion_file, blast_output)

            else:
                n= 120 *auto_proportion /min_reads_depth
                proportion_file = get_porportion(process_file, seq_number * 2, n, output_dir)
                min_reads_depth=determining_mitogenome_depth(args.species, config_path,args.thread,proportion_file, blast_output)
        with Pool(processes=args.thread) as pool:
            results = pool.starmap(process_header,
                                   [(h, min_reads_depth * args.filter_percentage, proportion_file) for h in header])
        for result in results:
            high_depth_kmer.update(result)

        if min_reads_depth <= 50:
            n = 1
        elif min_reads_depth <= 150:
            n = 20 / min_reads_depth
        elif min_reads_depth <= 500:
            n = 0.1
        else:
            n = 30 / min_reads_depth
        add_obtain_extract_file(proportion_file,output_dir,header,high_depth_kmer,kmer_length,args.accuracy,n)
    return seq_number, sum_len, min_seq, max_seq,estimate_mitogenome_depth
