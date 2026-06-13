#!/usr/bin/env python
# -*- coding: utf-8 -*-
import gzip
import os
import datetime
import random
from multiprocessing import Pool
import sys
import subprocess

def is_empty_file(path):
    try:
        return os.path.exists(path) and os.path.getsize(path)== 0
    except OSError:
        return None

def read_file_by_chunk(path, chunk_size=1024*1024*8, encoding="utf-8"):
    def _iter_chunks(f):
        carry = ""
        while True:
            chunk = f.read(chunk_size)
            if not chunk:
                break
            text = carry + chunk
            lines = text.split('\n')
            carry = lines.pop()
            if lines:
                yield '\n'.join(lines) + '\n'
        if carry:
            yield carry
    try:
        with gzip.open(path, 'rt', encoding=encoding) as f:
            yield from _iter_chunks(f)
    except (OSError, gzip.BadGzipFile):
        with open(path, 'r', encoding=encoding) as f:
            yield from _iter_chunks(f)

def open_input_text(path):
    try:
        f = gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="")
        f.read(1)      # 触发解压校验
        f.seek(0)
        return f
    except (OSError, gzip.BadGzipFile):
        f = open(path, "rt", encoding="utf-8", errors="replace", newline="")
        return f

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

def process_input_file(input_file,output_dir):
    seq_number = 0
    min_seq = float('inf')
    max_seq = 0
    sum_len = 0
    f_in = open_input_text(input_file)

    first_char = f_in.read(1)
    f_in.seek(0)
    if first_char == ">":  # FASTA
        unstardan_format = False
        for i, line in enumerate(f_in):
            if i in {2, 4, 6, 8, 10}:
                if not line.startswith(">"):
                    unstardan_format = True
                    print(f"({datetime.datetime.now()}) Modifying the file format")
                    break
            elif i == 11:
                break
        if unstardan_format:
            with open(os.path.join(output_dir,"process.fa"), "w",buffering=1024*1024*8) as f_w:
                current_seq_id = ""
                current_seq = ""
                for chunk in read_file_by_chunk(input_file):
                    for line in chunk.splitlines():
                        if not line:
                            continue
                        if line.startswith('>'):
                            if current_seq_id and current_seq:
                                seq_len = len(current_seq)
                                seq_number += 1
                                sum_len += seq_len
                                max_seq = max(max_seq, seq_len)
                                min_seq = min(min_seq, seq_len)
                                f_w.write(f"{current_seq_id}\n")
                                f_w.write(f"{current_seq}\n")
                                current_seq = ""
                            current_seq_id = line
                        else:
                            current_seq += line
                if current_seq_id and current_seq:
                    seq_len = len(current_seq)
                    seq_number += 1
                    sum_len += seq_len
                    max_seq = max(max_seq, seq_len)
                    min_seq = min(min_seq, seq_len)
                    f_w.write(f"{current_seq_id}\n")
                    f_w.write(f"{current_seq}\n")

            f_in.close()
            return seq_number, sum_len, min_seq, max_seq, os.path.join(output_dir,"process.fa")

        else:
            for chunk in read_file_by_chunk(input_file):
                for line in chunk.splitlines():
                    if not line:
                        continue
                    if line.startswith(">"):
                        seq_number += 1
                    else:
                        max_seq = max(max_seq, len(line.strip()))
                        min_seq = min(min_seq, len(line.strip()))
                        sum_len += len(line.strip())
            f_in.close()
            return seq_number, sum_len, min_seq, max_seq, input_file

    elif first_char == '@':
        with open(os.path.join(output_dir, "process.fa"), "w", buffering=1024 * 1024 * 8) as f_w:
            state=0
            print(f"({datetime.datetime.now()}) converting fq to fa")
            for chunk in read_file_by_chunk(input_file):
                for line in chunk.splitlines():
                    if not line:
                        continue
                    if state == 0:
                        if not line.startswith("@"):
                            raise ValueError(f"Wrong fastq format: {line}")
                        f_w.write(f'>{line[1:]}\n')
                        state = 1
                    elif state == 1:
                        seq_number += 1
                        sum_len += len(line.strip())
                        max_seq = max(max_seq, len(line.strip()))
                        min_seq = min(min_seq, len(line.strip()))
                        f_w.write(f'{line}\n')
                        state = 2
                    elif state == 2:
                        state=3
                    elif state == 3:
                        state=0
        f_in.close()
        return seq_number, sum_len, min_seq, max_seq, os.path.join(output_dir, "process.fa")
    else:
        print("file type unknow")
        f_in.close()
        sys.exit()

def determining_mitogenome_depth(args,process_file):
    config_path = os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), "config_database")
    os.makedirs(os.path.join(args.output_dir, "blast_output"), exist_ok=True)

    if args.species == "animal":
        prot_sequence = os.path.join(config_path, "Drosophila_gunungcola_CM045947.fasta")
    else:
        prot_sequence = os.path.join(config_path, "Arabidopsis_protein.fasta")
    if args.data_type == "HiFi":
        e_value="1e-10"
    else:
        e_value="1e-5"

    out_blast_db = os.path.join(args.output_dir, "blast_output","database")
    blast_result = os.path.join(args.output_dir,"blast_output","blast_result")

    command1 = ["makeblastdb", "-in", process_file, "-dbtype", "nucl", "-out", out_blast_db]
    result1=subprocess.run(command1,stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                       text=True,check=True)
    if result1.stdout:
        print(result1.stdout, end="")

    command2 = ["tblastn", "-num_descriptions","100000", "-num_alignments","100", "-num_threads",str(args.thread),\
       "-db", out_blast_db, "-query" ,prot_sequence, "-evalue", e_value, "-out", blast_result]
    subprocess.run(command2,text=True,check=True)

    if is_empty_file(blast_result):
        print("blast_result is an empty file. Please ensure your file is in standard FASTA or FASTQ format.")
        sys.exit()
    elif is_empty_file(blast_result) is None:
        print("Couldn't find file blast_result")
        sys.exit()

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

def process_header(head, threthold,file,kmer_length):
    dic = {}
    i=-1
    for chunk in read_file_by_chunk(file):
        for line in chunk.splitlines():
            i+=1
            if not line:
                continue
            if i % 2 == 1:
                start = 0
                while True:
                    start = line.strip().find(head, start)
                    if start == -1 or start + kmer_length > len(line.strip()):
                        break
                    kmer = line.strip()[start:start + kmer_length]
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

def get_porportion(path,line_number,proportion,output_dir):
    base_line = list(range(1, line_number, 2))
    select_line = set(random.sample(base_line, int(len(base_line) * proportion)))
    line_id=-1
    with open(os.path.join(output_dir, "proportion.fa"), "w", buffering=1024 * 1024 * 8) as f_w:
        for chunk in read_file_by_chunk(path):
            for line in chunk.splitlines():
                line_id += 1
                if line_id+1 in select_line or line_id in select_line:
                    f_w.write(f"{line}\n")
    return os.path.join(output_dir, "proportion.fa")

def obtain_extract_file(args,in_file,proportion,header,high_depth_kmer):
    with open(os.path.join(args.output_dir, "extract.fa"), "w",buffering=1024*1024*8) as w:
        base_line_number = 0
        for chunk in read_file_by_chunk(in_file):
            for line in chunk.splitlines():
                if line and line.startswith(">"):
                    base_line_number += 1
        line_id=-1
        need_line =set(random.sample(list(range(1, base_line_number*2, 2)),int(base_line_number*proportion)))
        for chunk in read_file_by_chunk(in_file):
            for line in chunk.splitlines():
                line_id += 1
                if line_id + 1 in need_line or line_id  in need_line:
                    if line_id % 2 == 0:
                        id_line = line
                    elif line_id % 2 == 1:
                        if len(line.strip()) > args.kmer_length:
                            kmers = set()
                            # for j in range(len(line.strip()) - kmer_length + 1):
                            #     kmer = line.strip()[j:j + kmer_length]
                            #     if kmer[0:3] in header:
                            #         kmers.add(kmer)
                            for head in header:
                                start = 0
                                while True:
                                    start = line.strip().find(head, start)
                                    if start == -1 or start + args.kmer_length > len(line.strip()):
                                        break
                                    kmer = line.strip()[start:start + args.kmer_length]
                                    start += 1
                                    kmers.add(kmer)
                            if len(kmers) > 0 and len(kmers & high_depth_kmer) / len(kmers) >= args.accuracy:
                                w.write(f"{id_line}\n")
                                w.write(f"{line}\n")
    print(f"({datetime.datetime.now()}) complete extracting")

def filter(args):
    print(f"({datetime.datetime.now()}) processing file")
    if args.__internal_seed:
        random_seed = args.__internal_seed
    else:
        random_seed = random.randint(1, 64)

    print(f"({datetime.datetime.now()}) random seed {random_seed}")

    if args.proportion >1 or args.proportion<0:
        print("The value of --proportion must be 0-1")
        sys.exit()
    if args.accuracy > 1 or args.accuracy < 0:
        print("The value of --accuracy must be 0-1")
        sys.exit()
    seq_number, sum_len, min_seq, max_seq, process_file = process_input_file(args.input_file, args.output_dir)

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
            proportion_file = get_porportion(process_file,seq_number*2,args.proportion,args.output_dir)

        if args.filter_depth:
            if args.normalize_depth:
                print("you have manually input an filter threshold ,cann't normalize mitogemome depth")
                sys.exit()
            else:
                min_reads_depth = args.filter_depth
                with Pool(processes=args.extract_parallel) as pool:
                    results = pool.starmap(process_header,
                                           [(h, min_reads_depth, proportion_file,args.kmer_length) for h in header])
                for result in results:
                    high_depth_kmer.update(result)

                obtain_extract_file(args,proportion_file,1,header,high_depth_kmer)
            reduction_radio=args.proportion
        else:
            min_reads_depth=determining_mitogenome_depth(args,proportion_file)

            with Pool(processes=args.extract_parallel) as pool:
                results = pool.starmap(process_header,
                [(h, min_reads_depth * args.filter_percentage, proportion_file,args.kmer_length) for h in header])
            for result in results:
                high_depth_kmer.update(result)
            ##output designated genome depth
            if args.normalize_depth :
                if args.normalize_depth >min_reads_depth or args.normalize_depth<0:
                    print("retain the maximum mitogenome depth")
                    obtain_extract_file(args,proportion_file,1,header,high_depth_kmer)
                    reduction_radio = args.proportion
                else:
                    n=args.normalize_depth/min_reads_depth
                    obtain_extract_file(args,proportion_file,n,header,high_depth_kmer)
                    reduction_radio = args.proportion*n

            else:
                if min_reads_depth <= 50:
                    n = 1
                elif min_reads_depth <= 150:
                    n = 20 / min_reads_depth
                elif min_reads_depth <= 500:
                    n = 0.1
                else:
                    n = 30 / min_reads_depth

                obtain_extract_file(args,proportion_file,n,header,high_depth_kmer)
                reduction_radio = args.proportion * n

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
        ###output maximum depth
        if args.normalize_depth <0:
            print("retain the maximum mitogenome depth")
            min_reads_depth = determining_mitogenome_depth(args,process_file)
            with Pool(processes=args.extract_parallel) as pool:
                results = pool.starmap(process_header,
                [(h, min_reads_depth * args.filter_percentage, process_file,args.kmer_length) for h in header])
            for result in results:
                high_depth_kmer.update(result)
            obtain_extract_file(args,process_file,1,header,high_depth_kmer)
            return seq_number, sum_len, min_seq, max_seq, 1

        if auto_proportion !=1:
            proportion_file = get_porportion(process_file,seq_number*2,auto_proportion,args.output_dir)
            print(f"({datetime.datetime.now()}) your file size is {file_size} GB,\
Automatically reduce the data set to {auto_proportion} times of its original size.")
        else:
            proportion_file = process_file
            print(f"({datetime.datetime.now()}) your file size is {file_size} GB")
        min_reads_depth=determining_mitogenome_depth(args,proportion_file)
        estimate_mitogenome_depth = min_reads_depth / auto_proportion
        temp_reduction_radio=auto_proportion

        ###normalize
        if args.normalize_depth > 0:
            ###
            if (abs(args.normalize_depth-min_reads_depth)<=5 and args.normalize_depth<=25) or \
                    (abs(args.normalize_depth-min_reads_depth)<=10 and 25<args.normalize_depth<=45) or\
                    (abs(args.normalize_depth-min_reads_depth)<=20 and 45<args.normalize_depth<=120) or\
                    (abs(args.normalize_depth-min_reads_depth)<=40 and args.normalize_depth>120):
                with Pool(processes=args.extract_parallel) as pool:
                    results = pool.starmap(process_header,
                    [(h, min_reads_depth * args.filter_percentage,proportion_file,args.kmer_length) for h in header])
                for result in results:
                    high_depth_kmer.update(result)
                obtain_extract_file(args,proportion_file,1,header,high_depth_kmer)
                reduction_radio = auto_proportion
            else:
                if args.normalize_depth < min_reads_depth:
                    with Pool(processes=args.extract_parallel) as pool:
                        results = pool.starmap(process_header,
                        [(h, min_reads_depth * args.filter_percentage, proportion_file,args.kmer_length) for h in header])
                    for result in results:
                        high_depth_kmer.update(result)
                    n=args.normalize_depth/min_reads_depth
                    obtain_extract_file(args,proportion_file,n,header,high_depth_kmer)
                    reduction_radio = auto_proportion*n
                else:
                    if args.normalize_depth < estimate_mitogenome_depth:
                        n = args.normalize_depth * auto_proportion / min_reads_depth
                        proportion_file = get_porportion(process_file, seq_number * 2, n, args.output_dir)
                        reduction_radio = n
                    else:
                        proportion_file=process_file
                        print("retain the maximum mitogenome depth")
                        reduction_radio = 1

                    min_reads_depth = determining_mitogenome_depth(args, proportion_file)

                    with Pool(processes=args.extract_parallel) as pool:
                        results = pool.starmap(process_header,
                        [(h, min_reads_depth * args.filter_percentage,proportion_file,args.kmer_length) for h in header])
                    for result in results:
                        high_depth_kmer.update(result)
                    obtain_extract_file(args,proportion_file,1,header,high_depth_kmer)

            return seq_number, sum_len, min_seq, max_seq,reduction_radio

        ####add data
        if min_reads_depth < 40 and auto_proportion != 1:
            ###if depth not exceed 120，tbalstn all
            if estimate_mitogenome_depth< 120 :
                proportion_file=process_file
                min_reads_depth=determining_mitogenome_depth(args,proportion_file)
                temp_reduction_radio = 1
            else:
                n= 120 *auto_proportion /min_reads_depth
                proportion_file = get_porportion(process_file, seq_number * 2, n, args.output_dir)
                min_reads_depth=determining_mitogenome_depth(args,proportion_file)
                temp_reduction_radio = n
        with Pool(processes=args.extract_parallel) as pool:
            results = pool.starmap(process_header,
            [(h, min_reads_depth * args.filter_percentage, proportion_file,args.kmer_length) for h in header])
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
        obtain_extract_file(args, proportion_file, n, header, high_depth_kmer)
        reduction_radio = temp_reduction_radio * n
    return seq_number, sum_len, min_seq, max_seq,reduction_radio