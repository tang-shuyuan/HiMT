import os
import sys
import datetime
import random
from multiprocessing import Pool
from fliter import process_input_file,generate_header,determining_mitogenome_depth,process_header
from assess import assese

def simple_gfa(assemble_species,input_file,output_dir,config_path):
    with open(input_file, "r") as f_in, \
            open(os.path.join(output_dir,"blast_output","out.fa"), "w") as f_w:
        for line in f_in:
            if line.startswith("S"):
                line_content = line.strip().split("\t")
                f_w.write(">" + line_content[1] + "\n")
                f_w.write(line_content[2] + "\n")
    fa_file = os.path.join(output_dir,"blast_output","out.fa")

    if assemble_species == "animal":
        prot_sequence = os.path.join(config_path, "Drosophila_gunungcola_CM045947.fasta")
    else:
        prot_sequence = os.path.join(config_path, "Arabidopsis_protein.fasta")
    out_blast_db = os.path.join(output_dir,"blast_output","temp_db")
    blast_result = os.path.join(output_dir, "blast_output","temp_blast_result")

    command4 = f"makeblastdb -in {fa_file} -dbtype nucl -out {out_blast_db}"
    os.system(command4)
    command5 = f"tblastn -num_descriptions 100000 -num_alignments 100 \
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
        id_new = id.copy()
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
                id = id_new.copy()
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
    print(f"({datetime.datetime.now()}) The size of the folder is: {total_size/(1024**3)} GB")


def assembly(args):
    output_dir = args.output_dir
    kmer_length = 21
    config_path = os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), "config_database")
    blast_output = os.path.join(output_dir, "blast_output")
    print(f"({datetime.datetime.now()}) processing file")
    process_file = process_input_file(args.input_file, output_dir)

    if args.proportion > 1 or args.proportion < 0:
        print("The value of --proportion must be 0-1")
        sys.exit()
    if args.accuracy > 1 or args.accuracy < 0:
        print("The value of --accuracy must be 0-1")
        sys.exit()
    if args.proportion !=1:
        with open(process_file, "r") as f:
            line_number = sum(1 for line in f)
            base_line = list(range(1, line_number, 2))
            select_line = set(random.sample(base_line, int(len(base_line) * args.proportion)))
            f.seek(0)
            with open(os.path.join(output_dir,"proportion.fa"), "w") as f_w:
                for i, line in enumerate(f):
                    if i + 1 in select_line or i in select_line:
                        f_w.write(line)
        process_file=os.path.join(output_dir,"proportion.fa")


    if args.base_number == 3:
        if args.head_number < 4 or args.head_number > 64:
            print("if base number is 3,head number only can be between 4 to 64 ")
            sys.exit()
    elif args.base_number == 4:
        if args.head_number == 4:
            args.head_number = 16
        if args.head_number < 16 or args.head_number > 256:
            print("if base number is 4,head number only can be be tween 16 to 256")
            sys.exit()
    else:
        print("title base number only can be 3 and 4")

    if args.base_number == 3 and args.head_number == 4:
        header = {x + y for x, y in zip(generate_header(1), random.sample(generate_header(2), 4))}
    else:
        header = set(random.sample(generate_header(args.base_number), args.head_number))
    print(f"({datetime.datetime.now()}) random head {header}")


    if args.fliter_depth:
        min_reads_depth=args.fliter_depth
        with Pool(processes=args.threat) as pool:
            results = pool.starmap(process_header, [(h,min_reads_depth,process_file) for h in header])
    else:
        min_reads_depth=determining_mitogenome_depth(args.species, config_path,args.thread,process_file, blast_output)

        with Pool(processes=args.thread) as pool:
            results = pool.starmap(process_header,[(h, min_reads_depth * args.fliter_percentage, process_file) for h in header])
    high_depth_kmer = set()
    for result in results:
        high_depth_kmer.update(result)
    ###提取reads
    if args.normalize_depth :
        if args.normalize_depth<0 or args.normalize_depth> min_reads_depth:
            n=1
        else:
            n= args.normalize_depth/min_reads_depth
    else:
        if min_reads_depth <= 50:
            n = 1
        elif min_reads_depth <= 150:
            n = 20 / min_depth
        elif min_reads_depth <= 500:
            n = 0.1
        else:
            n = 30 / min_reads_depth
    with open(process_file, "r") as f:
        line_number = sum(1 for line in f)
        base_line = list(range(1, line_number, 2))
        select_line = set(random.sample(base_line, int(len(base_line) * n)))
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
                            if len(kmers) > 0 and len(kmers & high_depth_kmer) / len(kmers) >= args.accuracy:
                                w.write(id_line)
                                w.write(line)
    print(f"({datetime.datetime.now()}) complete extracting")
    if args.no_flye_meta:
        os.system(f"flye --pacbio-hifi {os.path.join(output_dir, 'extract.fa')} \
            -o {os.path.join(output_dir, 'flye_output')} -t {args.thread}")
    else:
        os.system(f"flye --pacbio-hifi {os.path.join(output_dir,'extract.fa')} \
    -o {os.path.join(output_dir,'flye_output')} -t {args.thread} --meta")

    simple_gfa(args.species,os.path.join(output_dir,"flye_output/assembly_graph.gfa"),output_dir,config_path)
    args.input_file=os.path.join(output_dir, "extract.gfa")
    args.func=assese
    assese(args)
    get_hard_disk(output_dir)