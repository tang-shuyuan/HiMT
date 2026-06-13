#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys
from assemble import simple_gfa
import datetime
from assess import assess
from filter import open_input_text

def refassemble(args):
    ###converting gfa to fa
    with open(args.ref_genome,'r') as f:
        print(f"({datetime.datetime.now()}) checking refgenome format ")
        if f.read(1)==">":
            print(f"({datetime.datetime.now()}) refgenome is format ")
            format="fasta"
        else:
            print(f"({datetime.datetime.now()}) converting gfa to fasta")
            format="gfa"
        f.seek(0)
        if format=="gfa":
            with open(os.path.join(args.output_dir,"refgenome.fa"),'w') as f_w:
                for line in f:
                    if line.startswith("S"):
                        line_content = line.strip().split("\t")
                        f_w.write(">" + line_content[1] + "\n")
                        f_w.write(line_content[2] + "\n")
            args.ref_genome =os.path.join(args.output_dir,"refgenome.fa")

    if args.data_type.lower()=="hifi":
        map_option="map-hifi"
    elif args.data_type.lower()=="clr":
        map_option="map-pb"
    elif args.data_type.lower()=="ont":
        map_option="map-ont"
    print(f"({datetime.datetime.now()}) runing minimap2")
    os.system(f"minimap2 -x {map_option} {args.ref_genome} {args.input_file} -o {os.path.join(args.output_dir,'minimap2.paf')}")

    candidata_read=set()
    ###get mito id
    print(f"({datetime.datetime.now()}) extracting reads")
    with open(os.path.join(args.output_dir,"minimap2.paf"),"r") as f:
        for line in f:
            line_content=line.strip().split("\t")
            ref_length=int(line_content[6])
            read_length=int(line_content[1])
            ####
            if abs(int(line_content[7])-int(line_content[8]))/ref_length>= 0.9:
                candidata_read.add(line_content[0])
            elif abs(int(line_content[3])-int(line_content[2]))/read_length>args.percent and read_length>args.length:
                candidata_read.add(line_content[0])
    ###extract reads

    with open(os.path.join(args.output_dir,"candidata.fa"),"w",newline='\n',buffering=1024*1024*8) as f_w:
        f_in=open_input_text(args.input_file)
        first_item=f_in.read(1)
        f_in.seek(0)
        if first_item=='>':
            current_header = None
            current_sequence_lines = []
            for line in f_in:
                if line.startswith('>'):
                    if current_header is not None and current_header.strip().split()[0][1:] in candidata_read:
                        f_w.write(current_header)
                        f_w.writelines(current_sequence_lines)
                    current_header = line
                    current_sequence_lines = []
                else:
                    current_sequence_lines.append(line)
            if current_header is not None and current_header.strip().split()[0][1:] in candidata_read:
                f_w.write(current_header)
                f_w.writelines(current_sequence_lines)
        elif first_item=='@':
            for line in f_in:
                if line.startswith('@') and line.strip().split()[0][1:] in candidata_read:
                    f_w.write(">"+line[1:]+next(f_in))
        else:
            print("The file format may be incorrect")
            sys.exit(1)
        f_w.close()

    read_type_param = {"HiFi": "--pacbio-hifi", "ONT": "--nano-corr", "CLR": "--pacbio-corr"}
    if args.no_flye_meta:
        os.system(f"flye {read_type_param[args.data_type]} {os.path.join(args.output_dir, 'candidata.fa')} \
                -o {os.path.join(args.output_dir, 'flye_output')} -t {args.thread} ")
    else:
        os.system(f"flye {read_type_param[args.data_type]} {os.path.join(args.output_dir, 'candidata.fa')} \
        -o {os.path.join(args.output_dir, 'flye_output')} -t {args.thread} --meta")
    if os.path.exists(os.path.join(args.output_dir, 'flye_output','assembly_graph.gfa')):
        simple_gfa(args)
        if args.species == "plant":
            args.input_file = os.path.join(args.output_dir, "himt_refassemble.gfa")
            assess(args)
    else:
        print("flye not work")
        sys.exit(1)
