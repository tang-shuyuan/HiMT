#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import gzip
from assemble import simple_gfa
import datetime
from assess import auto_judge_category
from assess import assess
def refassemble(args):


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
    try:
        with gzip.open(args.input_file,'rb') as f,open(os.path.join(args.output_dir,"candidata.fa"),"w") as f_w:
            first_item=f.read(1).encode('utf-8')
            f.seek(0)
            if first_item=='>':
                for line in f:
                    line=line.decode('utf-8')
                    if line.startswith('>'):
                        if line.strip().split()[0][1:] in candidata_read:
                            f_w.write(line)
                            wirte_sequence = True
                        else:
                            wirte_sequence = False
                    elif wirte_sequence:
                        f_w.write(line)
            elif first_item=='@':
                if line.strip().split()[0][1:] in candidata_read:
                    f_w.write(">"+line[1:])
                    next_line = next(f).decode('utf-8')
                    f_w.write(next_line)
    except:
        with open(args.input_file, 'r') as f, open(os.path.join(args.output_dir, "candidata.fa"), "w") as f_w:
            first_item = f.read(1)
            f.seek(0)
            if first_item == '>':
                for line in f:
                    if line.startswith('>'):
                        if line.strip().split()[0][1:] in candidata_read:
                            f_w.write(line)
                            wirte_sequence = True
                        else:
                            wirte_sequence = False
                    elif wirte_sequence:
                        f_w.write(line)
            elif first_item == '@':
                if line.strip().split()[0][1:] in candidata_read:
                    f_w.write(">"+line[1:])
                    next_line = next(f).decode('utf-8')
                    f_w.write(next_line)
    read_type_param = {"HiFi": "--pacbio-hifi", "ONT": "--nano-corr", "CLR": "--pacbio-corr"}
    if args.no_flye_meta:
        os.system(f"flye {read_type_param[args.data_type]} {os.path.join(args.output_dir, 'candidata.fa')} \
                -o {os.path.join(args.output_dir, 'flye_output')} -t {args.thread}")
    else:
        os.system(f"flye {read_type_param[args.data_type]} {os.path.join(args.output_dir, 'candidata.fa')} \
        -o {os.path.join(args.output_dir, 'flye_output')} -t {args.thread} --meta")

    simple_gfa(args)
    if args.species=="plant":
        args.input_file=os.path.join(args.output_dir,"himt_refassemble.gfa")
        assess(args)
