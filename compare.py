#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import numpy as np
import random
import plotly.graph_objects as go
from assess import info
from assess import draw_line


def random_choice_color():
    color = random.choice(["rgba(49,124,183,0.3)", "rgba(183,34,48,0.3)", "rgba(246,178,147,0.6)"])
    return color
def compare(args):
    fig = go.Figure()
    args.input_file = args.query
    query = info(args)
    query.gap = np.pi / 45
    query.int()
    if query.gfa:
        query_fa = os.path.join(args.output_dir, 'assessment_out', "query.fa")
        os.rename(os.path.join(args.output_dir, 'assessment_out', "assess.fa"), query_fa)
    else:
        query_fa=args.query

    args.input_file=args.reference
    reference=info(args)
    reference.gap=np.pi/45
    reference.int()
    if reference.gfa:
        reference_fa=os.path.join(args.output_dir, 'assessment_out', "reference.fa")
        os.rename(os.path.join(args.output_dir, 'assessment_out', "assess.fa"), reference_fa)
    else:
        reference_fa=args.reference

    genome_gap = query.gap
    query.display = query.sum_len / (query.sum_len + reference.sum_len) * (2 * np.pi - 2 * genome_gap)
    reference.display = reference.sum_len / (query.sum_len + reference.sum_len) * (2 * np.pi - 2 * genome_gap)
    query.shift = np.pi / 2
    reference.shift = genome_gap + query.display + np.pi / 2

    ###draw synteny bar
    # query.draw_chromosome(fig,450,470)
    # reference.draw_chromosome(fig,450,470)

    ###draw GC content
    query.draw_GC_content(fig,550,620)
    reference.draw_GC_content(fig,550,620)
    ###draw chromosome
    query_color=random_choice_color()
    query.draw_chromosome(fig,690,760,same_color=query_color)
    reference_color=random_choice_color()
    while query_color==reference_color:
        reference_color=random_choice_color()

    reference.draw_chromosome(fig,690,760,same_color=reference_color)

    query.gfa = False
    reference.gfa = False
    query.draw_contig_depth(fig,690,760)

    reference.draw_contig_depth(fig,690,760)
    ####draw outer gene
    if args.category != "other":

        query.args.input_file = query_fa
        query.run_out_program()
        po_info, po_cds, ne_info, ne_cds,array, gene_list, max_copy=query.process_output_file()
        query.draw_gene(fig, po_info, po_cds, 1000, 20, 'red')
        query.draw_gene(fig, ne_info, ne_cds, 1000, -20, 'blue')

        reference.args.input_file = reference_fa
        reference.run_out_program()
        po_info, po_cds, ne_info, ne_cds, array, gene_list, max_copy =reference.process_output_file()
        reference.draw_gene(fig, po_info, po_cds, 1000, 20, 'red')
        reference.draw_gene(fig, ne_info, ne_cds, 1000, -20, 'blue')

    ###draw synteny
    synteny_db=os.path.join(args.output_dir,"synteny_db")
    synteny_out=os.path.join(args.output_dir,'synteny_blastn_out')
    os.system(f"makeblastdb -in {reference_fa} -dbtype nucl -out {synteny_db}")
    os.system(f"blastn  -db {synteny_db} -query {query_fa} -evalue 1e-10 -out {synteny_out} -outfmt 6")
    with open(synteny_out,'r') as f:
        region=[]
        for line in f:
            line_content = line.strip().split('\t')
            query_name = line_content[0]
            ref_name = line_content[1]
            query_start = int(line_content[6])
            query_end = int(line_content[7])
            ref_start = int(line_content[8])
            ref_end = int(line_content[9])
            if ref_start<ref_end:
                region.append((query_name,query_start,query_end,ref_name,ref_start,ref_end))
            else:
                region.append((query_name, query_start, query_end, ref_name,ref_end,ref_start))
        if len(region)>0:
            sorted_region=sorted(region,key=lambda x: (x[0], x[1], x[2]))
            query_need=[sorted_region[0]]
            for item in sorted_region:
                query_name, start, end = item[:3]
                chr,prev_start,prev_end=query_need[-1][:3]
                if query_name==chr:
                    if prev_start <= start and prev_end >= end :
                        continue
                    elif start <= prev_start and end >= prev_end:
                        query_need[-1] = item
                    elif (prev_end-start)/(end-start)>=0.8 and prev_end>start:
                        continue
                    else:
                        query_need.append(item)
                else:
                    query_need.append(item)
            sorted_region=sorted(region,key=lambda x: (x[3],x[4],x[5]))
            ref_need=[sorted_region[0]]
            for item in sorted_region:
                ref_name,start,end=item[3:]
                chr, prev_start, prev_end = ref_need[-1][3:]
                if ref_name==chr:
                    if prev_start <= start and prev_end >= end :
                        continue
                    elif start <= prev_start and end >= prev_end:
                        ref_need[-1] = item
                    elif (prev_end-start)/(end-start)>=0.8 and prev_end>start:
                        continue
                    else:
                        ref_need.append(item)
                else:
                    ref_need.append(item)
        draw_region=set(query_need)|set(ref_need)
        for item in draw_region:
            query_start_radian,query_end_radian=query.return_scope(item[0],item[1],item[2])
            ref_start_radian,ref_end_radian=reference.return_scope(item[3],item[4],item[5])
            radian_scope = max(abs(query_end_radian - query_start_radian), abs(ref_end_radian - ref_start_radian))
            line_number = int(radian_scope * 600 / np.pi + 1)
            query_theta = np.linspace(query_start_radian, query_end_radian, line_number)
            ref_theta = np.linspace(ref_end_radian, ref_start_radian, line_number)
            synteny_color = random_choice_color()
            ##synteny bar
            uniform_numbers = np.linspace(450, 470, 50)
            for r in uniform_numbers:
                draw_line(fig, query_start_radian, query_end_radian, r, query_name, width=2, color=synteny_color)
                draw_line(fig, ref_start_radian, ref_end_radian, r, ref_name, width=2, color=synteny_color)

            for i, value in enumerate(query_theta):
                if i==0 or i==len(query_theta)-1:
                    color="black"
                else:
                    color=synteny_color
                x1 = np.cos(value) * 400
                y1 = np.sin(value) * 400
                x2 = np.cos(ref_theta[i]) * 400
                y2 = np.sin(ref_theta[i]) * 400
                cx, cy = 0, 0
                t_values = np.linspace(0, 1, 50)
                x = (1 - t_values) ** 2 * x1 + 2 * (1 - t_values) * t_values * cx + t_values ** 2 * x2
                y = (1 - t_values) ** 2 * y1 + 2 * (1 - t_values) * t_values * cy + t_values ** 2 * y2
                fig.add_trace(
                    go.Scatter(
                        x=x,
                        y=y,
                        mode='lines',
                        name=None,
                        hoverinfo='name',
                        line=dict(color=color, width=1)
                    )
                )

    fig.update_layout(
        width=900,
        height=900,
        showlegend=False,
        plot_bgcolor='rgba(255, 255,255,0)',
        paper_bgcolor='rgba(255, 255,255,0)',
        xaxis=dict(title_text=None, range=[-1150, 1150],
                   showgrid=False, scaleanchor='y',
                   showticklabels=False,
                   zeroline=False, ),
        yaxis=dict(title_text=None, range=[-1100, 1100],
                   showgrid=False,
                   showticklabels=False,
                   zeroline=False, )
    )
    config = {
        'toImageButtonOptions': {
            'format': 'svg'
        }
    }
    fig.write_html(os.path.join(args.output_dir, "himt_compara.html"), config=config)
