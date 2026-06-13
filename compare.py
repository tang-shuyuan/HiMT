#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import numpy as np
import random
import plotly.graph_objects as go
from assess import info
from assess import draw_line
from assess import pre_process_file
def random_choice_color():
    color = random.choice(["rgba(49,124,183,0.3)", "rgba(183,34,48,0.3)", "rgba(246,178,147,0.6)"])
    return color
def compare(args):
    os.makedirs(os.path.join(args.output_dir, 'assessment_out'), exist_ok=True)
    fig = go.Figure()
    args.input_file = args.query
    query_fa = pre_process_file(args.query, os.path.join(args.output_dir, 'assessment_out', "query.fa"))
    args.processed_input_file = query_fa

    query = info(args)
    query.gap = np.pi / 45
    query.int()



    reference_fa = pre_process_file(args.reference, os.path.join(args.output_dir, 'assessment_out', "reference.fa"))
    args.processed_input_file = reference_fa
    args.input_file = args.reference
    reference = info(args)
    reference.gap = np.pi/45
    reference.int()



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


    ####draw outer gene
    if args.category != "other":

        query.run_out_program()
        po_info, po_cds, ne_info, ne_cds,array, gene_list, max_copy=query.process_output_file()
        query.draw_gene(fig, po_info, po_cds, 1000, 20, 'red')
        query.draw_gene(fig, ne_info, ne_cds, 1000, -20, 'blue')

        reference.run_out_program()
        po_info, po_cds, ne_info, ne_cds, array, gene_list, max_copy =reference.process_output_file()
        reference.draw_gene(fig, po_info, po_cds, 1000, 20, 'red')
        reference.draw_gene(fig, ne_info, ne_cds, 1000, -20, 'blue')

    ###draw synteny

    synteny_out = os.path.join(args.output_dir, 'synteny_minimap2.paf')
    os.system(f"minimap2 {reference_fa} {query_fa} > {synteny_out}")

    with open(synteny_out,'r') as f:
        for line in f:
            line_content=line.split("\t")
            query_name = line_content[0]
            query_start = int(line_content[2])
            query_end = int(line_content[3])

            ref_name=line_content[5]
            ref_start = int(line_content[7])
            ref_end = int(line_content[8])

            query_start_radian, query_end_radian = query.return_scope(query_name, query_start, query_end)

            ref_start_radian, ref_end_radian = reference.return_scope(ref_name, ref_start, ref_end)

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
    fig.write_html(os.path.join(args.output_dir, "himt_compare.html"), config=config)
