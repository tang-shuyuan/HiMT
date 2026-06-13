#!/usr/bin/env python
# -*- coding: utf-8 -*-
import json
import os
import sys
import random
import argparse
import plotly.graph_objects as go
# import plotly.io as pio
import numpy as np
import base64
import copy
from pathlib import Path
from view import build_render_options, validate_args, run as render_gfa_view
def precision_to_json(fig, precision=2):
    for trace in fig.data:
        if hasattr(trace, "x") and isinstance(trace.x, (list, np.ndarray)):
            trace.x = np.round(trace.x, precision).tolist()
        if hasattr(trace, "y") and isinstance(trace.y, (list, np.ndarray)):
            trace.y = np.round(trace.y, precision).tolist()
    return fig.to_json()

def fig_to_json(fig,name="himt"):
    config_dict = {
        'toImageButtonOptions': {
            'format': 'svg',
            'filename': f'{name}',
            'scale': 1
        }
    }
    json_str = fig.to_json()
    data_dict = json.loads(json_str)
    data_dict['config'] = config_dict
    full_json = json.dumps(data_dict)
    return full_json

def return_content(file):
    with open(file,'r',encoding='utf-8') as f:
        content=f.read()
    return content

def pre_process_file(file,create_file):
    with open(file,'r') as f:
        rewrite = False
        gfa=True
        for line in f:
            if line.startswith(">"):
                gfa=False
                break
        f.seek(0)
        for line in f:
            if line.startswith(">"):
                content=line.split()
                if line.strip()!=content[0]:
                    rewrite = True
                    break

        if rewrite:
            with open(create_file,'w') as f_w:
                f.seek(0)
                for line in f:
                    if line.startswith(">"):
                        content = line.split()
                        f_w.write(content[0]+"\n")
                    else:
                        f_w.write(line)
            return create_file
        elif gfa:
            with open(create_file, 'w') as f_w:
                f.seek(0)
                for line in f:
                    if line.startswith("S"):
                        line_content = line.strip().split("\t")
                        f_w.write(">" + line_content[1] + "\n" + line_content[2] + "\n")
            return create_file
        else:
            return file

def auto_judge_category(args):
    outdb = os.path.join(args.output_dir, 'assessment_out', 'blastdb')

    config_file=os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), "config_database",
                            "Arabidopsis_chr_protein.fasta")

    blast_out=os.path.join(args.output_dir, 'assessment_out', 'blast.out')

    os.system(f"makeblastdb -in {args.processed_input_file} -dbtype nucl -out {outdb}")

    os.system(f"tblastn  -db {outdb} -query {config_file} -evalue 1e-10 -out {blast_out}  -outfmt 6")

    total_chloroplast_gene_set=set()
    match_gene_set=set()
    with open(config_file,'r') as f:
        for line in f:
            if line.startswith(">"):
                total_chloroplast_gene_set.add(line.strip()[1:])
    with open(blast_out,'r') as f:
        for line in f:
            match_gene_set.add(line.split("\t")[0])
    if len(match_gene_set)/len(total_chloroplast_gene_set)>=0.9:
        return 'chloroplast'
    else:
        return 'mitochondrial'


def draw_circos_ring(fig, start, end, r1, r2, name, width=1, color="black"):
    theta = np.linspace(start, end, 100)
    x_1 = r1 * np.cos(theta)
    y_1 = r1 * np.sin(theta)
    fig.add_trace(go.Scatter(x=x_1, y=y_1, mode='lines', marker=dict(size=0), hoverinfo='name', name=name,
                             line=dict(color=color, width=width)))
    x_2 = r2 * np.cos(theta)
    y_2 = r2 * np.sin(theta)
    fig.add_trace(go.Scatter(x=x_2, y=y_2, mode='lines', marker=dict(size=0), hoverinfo='name', name=name,
                             line=dict(color=color, width=width)))

    x1, y1 = x_1[0], y_1[0]
    x2, y2 = x_2[0], y_2[0]
    fig.add_trace(go.Scatter(x=[x1, x2], y=[y1, y2], mode='lines', line=dict(color=color, width=width)))

    x1, y1 = x_1[-1], y_1[-1]
    x2, y2 = x_2[-1], y_2[-1]
    fig.add_trace(go.Scatter(x=[x1, x2], y=[y1, y2], mode='lines', line=dict(color=color, width=width)))


def draw_line(fig, start, end, r, name, width=1,color="black"):
    theta = np.linspace(start, end, 100)
    x_1 = r * np.cos(theta)
    y_1 = r * np.sin(theta)
    fig.add_trace(go.Scatter(x=x_1, y=y_1, mode='lines',
                             marker=dict(size=0), hoverinfo='name',
                             name=name, line=dict(color=color, width=width)))

def random_choice_color():
    color=random.choice(["rgba(49,124,183,0.3)", "rgba(183,34,48,0.3)","rgba(246,178,147,0.6)"])
    return color
def intersection(region1, region2):
    return max(region1[0], region2[0]) < min(region1[1], region2[1])


def render_local_gfa_svg(input_file, output_file):
    render_args = argparse.Namespace(
        width=1600,
        height=1200,
        node_label="name-depth-length",
        edge_label="name",
        font_size=None,
        node_font_size=16,
        edge_font_size=12,
        font_color=None,
        seed=1,
    )
    validate_args(render_args)
    options = build_render_options(render_args)
    render_gfa_view(Path(input_file).expanduser().resolve(),
                    Path(output_file).expanduser().resolve(),
                    render_args.width,
                    render_args.height,
                    options,
                    render_args.seed)


def write_simple_animal_mito_report(args, process):
    os.makedirs(os.path.join(args.output_dir, "bandage_graph"), exist_ok=True)
    bandage_svg = os.path.join(args.output_dir, "bandage_graph", "mitochondrial_bandage.svg")
    render_local_gfa_svg(args.input_file, bandage_svg)
    with open(bandage_svg, "rb") as f:
        bandage_image = base64.b64encode(f.read()).decode("utf-8")

    sum_len = 0
    n50 = 0
    for contig_len in process.contig_length:
        sum_len += contig_len
        if sum_len >= process.sum_len / 2:
            n50 = contig_len
            break

    circular_contig = process.draw_sequence_connection(go.Figure()) if process.gfa else "No information"
    gc_content = round(sum(process.contig_GC_number) / process.sum_len, 4) if process.sum_len else 0
    depth_values = [x for x in process.contig_depth if x is not None]
    mean_depth = round(sum(depth_values) / len(depth_values), 2) if depth_values else "No information"

    basic_rows = [
        ("Category", "animal mitochondrial"),
        ("Input type", "gfa" if process.gfa else "fasta"),
        ("Total length", process.sum_len),
        ("Total contig number", len(process.contig_name)),
        ("Number of closed circular contigs", circular_contig),
        ("Total GC content", gc_content),
        ("N count", process.N_number),
        ("Minimum length", min(process.contig_length) if process.contig_length else 0),
        ("Maximum length", max(process.contig_length) if process.contig_length else 0),
        ("N50", n50),
        ("Mean depth", mean_depth),
    ]

    contig_headers = ["Index", "Contig name", "GC content", "Length"]
    contig_values = [
        list(range(1, len(process.contig_name) + 1)),
        process.contig_name,
        [round(process.contig_GC_number[i] / process.contig_length[i], 4) if process.contig_length[i] else 0 for i in range(len(process.contig_name))],
        process.contig_length,
    ]
    if process.gfa and len(process.contig_name) == len(process.contig_depth):
        contig_headers.append("Depth")
        contig_values.append(process.contig_depth)

    basic_table = go.Figure()
    basic_table.add_trace(go.Table(
        header=dict(values=["Animal mitochondrial basic information", "Value"],
                    fill_color='rgba(49,124,183, 0.3)', line_color='black'),
        cells=dict(values=[[row[0] for row in basic_rows], [row[1] for row in basic_rows]],
                   fill_color='rgba(255, 255,255, 0)', line_color='black')
    ))
    basic_table.update_layout(height=320, width=820, paper_bgcolor='rgba(255, 255,255,0)',
                              margin=dict(l=1, t=1, b=1), font=dict(family="Arial", color="black"))

    contig_table = go.Figure()
    contig_table.add_trace(go.Table(
        header=dict(values=contig_headers, fill_color='rgba(49,124,183, 0.3)', line_color='black'),
        cells=dict(values=contig_values, fill_color='rgba(255, 255,255, 0)', line_color='black')
    ))
    contig_table.update_layout(height=220, width=820, paper_bgcolor='rgba(255, 255,255,0)',
                               margin=dict(l=1, t=1, b=1), font=dict(family="Arial", color="black"))

    plotly_content = return_content(os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), "config_database", "Plotly.js"))
    with open(os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), "config_database", "himt_logo.svg"), 'rb') as image_file:
        logo_image = base64.b64encode(image_file.read()).decode("utf-8")

    html_content = f"""
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <title>HiMT animal mitochondrial report</title>
        <script>{plotly_content}</script>
        <style>
            body {{
                font-family: 'Segoe UI', Arial, sans-serif;
                line-height: 1.6;
                max-width: 1000px;
                margin: 0 auto;
                padding: 20px;
                background-color: #f9f9f9;
            }}
            .header {{
                width: 100vw;
                position: relative;
                left: 50%;
                transform: translateX(-50%);
                padding: 20px 0;
                background: rgb(151, 25, 23);
                margin-bottom: 20px;
                box-shadow: 0 2px 15px rgba(0,0,0,0.1);
            }}
            .header-content {{
                max-width: 1000px;
                margin: 0 auto;
                padding: 0 20px;
                justify-content: center;
                display: flex;
                align-items: center;
                gap: 20px;
            }}
            .section {{
                background: white;
                padding: 25px;
                border-radius: 8px;
                box-shadow: 0 2px 5px rgba(0,0,0,0.1);
                margin-bottom: 25px;
            }}
            .chart-container {{
                display: flex;
                justify-content: center;
                align-items: center;
                padding: 0 0 10px 0;
                border: 1px solid #e0e0e0;
                border-radius: 6px;
            }}
            .figure-caption {{
                font-size: 0.9em;
                color: #666;
                text-align: center;
                margin-top: 10px;
            }}
        </style>
    </head>
    <body>
        <div class="header">
            <div class="header-content">
                <img src="data:image/svg+xml;base64,{logo_image}" alt="Toolkit Logo" style="width: 180px; height: auto;">
                <div style="flex: 1;">
                    <h1 style="margin: 0 0 20px 0;font-size:40px;font-family: Arial,sans-serif">
                        <a href="https://github.com/tang-shuyuan/HiMT" style="text-decoration: none;color:white;" target="_blank">HiMT</a>
                    </h1>
                    <p style="margin: 0; line-height: 1.6;font-size: 20px;color:white">
                        Simplified assessment report for animal mitochondrial assemblies
                    </p>
                </div>
            </div>
        </div>

        <div class="section">
            <h2>The animal mitochondrial genome assembly graph</h2>
            <div class="chart-container">
                <img src="data:image/svg+xml;base64,{bandage_image}" alt="Assembly Graph" style="width: 760px; height: auto;">
            </div>
            <p class="figure-caption">Visualization of the input GFA assembly graph.</p>
        </div>

        <div class="section">
            <h2>Basic information</h2>
            <div class="chart-container" style="padding: 20px 0"><div id="table1-chart"></div></div>
        </div>

        <div class="section">
            <h2>Contig information</h2>
            <div class="chart-container" style="padding: 20px 0"><div id="table2-chart"></div></div>
        </div>

        <script>
            var tableLayout1 = {fig_to_json(basic_table, name="animal_mito_basic")};
            Plotly.newPlot('table1-chart', tableLayout1.data, tableLayout1.layout, tableLayout1.config);
            var tableLayout2 = {fig_to_json(contig_table, name="animal_mito_contigs")};
            Plotly.newPlot('table2-chart', tableLayout2.data, tableLayout2.layout, tableLayout2.config);
        </script>
    </body>
    </html>"""

    with open(os.path.join(args.output_dir, f"himt_mitochondrial.html"), "w", encoding="utf-8") as f:
        f.write(html_content)


class info:
    def __init__(self, args):
        self.args = copy.deepcopy(args)
        self.contig_name = []
        self.contig_length = []
        self.contig_info = []
        self.contig_depth = []
        self.contig_GC_number = []
        self.N_number = 0
        self.shift=np.pi/2
        self.display = 3 * np.pi / 2
        self.gap = np.pi / 90


    def int(self):

        self.config_file = self._get_config_file()
        self.name_prefix = self._get_name_prefix()
        self.get_info()


    def _get_config_file(self):
        if self.args.category == 'mitochondrial':
            return os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), "config_database",
                                "Arabidopsis_protein.fasta")
        elif self.args.category == 'chloroplast':
            return os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), "config_database",
                                "Arabidopsis_chr_protein.fasta")
        else:
            return None

    def _get_name_prefix(self):
        if self.args.category == 'mitochondrial':
            self.picture_height = 800
            return "mitochondrial"
        elif self.args.category == 'chloroplast':
            self.picture_height = 1800
            return "chloroplast"
        else:
            return None

    def get_info(self):
        self.gfa=True
        with open(self.args.input_file, "r") as f:
            for line in f:
                if line.startswith(">"):
                    self.gfa=False
                    break
        if not self.gfa:
            with open(self.args.processed_input_file,'r') as f:
                first_line = f.readline()
                contig = first_line.strip()[1:]
                length = 0
                sequence = ''
                GC_content = 0
                for line in f:
                    if line.startswith(">"):
                        self.contig_name.append(contig)
                        self.contig_length.append(length)
                        self.contig_info.append(sequence)
                        self.contig_GC_number.append(GC_content)
                        contig = line.strip()[1:]
                        sequence = ''
                        length = 0
                        GC_content = 0
                    else:
                        length += len(line.strip())
                        sequence += line.strip()
                        G_number = line.strip().upper().count('G')
                        C_number = line.strip().upper().count('C')
                        GC_content += G_number + C_number
                        self.N_number += line.strip().upper().count('N')
                self.contig_name.append(contig)
                self.contig_info.append(sequence)
                self.contig_length.append(length)
                self.contig_GC_number.append(GC_content)
                self.contig_depth = [0] * len(self.contig_name)

        else:
            with open(self.args.input_file,'r') as f:
                for line in f:
                    if line.startswith("S"):
                        line_content = line.strip().split("\t")
                        G_number = line_content[2].strip().upper().count("G")
                        C_number = line_content[2].strip().upper().count("C")
                        self.contig_GC_number.append(G_number + C_number)
                        self.N_number += line_content[2].strip().upper().count("N")
                        self.contig_name.append(line_content[1].strip())
                        self.contig_length.append(len(line_content[2].strip()))
                        self.contig_info.append(line_content[2].strip())
                        s_option = line_content[3:]
                        s_option_split = []
                        if len(s_option) > 0:
                            for item in s_option:
                                s_option_split.extend(item.upper().split(":"))
                            if "DP" in s_option_split:
                                self.contig_depth.append(int(float(s_option_split[s_option_split.index("DP") + 2])))
                            elif "RC" in s_option_split and "LN" in s_option_split:
                                rc = float(s_option_split[s_option_split.index("RC") + 2])
                                ln = float(s_option_split[s_option_split.index("LN") + 2])
                                self.contig_depth.append(int(round(rc / ln, 0)))
        combined_lists = list(
            zip(self.contig_length, self.contig_name, self.contig_info, self.contig_depth, self.contig_GC_number))
        sorted_combined = sorted(combined_lists, key=lambda x: x[0], reverse=True)
        self.contig_length, self.contig_name, self.contig_info, self.contig_depth, self.contig_GC_number = zip(
            *sorted_combined)
        self.contig_length, self.contig_name, self.contig_info, self.contig_depth, self.contig_GC_number = list(
            self.contig_length), \
            list(self.contig_name), list(self.contig_info), list(self.contig_depth), list(self.contig_GC_number)
        self.sum_len = sum(self.contig_length)
        while self.gap * len(self.contig_name) > self.display / 2:
            self.gap = self.gap / 2



    def run_out_program(self):
        outdb = os.path.join(self.args.output_dir, 'assessment_out', 'blastdb')

        os.system(f"makeblastdb -in {self.args.processed_input_file} -dbtype nucl -out {outdb}")

        os.system(f"tblastn  -db {outdb} -query {self.config_file} \
                    -evalue 1e-10 -out {os.path.join(self.args.output_dir, 'assessment_out', f'{self.name_prefix}_blast.out')} \
                    -outfmt 6")

        os.system(f"miniprot {self.args.processed_input_file} {self.config_file} \
        -L 0 -p 0 --outs=0 --outc=0 -u --gff >{os.path.join(self.args.output_dir, f'{self.name_prefix}_raw_miniprot.gff3')} ")


        ####remove same element
        uniq_element=set()
        with open(os.path.join(self.args.output_dir, f'{self.name_prefix}_raw_miniprot.gff3'),'r') as f,\
                open(os.path.join(self.args.output_dir, f'{self.name_prefix}_miniprot.gff3'),'w') as f_w:
            write_content = True
            for line in f:
                if line.startswith("##") :
                    line_content = line.split()[:12]
                    line_content = "&".join(map(str,line_content))
                    if line_content not in uniq_element:
                        uniq_element.add(line_content)
                        write_content=True
                    else:
                        write_content=False

                if write_content:
                    f_w.write(line)
        os.remove(os.path.join(self.args.output_dir, f'{self.name_prefix}_raw_miniprot.gff3'))

    def process_output_file(self):
        ##get gene completeness
        with open(os.path.join(self.args.output_dir, f'{self.name_prefix}_miniprot.gff3'), 'r') as f:
            dic = {}
            query_length = {}
            self.conserved_gene_number = {key: 0 for key in self.contig_name}

            for line in f:
                if line.startswith('##PAF'):
                    line_content = line.strip().split("\t")
                    conserved_gene = line_content[1]
                    alignment_length = int(line_content[4]) - int(line_content[3])
                    gene_length = int(line_content[2])
                    completeness = int(round(alignment_length / gene_length * 100, 0))
                    number_contig_name = line_content[6]
                    if conserved_gene not in query_length:
                        query_length[conserved_gene] = gene_length
                    if line_content[5] != '*':
                        if conserved_gene in dic:
                            dic[conserved_gene].append(completeness)
                        else:
                            dic[conserved_gene] = [completeness]
                        self.conserved_gene_number[number_contig_name] += 1
            ##extract gene info
        positive_gene_info = {}
        positive_cds_info = {}
        negative_gene_info = {}
        negative_cds_info = {}
        with open(os.path.join(self.args.output_dir, f'{self.name_prefix}_miniprot.gff3'), 'r') as f:
            for line in f:
                if not line.startswith("#"):
                    line_content = line.strip().split('\t')
                    gene_start = int(line_content[3])
                    gene_end = int(line_content[4])
                    contig = line_content[0]
                    gene_name = line_content[-1].split(';')[-1].split(' ')[0].split("=")[1]
                    parent = line_content[8].split(";")[0].split("=")[1]
                    strand = line_content[6]

                    radian_start = self.return_radian(contig, gene_start)
                    radian_end = self.return_radian(contig, gene_end)

                    if line_content[2].upper() == "MRNA":
                        if strand == "+":
                            positive_gene_info[parent + "&" + gene_name] = (radian_start, radian_end)
                        else:
                            negative_gene_info[parent + "&" + gene_name] = (radian_start, radian_end)
                    elif line_content[2].upper() == "CDS":
                        if strand == "+":
                            positive_cds_info[(radian_start, radian_end, parent)] = parent
                        else:
                            negative_cds_info[(radian_start, radian_end, parent)] = parent

        with open(os.path.join(self.args.output_dir, 'assessment_out', f'{self.name_prefix}_blast.out'), 'r') as f:
            for line in f:
                line_content = line.strip().split("\t")
                conserved_gene = line_content[0]
                number_contig_name = line_content[1]
                if conserved_gene not in dic:
                    print(f'add {conserved_gene} via blast')
                    alignment_length = int(line_content[7]) - int(line_content[6])
                    completeness = int(round(alignment_length / query_length[conserved_gene] * 100, 0))
                    dic[conserved_gene] = [completeness]
                    self.conserved_gene_number[number_contig_name] += 1
                    if int(line_content[9]) > int(line_content[8]):

                        radian_start, radian_end = self.return_scope(number_contig_name, int(line_content[8]),
                                                                     int(line_content[9]))
                        positive_gene_info[conserved_gene + "&" + conserved_gene] =(radian_start, radian_end)

                        positive_cds_info[(radian_start, radian_end, conserved_gene)] = conserved_gene
                        # self.draw_exon(fig, radian_start, radian_end, 900, 20, conserved_gene, 'red')
                    elif int(line_content[9]) < int(line_content[8]):
                        radian_start, radian_end = self.return_scope(number_contig_name, int(line_content[9]),
                                                                     int(line_content[8]))
                        negative_gene_info[conserved_gene + "&" + conserved_gene] = (radian_start, radian_end)
                        negative_cds_info[(radian_start, radian_end, conserved_gene)] = conserved_gene
                        # self.draw_exon(fig, radian_start, radian_end, 900, -20, conserved_gene, 'blue')
        array = []
        gene_list = []
        max_copy = 0
        for gene, compeleteness_list in dic.items():
            gene_list.append(gene)
            max_copy = max(max_copy, len(compeleteness_list))
            array.append(sorted(compeleteness_list, reverse=True))
        for i, li in enumerate(array):
            if len(li) != max_copy:
                li.extend([0] * (max_copy - len(li)))
                array[i] = li
        with open(self.config_file, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    if line.strip()[1:] not in gene_list:
                        gene_list.append(line.strip()[1:])
                        array.append([0] * max_copy)
        return positive_gene_info, positive_cds_info, negative_gene_info, negative_cds_info,array, gene_list, max_copy

    # def draw_heatmap(self):
    def draw_exon(self, fig, radian_start, radian_end, r, strech_r, gene_name, color):
        line_number = int(abs(radian_end - radian_start) * 1500 / np.pi + 20)
        theta = np.linspace(radian_start, radian_end, line_number)
        x_1 = r * np.cos(theta)
        y_1 = r * np.sin(theta)
        x_2 = (r + strech_r) * np.cos(theta)
        y_2 = (r + strech_r) * np.sin(theta)

        scatter_data = [go.Scatter(x=[x1, x2], y=[y1, y2], mode='lines', name=gene_name, hoverinfo='name',
                                   line=dict(color=color, width=1))
                        for x1, y1, x2, y2 in zip(x_1, y_1, x_2, y_2)]

        fig.add_traces(scatter_data)

    def draw_gene(self, fig, gene_info, cds_info, r, strech_r, color):
        display = self.display - ((len(self.contig_name) - 1) * self.gap)
        last_muti_region = (0, 0)
        single_exon = set()
        muti_exon = set()
        gene_annotation = []
        origin_r = r
        for value in cds_info.values():
            if value in single_exon:
                muti_exon.add(value)
            else:
                single_exon.add(value)
        single_exon = single_exon - muti_exon
        gene_info = {key: value for key, value in sorted(gene_info.items(), key=lambda x: x[1])}

        for key, value in gene_info.items():
            parent = key.split("&")[0]
            gene_name = key.split("&")[1]
            radian_start = value[0]
            radian_end = value[1]

            gene_center = (radian_start + radian_end) / 2

            if parent in single_exon:
                self.draw_exon(fig, radian_start, radian_end, origin_r, strech_r, gene_name, color)

                x1 = (origin_r + strech_r) * np.cos(gene_center)
                x2 = (origin_r + strech_r * 3 / 2) * np.cos(gene_center)
                y1 = (origin_r + strech_r) * np.sin(gene_center)
                y2 = (origin_r + strech_r * 3 / 2) * np.sin(gene_center)

                fig.add_trace(go.Scatter(x=[x1, x2], y=[y1, y2], mode='lines', name=gene_name, hoverinfo='name',
                                         line=dict(color="black", width=1)))

                angel_deg = gene_center * (180 / np.pi)
                gene_annotation.append([gene_center, gene_name, origin_r + strech_r * 4, -angel_deg])

            else:
                all_region = sorted(set(value for value in gene_info.values()))
                all_region.append((0, 0))
                index = all_region.index(value)
                if intersection(all_region[index - 1], value) or intersection(value, all_region[index + 1]):
                    if intersection(last_muti_region, value):
                        r = r + strech_r + strech_r / 2
                    else:
                        r = origin_r + strech_r + strech_r / 2
                else:
                    r = origin_r
                xian = []
                for re, pa in cds_info.items():
                    if pa == parent:
                        xian.append((re[0] + re[1]) / 2)
                        self.draw_exon(fig, re[0], re[1], r, strech_r, gene_name, color)
                xian = sorted(xian)
                for i in range(1, len(xian)):
                    draw_line(fig, xian[i - 1], xian[i], r + strech_r / 2, gene_name)
                last_muti_region = value

                x1 = (r + strech_r / 2) * np.cos(gene_center)
                x2 = (r + strech_r) * np.cos(gene_center)
                y1 = (r + strech_r / 2) * np.sin(gene_center)
                y2 = (r + strech_r) * np.sin(gene_center)
                fig.add_trace(go.Scatter(x=[x1, x2], y=[y1, y2], mode='lines', name=gene_name, hoverinfo='name',
                                         line=dict(color="black", width=1)))

                angel_deg = gene_center * (180 / np.pi)
                gene_annotation.append([gene_center, gene_name, r + strech_r * 4, -angel_deg])

        gene_annotation = sorted(gene_annotation, key=lambda x: x[0])
        temp_position = -10
        for index in range(len(gene_annotation)):
            content = gene_annotation[index]
            gene_center = content[0]
            gene_name = content[1]
            x2 = content[2] * np.cos(gene_center)
            y2 = content[2] * np.sin(gene_center)
            spin_angel = content[3]
            gap_angel=np.pi/120
            if gene_center - temp_position > gap_angel:
                fig.add_annotation(x=x2, y=y2, text=gene_name, showarrow=False, font=dict(size=12, color='black',family="Arial"),
                                   textangle=spin_angel)
                temp_position = gene_center
            else:
                x2 = content[2] * np.cos(temp_position + gap_angel)
                y2 = content[2] * np.sin(temp_position + gap_angel)
                fig.add_annotation(x=x2, y=y2, text=gene_name, showarrow=False, font=dict(size=12, color='black',family="Arial"),
                                   textangle=spin_angel)
                temp_position = temp_position + gap_angel
        for index, name in enumerate(self.contig_name):
            start = display * (
                sum(self.contig_length[:index])) / self.sum_len + index * self.gap + self.shift
            end = display * (
                sum(self.contig_length[:index + 1])) / self.sum_len + index * self.gap + self.shift
            draw_line(fig, start, end, origin_r, '', 1)


    ###draw chromosome
    def draw_chromosome(self,fig,r1,r2,same_color=None):
        display = self.display - ((len(self.contig_name) - 1) * self.gap)
        for index, name in enumerate(self.contig_name):
            start = display * (
                sum(self.contig_length[:index])) / self.sum_len + index * self.gap + self.shift
            end = display * (
                sum(self.contig_length[:index + 1])) / self.sum_len + index * self.gap + self.shift
            if end - start < 0.04 * len(name):
                txet = ''
            else:
                txet = name

            theta = np.linspace(start, end, 100)
            x = (r1+r2)/2 * np.cos(theta)
            y = (r1+r2)/2 * np.sin(theta)
            h = 1e-6
            dx_dtheta_approx = (x[51] - x[50]) / h
            dy_dtheta_approx = (y[51] - y[50]) / h
            slope = dy_dtheta_approx / dx_dtheta_approx
            angel_rad = np.arctan(slope)
            angel_deg = angel_rad * (180 / np.pi)
            uniform_numbers = np.linspace(r1, r2, 50)
            if same_color:
                color=same_color
            else:
                color = random_choice_color()
            for item in uniform_numbers:
                draw_line(fig, start, end, item, name, width=2, color=color)


            fig.add_annotation(x=x[50], y=y[50], text=txet, showarrow=False,
                                font=dict(size=12, color='black'), textangle=-angel_deg)


    ###draw GC content
    def draw_GC_content(self,fig,r1,r2):
        display = self.display - ((len(self.contig_name) - 1) * self.gap)
        for index, name in enumerate(self.contig_name):
            start = display * (
                sum(self.contig_length[:index])) / self.sum_len + index * self.gap + self.shift
            end = display * (
                sum(self.contig_length[:index + 1])) / self.sum_len + index * self.gap + self.shift
            sequence = self.contig_info[index]
            step = int(self.sum_len / 1000)
            mid = int(step / 2)
            for i in range(mid, len(sequence), step):
                G_content = sequence[i - mid:i + mid].upper().count("G")
                C_content = sequence[i - mid:i + mid].upper().count("C")
                GC_content = round((G_content + C_content) / step, 2)
                radian = display * (
                        sum(self.contig_length[:index]) + i) / self.sum_len + index * self.gap + self.shift

                x = np.cos(radian) * r1
                y = np.sin(radian) * r1
                x1 = np.cos(radian) * (r1 + (r2-r1) * GC_content)
                y1 = np.sin(radian) * (r1+ (r2-r1) * GC_content)
                fig.add_trace(go.Scatter(x=[x, x1], y=[y, y1], mode='lines', name=GC_content, hoverinfo='name',
                                         line=dict(color="rgba(151, 25, 23,0.8)", width=1)))
            draw_circos_ring(fig, start, end, r1, r2, name)


    ###draw contig depth
    def draw_contig_depth(self,fig,r1,r2):
        display = self.display - ((len(self.contig_name) - 1) * self.gap)
        color = "rgba(49,124,183, 0.6)"
        for index, name in enumerate(self.contig_name):
            start = display * (
                sum(self.contig_length[:index])) / self.sum_len + index * self.gap + self.shift
            end = display * (
                sum(self.contig_length[:index + 1])) / self.sum_len + index * self.gap + self.shift
            if self.gfa:
                max_depth = max(self.contig_depth)
                depth = self.contig_depth[index]
                length=(r2-r1)*depth/max_depth
                uniform_numbers = np.linspace(r1, r1+length, 100)

                for item in uniform_numbers:
                    draw_line(fig, start, end, item, depth, width=1, color=color)

                # theta = np.linspace(start, end, int((end - start) / display * 1500 + 5))
                # x = np.cos(theta) * r1
                # y = np.sin(theta) * r1
                # x1 = np.cos(theta) * (r1 + (r2-r1) * depth / max_depth)
                # y1 = np.sin(theta) * (r1 + (r2-r1) * depth / max_depth)
                # scatter_data = [go.Scatter(x=[x1, x2], y=[y1, y2], mode='lines', name=depth, hoverinfo='name',
                #                            line=dict(color='rgba(128,128,128, 0.5)', width=2))
                #                 for x1, y1, x2, y2 in zip(x, y, x1, y1)]
                #
                # fig.add_traces(scatter_data)
            draw_circos_ring(fig, start, end, r1, r2, name=name)

    def return_radian(self, contig, coordinate):
        display = self.display - ((len(self.contig_name) - 1) * self.gap)
        index = self.contig_name.index(contig)
        radian = display * (sum(self.contig_length[:index])+coordinate)/self.sum_len +index*self.gap+ self.shift
        return radian

    def return_scope(self, contig, start, end):
        start_radian = self.return_radian(contig, start)
        end_radian = self.return_radian(contig, end)
        return start_radian, end_radian

    def add_annotation(self, fig):
        fig.add_annotation(x=0, y=460, showarrow=False, xanchor='left', yanchor='middle',
                           font=dict(size=16, color='black',family="Arial"),
                           text='① Sequence ID')

        fig.add_annotation(x=0, y=585, showarrow=False, xanchor='left', yanchor='middle',
                           font=dict(size=16, color='black',family="Arial"),
                           text='② GC content')
        fig.add_annotation(x=0, y=725, showarrow=False, xanchor='left', yanchor='middle',
                           font=dict(size=16, color='black',family="Arial"),
                           text='③ Depth')
        fig.add_annotation(x=0, y=1000, showarrow=False, xanchor='left', yanchor='middle',
                           font=dict(size=16, color='black',family="Arial"),
                           text='④ Conserved genes')
        # fig.add_annotation(x=400, y=300, showarrow=False, xanchor='left', yanchor='middle',
        #                    font=dict(size=12, color='black'),
        #                    text='''total length: 100<br> total GC_content: 0.45<br>contig number: 3''',
        #                    row=2, col=1)
    def draw_sequence_connection(self,fig,color="gray"):
        circular_contig=set()
        with open(self.args.input_file,'r') as f:
            for line in f:
                if not line.startswith("#"):
                    if  line.startswith("L"):
                        line_content=line.strip().split("\t")
                        contig1=line_content[1]
                        contig1_length=self.contig_length[self.contig_name.index(contig1)]
                        contig2=line_content[3]
                        contig2_length=self.contig_length[self.contig_name.index(contig2)]
                        if contig1==contig2:
                            circular_contig.add(contig2)
                            if line_content[2]=="+":
                                position1_radian=self.return_radian(contig1,contig1_length)
                            else:
                                position1_radian=self.return_radian(contig1,0)

                            if line_content[4]=="+":
                                position2_radian=self.return_radian(contig2,0)
                            else:
                                position2_radian=self.return_radian(contig2,contig2_length)
                            x1 = np.cos(position1_radian) * 400
                            y1 = np.sin(position1_radian) * 400
                            x2 = np.cos(position2_radian) * 400
                            y2 = np.sin(position2_radian) * 400
                            cx, cy = 0, 0
                            t_values = np.linspace(0, 1, 50)
                            x = (1 - t_values) ** 2 * x1 + 2 * (1 - t_values) * t_values * cx + t_values ** 2 * x2
                            y = (1 - t_values) ** 2 * y1 + 2 * (1 - t_values) * t_values * cy + t_values ** 2 * y2
                            fig.add_trace(
                                go.Scatter(
                                    x=x,
                                    y=y,
                                    mode='lines',
                                    name=contig1+line_content[1]+contig2+line_content[4],
                                    hoverinfo='name',
                                    line=dict(color=color, width=2)
                                )
                            )
        return len(circular_contig)

###heatmap
def assess(args):
    os.makedirs(os.path.join(args.output_dir, 'assessment_out'), exist_ok=True)

    args.processed_input_file=pre_process_file(args.input_file,os.path.join(args.output_dir, 'assessment_out','assess.fa'))

    if not args.category:
        args.category = auto_judge_category(args)
    process = info(args)
    process.int()
    if getattr(args, "species", "plant") == "animal" :
        write_simple_animal_mito_report(args, process)
        return
    process.run_out_program()
    fig1 = go.Figure()
    fig2 = go.Figure()
    po_info, po_cds, ne_info, ne_cds,array, gene_list, max_copy = process.process_output_file()

    ###draw heatmap

    colorscale = [[0, 'rgb(255,255,255)'], [1, 'rgb(151, 25, 23)']]
    fig1.add_trace(go.Heatmap(
        z=array,
        x=[f'copy{i + 1}' for i in range(max_copy)],
        y=gene_list,
        colorscale=colorscale,
        coloraxis="coloraxis1",
        zsmooth=None))

    for i in range(len(array)):
        for j in range(len(array[i])):
            val = array[i][j]
            if val > 0:
                fig1.add_annotation(
                    x=f'copy{j + 1}',
                    y=gene_list[i],
                    text=val,
                    showarrow=False,
                    xref='x',
                    yref='y',
                    font=dict(color='black',family="Arial"),
                )
    fig1.update_layout(
        width=550,
        height=process.picture_height,
        showlegend=False,
        plot_bgcolor='rgba(255, 255,255,0)',
        paper_bgcolor='rgba(255, 255,255,0)',
        font=dict(family="Arial",color="black"),
        xaxis=dict(title_font=dict(size=16, color='black',family="Arial"),
                   showgrid=False,
                   zeroline=False,
                   ),
        yaxis=dict(title_font=dict(size=16, color='black',family="Arial"), title_text='Conserved protein coding genes',
                   showgrid=False,
                   zeroline=False),
        coloraxis1=dict(colorscale=colorscale,
                        colorbar=dict(
                            len=0.25,
                            x=1.05,
                            y=0.8,
                            title=dict(text="Gene integrity", side="right", font=dict(family="Arial",size=16, color='black')),
                        )),
        shapes=[
            dict(
                type="rect",
                xref="x", yref="y",
                x0=-0.5, y0=-0.5, x1=max_copy - 0.5, y1=len(gene_list) - 0.5,
                line=dict(width=1, color="black")
            )
        ]
    )




    ##draw circos
    process.draw_chromosome(fig2,420,480)
    process.draw_GC_content(fig2,550,620)
    process.draw_contig_depth(fig2,690,760)

    #### draw gene
    process.draw_gene(fig2, po_info, po_cds, 1000, 25, 'red')
    process.draw_gene(fig2, ne_info, ne_cds, 1000, -25, 'rgb(13,168,236)')
    ###add annotation information
    process.add_annotation(fig2)


    ####draw contig connection
    if process.gfa:
        circular_contig=process.draw_sequence_connection(fig2)
    else:
        circular_contig="No information "


    fig2.update_layout(
        width=800,
        height=800,
        showlegend=False,
        plot_bgcolor='rgba(255, 255,255,0)',
        paper_bgcolor='rgba(255, 255,255,0)',
        font=dict(family="Arial"),
        xaxis=dict(title_text=None, range=[-1150, 1150],
                   showgrid=False, scaleanchor='y',
                   showticklabels=False,
                   zeroline=False, ),
        yaxis=dict(title_text=None, range=[-1150, 1150],
                   showgrid=False,
                   showticklabels=False,
                   zeroline=False, )
    )

    precision_to_json(fig2)
    ###draw table 1
    sum_len = 0
    for N50 in process.contig_length:
        sum_len += N50
        if sum_len >= process.sum_len / 2:
            break
    information = ["Total length", "Total contig number", "N count", "Number of closed circular contigs","Total GC content", "Minimum length",
                   "Maximum length", "N50"]
    value = [process.sum_len, len(process.contig_name), process.N_number,circular_contig,
             round(sum(process.contig_GC_number) / process.sum_len, 2),
             min(process.contig_length), max(process.contig_length), N50]
    table1 = go.Figure()
    table1.add_trace(go.Table(header=dict(values=[f"{process.name_prefix.title()} basic information", "Value"],
                                          fill_color='rgba(49,124,183, 0.3)', line_color='black'),
                              cells=dict(values=[information, value], fill_color='rgba(255, 255,255, 0)',
                                         line_color='black')))
    table1.update_layout(height=200,width=800, paper_bgcolor='rgba(255, 255,255,0)', margin=dict(l=1, t=1, b=1),
                         font=dict(family="Arial", color="black"))

    # table_html = pio.to_html(table1, full_html=False)
    ###draw table 2
    contig_GC_content = [round(process.contig_GC_number[i] / process.contig_length[i], 2) for i in
                         range(len(process.contig_name))]
    conserved_gene_number_list = []

    for name in process.contig_name:
        conserved_gene_number_list.append(process.conserved_gene_number[name])
    if process.gfa and len(process.contig_name) == len(process.contig_depth):
        header = ["Index", "Contig name", "Conserved gene number", "GC content", "Length", "Depth"]
        information = [list(range(1, len(process.contig_name) + 1)), process.contig_name, conserved_gene_number_list,
                       contig_GC_content, process.contig_length, process.contig_depth]
    else:
        header = ["Index", "Contig name", "Conserved gene number", "GC content", "Length"]
        information = [list(range(1, len(process.contig_name) + 1)), process.contig_name, conserved_gene_number_list,
                       contig_GC_content, process.contig_length]
    table2 = go.Figure()
    table2.add_trace(go.Table(header=dict(values=header, fill_color='rgba(49,124,183, 0.3)', line_color='black')
                              , cells=dict(values=information, fill_color='rgba(255, 255,255, 0)', line_color='black')))
    table2.update_layout(height=200,width=800, paper_bgcolor='rgba(255, 255,255,0)', margin=dict(l=1, t=1, b=1),
                         font=dict(family="Arial", color="black"))

    plotly_content=return_content(os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), "config_database",
                                "Plotly.js"))
    with open(os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), "config_database",
                                "himt_logo.svg"),'rb') as image_file:
        logo_image = base64.b64encode(image_file.read()).decode("utf-8")

    ### Obtain bandage figure
    os.makedirs(os.path.join(args.output_dir, "bandage_graph"), exist_ok=True)
    bandage_svg = os.path.join(args.output_dir, 'bandage_graph', f'{args.category}_bandage.svg')
    render_local_gfa_svg(args.input_file, bandage_svg)


    with open(bandage_svg,'rb') as f:
        bandage_image=base64.b64encode(f.read()).decode("utf-8")

    html_content = f"""
        <!DOCTYPE html>
        <html lang="en">
        <head>
            <meta charset="UTF-8">
            <title>HiMT report</title>
            <!-- 引入 Plotly.js -->
            <script >
        {plotly_content}
        </script>
            <style>
                 body {{
                font-family: 'Segoe UI', Arial, sans-serif;
                line-height: 1.6;
                 max-width: 1000px;
                margin: 0 auto;
                padding: 20px;
                background-color: #f9f9f9;}}
         .header {{
                width: 100vw;          
                position: relative;    
                left: 50%;             
                transform: translateX(-50%);
                padding: 20px 0;     
                background: rgb(151, 25, 23);
                margin-bottom: 20px;
                box-shadow: 0 2px 15px rgba(0,0,0,0.1);
            }}  
        .header-content {{
               max-width: 1000px;    
                margin: 0 auto;       
                padding: 0 20px;   
                justify-content: center;   
                display: flex;
                align-items: center;
                gap: 20px;
           }} 


                .section {{
                background: white;
                padding: 25px;
                border-radius: 8px;
                box-shadow: 0 2px 5px rgba(0,0,0,0.1);
                margin-bottom: 25px;
            }}
        .chart-container {{
               display: flex;
                justify-content: center;
                align-items: center;
                padding: 0px 0 10px 0;
                border: 1px solid #e0e0e0;  /* 添加边框 */
                border-radius: 6px;
            }}
             .figure-caption {{
                font-size: 0.9em;
                color: #666;
                text-align: center;
                margin-top: 10px;
            }}

            </style>
        </head>

        <body>
        <div class="header">

        <div class="header-content">
        <img src="data:image/svg+xml;base64,{logo_image}" alt="Toolkit Logo" style="width: 180px; height: auto;">
        <div style="flex: 1;">
            <h1 style="margin: 0 0 20px 0;font-size:40px;font-family: Arial,sans-serif">
            <a href="https://github.com/tang-shuyuan/HiMT" style="text-decoration: none;color:white;" target="_blank">HiMT</a> </h1>
            <p style="margin: 0; line-height: 1.6;font-size: 20px;color:white">
                An Integrative Toolkit for Assembling Organelle Genomes Using HiFi Reads</p></div>
                </div>
        </div>

        <div class="section">
            <h2>The {args.category} genome assembly graph</h2>
            <div class="chart-container">
            <img src="data:image/svg+xml;base64,{bandage_image}" alt="Bandage Image" style="margin-left: 20px;
            width: 600px; height: auto;"></div>
            <p class="figure-caption">Visualization of {args.category} structure .</p>
            
        </div>

         <div class="section">
                 <h2>The completeness of core genes</h2>
               <div class="chart-container"> <div id="figure1"></div></div>
                <p class="figure-caption">With conserved genes from <i>Arabidopsis thaliana</i> {process.name_prefix} genome as a reference, 
                 assess whether the assembly result encompasses all core genes and the integrity of each one.</p>
         </div>

         <div class="section">
                <h2>Basic information visualization </h2>
               <div class="chart-container"> <div id="figure2"></div></div>
                <p class="figure-caption">
                 The inner lines represent the circular sequence.
                 Starting from the innermost part and moving outwards, the elements represent the sequence ID, 
                 GC content, sequencing depth, and conserved genes respectively.</p>
         </div>

         <div class="section">
             <h2 >Basic information of {process.name_prefix} genome </h2>
                <div class="chart-container" style="padding: 20px 0"> <div id="table1-chart" ></div></div>
         </div>

         <div class="section">       
             <h2>Contig information </h2>
                 <div class="chart-container"  style="padding: 20px 0"><div id="table2-chart"> </div></div>
         </div>   

            <script>

                var layout1 = {fig_to_json(fig1)};
                Plotly.newPlot('figure1', layout1.data, layout1.layout,layout1.config);
                var layout2 = {fig_to_json(fig2)};
                Plotly.newPlot('figure2', layout2.data, layout2.layout,layout2.config);
                //table1
                var tableLayout = {fig_to_json(table1)};
                Plotly.newPlot('table1-chart', tableLayout.data, tableLayout.layout,tableLayout.config);
                //table2
                 var tableLayout = {fig_to_json(table2)};
                Plotly.newPlot('table2-chart', tableLayout.data, tableLayout.layout,tableLayout.config);      
            </script>
        </body>
        </html>"""
    with open(os.path.join(process.args.output_dir, f"himt_{args.category}.html"), 'w') as f:
        f.write(html_content)
        if hasattr(args, "table3_value"):
            information = [["Num_seqs", "Sum_len", "Min_len", "Avg_len", "Max_len",
                            f"Estimated {process.name_prefix} genome depth",
                            "Num_seqs(extracted)", "Sum_len(extracted)", "Min_len(extracted)", "Avg_len(extracted)",
                            "Max_len(extracted)",
                            f"HiMT output {process.name_prefix} genome depth"], args.table3_value]

            table3 = go.Figure(data=[
                go.Table(header=dict(values=["Data statistics information", "value"],
                                     fill_color='rgba(49,124,183, 0.3)', line_color='black'),
                         cells=dict(values=information, fill_color='rgba(255, 255,255, 0)', line_color='black'))])
            table3.update_layout(height=300,width=800, paper_bgcolor='rgba(255, 255,255,0)', margin=dict(l=1, t=1, b=1),
                                 font=dict(family="Arial", color="black"))
            add = f"""
            <html>
            <body>
            <div class="section">   
            <h2>Data information </h2>
             <div class="chart-container"  style="padding: 20px 0"><div id="table3-chart"></div></div>
            </div>
            <script>
            var tableLayout = {fig_to_json(table3)};
            Plotly.newPlot('table3-chart', tableLayout.data, tableLayout.layout,tableLayout.config);      
            </script>
            </body>
            </html>"""
            f.write(add)

    # fig.write_html(os.path.join(process.args.output_dir, f"himt_{args.category}.html"))
