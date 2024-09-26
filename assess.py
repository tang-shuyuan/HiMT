import os
import sys
import pandas as pd
import base64
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
def assese(args):
    config_path = os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), "config_database")
    basic_information_list=["Total lenght","Total contig number","N_count","Total GC content(%)",\
                            "Minimum length","Maximum lenght","N50"]
    contig_name=[]
    contig_GC_number=[]
    contig_length=[]
    with open(args.input_file,"r") as f:
        first_line=f.readline()
        N_count = 0
        if first_line.startswith(">"):
            contig = first_line.strip()[1:]
            GC_number = 0
            length = 0

            for line in f:
                if line.startswith(">"):
                    contig_name.append(contig)
                    contig_GC_number.append(GC_number)
                    contig_length.append(length)
                    contig=line.strip()[1:]
                    GC_number=0
                    length=0
                else:
                    length+=len(line.strip())
                    for item in line.strip():
                        if item.upper() =="G" or item.upper()=="C":
                            GC_number+=1
                        elif item.upper()=="N":
                            N_count+=1
            contig_name.append(contig)
            contig_GC_number.append(GC_number)
            contig_length.append(length)

            os.system(f"miniprot {args.input_file} {os.path.join(config_path,'Arabidopsis_protein.fasta')} \
-p 0.05 --outs=0.05 --outc=0.01 -u >{os.path.join(args.output_dir,'miniprot.paf')} ")

            os.system(f"miniprot {args.input_file} {os.path.join(config_path, 'Arabidopsis_protein.fasta')} \
            -p 0.05 --outs=0.05 --outc=0.01 --gff >{os.path.join(args.output_dir, 'miniprot.gff')} ")



            os.makedirs(os.path.join(args.output_dir,'assesement_out'), exist_ok=True)
            outdb=os.path.join(args.output_dir,'assesement_out','blastdb')
            os.system(f"makeblastdb -in {args.input_file} -dbtype nucl -out {outdb}")
            os.system(f"tblastn  -db {outdb} -query {os.path.join(config_path,'Arabidopsis_protein.fasta')} \
            -evalue 1e-10 -out {os.path.join(args.output_dir,'assesement_out','blast.out')} -outfmt 6")


        else:
            f.seek(0)
            contig_depth=[]
            os.makedirs(os.path.join(args.output_dir, 'assesement_out'), exist_ok=True)
            with open(os.path.join(args.output_dir,'assesement_out',"assese.fa"),"w") as f_w:
                for line in f:
                    if line.startswith("S"):
                        line_content = line.strip().split("\t")
                        f_w.write(">" + line_content[1] + "\n"+line_content[2] + "\n")
                        GC_number=0
                        for base in line_content[2].strip():
                            if base.upper()=="G" or base.upper()=="C":
                                GC_number+=1
                            elif base.upper()=="N":
                                N_count+=1
                        contig_name.append(line_content[1].strip())
                        contig_GC_number.append(GC_number)
                        contig_length.append(len(line_content[2].strip()))

                        s_option=line_content[3:]
                        s_option_split=[]
                        if len(s_option)>0:
                            for item in s_option:
                                s_option_split.extend(item.upper().split(":"))

                            if "DP" in s_option_split:
                                contig_depth.append(s_option_split[s_option_split.index("DP")+2])
                            elif "RC" in s_option_split and "LN" in s_option_split :
                                rc=float(s_option_split[s_option_split.index("RC")+2])
                                ln=float(s_option_split[s_option_split.index("LN")+2])
                                contig_depth.append(int(round(rc/ln,0)))

            os.system(f"miniprot {os.path.join(args.output_dir,'assesement_out','assese.fa')} \
{os.path.join(config_path,'Arabidopsis_protein.fasta')} \
-p 0.05 --outs=0.05 --outc=0.01 -u >{os.path.join(args.output_dir,'miniprot.paf')} ")
            os.system(f"miniprot {os.path.join(args.output_dir, 'assesement_out', 'assese.fa')} \
            {os.path.join(config_path, 'Arabidopsis_protein.fasta')} \
            -p 0.05 --outs=0.05 --outc=0.01 --gff > {os.path.join(args.output_dir, 'miniprot.gff')} ")

            outdb = os.path.join(args.output_dir, 'assesement_out', 'blastdb')
            os.system(f"makeblastdb -in {os.path.join(args.output_dir,'assesement_out','assese.fa')} -dbtype nucl -out {outdb}")
            os.system(f"tblastn  -db {outdb} -query {os.path.join(config_path, 'Arabidopsis_protein.fasta')} \
                       -evalue 1e-10 -out {os.path.join(args.output_dir, 'assesement_out','blast.out')} -outfmt 6")

    with open(os.path.join(args.output_dir,'miniprot.paf'),'r') as f:
        dic={}
        query_length={}
        for i, line in enumerate(f):
            line_content = line.strip().split("\t")
            conserved_gene = line_content[0]
            alignment_length = int(line_content[3]) - int(line_content[2])
            gene_length = int(line_content[1])
            completeness = int(round(alignment_length / gene_length * 100, 0))
            if conserved_gene not in query_length:
                query_length[conserved_gene]=gene_length
            if line_content[4] !='*':
                if conserved_gene in dic:
                    dic[conserved_gene].append(completeness)
                else:
                    dic[conserved_gene] = [completeness]
    with open(os.path.join(args.output_dir, 'assesement_out','blast.out'),'r') as f:
        for line in f:
            line_content=line.strip().split("\t")
            conserved_gene=line_content[0]
            if conserved_gene not in dic:
                print(f'add the {conserved_gene} via blast')
                alignment_length=int(line_content[7])-int(line_content[6])
                completeness=int(round(alignment_length/query_length[conserved_gene]*100,0))
                dic[conserved_gene]=[completeness]

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
    with open(os.path.join(config_path,'Arabidopsis_protein.fasta'),'r') as f:
        for line in f:
            if line.startswith('>'):
                if line.strip()[1:] not in gene_list:
                    gene_list.append(line.strip()[1:])
                    array.append([0]*max_copy)

###heatmap
    colors = [(1, 1, 1), (0, 0.5, 0)]
    n_bins = 100  # 颜色渐变的步数
    cm = mcolors.LinearSegmentedColormap.from_list('custom_cmap', colors, N=n_bins)

    fig, ax = plt.subplots(figsize=(8, 10))
    cax = ax.matshow(array, cmap=cm)
    for i in range(len(array)):
        for j in range(len(array[i])):
            val = array[i][j]
            if val > 0:
                ax.text(j, i, f'{val}', ha='center', va='center', color='black')

    # 添加颜色渐变条并设置标签
    cbar = fig.colorbar(cax)
    cbar.set_label('Gene Completeness')

    ax.set_xticks(range(max_copy))
    ax.set_yticks(range(len(gene_list)))
    ax.set_xticklabels([f'copy{i + 1}' for i in range(max_copy)])
    ax.set_yticklabels(gene_list)
    ax.set_aspect(aspect='auto')
    # plt.xticks(rotation=90)
    # plt.title('fig.1 Statistics of completeness of mitochondrial genome protein coding genes')
    plt.xlabel('gene copy')
    plt.ylabel('conserved protein coding genes')
    plt.savefig(os.path.join(args.output_dir, "hetmap.png"), dpi=300, bbox_inches='tight')

    # 将图片转换为Base64编码
    with open(os.path.join(args.output_dir,"hetmap.png"), "rb") as img_file:
        b64_string = base64.b64encode(img_file.read()).decode('utf-8')



    # ###table1基本统计信息表格
    basic_information_corresponding=[sum(contig_length),len(contig_name),\
    N_count,int(round(sum(contig_GC_number)/sum(contig_length)*100,0)),min(contig_length),max(contig_length),]
    value=0
    for N50 in sorted(contig_length, reverse=True):
        value+=N50
        if value >=sum(contig_length)/2:
            basic_information_corresponding.append(N50)
            break
    table1={'statistic information':basic_information_list,'value':basic_information_corresponding}
    df1=pd.DataFrame(table1)

    ###table2基本统计信息
    contig_GC_content=[int(round(contig_GC_number[i]/contig_length[i]*100,0)) for i in range(len(contig_name))]
    if not first_line.startswith(">") and len(contig_name)==len(contig_depth):
        table2={'contig name':contig_name,'contig GC content(%)':contig_GC_content,'contig length':contig_length,'contig depth':contig_depth}
    else:
        table2={'contig name':contig_name,'contig GC content(%)':contig_GC_content,'contig length':contig_length}



    df2=pd.DataFrame(table2)
    df2.insert(0, 'Index', range(1, 1 + len(df2)))
    table1_html = df1.to_html(index=False)
    table2_html = df2.to_html(index=False)

    table1_html = table1_html.replace('<table border="1" class="dataframe">',
                                      '<table border="1" class="dataframe" style="font-size: 44px;">')
    table2_html = table2_html.replace('<table border="1" class="dataframe">',
                                      '<table border="1" class="dataframe" style="font-size: 44px;">')

    html_content = f'''
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>HiMT Report</title>
    </head>
    <body>
        <h1 style="font-size: 54px;" >fig.1 Statistics of completeness of conserved mitochondrial protein coding genes</h1>
        <img src="data:image/png;base64,{b64_string}" alt="Gene Expression Heatmap">
        <h2 style="font-size: 54px;">Table 1 Basic information of mitochondrial genome</h2>
        {table1_html}
        <h2 style="font-size: 54px;">Table 2 Contig information</h2>
        {table2_html}
    </body>
    </html>
    '''
    with open(os.path.join(args.output_dir, "heatmap_and_tables.html"), "w") as html_file:
        html_file.write(html_content)