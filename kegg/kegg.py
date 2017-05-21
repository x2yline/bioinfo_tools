#coding: utf-8
# 参考:https://wenku.baidu.com/view/a2539cb2fd0a79563c1e729c.html
def read_trans_file(file, target_column = [1,2]):
    target_dict = {}
    if file.split('\\')[-1][0:3] == 'GPL':
        with open(file) as f:
            for line in f:
                if line:
                    if line[0] not in ['#','!' '^']:
                        break
        target_column[0] = line.strip().split('\t').index('ID')
        target_column[1] = line.strip().split('\t').index('Entrez_Gene_ID')
    with open(file) as f:
        for line in f:
            if line[0] not in ['#','!' '^']:
                line_list = line.strip().split('\t')
                geneid = line_list[target_column[0]]
                gene_trans = line_list[target_column[1]]
                if geneid not in target_dict.keys():
                    target_dict[geneid] = gene_trans
                else:
                    target_dict[geneid] = [gene_trans] + list(target_dict[geneid])
    return(target_dict)
    
def parse_kegg(file_path):
    '''file_path 为本地文件 hsa00001.keg
    kegg中的基因为gene symbol'''
    kegg_dict={}
    with open(file_path,'r') as f:
        for line in f:
            if line.startswith('C'):
                try:
                    path_num = line[1:].strip()
                    kegg_dict[path_num]= []
                except:
                    print('Warning: '+line)
            elif line.startswith('D'):
                try:
                    gene_name = line[:line.find(';')].split()[-1]
                    kegg_dict[path_num].append(gene_name)
                except:
                    print('Warning: '+line)
    non_empty_path = []
    all_gene = []
    for keys in kegg_dict.keys():
        if kegg_dict[keys]:
           non_empty_path.append(keys)
           all_gene += kegg_dict[keys]
    print('The number of pathways in kegg is {}\nThe number of gene in kegg is {}'.format(len(set(non_empty_path)),len(set(all_gene))))
    return(kegg_dict, list(set(all_gene)))
    
def fish_test(x, M, n, N):
    '''x is the hitnums
    n is the genes nums in the specific pathway
    M is the total background genes
    N is the diff genes'''
    from scipy.stats import hypergeom
    pVal = hypergeom.sf(x, M, n, N, loc=0)
    return(pVal)

def get_N_in_kegg(gene_list, kegg_genes):
    gene_list_kegg = []
    for i in gene_list:
        if i in kegg_genes:
            gene_list_kegg.append(i)
    return(len(gene_list_kegg))

        
def kegg_enrich(gene_list, gene_symbol2id, file_path):
    '''gene_list为差异基因的列表,
    gene_symbol2id为字典,
    file_path为kegg文件数据的路径'''
    kegg_dict, kegg_genes = parse_kegg(file_path)
    N = get_N_in_kegg(gene_list, kegg_genes)
    M = len(kegg_genes)
    enrich_result = {}
    for key, vals in kegg_dict.items():
        if not vals:
            continue
        hit_special_path = []
        if ':' in key:
            kegg_map_url = 'http://www.kegg.jp/kegg-bin/show_pathway?map=' + key.split('PATH:')[-1].strip().strip(']') + '&multi_query='
        for gene in vals:
            if gene in gene_list:
                hit_special_path.append(gene)
                kegg_map_url += gene_symbol2id[gene]+'+red%2Cblue%0D%0A'
        if not hit_special_path:
            continue
        pval = fish_test(len(hit_special_path), M, len(vals), N)
        rich_factor = len(hit_special_path)/len(vals)
        enrich_result[key] = [len(hit_special_path), len(vals), pval, rich_factor, hit_special_path, kegg_map_url]
    
    import pandas as pd
    result = pd.DataFrame.from_dict(enrich_result)
    result = pd.DataFrame.transpose(result)
    result.columns = ['counts', 'path_genes', 'p_val', 'rich_factor', 'hit_genes', 'url']
    result = result.sort_values(by='p_val',axis=0)
    return(result)
def enrichment_plot(data, item_list=['pathway','gene_number', 'q_val', 'rich_factor']):
    '''data为字典
    其键为item_list
    item_list的顺序为[pathway_discription,
    hit_gene_number, q_val, rich_factor]'''
    import matplotlib.pyplot as plt
    import numpy as np
    import matplotlib as mpl
    x = [float(i) for i in data[item_list[3]]]
    y = [float(j) for j in range(len(data[item_list[0]]))]
    qvals = [float(i) for i in data[item_list[2]]]

    parameters =  np.linspace(np.min(qvals), np.max(qvals), len(y))
    norm = mpl.colors.Normalize(
        vmin=np.min(parameters),
        vmax=np.max(parameters))
    
    c_m = mpl.cm.spring#autumn
    s_m = mpl.cm.ScalarMappable(cmap=c_m, norm=norm)
    s_m.set_array([])
    fig = plt.figure(figsize=(11,12))
    fig.patch.set_facecolor('w')
    fig.suptitle('KEGG Enrichment', fontsize=24)
    ax = fig.add_axes([0.43, 0.1, 0.5, 0.8])

    ax.set_ylim([np.min(y), np.max(y)])

    ax.set_yticks([ j for j in range(len(data[item_list[0]]))])
    ax.set_yticklabels(data[item_list[0]], fontsize=18, color='k')

    for i in range(len(data[item_list[3]])):
        ax.plot(float(data[item_list[3]][i]), i, 'bo',
                markersize=float(data[item_list[1]][i])*1.2+2 ,clip_on=False,
                color=s_m.to_rgba(qvals[i]),
                 markeredgewidth=0.0)
    ax.set_xlim([np.min(x) - (np.max(x) - np.min(x))/20, np.max(x) + (np.max(x) - np.min(x))/20])
    ax.set_xticks(np.linspace(np.min(x) - (np.max(x) - np.min(x))/20, np.max(x) + (np.max(x) - np.min(x))/20, 5))
    ax.set_xticklabels([round(float(i), 2) for i in np.linspace(np.min(x) - (np.max(x) - np.min(x))/20, np.max(x) + (np.max(x) - np.min(x))/20, 5)])
    ax.get_xaxis().tick_bottom()
    ax.set_xlabel('Rich factor', fontsize=20)
    ax.set_ylabel('Kegg Pathways', fontsize=20)
    ax.get_yaxis().tick_left()
    ax.get_xaxis().set_tick_params(direction='out')
    ax.get_yaxis().set_tick_params(direction='out')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_position(('outward', 10))
    ax.spines['bottom'].set_position(('outward', 30))
    axbar = fig.add_axes([0.85, 0.4, 0.1, 0.3])
    axbar.set_axis_off()
    cb = plt.colorbar(s_m)
    cb.outline.set_visible(False)
    cb.set_label('Q value', labelpad=-20, y=1.07, rotation=0)
    axmarker = fig.add_axes([0.89, 0.75, 0.1, 0.1])
    axmarker.set_axis_off()
    axmarker.text(0, 1.2, 'Gene number', ha='center')
    for i in range(3):
        axmarker.plot(0, i*0.5, 'ko', markersize=(i*5+5)*1.2+2, clip_on=False)
        axmarker.text(0.02, i*0.5, str(i*5+5), clip_on=False, color='k', ha='left', va='center')
    plt.savefig('./kegg_enrichment_result/enrichment.png', dpi=100)
    plt.show()


# 需要文件geneid2symbol第一列为taxid, 第二列为geneid第三列为genesymbol
# 需要hsa00001.keg文件
# 需要diff_gene.txt作为富集的原始材料, 可以为symbol也可以为geneid
# 其余均不用做修改

gene_symbol2id = read_trans_file('geneid2symbol', [2,1])
geneid2symbol = read_trans_file('geneid2symbol', [1,2])
file_path = 'hsa00001.keg'

gene_list = []
with open("diff_gene.txt", 'r') as f:
    for line in f:
        gene_list.append(line.strip().split()[0].strip())
gene_list_new = []
if gene_list[0] in gene_symbol2id.values():
    for i in gene_list:
        gene_list_new.append(geneid2symbol[i])
gene_list = gene_list_new

enrichment_result = kegg_enrich(gene_list, gene_symbol2id, file_path)

import os
try:
	os.mkdir('kegg_enrichment_result')
except:
	pass
enrichment_result.to_csv('./kegg_enrichment_result/enrichment_result.csv')
def data_extract(enrichment_result, topnum=20):  
    enrichment_result = enrichment_result.iloc[:topnum,]
    a = enrichment_result.transpose()
    data = a.to_dict(orient="index")
    old_keys = list(data.keys())
    data['pathway'] = list(data['counts'].keys())
    for i in old_keys:
        data[i] = list(data[i].values())  
    data_new = {}
    data_new['pathway'] = [i.split('[')[0][i.find(' '):].strip() for i in data['pathway']]
    data_new['gene_number'] = data['counts']
    data_new['q_val'] = data['p_val']
    data_new['rich_factor'] = data['rich_factor']
    return(data_new)

data_new = data_extract(enrichment_result, topnum=30)
enrichment_plot(data_new, item_list=['pathway','gene_number', 'q_val', 'rich_factor'])

