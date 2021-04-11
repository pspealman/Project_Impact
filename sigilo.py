# -*- coding: utf-8 -*-
"""
Created on Sun Feb  9 23:03:48 2020

sigilo - for the generation and analysis of PICRUSt2 output

python sigilo.py --generate_heatmap -i pred_metagenome_unstrat.tsv -sig significant_objects_file -o metagenome

python sigilo.py --ko_enrichment -c pred_metagenome_contrib.tsv -t taxonomy.tsv -k ko00001.keg -pct 0.1 -pval 0.05 -o metagenome_contrib

python sigilo.py --versatility -c pred_metagenome_contrib.tsv -t taxonomy.tsv -k ko00001.keg -pct 0.1 -pval 0.05 -o metagenome_contrib

python sigilo.py --pathway_enrichment -c pred_metagenome_contrib.tsv -t taxonomy.tsv -k ko00001.keg -pct 0.05 -pval 0.05 -o metagenome_contrib_CPP
python scripts\sigilo_v8.py --pathway_enrichment -c picrust2_out_pipeline\KO_metagenome_out\pred_metagenome_contrib.tsv -t qiime_results\taxonomy.tsv -k C:\Gresham\Project_Gravimondo\ko00001.keg -pct 0.05 -pval 0.05 -o metagenome_contrib

python scripts/sigilo.py --asv2taxa -f qiime_results/feature-table.biom.txt -t qiime_results/taxonomy.tsv -o qiime_results/asv2taxa_feature_taxonomy.tab

python scripts/sigilo.py --bioindicator -f qiime_results/feature-table.biom.txt -t qiime_results/taxonomy.tsv -o qiime_results/bioi_


Version: Public 0.1 (I refuse to give up because I haven't tried all possible ways.)
Version: Beta 0.4 (Error Construct): added pathway_enrichment
Version: Beta 0.5 (Restaurant Dismiss): added versatility
Version: Beta 0.6 (Obscure Threat): added summative taxonomic levels
Version: Beta 0.7 (Leader Route):
    _x_ added asv2fa function to directly calculate the functional abundance of any taxa
Version: Beta 0.8 (Exhibition Brother)
    _x_ recode for Project Impact (2 sample sites with 3 replicates each)
Version: Beta 0.9 (Tick Ghostwriter)
    _x_ added asv2taxa function to generate unnormalized taxa abundances 
    _x_ added bioindicator function (2 sample sites with 3 replicates each) using kw

python scripts/sigilo.py --generate_heatmap -i picrust2_out_pipeline/KO_metagenome_out/pred_metagenome_unstrat.tsv -sig Aldex_results/aldex_significant_all_glm.ep_KO.csv -k metadata/ko00001.keg -o sigilo_results/
python scripts/sigilo.py --pathway_enrichment -c picrust2_out_pipeline\KO_metagenome_out\pred_metagenome_contrib.tsv -t qiime_results\taxonomy.tsv -k C:\Gresham\Project_Gravimondo\ko00001.keg -pct 0.05 -pval 0.05 -o metagenome_contrib

    
python scripts/sigilo.py --asv2fa -level 4 -select novembro_results_2/site_specific_family_enrichment.tab -c picrust2_out_pipeline/KO_metagenome_out/pred_metagenome_contrib.tsv -t qiime_results/taxonomy.tsv -k metadata/ko00001.keg -o sigilo_results_2/asv2fa.tab
python scripts/sigilo.py --asv2nsti -level 4 -select novembro_results_2/site_specific_family_enrichment.tab -n picrust2_out_pipeline/marker_predicted_and_nsti.tsv -t qiime_results/taxonomy.tsv -o sigilo_results_2/asv2nsti.tab



@author: pspea
"""
import plotly.graph_objects as go
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-hm',"--generate_heatmap", action='store_true')
parser.add_argument('-i',"--input_abundance_file")
parser.add_argument('-sig',"--significant_objects_file")

parser.add_argument('-ko',"--ko_enrichment", action='store_true')
parser.add_argument('-c',"--contrib_file")
parser.add_argument('-n',"--nsti_file")
parser.add_argument('-t',"--taxonomy")
parser.add_argument('-select', '--select_taxa_list')
parser.add_argument('-pct', '--pct_threshold')
parser.add_argument('-pval', '--pval_threshold')
parser.add_argument('-k',"--kegg_file")

parser.add_argument('-plot_ko',"--plot_ko_enrichment", action='store_true')
parser.add_argument('-s',"--taxa_source")

parser.add_argument('-a2f',"--asv2fa", action="store_true")
parser.add_argument('-level',"--taxonomic_level")

parser.add_argument('-a2n',"--asv2nsti", action="store_true")

parser.add_argument('-pathway',"--pathway_enrichment", action='store_true')
parser.add_argument('-plot_pathway',"--plot_pathway_enrichment", action='store_true')

parser.add_argument('-a2t',"--asv2taxa", action='store_true')
parser.add_argument('-f',"--feature_file")

parser.add_argument('-bioi',"--bioindicator", action='store_true')



parser.add_argument('-o',"--output_file")

args = parser.parse_args()

#import packages:
import numpy as np
import scipy.stats as stats
import pickle

if args.taxa_source:
    taxa_source = args.taxa_source
else:
    taxa_source = 'silva'

if taxa_source.lower() == 'silva' or taxa_source.lower() == 's':
    prefixe = ['D_0__','D_1__','D_2__','D_3__','D_4__','D_5__','D_6__','D_7__','D_8__','D_9__','D_10__','D_11__','D_12__','D_13__','D_14__']
   
if 'green' in taxa_source.lower() or taxa_source.lower() == 'gg':
    prefixe = ['k__','p__','c__','o__','f__','g__','s__']
    
rank_order = ['species','genus','family','order','class','phylum','kingdom']
convert_taxa_to_rank = {'kingdom':0, 'phylum':1, 'class':2, 'order':3, 'family':4, 'genus':5, 'species': 6}

#globals
ko_dict = {}
universal_ko_lookup = {}

def load_sig_obj(so_file_name = args.significant_objects_file):
    so_file = open(so_file_name)
    sig_object_set = set()
    
    for line in so_file:
        ko = line.split(',')[0]
        if len(ko) > 0:
            ko = ko.replace('"','')
            sig_object_set.add(ko)
    so_file.close()
    
    return(sig_object_set)
    
def load_ko(kegg_file_name):
    universal_ko_dict = {}
    
    keg_file = open(kegg_file_name)
    
    for line in keg_file:
        if line[0] =='A':
            #A09100 Metabolism
            a_ko = line[1:6]
            universal_ko_lookup[a_ko]=line[7:].strip()
    
        if line[0] == 'B':
            #B  09102 Energy metabolism
            if len(line.strip()) > 1:
                b_ko = line[3:8]
                universal_ko_lookup[b_ko]=line[8:].strip()
    
        if line[0] == 'C':
            #C    00010 Glycolysis / Gluconeogenesis [PATH:ko00010]
            c_ko = line[5:10]
            universal_ko_lookup[c_ko]=line[10:].strip()
    
        if line[0] == 'D':
            #D      K00844  HK; hexokinase [EC:2.7.1.1]
            d_ko = line[7:13]
            universal_ko_lookup[d_ko]=line[13:].strip()
            
            if d_ko not in universal_ko_dict:
                universal_ko_dict[d_ko] = {'A':set(), 'B':set(), 'C':set()}
                
            universal_ko_dict[d_ko]['A'].add(a_ko)
            universal_ko_dict[d_ko]['B'].add(b_ko)
            universal_ko_dict[d_ko]['C'].add(c_ko)
    
    keg_file.close()  
    
    return(universal_ko_dict, universal_ko_lookup)
    
def ko_to_pathway(universal_ko_dict, universal_ko_lookup):
        
    ko_super_set_type = {'Photosynthesis [PATH:ko00195]':set(), 
                         'Sulfur metabolism [PATH:ko00920]':set(), 
                         'Nitrogen metabolism [PATH:ko00910]': set(), 
                         'Methane metabolism [PATH:ko00680]': set(), 
                         'Pentose phosphate pathway [PATH:ko00030]': set(),
                         'Carbon fixation in photosynthetic organisms [PATH:ko00710]':set(),
                         'Carbon fixation pathways in prokaryotes [PATH:ko00720]':set()
                         }
       
    ko_to_path_lookup = {}
    
    write_out = set()
    
    for ko in universal_ko_dict:
        for each in universal_ko_dict[ko]['C']:
            
            if universal_ko_lookup[each] == 'Carbon fixation pathways in prokaryotes [PATH:ko00720]':
                ko_super_set_type['Carbon fixation pathways in prokaryotes [PATH:ko00720]'].add(ko)
                ko_to_path_lookup[ko]='Carbon fixation pathways in prokaryotes'
                write_out.add(ko)
                
            if universal_ko_lookup[each] == 'Carbon fixation in photosynthetic organisms [PATH:ko00710]':
                ko_super_set_type['Carbon fixation in photosynthetic organisms [PATH:ko00710]'].add(ko)
                ko_to_path_lookup[ko]='Carbon fixation in photosynthetic organisms'
                write_out.add(ko)
                
            if universal_ko_lookup[each] == 'Photosynthesis [PATH:ko00195]':
                ko_super_set_type['Photosynthesis [PATH:ko00195]'].add(ko)
                ko_to_path_lookup[ko]='Photosynthesis'
                write_out.add(ko)
            
            if universal_ko_lookup[each] == 'Methane metabolism [PATH:ko00680]':
                ko_super_set_type['Methane metabolism [PATH:ko00680]'].add(ko)
                ko_to_path_lookup[ko]='Methane'
                write_out.add(ko)

            if universal_ko_lookup[each] == 'Sulfur metabolism [PATH:ko00920]':
                ko_super_set_type['Sulfur metabolism [PATH:ko00920]'].add(ko)
                ko_to_path_lookup[ko]='Sulfur metabolism'
                write_out.add(ko)
                
            if universal_ko_lookup[each] == 'Nitrogen metabolism [PATH:ko00910]':
                ko_super_set_type['Nitrogen metabolism [PATH:ko00910]'].add(ko)
                ko_to_path_lookup[ko]='Nitrogen metabolism'
                write_out.add(ko)
                
            if universal_ko_lookup[each] == 'Pentose phosphate pathway [PATH:ko00030]':
                ko_super_set_type['Pentose phosphate pathway [PATH:ko00030]'].add(ko)
                ko_to_path_lookup[ko]='Pentose phosphate pathway'
                write_out.add(ko)

    ko_outfile_name = args.output_file + '_ko_list.log'
    
    ko_outfile = open(ko_outfile_name, 'w')
                
    for ko in write_out:
        outline = ('{}\n').format(ko)
        ko_outfile.write(outline)
                    
    ko_outfile.close()
                
    return(ko_super_set_type, ko_to_path_lookup)

def parse_line(line, runmode='log'):
    global universal_ko_lookup
    
    print(line)
    if runmode == 'log':
        KO = line.split('\t')[0].strip()
        #
        P1 = np.log(max(float(line.split('\t')[1]),1))
        P2 = np.log(max(float(line.split('\t')[2]),1))
        P3 = np.log(max(float(line.split('\t')[3]),1))
        #
        S1 = np.log(max(float(line.split('\t')[4]),1))
        S2 = np.log(max(float(line.split('\t')[5]),1))
        S3 = np.log(max(float(line.split('\t')[6]),1))        
        #
        #
        if KO in universal_ko_lookup:
            description = universal_ko_lookup[KO]
        else:
            description = KO
    else:
        KO = line.split('\t')[0].strip()
        #
        P1 = float(line.split('\t')[1])
        P2 = float(line.split('\t')[2])
        P3 = float(line.split('\t')[3])
        #
        S1 = float(line.split('\t')[4])
        S2 = float(line.split('\t')[5])
        S3 = float(line.split('\t')[6])        
        #
#        V1 = float(line.split('\t')[7])
#        V2 = float(line.split('\t')[8])
#        V3 = float(line.split('\t')[9])
        #
        if KO in universal_ko_lookup:
            description = universal_ko_lookup[KO]
        else:
            description = KO
    
    #return(KO, P1, P2, P3, S1, S2, S3, V1, V2, V3, description)
    return(KO, P1, P2, P3, S1, S2, S3, description)
   
def go_rep_heatmap(ko_master_dict, output_filename, max_value):
    global ko_dict
        
    title_is = ('Heatmap of Significant KOs').format()

    site_list = ["Cpn_1", "Con_2", "Con_3", "Imp_1", "Imp_2", "Imp_3"]
    ko_list = []
    value_array = []
    
    for ko, data in ko_master_dict.items():
        print(ko, data)
        value_array.append(data)
        ko_line = ('{}: {}').format(ko, ko_dict[ko])
        ko_list.append(ko_line)
        
    ko_list.sort()

    fig = go.Figure(data=go.Heatmap(
                    z=value_array,
                    x=site_list,
                    y=ko_list,
                    hoverongaps = False,
                    zmin=0, zmax=max_value))
    
    fig.update_layout(
            title=title_is,
            autosize=False,
            width=3600,
            height=3600,
            margin=go.layout.Margin(
                    l=50,
                    r=50,
                    b=100,
                    t=100,
                    pad=4
                    )
            )

    #fig.show()
    fig.write_image(output_filename)
        
def go_median_heatmap(ko_master_dict, output_filename, max_value):
    global ko_dict
          
    title_is = ('Heatmap of Significant KOs').format()

    site_list = ["Con", "Imp"]
    ko_list = []
    value_array = []
    
    for ko, data in ko_master_dict.items():
        data=ko_master_dict[ko]
        print(ko, data)
        value_array.append(data)
        ko_line = ('{}: {}').format(ko, ko_dict[ko])
        ko_list.append(ko_line)
        
    ko_list.sort()

    fig = go.Figure(data=go.Heatmap(
                    z=value_array,
                    x=site_list,
                    y=ko_list,
                    hoverongaps = False,
                    zmin=0, zmax=max_value))
    
    fig.update_layout(
            title=title_is,
            autosize=False,
            width=3600,
            height=3600,
            margin=go.layout.Margin(
                    l=50,
                    r=50,
                    b=100,
                    t=100,
                    pad=4
                    )
            )
    
    #fig.show()
    fig.write_image(output_filename)
    
def return_log10(each_set):
    new_set = []
    for each_obs in each_set:
        if each_obs == 0:
            each_obs = 1
        else:
            each_obs = np.log10(each_obs)

        new_set.append(each_obs)
    
    return(new_set)
    
def min_max(is_05_min, is_05_max, is_all_min, is_all_max, abun_val_file_name, ko_set, ko_list):
    abun_val_file = open(abun_val_file_name)
        
    for line in abun_val_file:
        #KO	P1.1	P1.2	P1.3	P2.1	P2.2	P2.3	P3.1	P3.2	P3.3	description
        if line.split('\t')[0] != 'function':
            KO, P1, P2, P3, S1, S2, S3, description = parse_line(line)
            
            if KO in ko_set and KO in ko_list:
                is_05_min = min([P1, P2, P3, S1, S2, S3, is_05_min])  
                is_05_max = max([P1, P2, P3, S1, S2, S3, is_05_max])     
                
            if KO in ko_list:
                is_all_min = min([P1, P2, P3, S1, S2, S3, is_all_min])
                is_all_max = max([P1, P2, P3, S1, S2, S3, is_all_max])
                
    abun_val_file.close()
    
    return(is_05_min, is_05_max, is_all_min, is_all_max)
        
def with_sig_object():
    output_file_name = args.output_file
    #universal_ko_dict, universal_ko_lookup = load_ko('C:/Gresham/Project_Gravimondo/Project_Osun/sigilo/ko00001.keg')
    universal_ko_dict, universal_ko_lookup = load_ko(args.kegg_file)
    ko_super_set_type, ko_to_path_lookup = ko_to_pathway(universal_ko_dict, universal_ko_lookup)
    so_file = (args.significant_objects_file)
    #so_file = ('C:/Gresham/Project_Gravimondo/quarantine/full_pred_metagenome_unstrat.tsv/Aldex_results_old/aldex_significant_all_glm.ep_KO.csv')
    ko_set = load_sig_obj(so_file)
    
    is_val_name = ('{}isval.tab').format(output_file_name)
    is_val_file = open(is_val_name,'w')
    
    ko_val_name = ('{}koval.tab').format(output_file_name)
    ko_val_file = open(ko_val_name,'w')
    
    is_05_min = 1000
    is_05_max = 0
    is_all_min = 1000
    is_all_max = 0
    
    for super_condition, ko_select in ko_super_set_type.items():
        
        abun_val_file_name = args.input_abundance_file
        ko_list = []
        for ko in ko_select:
            ko_list.append(ko)
                
        is_05_min, is_05_max, is_all_min, is_all_max = min_max(is_05_min, is_05_max, is_all_min, is_all_max, abun_val_file_name, ko_set, ko_list)
#        print(is_05_min, is_05_max, is_all_min, is_all_max)
#        1/0
        
    for super_condition, ko_select in ko_super_set_type.items():
        outline = ('{}\t{}\t{}\t{}\n').format('KO', 'Con', 'Imp', super_condition)
        ko_val_file.write(outline)
        print(super_condition)
        
        if len(ko_select)>0:
            condition = ko_to_path_lookup[list(ko_select)[0]]
            #ko_dict = {}
            rep_sig_ko_dict = {}
            median_sig_ko_dict = {}
            rep_all_ko_dict = {}
            median_all_ko_dict = {}
            sig_max_value = 0
            all_max_value = 0
            
            all_rep_heatmap_name = ('{}{}_all_replicates_heatmap.pdf').format(output_file_name, condition)
            sig_rep_heatmap_name = ('{}{}_0.05_replicates_heatmap.pdf').format(output_file_name, condition)
            all_median_heatmap_name = ('{}{}_all_median_heatmap.pdf').format(output_file_name, condition)
            sig_median_heatmap_name = ('{}{}_0.05_median_heatmap.pdf').format(output_file_name, condition)
            
            ko_list = []
            for ko in ko_select:
                ko_list.append(ko)
                
            abun_val_file_name = args.input_abundance_file
            
            abun_val_file = open(abun_val_file_name)
            
            for line in abun_val_file:
                #KO	P1.1	P1.2	P1.3	P2.1	P2.2	P2.3	P3.1	P3.2	P3.3	description
                if line.split('\t')[0] != 'function':
                    KO, P1, P2, P3, S1, S2, S3, description = parse_line(line)
                    
                    if KO in ko_set and KO in ko_list:
                        sig_max_value = max([P1, P2, P3, S1, S2, S3, sig_max_value])
                                
                        if KO not in rep_sig_ko_dict:
                            ko_dict[KO]=description.strip()
                            rep_sig_ko_dict[KO]=[P1, P2, P3, S1, S2, S3]
                            median_sig_ko_dict[KO]=[np.median([P1, P2, P3]), np.median([S1, S2, S3])]
                            
                            outline = ('{}\t{}\t{}\t{}\n').format(KO, np.median([P1, P2, P3]), np.median([S1, S2, S3]), description)
                            ko_val_file.write(outline)
                        else:
                            print('Error: Duplicate KO Identified')
                            quit()
                            
                    if KO in ko_list:
                        all_max_value  = max([P1, P2, P3, S1, S2, S3, all_max_value])
                                
                        if KO not in rep_all_ko_dict:
                            ko_dict[KO]=description.strip()
                            rep_all_ko_dict[KO]=[P1, P2, P3, S1, S2, S3]
                            median_all_ko_dict[KO]=[np.median([P1, P2, P3]), np.median([S1, S2, S3])]
                        else:
                            print('Error: Duplicate KO Identified')
                            quit()
            
            abun_val_file.close()
            
            go_rep_heatmap(rep_sig_ko_dict, sig_rep_heatmap_name, is_05_max)
        
            go_median_heatmap(median_sig_ko_dict, sig_median_heatmap_name, is_05_max)
            
            go_rep_heatmap(rep_all_ko_dict, all_rep_heatmap_name, is_all_min)
        
            go_median_heatmap(median_all_ko_dict, all_median_heatmap_name, is_all_min)
            
            outline = ('For {}, {} significant of {} detected, out of {} total.\n').format(condition, len(rep_sig_ko_dict), len(rep_all_ko_dict), len(ko_list))
            is_val_file.write(outline)
        
    is_val_file.close()
    ko_val_file.close()

### TODO parse fuck
def parse_taxonomy(taxa_file_name, taxa_cutoff_num):
    taxa_file = open(taxa_file_name)
    
    taxa_set = set()
    taxa_to_asv_dict = {}
    asv_to_taxa_dict = {}
        
    for line in taxa_file:
        if line[0]!='#':
            line = line.replace('"','')
            line = line.strip()
            asv = line.split('\t')[0]
            taxa = line.split('\t')[1]
            for each in prefixe:
                taxa = taxa.replace(each,'')

            taxa = taxa.replace(';','_')
            while taxa[-1] == '_':
                taxa = taxa[:-1]
                                
            if taxa.count('_') >= taxa_cutoff_num:
                taxa_set.add(taxa)
                
                if taxa.count('_') > taxa_cutoff_num:
                    taxa_list = taxa.split('_')[:taxa_cutoff_num+1]
                    taxa = ''
                    for each in taxa_list:
                        taxa+=str(each)+'_'
                    
                    if taxa[-1] == '_':
                        taxa = taxa[:-1]
                               
                if taxa not in taxa_to_asv_dict:
                    taxa_to_asv_dict[taxa] = [asv]
                else:
                    taxa_to_asv_dict[taxa].append(asv)
                    
                if asv not in asv_to_taxa_dict:
                    asv_to_taxa_dict[asv] = taxa
    
                else:
                    print('err')
                
    taxa_file.close()
        
    return(taxa_to_asv_dict, asv_to_taxa_dict, taxa_set)
    
def convert_sample_to_site(sample):
    if sample in ['P1', 'P2', 'P3']:
        return('P')
    if sample in ['M1', 'M2', 'M3']:
        return('S')
        
def parse_pathway_contrib(contrib_name,  taxa_dict, level, convert_to_site=True):
    global ordered_list, ko_to_path_lookup
    contrib = open(contrib_name)
    ko_dict = {'total':{}}
#    taxonomic_level_total = {}
    ko_tfa_dict = {}
    pathway_totals_dict = {}
    
    for line in contrib:
        if 'function' not in line:
            if convert_to_site:
                site = convert_sample_to_site(line.split('\t')[0])
            else:
                site = line.split('\t')[0]
                                
            ko = line.split('\t')[1]

            if ko in ordered_list:
                asv = line.split('\t')[2]
                
                taxon_function_abun = float(line.split('\t')[6])
                
                if asv in taxa_dict:
                    taxon = taxa_dict[asv]
                                                                                                                
                    if site not in ko_dict:
                        ko_dict[site] = {}
                    if level not in ko_dict[site]:
                        ko_dict[site][level] = {}    
                    if ko not in ko_dict[site][level]:
                        ko_dict[site][level][ko] = {}                        
                    if taxon not in ko_dict[site][level][ko]:
                        ko_dict[site][level][ko][taxon] = 0
                        
                    ko_dict[site][level][ko][taxon] += taxon_function_abun
                        
                    pathway = ko_to_path_lookup[ko]
                    
                    #pathway_totals_dict[level][site][pathway]
                    if level not in pathway_totals_dict:
                        pathway_totals_dict[level] = {} 
                    if site not in pathway_totals_dict[level]:
                        pathway_totals_dict[level][site] = {} 
                    if pathway not in pathway_totals_dict[level][site]:
                        pathway_totals_dict[level][site][pathway] = 0
                                                                                          
                    pathway_totals_dict[level][site][pathway] += taxon_function_abun
                            
                    if ko not in ko_tfa_dict:
                        ko_tfa_dict[ko] = set()
                    else:
                        ko_tfa_dict[ko].add(site)
                    
                    taxon = taxon.rsplit(';',1)[0]
                         
    contrib.close()
    
    return(ko_dict, ko_tfa_dict, pathway_totals_dict)
            
def run_mwu(x_set, y_set):
    if sum(x_set) > 30 or sum(y_set) > 30:
        _w, p_xvy = stats.mannwhitneyu(x_set, y_set)
        return(p_xvy)
    else:
        return(1)
        
def fill_out(inlist, length):
    
    while len(inlist) < length:
        inlist.append(0)
        
    return(inlist)
    
def try_kruskal_3(p_is, s_is, v_is):
    if (sum(p_is) >= 5) and (sum(s_is) >= 5) and (sum(v_is) >= 5):
        if (p_is != s_is) and (s_is != v_is):
            is_list = [p_is, s_is, v_is]
            for x in range(len(is_list)):
                for y in range(len(is_list)):
                    if x != y:
                        if set(is_list[x]) == set(is_list[y]):
                            return(False)
            return(True)
            
    return(False)
    
def try_kruskal_2(p_is, s_is):
    if (sum(p_is) >= 5) and (sum(s_is) >= 5):
        if (p_is != s_is):
            is_list = [p_is, s_is]
            for x in range(len(is_list)):
                for y in range(len(is_list)):
                    if x != y:
                        if set(is_list[x]) == set(is_list[y]):
                            return(False)
            return(True)
            
    return(False)
    
def run_chi2(x_num, x_den, y_num, y_den):
    if (x_num) > 5 or (y_num) > 5:
        obs = np.array([[max(x_num,1), x_den], [max(y_num,1), y_den]])
        chi2, pval, dof, expected = stats.chi2_contingency(obs, correction=True)
        return(pval)
    else:
        return(1)
        
def run_chi2x3(x_num, x_den, y_num, y_den, z_num, z_den):
    if (x_num) > 5 or (y_num) > 5 or (z_num) > 5:
        obs = np.array([[max(x_num,1), x_den], [max(y_num,1), y_den], [max(z_num,1), z_den]])
        chi2, pval, dof, expected = stats.chi2_contingency(obs, correction=True)
        return(pval)
    else:
        return(1)
        
def clean_up_name(name):
    name = name.strip()
    
    while name[-1] == ',':
        name=name[:-1]
        name = name.strip()
    
    return(name)
#    
def make_plot_figures(file_name):
    #for rank in ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']:
        #file_name = ('plot_{}_{}.p').format(args.output_file, rank)
    
        file = open(file_name,'rb')
        plot_round = pickle.load(file)
        file.close()
        
        for uname in plot_round:
            p_tag = ('Conserved')
            s_tag = ('Impacted')
            #v_tag = ('Valley')
            x_data = p_tag, s_tag
            
            p_list = return_log10(plot_round[uname]['P'])
            s_list = return_log10(plot_round[uname]['S'])
            #v_list = return_log10(plot_round[uname]['V'])
            y_data = p_list, s_list

            print(uname)
            uname = clean_up_name(uname)
          
            fig = go.Figure()
            
            colors = 'rgba(44, 160, 101, 0.5)', 'rgba(128, 0, 128, 0.5)'
            #, 'rgba(93, 164, 214, 0.5)'
            
            outfile_name = ('{}_enrichment.pdf').format(uname)
            print(outfile_name)
            
            for xd, yd, clrs in zip(x_data, y_data, colors):
                    fig.add_trace(go.Box(
                        #,
                        y=yd,
                        name=xd,
                        
                        boxpoints='all',
                        notched=True,
                        jitter=0.5,
                        whiskerwidth=0.2,
                        fillcolor=clrs,
                        line_color=clrs,
                        marker_size=5,
                        line_width=1,
                        showlegend=False)
                    )
                    
            fig.update_layout(
                title=uname.replace(',', '\n'),
                xaxis_title="Sample Site",
                yaxis_title="Log10(Relative Functional Abundance)",
                        font=dict(
                            family="Courier New, monospace",
                            size=10,
                            color="#7f7f7f"
                        )
            )
            
                    
            #fig.show()
            fig.write_image(outfile_name)
        
def pathway_build(taxa_cutoff_name, pct_threshold=0.05, pval_threshold=0.05):    
    outfile_name = ('{}_{}.log').format(args.output_file, taxa_cutoff_name)
    outfile = open(outfile_name,'w')
    
    header = ('Level\tTaxa\tPathway\tKW-pval\tsubl_taxa\tsubl_total\tinter_taxa\tinter_total\tsupra_taxa\tsupra_total\tko\n')
    outfile.write(header)
    
    correction = 1
    pval_ct = 0
    global ko_to_path_lookup, ko_dict, pathway_totals_dict
    
    #correction = len(ko_to_path_lookup)

    pathway_round = {}
    plot_round = {}
    pct_round = {}
    
    for site in ko_dict:
        if site != 'totals':
            for level in ko_dict[site]:
                for ko in ko_dict[site][level]:                      
                    for taxon, taxa_fa in ko_dict[site][level][ko].items():
                        if ko in ko_to_path_lookup:
                            pathway = ko_to_path_lookup[ko]
                           
                            if level not in pathway_round:
                                pathway_round[level] = {}
                            if taxon not in pathway_round[level]:
                                pathway_round[level][taxon] = {}
                            if pathway not in pathway_round[level][taxon]:
                                pathway_round[level][taxon][pathway] = {}
                            if site not in pathway_round[level][taxon][pathway]:
                                pathway_round[level][taxon][pathway][site] = {}
                            if ko not in pathway_round[level][taxon][pathway][site]:
                                pathway_round[level][taxon][pathway][site][ko] = []
                            
                            pathway_round[level][taxon][pathway][site][ko].append(taxa_fa) 
                            
                                                    
    for level in pathway_round:
        for taxon in pathway_round[level]:
            for pathway in pathway_round[level][taxon]:
                temp_path = {}
                temp_total = {}
                for site in pathway_round[level][taxon][pathway]:
                    temp_path[site]=[]
                    for ko in pathway_round[level][taxon][pathway][site]:
                        for taxa_fa in pathway_round[level][taxon][pathway][site][ko]:
                            temp_path[site].append(taxa_fa)   
                    
                    #print(pathway_totals_dict)
                    temp_total[site] = pathway_totals_dict[level][site][pathway]
                
                uname = ('{}_{}_{}').format(level, taxon, pathway)
                p_is = [0]
                s_is = [0]
                #v_is = [0]
                
                p_den = 1
                s_den = 1
                #v_den = 1
                
                if 'P' in temp_path:
                    p_is = temp_path['P']
                    p_den = temp_total['P']
                if 'S' in temp_path:
                    s_is = temp_path['S']
                    s_den = temp_total['S']
                    
                p_is = fill_out(p_is, max([len(p_is), len(s_is)]))
                s_is = fill_out(s_is, max([len(p_is), len(s_is)]))
                #v_is = fill_out(v_is, max([len(p_is), len(s_is), len(v_is)]))
                                
                runmode = 'kruskal'
                if runmode == 'kruskal':
                    if try_kruskal_2(p_is, s_is):
                        _stat, pval = stats.kruskal(p_is, s_is)                    
                    else:
                        pval = 1
                else:
                    #pval = run_chi2x3(sum(p_is), p_den, sum(s_is), s_den, sum(v_is), v_den)
                    pval = run_chi2(sum(p_is), p_den, sum(s_is), s_den)

                if pval*correction <= pval_threshold:
                    valcheck = set()
                    site_set = set()
                    for site in temp_path:
                        if len(temp_path[site]) >= 2:
                            site_set.add(site)
                            
                    for x_site in site_set:
                        for y_site in site_set:
                            if x_site != y_site:
                                x_set = temp_path[x_site]
                                y_set = temp_path[y_site]
                                
                                x_den = temp_total[x_site]
                                y_den = temp_total[y_site]
                                
                                if (sum(x_set) >= pct_threshold*x_den) or (sum(y_set) >= pct_threshold*y_den):
                                
                                    if runmode == 'kruskal':
                                        pval_2 = run_mwu(x_set, y_set)
                                    else:
                                        pval_2 = run_chi2(sum(x_set), x_den, sum(y_set), y_den)
                                    
                                    if (np.median(y_set) == 0):
                                        x_ratio = np.median(x_set)
                                    else:
                                        x_ratio = np.median(x_set)/np.median(y_set)
                                    
                                    if (np.median(x_set) == 0):
                                        y_ratio = np.median(y_set)
                                    else:
                                        y_ratio = np.median(y_set)/np.median(x_set)
                                        
                                    if ((pval_2*correction) <= pval_threshold) and ((x_ratio >= 1+pct_threshold) or (y_ratio >= 1+pct_threshold)):
                                        print(pval, pval_2, x_ratio, y_ratio)
                                        checkname = ('{}v{}').format(min(x_site,y_site),max(x_site,y_site))
                                        valcheck.add(checkname)
                        
                    if len(valcheck) >= 1:
                        p_pct = sum(p_is)/p_den
                        s_pct = sum(s_is)/s_den
                        #v_pct = sum(v_is)/v_den                                    
                        #outline = ('{level}\t{taxon}\t{pathway}\t{pval}\t{p_is}\t{p_den}\t{s_is}\t{s_den}\t{v_is}\t{v_den}\t{ko}\n').format(level=level, taxon=taxon, pathway=pathway, pval=(pval*correction), p_is=sum(p_is), p_den=p_den, s_is=sum(s_is), s_den=s_den, v_is=sum(v_is), v_den=v_den, ko=list(pathway_round[level][taxon][pathway][site].keys()))
                        outline = ('{level}\t{taxon}\t{pathway}\t{pval}\t{p_is}\t{p_den}\t{s_is}\t{s_den}\t{ko}\n').format(level=level, taxon=taxon, pathway=pathway, pval=(pval*correction), p_is=sum(p_is), p_den=p_den, s_is=sum(s_is), s_den=s_den, ko=list(pathway_round[level][taxon][pathway][site].keys())) 
                        outfile.write(outline)

                        pval_ct+=1
                        
                        if uname in plot_round:
                            print('uname duplicate', plot_round[uname])
                        
                        plot_round[uname]={'P':p_is, 'S': s_is, 'ko':list(pathway_round[level][taxon][pathway][site].keys())}
                        pct_round[uname]={'P':p_pct, 'S': s_pct, 'ko':list(pathway_round[level][taxon][pathway][site].keys())}
                        
    
    pickle_name = ('{}_pct_{}.p').format(args.output_file, taxa_cutoff_name)
    pickle.dump(pct_round, open(pickle_name, 'wb'))
    
    pickle_name = ('{}_plot_{}.p').format(args.output_file, taxa_cutoff_name)
    pickle.dump(plot_round, open(pickle_name, 'wb'))
    
    if len(plot_round) > 0:
        make_plot_figures(pickle_name)

    outfile.close()
    
def parse_select_taxa():
    select_taxa = set()
    select_file = args.select_taxa_list
    
    for line in select_file:
        select_taxa.add(line.strip())
    
    select_file.close()
    return(select_taxa)  
    
def check_for_null(taxa):
    taxa_list = taxa.split(';')
    
    taxa_string = ''
    
    for each in ['D_0__', 'D_1__', 'D_2__', 'D_3__', 'D_4__', 'D_5__', 'D_6__']:
        for every in taxa_list:
            if each in every:
                new_taxa, test = every.split(each)
                print(test)
                
                if len(test) > 0:
                    taxa_string += every + ';'
                else:
                    return(taxa_string[:-1])
                    
    return(taxa_string[:-1])
    
def clean_up_taxa(taxa):
    if 'Unassigned' in taxa:
        return(False)
    
    if 'D_7__' in taxa:
        taxa = taxa.split('D_7__')[0]
        
    taxa = check_for_null(taxa)
    
    taxa = taxa.replace('D_0__','k__')
    taxa = taxa.replace('D_1__','p__')
    taxa = taxa.replace('D_2__','c__')    
    taxa = taxa.replace('D_3__','o__')
    taxa = taxa.replace('D_4__','f__')
    taxa = taxa.replace('D_5__','g__')
    taxa = taxa.replace('D_6__','s__')
    
    return(taxa)
    
def check_sum(abundance_list):
    total = 0

    for each in abundance_list:
        total += float(each)
        
    return(total)
    
if args.generate_heatmap:
    with_sig_object()
        
if args.pathway_enrichment:
    if args.pct_threshold:
        pct_threshold = float(args.pct_threshold)
    else:
        pct_threshold = 1
        
    if args.pval_threshold:
        pval_threshold = float(args.pval_threshold)
    else:
        pval_threshold = 0.05
    
    ordered_list = []
    
    universal_ko_dict, universal_ko_lookup = load_ko(args.kegg_file)
    ko_super_set_type, ko_to_path_lookup = ko_to_pathway(universal_ko_dict, universal_ko_lookup)
    
    taxonomy_name = args.taxonomy
    contrib_name = args.contrib_file
    
    for superset, ko_list in ko_super_set_type.items():
        for ko in ko_list:
            ordered_list.append(ko)
    
    for taxa_cutoff_name in rank_order:
        taxa_cutoff_num = convert_taxa_to_rank[taxa_cutoff_name]
        
        print('For ', taxa_cutoff_name)
        print('... running parse_taxonomy')
        taxa_to_asv_dict, taxa_dict, taxa_set = parse_taxonomy(taxonomy_name, taxa_cutoff_num)
        print('... running parse_contrib')
        ko_dict, ko_tfa_dict, pathway_totals_dict = parse_pathway_contrib(contrib_name, taxa_dict, taxa_cutoff_name)
        print('... running pathway_build')
        pathway_build(taxa_cutoff_name, pct_threshold, pval_threshold)
        

def normalize_counts(asv_counts, taxa_to_asv_dict):
    impacted_1, impacted_2, impacted_3 = 0, 0, 0
    pristine_1, pristine_2, pristine_3 = 0, 0, 0
    
    for asv in asv_counts:
        i1, i2, i3, p1, p2, p3 = asv_counts[asv]

        impacted_1 += i1
        impacted_2 += i2
        impacted_3 += i3
        #
        pristine_1 += p1
        pristine_2 += p2
        pristine_3 += p3
        #
    
    global_min = min(min([pristine_1, pristine_2, pristine_3]),
                     min([impacted_1, impacted_2, impacted_3])
                     )
    
    #Define number of observations per site for the purposes of downsampling
    pristine_cor_1, pristine_cor_2, pristine_cor_3 = global_min/pristine_1, global_min/pristine_2, global_min/pristine_3
    impacted_cor_1, impacted_cor_2, impacted_cor_3 = global_min/impacted_1, global_min/impacted_2, global_min/impacted_3  

    taxa_to_counts = {}
    for taxa, asvs in taxa_to_asv_dict.items():
        i1, i2, i3 = 0, 0, 0
        p1, p2, p3 = 0, 0, 0
        
        for asv in asvs:
            if asv in asv_counts:
               impacted_1, impacted_2, impacted_3, pristine_1, pristine_2, pristine_3 = asv_counts[asv]
               i1 += impacted_1
               i2 += impacted_2
               i3 += impacted_3
               
               p1 += pristine_1
               p2 += pristine_2
               p3 += pristine_3
           
        taxa_to_counts[taxa] = [i1, i2, i3, p1, p2, p3]
        
        #taxa_dict = {'total':[0,0], 'observed':0}
    
    taxa_dict = {}  
    for taxa, counts in taxa_to_counts.items():
        #
        impacted_1 = counts[0] * impacted_cor_1
        impacted_2 = counts[1] * impacted_cor_2
        impacted_3 = counts[2] * impacted_cor_3
        impacted = [impacted_1, impacted_2, impacted_3]
        #
        pristine_1 = counts[3] * pristine_cor_1
        pristine_2 = counts[4] * pristine_cor_2
        pristine_3 = counts[5] * pristine_cor_3
        pristine = [pristine_1, pristine_2, pristine_3]
        #

        #           
        if taxa not in taxa_dict:
            taxa_dict[taxa] = [impacted, pristine]
        else:
            for index in range(3):
                taxa_dict[taxa][0][index] += impacted[index]
                taxa_dict[taxa][1][index] += pristine[index]
            
#        taxa_dict['total'][0] += sum(impacted)
#        taxa_dict['total'][1] += sum(pristine)
#        
#        taxa_dict['observed'] += len([s for s in pristine if s != 1]) + len([s for s in impacted if s != 1])

    return(taxa_dict)
    
def build_asv_counts(feature_table_name):
    #this is to build the counts table
    asv_file = open(feature_table_name)
    
    asv_counts = {}
    total_count = 0
    
    
    for line in asv_file:
        if line[0]!='#' and 'Feature ID' not in line:
            #asv ID	M1	M2	M3	P1	P2	P3
            line = line.strip()
            asv = line.split('\t')[0].strip()# 
            #
            impacted_1 = float(line.split('\t')[4])
            impacted_2 = float(line.split('\t')[5])  
            impacted_3 = float(line.split('\t')[6])
            #
            pristine_1 = float(line.split('\t')[1])
            pristine_2 = float(line.split('\t')[2])  
            pristine_3 = float(line.split('\t')[3])
            #
            asv_counts[asv] = [impacted_1, impacted_2, impacted_3, pristine_1, pristine_2, pristine_3]
            total_count += sum([impacted_1, impacted_2, impacted_3, pristine_1, pristine_2, pristine_3])
    #        
    asv_file.close()
    return(asv_counts, total_count)
    
def build_taxa_asv_association(taxa_file_name):
    taxa_set = set()
    taxa_to_asv_dict = {}
    asv_to_taxa_dict = {}
    
    taxa_file = open(taxa_file_name)
    
    for line in taxa_file:
        if line[0]!='#':
            line = line.replace('"','')
            line = line.strip()
            asv = line.split('\t')[0]
            taxa = line.split('\t')[1]

            for each in prefixe:
                taxa = taxa.replace(each,'')

            taxa = taxa.replace(';','_')
            while taxa[-1] == '_':
                taxa = taxa[:-1]
                                   
            if taxa.count('_') >= taxa_cutoff_num:
                taxa_set.add(taxa)
                
                if taxa.count('_') > taxa_cutoff_num:
                    taxa_list = taxa.split('_')[:taxa_cutoff_num+1]
                    taxa = ''
                    for each in taxa_list:
                        taxa+=str(each)+'_'
                    
                    if taxa[-1] == '_':
                        taxa = taxa[:-1]
                                
                if taxa not in taxa_to_asv_dict:
                    taxa_to_asv_dict[taxa] = [asv]
                else:
                    taxa_to_asv_dict[taxa].append(asv)
                    
                if asv not in asv_to_taxa_dict:
                    asv_to_taxa_dict[asv] = taxa
    
                else:
                    print('err')
                
    taxa_file.close()
    
    return(taxa_set, taxa_to_asv_dict, asv_to_taxa_dict)

def asv_taxa_clean_up(taxa, level):    
    if taxa.count(';') >= level:
        taxa_temp = ''
        taxa_list = taxa.split(';')[:level+1]
        for each in taxa_list:
            if '__' in each:
                taxa_temp += each.split('__')[1]+'_'
            else:
                taxa_temp += each+'_'

        if taxa_temp[-1] == '_':
            taxa_temp = taxa_temp[:-1]
            
        taxa = taxa_temp
        return(taxa)

    else:
        taxa_temp = ''
        taxa_list = taxa.split(';')[:level]
        
        for each in taxa_list:
            print(each, taxa_list)
            if '__' in each:
                taxa_temp += each.split('__')[1]+'_'
            else:
                taxa_temp += each+'_'

        if taxa_temp[-1] == '_':
            taxa_temp = taxa_temp[:-1]
            
        taxa = taxa_temp
        
        return(taxa)

def bioi_stats(taxa_dict, taxa_cutoff_name):    
    header = ('taxa\t'
              'Normalized_Impacted_1\tNormalized_Pristine_1\t'
              'Normalized_Impacted_2\tNormalized_Pristine_2\t'
              'Normalized_Impacted_3\tNormalized_Pristine_3\t'
              'KW_stat\tKW_pval\tPristine_to_Impacted_ratio\n')
    
    all_outfile_name = ('{}_{}.tab').format(args.output_file, taxa_cutoff_name)
    all_outfile = open(all_outfile_name, 'w')
    all_outfile.write(header)
    
    for taxa, taxa_array in taxa_dict.items():
        imp_s = taxa_array[0]
        pri_s = taxa_array[1]
        
        if try_kruskal_2(imp_s, pri_s):
            kw_stat, kw_pval = stats.kruskal(imp_s, pri_s)                    
        else:
            kw_stat = 0
            kw_pval = 1
                                        
        outline = ('{taxa}\t'
                   '{imp1}\t{imp2}\t{imp3}\t'
                   '{pri1}\t{pri2}\t{pri3}\t'
                   '{kw_stat}\t{kw_pval}\t{ratio}\n').format(taxa=taxa,
                    imp1=imp_s[0], imp2=imp_s[1], imp3=imp_s[2],
                    pri1=pri_s[0], pri2=pri_s[1], pri3=pri_s[2],
                    kw_stat=kw_stat, kw_pval=kw_pval, ratio=sum(taxa_array[1])/max(sum(taxa_array[0]),1))
        
        all_outfile.write(outline)

    all_outfile.close()
        
if args.asv2taxa:
    if args.taxonomic_level: 
        level = int(args.taxonomic_level)         
    else:
        level = 7
        
    taxonomy = open(args.taxonomy)
        
    asv_dict = {}
    for line in taxonomy:
        if 'Feature' not in line:
            asv, taxa, _confidence = line.split('\t')
            taxa = asv_taxa_clean_up(taxa, level)
            
            if taxa:
                asv_dict[asv] = asv_taxa_clean_up(taxa, level)
            else:
                outline = ('Unassigned at {}').format(level)
                asv_dict[asv] = outline
                
    taxonomy.close()
        
    feature = open(args.feature_file)
    
    taxa_count_check = {}
    
    for line in feature:
        if line[0]!='#':
            line = line.strip()
            #asv ID	P1	P2	P3	M1	M2	M3
            asv, P1, P2, P3, M1, M2, M3 = line.split('\t')
            taxa = asv_dict[asv.strip()]
            
            if taxa not in taxa_count_check:
                taxa_count_check[taxa] = check_sum([P1, P2, P3, M1, M2, M3])
            else:
                taxa_count_check[taxa] += check_sum([P1, P2, P3, M1, M2, M3])
    
    feature.close()
    
    feature = open(args.feature_file)
    outfile = open(args.output_file, 'w')
    total_of_each_taxa = {}
    
    header = ('asv\ttaxa\tP1\tP2\tP3\tM1\tM2\tM3\n')
    
    outfile.write(header)
    for line in feature:
        if line[0]!='#':
            line = line.strip()
            #asv ID	P1	P2	P3	M1	M2	M3
            asv, P1, P2, P3, M1, M2, M3 = line.split('\t')
            taxa = asv_dict[asv.strip()]

            if taxa not in total_of_each_taxa:
                total_of_each_taxa[taxa] = {'P1':0, 'P2':0, 'P3':0, 'M1':0, 'M2':0, 'M3':0}
                
            total_of_each_taxa[taxa]['P1']+=float(P1)
            total_of_each_taxa[taxa]['P2']+=float(P2)
            total_of_each_taxa[taxa]['P3']+=float(P3)
            
            total_of_each_taxa[taxa]['M1']+=float(M1)
            total_of_each_taxa[taxa]['M2']+=float(M2)
            total_of_each_taxa[taxa]['M3']+=float(M3)
                        
            if taxa_count_check[taxa]>=5:
                if taxa:
                    outline = ('{asv}\t{taxa}\t{P1}\t{P2}\t{P3}\t{M1}\t{M2}\t{M3}\n').format(
                            asv=asv, taxa=taxa, 
                            P1=P1, P2=P2, P3=P3, 
                            M1=M1, M2=M2, M3=M3)
                    
                    outfile.write(outline)
                    
    outfile.close()
    feature.close()

    outfile = open(args.output_file, 'w')
    header = ('taxa\tP1\tP2\tP3\tM1\tM2\tM3\n')
    
    outfile.write(header)
    
    for taxa in total_of_each_taxa:
        outline = ('{taxa}\t{P1}\t{P2}\t{P3}\t{M1}\t{M2}\t{M3}\n').format(
            taxa=taxa, 
            P1=total_of_each_taxa[taxa]['P1'], P2=total_of_each_taxa[taxa]['P2'], P3=total_of_each_taxa[taxa]['P3'], 
            M1=total_of_each_taxa[taxa]['M1'], M2=total_of_each_taxa[taxa]['M2'], M3=total_of_each_taxa[taxa]['M3'])
        
        outfile.write(outline)
    outfile.close()
    
if args.bioindicator:
    for taxa_cutoff_name in rank_order:
        taxa_cutoff_num = convert_taxa_to_rank[taxa_cutoff_name]
        
        feature_file_name = args.feature_file
        taxa_file_name = args.taxonomy
        
        taxa_set, taxa_to_asv_dict, asv_to_taxa_dict = build_taxa_asv_association(taxa_file_name)
            
        asv_counts, _total_count = build_asv_counts(feature_file_name)
        taxa_dict = normalize_counts(asv_counts, taxa_to_asv_dict)
      
        bioi_stats(taxa_dict, taxa_cutoff_name)

if args.asv2fa: 
    level = int(args.taxonomic_level)         
    taxonomy = open(args.taxonomy)
    
    asv_to_taxa_dict = {}
    taxa_to_asv_dict = {}
    
    for line in taxonomy:
        if 'Feature' not in line:
            asv, taxa, _confidence = line.split('\t')
            taxa = asv_taxa_clean_up(taxa, level)
            if taxa:
                if asv not in asv_to_taxa_dict:
                    asv_to_taxa_dict[asv] = taxa
                    
                else:
                    print('something wrong')
                    
                if taxa not in taxa_to_asv_dict:
                    taxa_to_asv_dict[taxa] = set()
                
                taxa_to_asv_dict[taxa].add(asv)
                                    
    taxonomy.close()
            
    assigned_asv = {}
    
    select_taxa_list = open(args.select_taxa_list)
    
    for line in select_taxa_list:
        if line[0] != '#':
            taxa = line.split('\t')[0]
            #taxa = asv_taxa_clean_up(taxa, level)
            if taxa:
                if taxa in taxa_to_asv_dict:
                    assigned_asv[taxa] = taxa_to_asv_dict[taxa]
    
    select_taxa_list.close()
        
    universal_ko_dict, universal_ko_lookup = load_ko(args.kegg_file)
    ko_super_set_type, ko_to_path_lookup = ko_to_pathway(universal_ko_dict, universal_ko_lookup)
                        
    contrib = open(args.contrib_file)
    
    taxa_fun = {}
    
    for line in contrib:
        if 'function' not in line:
            #P1-1	K00001	01edd66886a699ad420ca0d8db401937	154.0	1.091424521615875	2	308.0	2.18284904323175
            ko = line.split('\t')[1]
            
            if ko in ko_to_path_lookup:                
                asv = line.split('\t')[2]
                
                if asv in asv_to_taxa_dict:
                    taxa = asv_to_taxa_dict[asv]
                    
                    #if taxa in assigned_asv:
                    site = convert_sample_to_site(line.split('\t')[0])
                    
                    taxon_rel_function_abun = float(line.split('\t')[7])
                    
                    if ko not in taxa_fun:
                        taxa_fun[ko] = {'total':{'P':0, 'S':0}}
                    
                    if taxa not in taxa_fun[ko]:
                        taxa_fun[ko][taxa] = {'P':0, 'S':0}
                        
                    taxa_fun[ko]['total'][site] += taxon_rel_function_abun
                    taxa_fun[ko][taxa][site] += taxon_rel_function_abun
    
    contrib.close()
    
    taxa_enrichment_by_site = {}
        
    outfile_base = args.output_file
    pathfile = open(outfile_base + '_pathways.log','w')
        
    for ko in taxa_fun:
        pathway = ko_to_path_lookup[ko]
        
        for taxa in taxa_fun[ko]:
            if (taxa != 'total') and (taxa in assigned_asv):
                for site in taxa_fun[ko][taxa]:
                    if (taxa_fun[ko][taxa][site]) > 0:                        
                        kts_ratio = (taxa_fun[ko][taxa][site]) / (taxa_fun[ko]['total'][site])
    
                        if kts_ratio >= 0.10:
                            if site not in taxa_enrichment_by_site:
                                taxa_enrichment_by_site[site] = {}
                            
                            if taxa not in taxa_enrichment_by_site[site]:
                                taxa_enrichment_by_site[site][taxa] = {}
                            
                            if pathway not in taxa_enrichment_by_site[site][taxa]:
                                taxa_enrichment_by_site[site][taxa][pathway] = {}
                                
                            if ko not in taxa_enrichment_by_site[site][taxa][pathway]:
                                taxa_enrichment_by_site[site][taxa][pathway][ko] = kts_ratio
                        
    
    for site in taxa_enrichment_by_site:
        for taxa in taxa_enrichment_by_site[site]:
            ko_name = ('{}_{}_{}_ko.log').format(outfile_base, taxa, site) 
            kofile = open(ko_name,'w')
            
            taxa_ko_set = set()
            for pathway in taxa_enrichment_by_site[site][taxa]:
                ko_set = taxa_enrichment_by_site[site][taxa][pathway]
                
                if len(ko_set) >= 3:
                    outline = ('{}\t{}\t{}\t{}\n').format(site,taxa,pathway,ko_set)
                    pathfile.write(outline)
                    
                for each in ko_set:
                    print(each)
                    taxa_ko_set.add(each)
            
            if len(taxa_ko_set) >= 3:
                ct = 0 
                for each in taxa_ko_set:
                    outline = ('gene_{}\t{}\n').format(ct, each)
                    ct+=1
                    kofile.write(outline)
            
            kofile.close()
    
    pathfile.close()         
      

if args.asv2nsti: 
    level = int(args.taxonomic_level)         
    taxonomy = open(args.taxonomy)
    
    asv_to_taxa_dict = {}
    taxa_to_asv_dict = {}
    
    for line in taxonomy:
        if 'Feature' not in line:
            asv, taxa, _confidence = line.split('\t')
            taxa = asv_taxa_clean_up(taxa, level)
            if taxa:
                if asv not in asv_to_taxa_dict:
                    asv_to_taxa_dict[asv] = taxa
                    
                else:
                    print('something wrong')
                    
                if taxa not in taxa_to_asv_dict:
                    taxa_to_asv_dict[taxa] = set()
                
                taxa_to_asv_dict[taxa].add(asv)
                                    
    taxonomy.close()
            
    assigned_asv = {}
    
    select_taxa_list = open(args.select_taxa_list)
    
    for line in select_taxa_list:
        if line[0] != '#':
            taxa = line.split('\t')[0]
            #taxa = asv_taxa_clean_up(taxa, level)
            if taxa:
                if taxa in taxa_to_asv_dict:
                    assigned_asv[taxa] = taxa_to_asv_dict[taxa]
    
    select_taxa_list.close()
        
    #universal_ko_dict, universal_ko_lookup = load_ko(args.kegg_file)
    #ko_super_set_type, ko_to_path_lookup = ko_to_pathway(universal_ko_dict, universal_ko_lookup)
                        
    contrib = open(args.nsti_file)
    
    nsti_fun = {}
    
    for line in contrib:
        if 'NSTI' not in line:
            #P1-1	K00001	01edd66886a699ad420ca0d8db401937	154.0	1.091424521615875	2	308.0	2.18284904323175
            asv = line.split('\t')[0]
            nsti = float(line.split('\t')[2])
            
            process = False
            
            if asv in asv_to_taxa_dict:
                taxa = asv_to_taxa_dict[asv]
                
                if taxa in assigned_asv:
                    if taxa not in nsti_fun:
                        nsti_fun[taxa] = []
                        
                    nsti_fun[taxa].append(nsti)
                
    contrib.close()
        
    outfile_name = args.output_file
    outfile = open(outfile_name,'w')
        
    for taxa in nsti_fun:
        nsti_list = nsti_fun[taxa]
        
        outline = ('{}\t{}\t{}\t{}\n').format(taxa, np.mean(nsti_list), np.median(nsti_list), np.std(nsti_list))
        print(outline)
        outfile.write(outline)
    
    outfile.close()

    
def define_type(s_taxa, p_taxa, raw_s_fa, raw_p_fa):
    coord_dict = {'S':[0,0], 'P':[0,0]} 
    s_fa = np.median(raw_s_fa)
    p_fa = np.median(raw_p_fa)
    
    coord_dict = {'S':[s_taxa, s_fa], 'P':[p_taxa, p_fa]}
    
    return(coord_dict)      

def plot_funbubbles():
    import plotly.graph_objects as go
    import pickle 

    file = open('C:/Gresham/Project_Gravimondo/Project_Impact_2/sigilo_results/pathway/pct_metagenome_contrib_family.p','rb')
    fun_dict = pickle.load(file)
    file.close()
    
    condition_dict = {}
    site_abundance_dict = {}
    
    taxon_set = set()
        
    for each in fun_dict:
        taxa = each.split('_')
        condition = taxa[-1]
        if condition not in condition_dict:
            condition_dict[condition] = set()
        
        taxa_list = taxa[1:-1]
        taxon = ''
        for taxa in taxa_list:
            taxon += taxa +'_'
        taxon = taxon[:-1]
        
        condition_dict[condition].add(taxon)
        taxon_set.add(taxon)
        
        for site in fun_dict[each]:
            if site != 'ko':
                if condition not in site_abundance_dict:
                    site_abundance_dict[condition] = {}
                if taxon not in site_abundance_dict[condition]:
                    site_abundance_dict[condition][taxon] = {'S':0, 'P':0}
                  
#                val = np.median(fun_dict[each][site])
#                
#                if val == 0:
#                    val = 1 
#                
#                val = np.log10(val)
#                
#                if val < 0:
#                    val = 0
                
                site_abundance_dict[condition][taxon][site] += fun_dict[each][site]
                        
    for condition, taxa_set in condition_dict.items():
        x_compound_dict = {'S':[], 'P':[]}
        y_compound_dict = {'S':[], 'P':[]}
    
        for taxon in taxa_set:
            for site in site_abundance_dict[condition][taxon]:
                val = site_abundance_dict[condition][taxon][site]

                x_compound_dict[site].append(val)
                y_compound_dict[site].append(taxon)
            
        outfile_name = ('family_{}.pdf').format(condition)
        
        import plotly.graph_objects as go
        
        fig = go.Figure()
        fig.add_trace(go.Scatter(
            x=x_compound_dict['S'],
            y=y_compound_dict['S'],
            marker=dict(color='rgba(128, 0, 128, 0.5)', size=6),
            mode="markers",
            name="Impact",
        ))
                
        fig.add_trace(go.Scatter(
            x=x_compound_dict['P'],
            y=y_compound_dict['P'],
            marker=dict(color='rgba(44, 160, 101, 0.5)', size=6),
            mode="markers",
            name="Conserved",
        ))

        fig.update_layout(title=condition,
                          xaxis_title="Percent Relative Abundance",
                          yaxis_title="Taxa",
                          font_size=10,
                          width=1500,
                          height=600)
        
        fig.show()

        fig.write_image(outfile_name)
        
        print(len(taxon_set))
        
def apply_versatility():
    global ko_dict
    
    versa_dict = {}

    for site in ko_dict:
        for level in ko_dict[site]:
            if level not in versa_dict:
                versa_dict[level] = {}
            if site not in versa_dict[level]:
                # a set of all taxa in a site at a level
                #   this is for normalization ot the KO specific numbers of taxa
                versa_dict[level][site] = {'total_taxa': set()}
                
            for ko in ko_dict[site][level]:
                if ko not in versa_dict[level][site]:
                    # the KO of every site should have:
                    #   how many taxa have the ko:
                    #       versa_dict[level][site][ko]['taxa_num'].add(taxon)
                    #   a taxa specific functional abundance:
                    #       versa_dict[level][site][ko][taxon] += (ko_dict[site][level][ko][taxon])
                    #   a distribution of all functional abundances
                    #       versa_dict[level][site][ko]['fa_dist'].append(versa_dict[level][site][ko][taxon])
                    versa_dict[level][site][ko] = {'taxa_num':set(), 'fa_dist':[]}
                
                for taxon in ko_dict[site][level][ko]:
                    versa_dict[level][site][ko]['taxa_num'].add(taxon)
                    versa_dict[level][site]['total_taxa'].add(taxon)
                    
                    if taxon not in versa_dict[level][site][ko]:
                        versa_dict[level][site][ko][taxon] = 0
                    
                    #The KO functional abundance of the a specific taxa
                    #   not normalised
                    #   ko_dict[site][level][ko]['total']
                    versa_dict[level][site][ko][taxon] += (ko_dict[site][level][ko][taxon])
                    versa_dict[level][site][ko]['fa_dist'].append(versa_dict[level][site][ko][taxon])

    fet_ct = 0
    versatile_dict = {}
    
    for level in versa_dict:
        temp_site_dict = {}
        site_list = {}
        
        for site in versa_dict[level]:
            # a set of all taxa in a site at a level
            total_taxa = len(versa_dict[level][site]['total_taxa'])

            for ko in versa_dict[level][site]:
                #filters 'total_taxa'
                if 'taxa_num' in versa_dict[level][site][ko]:
                    #only consider those sites with ko in each replicate
                    if ko not in site_list:
                        site_list[ko] = set()
                        
                    site_list[ko].add(site)
                    
                    #   how many taxa have the ko:
                    taxa_num = len(versa_dict[level][site][ko]['taxa_num'])
                    fa_dist = versa_dict[level][site][ko]['fa_dist']
                    
                    if ko not in temp_site_dict:
                        temp_site_dict[ko]={}
                        
                    if site not in temp_site_dict[ko]:
                        temp_site_dict[ko][site]={'taxa_num': 0, 'total_taxa': total_taxa, 'fa_dist': []}
                        
                    temp_site_dict[ko][site]['taxa_num'] += taxa_num
                    if len(fa_dist) > 1:
                        for each in fa_dist:
                            temp_site_dict[ko][site]['fa_dist'].append(each)
                    else:
                        temp_site_dict[ko][site]['fa_dist'].append(fa_dist[0])                                       
                    #temp_site_dict[ko][site]['total_taxa'] = total_taxa
                        
    
        for ko in temp_site_dict:
            process =  True
            #ko_sites = site_list[ko]
#            print(ko_sites)
            #for each in ['S', 'V', 'P']:
#            for each in ['S1', 'S2', 'S3', 'V1', 'V2', 'V3', 'P1', 'P2', 'P3']:
#                if each not in ko_sites:
#                    process = False
                    
            if process:
                uname = ('{}_{}').format(level, ko)
                print(uname)
#                for x_sample_set in [['M1', 'M2', 'M3'], ['P1', 'P2', 'P3']]:
#                    for y_sample_set in [['M1', 'M2', 'M3'], ['P1', 'P2', 'P3']]:
#                        if x_sample_set != y_sample_set:
                x_set = []
                x_t_set = []
                x_fa = []
                x_ratio = []
                
                y_set = []
                y_t_set = []
                y_fa = []
                y_ratio = []
                
                z_set = []
                z_t_set = []
                z_fa = []
                z_ratio = []
                
                for x in ['S1', 'S2', 'S3']:
                    if x in temp_site_dict[ko]:
                        x_num = temp_site_dict[ko][x]['taxa_num']
                        x_total = temp_site_dict[ko][x]['total_taxa']
                        x_fa_dist = temp_site_dict[ko][x]['fa_dist']
                        
                    else:
                        x_num = 0
                        x_total = 1
                        x_fa_dist = [1]
                    
                    x_set.append(x_num)
                    x_t_set.append(x_total)
                    
                    if len(x_fa_dist) > 1:
                        for each in x_fa_dist:
                            x_fa.append(each)
                            x_ratio.append((x_num/x_total)/each)
                    else:
                        x_ratio.append((x_num/x_total)/fa_dist[0])
                        x_fa.append(fa_dist[0])
                        
                    
#                for y in ['V1', 'V2']:
#                    if y in temp_site_dict[ko]:
#                        y_num = temp_site_dict[ko][y]['taxa_num']
#                        y_total = temp_site_dict[ko][y]['total_taxa']
#                        y_fa_dist = temp_site_dict[ko][y]['fa_dist']
#                    else:
#                        y_num = 0
#                        y_total = 1
#                        y_fa_dist = [1]
#                        
#                    y_set.append(y_num)
#                    y_t_set.append(y_total)
#                    
#                    if len(y_fa_dist) > 1:
#                        for each in y_fa_dist:
#                            y_fa.append(each)
#                            y_ratio.append((y_num/y_total)/each)
#                    else:
#                        y_ratio.append((y_num/y_total)/fa_dist[0])
#                        y_fa.append(fa_dist[0])
                        
                for z in ['P1', 'P2', 'P3']:
                    if z in temp_site_dict[ko]:
                        z_num = temp_site_dict[ko][z]['taxa_num']
                        z_total = temp_site_dict[ko][z]['total_taxa']
                        z_fa_dist = temp_site_dict[ko][z]['fa_dist']
                    else:
                        z_num = 0
                        z_total = 1
                        z_fa_dist = [1]
                    
                    z_set.append(z_num)
                    z_t_set.append(z_total)
                    
                    if len(z_fa_dist) > 1:
                        for each in z_fa_dist:
                            z_fa.append(each)
                            z_ratio.append((z_num/z_total)/each)
                    else:
                        z_ratio.append((z_num/z_total)/fa_dist[0])
                        z_fa.append(fa_dist[0])
                        
                print(x_ratio)
                #print(y_ratio)
                print(z_ratio)
                
                #_u, mwu_pval_1 = run_mwu(x_ratio, y_ratio) 
                #_u, mwu_pval_2 = run_mwu(y_ratio, z_ratio)
                _u, mwu_pval_3 = run_mwu(x_ratio, z_ratio)

                if (mwu_pval_3) <= 0.05:
#                obs = np.array([x_set, y_set, z_set])
#                chi2, pval, dof, expected = stats.chi2_contingency(obs, correction=True)
                
#                if pval <= 0.05:
#                    if uname not in versatile_dict:
#                        versatile_dict[uname] = {}
#                        
#                    versatile_dict[uname]['M_fa'] = x_fa
#                    versatile_dict[uname]['P_fa'] = y_fa
#                    versatile_dict[uname]['M_num'] = sum(x_set)
#                    versatile_dict[uname]['P_num'] = sum(y_set)
#                    versatile_dict[uname]['M_norm'] = sum(x_t_set)
#                    versatile_dict[uname]['P_norm'] = sum(y_t_set)
                                                    
                #_od, fet_pval= stats.fisher_exact([[sum(x_set), sum(x_t_set)-sum(x_set)],[sum(y_set), sum(y_t_set)-sum(y_set)]])
                #if pval <= 0.05:
                    print(ko)
                    if uname not in versatile_dict:
                        versatile_dict[uname] = {}
                        
                    #coord_dict = define_type((sum(x_set)/sum(x_t_set)), (sum(y_set)/sum(y_t_set)), x_fa, y_fa)
                        
                    versatile_dict[uname]['S_fa'] = x_fa
                    #versatile_dict[uname]['V_fa'] = y_fa
                    versatile_dict[uname]['P_fa'] = z_fa
                    versatile_dict[uname]['S_num'] = sum(x_set)
                    #versatile_dict[uname]['V_num'] = sum(y_set)
                    versatile_dict[uname]['P_num'] = sum(z_set)
                    versatile_dict[uname]['S_norm'] = sum(x_t_set)
                    #versatile_dict[uname]['V_norm'] = sum(y_t_set)
                    versatile_dict[uname]['P_norm'] = sum(z_t_set)
                    #versatile_dict[uname]['istype'] = istype
                    #versatile_dict[uname]['ratio_taxa'] = ratio_taxa
                    #versatile_dict[uname]['ratio_fa'] = ratio_fa
                    versatile_dict[uname]['S_ratio'] = x_ratio
                    #versatile_dict[uname]['V_ratio'] = y_ratio
                    versatile_dict[uname]['P_ratio'] = z_ratio                    
#                print(x_fa)
#                print(y_fa)

                    fet_ct+=1
#                    if level not in versatile_dict:
#                        versatile_dict[level] = set()
#                    versatile_dict[level].add(ko)
#                    
                #print(obs, pval, chi_ct, fet_pval, fet_ct)
                    
    pickle_name = ('plot_versatility.p').format()
    pickle.dump(versatile_dict, open(pickle_name, 'wb'))
                     
    return(versa_dict)

def make_versatility_figures():
    import plotly.express as px
    convert_rank_to_taxa = {0:'kingdom', 1:'phylum', 2:'class', 3:'order', 4:'family', 5:'genus', 6:'species'}
    
    file = open('plot_versatility.p','rb')
    versatile_dict = pickle.load(file)
    file.close()
    
    for uname in versatile_dict:
        i_tag = ('Impacted, n = {} of {}: {}. Median = {}').format(versatile_dict[uname]['M_num'], versatile_dict[uname]['M_norm'], round(versatile_dict[uname]['M_num']/versatile_dict[uname]['M_norm'],2), np.median(versatile_dict[uname]['M_fa']))
        p_tag = ('Pristine, n = {} of {}: {}. Median = {}').format(versatile_dict[uname]['P_num'], versatile_dict[uname]['P_norm'], round(versatile_dict[uname]['P_num']/versatile_dict[uname]['P_norm'],2), np.median(versatile_dict[uname]['P_fa']))
        #s_tag = ('Source').format(path)
        #v_tag = ('Valley').format(path)
        x_data = i_tag, p_tag
        
        i_list = return_log10(versatile_dict[uname]['M_fa'])
        p_list = return_log10(versatile_dict[uname]['P_fa'])
        #s_list = return_log10(plot_round[uname]['S'])
        #v_list = return_log10(plot_round[uname]['V'])
        y_data = i_list, p_list
    
        ko = uname.split('_')[1]
        
        taxa_level = convert_rank_to_taxa[int(uname.split('_')[0])]
                        
        fig = go.Figure()
        colors = 'rgba(255, 144, 14, 0.5)', 'rgba(44, 160, 101, 0.5)'

        outfile_name = ('{}_{}_versatility.pdf').format(taxa_level, ko)
        print(outfile_name)
        
        for xd, yd, cls in zip(x_data, y_data, colors):
                fig.add_trace(go.Box(
                    #,
                    y=yd,
                    name=xd,
                    boxpoints='all',
                    notched=True,
                    jitter=0.5,
                    whiskerwidth=0.2,
                    fillcolor=cls,
                    line_color=cls,
                    marker_size=5,
                    line_width=1,
                    showlegend=False)
                )
                
        fig.update_layout(
            title=ko,
            xaxis_title="Sample Site",
            yaxis_title="Log10(Relative Functional Abundance)",
                    font=dict(
                        family="Courier New, monospace",
                        size=10,
                        color="#7f7f7f"
                    )
        )
        
                
        fig.show()
        fig.write_image(outfile_name)

    
    x_list = []
    y_list = []
    size_list = []
    istype_list = []
    uname_list = []
    type_dict = {}    
    for uname in versatile_dict:
        if 'K02591' in uname:
            size_list.append(int(uname.split('_')[0])**3)
            x_list.append(np.log10(versatile_dict[uname]['ratio_taxa']))
            y_list.append(np.log10(versatile_dict[uname]['ratio_fa']))
            
            istype = versatile_dict[uname]['istype']
            istype_list.append(istype)
            uname_list.append(uname)
            
            if istype not in type_dict:
                type_dict[istype] = []
            type_dict[istype].append(uname)
        
    fig = px.scatter(x=x_list, y=y_list, color=istype_list, size=size_list)
    fig.show()
    fig.write_image('Nitrogen_example.pdf')
    
    print(type_dict)

#if args.versatility:
#    #TODO fix 'select_taxa'
#    #select_taxa_set = parse_select_taxa()
#    taxonomy_name = args.taxonomy        
#    universal_ko_dict, universal_ko_lookup = load_ko(args.kegg_file)
#    contrib_name = args.contrib_file
#    
#    pct_threshold = float(args.pct_threshold)
#    pval_threshold = float(args.pval_threshold)
#    
#    ko_super_set_type, ko_to_path_lookup = ko_to_pathway(universal_ko_dict, universal_ko_lookup)
#    
#    ordered_list = []
#    for superset, ko_list in ko_super_set_type.items():
#        for ko in ko_list:
#            ordered_list.append(ko)
#                                    
#    for taxa_cutoff_name in rank_order:
#            taxa_cutoff_num = convert_taxa_to_rank[taxa_cutoff_name]
#            
#            print('For ', taxa_cutoff_name)
#            print('... running parse_taxonomy')
#            taxa_to_asv_dict, taxa_dict, taxa_set = parse_taxonomy(taxonomy_name, taxa_cutoff_num)
#            print('... running parse_contrib')
#            ko_dict, ko_tfa_dict, pathway_totals_dict = parse_pathway_contrib(contrib_name, taxa_dict, taxa_cutoff_name, False)
#    
##    print('running parse_contrib')
##    ko_dict, ko_tfa_dict = parse_pathway_contrib(contrib_name, False)
#
#            print('running apply_versatility')
#            versa_dict = apply_versatility()
#       

if False:
    min_list = [0, 0, 0]
    max_list = [0, 0, 0]
    total_ct = 0
    min_ct = 0
    max_ct = 0
    
    infile = open('koval.tab')
    
    for line in infile:
        if 'KO' not in line:
            total_ct +=1
            ko = line.split('\t')[0]
            subl = float(line.split('\t')[1])
            inte = float(line.split('\t')[2])
            supl = float(line.split('\t')[3])

            if min([subl, inte, supl]) == subl:
                min_list[0]+=1
            
            if max([subl, inte, supl]) == subl:
                max_list[0]+=1
            
            if min([subl, inte, supl]) == inte:
                min_list[1]+=1
            
            if max([subl, inte, supl]) == inte:
                max_list[1]+=1

            if min([subl, inte, supl]) == supl:
                min_list[2]+=1
            
            if max([subl, inte, supl]) == supl:
                max_list[2]+=1
        else:
            print(line)
            
            
            