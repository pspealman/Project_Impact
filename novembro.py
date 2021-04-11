# -*- coding: utf-8 -*-
"""
Created on Tue Oct 22 09:29:45 2019

Novembro - for identifying taxa enrichments from Qiime generated dat

python novembro.py -f feature-table.biom.txt -t taxonomy.tsv -s silva -o taxa_counts.tab

python /scratch/ps163/Dr_Carolina/scripts/novembro.py \
    -f /scratch/ps163/Dr_Carolina/Project_Impact/qiime_results/feature-table.biom.txt \
    -t /scratch/ps163/Dr_Carolina/Project_Impact/qiime_results/taxonomy.tsv \
    -s silva \
    -o /scratch/ps163/Dr_Carolina/Project_Impact/qiime_results/taxa_counts.tab

NB: feature-table.biom.txt should be generated from converting the qiime feature_table.biom file to a tsv
#   biom convert -i feature-table.biom -o feature-table.biom.txt --to-tsv

Version: Public 1.0 (I fill my day with hope and face it with joy.)
Version: Public 1.1 (Piano Harmony)
    _x_ added Kruskal-Wallis test for sets with triplicate
    _x_ retained Chi2 support for sets without triplicate
    _x_ kludge taxa_raw_dict iteration to handle non-triplicate sets
Version: Public 1.2 (Diamond Retiree)
    _x_ added unique_function
    

@author: ps163@nyu.edu
"""

import numpy as np
import scipy.stats as stats

import plotly.graph_objects as go
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-f',"--feature_table")
parser.add_argument('-t',"--taxonomy")
parser.add_argument('-s',"--taxa_source")
parser.add_argument('-o',"--output_file")
args = parser.parse_args()

parser.add_argument('-u',"--unique_sets")


feature_table_name = args.feature_table 
taxa_file_name = args.taxonomy

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

#filter_previous = []
pass_dict = {'criteria':0, 'pval':0, 'p_pval':0, 'pass_set':0, 'figure_dict':0}

simplified_enrichment = {}

def increment_pass_dict(p_1v2, p_1v3=1, p_2v3=1):
    xlist = [p_1v2, p_1v3, p_2v3]
    
    for each_p in xlist:
        if each_p <= 0.05:
            return(1)
    return(0)
    
def evaluate_the_mwu(p_val, w_val, testname, max_effect_size, min_effect_size, pval_threshold, pass_set):
    if (p_val <= pval_threshold) and (w_val >= max_effect_size or w_val <= min_effect_size):
        pass_set.append(testname)
    return(pass_set)
    
def evaluate_the_chi2(p_val, w_val, testname, max_effect_size, min_effect_size, pval_threshold, pass_set):
    if (p_val <= pval_threshold) and (w_val >= max_effect_size or w_val <= min_effect_size):
        pass_set.append(testname)
    return(pass_set)
    
def obs_counter(list_obs_ct):
    obs = 0
    for obs_ct in list_obs_ct:
        if obs_ct > 0:
            obs+=1
    
    if obs > 1:
        return(True)
    else:
        return(False)

def find_correction_value(otu_counts):
    impacted_1, impacted_2, impacted_3 = 0, 0, 0
    pristine_1, pristine_2, pristine_3 = 0, 0, 0
    
    for otu in otu_counts:
        i1, i2, i3, p1, p2, p3 = otu_counts[otu]

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

    return(impacted_cor_1, impacted_cor_2, impacted_cor_3, pristine_cor_1, pristine_cor_2, pristine_cor_3)

def build_otu_counts(feature_table_name):
#this is to build the 
    otu_file = open(feature_table_name)
    
    otu_counts = {}
    
    for line in otu_file:
        if line[0]!='#' and 'Feature ID' not in line:
            #OTU ID	M1	M2	M3	P1	P2	P3
            line = line.strip()
            otu = line.split('\t')[0].strip()# 
            #
            impacted_1 = float(line.split('\t')[1])
            impacted_2 = float(line.split('\t')[2])  
            impacted_3 = float(line.split('\t')[3])
            #
            pristine_1 = float(line.split('\t')[4])
            pristine_2 = float(line.split('\t')[5])  
            pristine_3 = float(line.split('\t')[6])
            #
            otu_counts[otu] = [impacted_1, impacted_2, impacted_3, pristine_1, pristine_2, pristine_3]
    #        
    otu_file.close()
    return(otu_counts)
    
def criteria(i_set, p_set, p_1v2, pval, taxa, pct_effect_size=0.05, pval_threshold=0.05):
    pass_set = []
    #log_fold_diff = []
    
    global pass_dict
    
    pass_dict['criteria']+=1
    #print(taxa)
    if (sum(p_set)+ sum(i_set)) >= 100:
        #if (pval <= pval_threshold):
        pass_dict['pval']+=1
        max_effect_size = (1+pct_effect_size)
        min_effect_size = (1-pct_effect_size)
        
        p_mean = np.mean(p_set)
        i_mean = np.mean(i_set)
        
        if p_mean == 0:
            p_mean = 1
        if i_mean == 0:
            i_mean = 1
                          
        w_1v2 = (p_mean/i_mean)
        
        #log_fold_diff = [w_1v2, w_1v3, w_1v4, w_2v3, w_2v4, w_3v4]
        
        pass_dict['p_pval'] += increment_pass_dict(p_1v2)
        #print((p_1v2, p_1v3, p_2v3))
        #print(w_1v2, w_1v3, w_2v3)
        #
        pass_set = evaluate_the_chi2(p_1v2, w_1v2, '1v2', max_effect_size, min_effect_size, pval_threshold, pass_set)
    
        if len(pass_set) >= 1:
            return(True)            

    return(False)
        
def test_max(set_list, max_obs):
    for each_set in set_list:
        if max(each_set) >= max_obs:
            max_obs = max(each_set)
    return(max_obs)

def return_log10(each_set, fraction_correction=False):
    new_set = []
    
    for each_obs in each_set:
        if fraction_correction:
            if each_obs < 1:
                each_obs = 0
            else:
                each_obs = np.log10(each_obs)
        else:
            if each_obs == 0:
                each_obs = 0
            else:
                each_obs = np.log10(each_obs)

        new_set.append(each_obs)
                
    return(new_set)
    
def return_deets(x_array, y_array):
    return(x_array, np.mean(x_array), np.std(x_array), y_array, obs_counter(y_array))
    
def run_kruskal(x_set, y_set):
    if sum(x_set) > 30 or sum(y_set) > 30:
        _w, p_xvy = stats.kruskal(x_set, y_set)
        return(p_xvy)
    else:
        return(1)
        
def run_chi2(x_num, x_den, y_num, y_den):
    if (x_num) > 5 or (y_num) > 5:
        obs = np.array([[max(x_num,1), x_den], [max(y_num,1), y_den]])
        chi2, pval, dof, expected = stats.chi2_contingency(obs, correction=True)
        return(pval)
    else:
        return(1)
        
def mod_null_set(isset):
    if len(isset) <1:
        isset.append(0)
    return(isset)
    
def simplify_enrichment(taxa, p12, i_med, p_med):
    global simplified_enrichment
    
    print(taxa, p12, i_med, p_med)
    simplified_enrichment[taxa] = 'complex'
    
    if p12<=0.05:
        if p_med > i_med:
            simplified_enrichment[taxa]='P_high'
        else:
            simplified_enrichment[taxa]='I_high'        
        
    return(simplified_enrichment[taxa])
    
    
def plot_top10_taxa(all_taxa_dict, prefix_name, pct_threshold=0.01):
    #taxa_dict[taxa][0] += sum(pristine)
    #taxa_dict['total'][0] += sum(pristine)
    
    logfile_name = ('taxonomic_abundance_{}_{}_{}.log').format(prefix_name, pct_threshold, taxa_cutoff_name)
    logfile = open(logfile_name, 'w')
    header = ('taxa\tsource_log10_abundance\tvalley_log10_abundance\tmangrove_log10_abundance\n')
    logfile.write(header)
    
    specific_taxa_dict = {'taxa':[],'P':[], 'I':[]}
    
    pct_i = taxa_dict['total'][0]*pct_threshold
    pct_p = taxa_dict['total'][1]*pct_threshold
    
    minimum_threshold = min([pct_i, pct_p])
    
    taxa_list = list(all_taxa_dict.keys())
    taxa_list.sort(reverse=True)
    
    for taxa in taxa_list:
        abundances = all_taxa_dict[taxa]
        if taxa != 'total' and taxa != 'observed':
            s_i, s_p = sum(abundances[0]), sum(abundances[1])
            if max([s_i, s_p]) >= minimum_threshold:
                l_i, l_p = sum(return_log10(abundances[0], True)), sum(return_log10(abundances[1], True))
                specific_taxa_dict['taxa'].append(taxa) 
                specific_taxa_dict['I'].append(l_i)
                specific_taxa_dict['P'].append(l_p)

                outline = ('{}\t{}\t{}\n').format(taxa, l_i, l_p)
                logfile.write(outline)
    
    logfile.close()        
    outfile_name = ('taxonomic_abundance_{}_{}_{}_.pdf').format(prefix_name, pct_threshold, taxa_cutoff_name)
        
    fig = go.Figure()
    fig.add_trace(go.Scatter(
        y=specific_taxa_dict['taxa'],
        x=specific_taxa_dict['P'],
        marker=dict(color='rgba(93, 164, 214, 0.5)', size=10),
        mode="markers",
        name="Pristine",
    ))
        
    fig.add_trace(go.Scatter(
        y=specific_taxa_dict['taxa'],
        x=specific_taxa_dict['I'],
        marker=dict(color='rgba(229,43,80, 0.5)', size=10),
        mode="markers",
        name="Impacted",
    ))

    fig.update_layout(title='Family Level Taxonomic Abundances',
                      xaxis_title="Log10 Taxa Abundance",
                      yaxis_title="Taxa",
                      font_size=10,
                      width=1500,
                      height=1500)
    
    fig.write_image(outfile_name)
        
for taxa_cutoff_name in rank_order:
    taxa_cutoff_num = convert_taxa_to_rank[taxa_cutoff_name]
    
    taxa_set = set()
    taxa_to_otu_dict = {}
    otu_to_taxa_dict = {}
    
    taxa_file = open(taxa_file_name)
    
    for line in taxa_file:
        if line[0]!='#':
            line = line.replace('"','')
            line = line.strip()
            otu = line.split('\t')[0]
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
                                
                if taxa not in taxa_to_otu_dict:
                    taxa_to_otu_dict[taxa] = [otu]
                else:
                    taxa_to_otu_dict[taxa].append(otu)
                    
                if otu not in otu_to_taxa_dict:
                    otu_to_taxa_dict[otu] = taxa
    
                else:
                    print('err')
                
    taxa_file.close()
    
    otu_counts = build_otu_counts(feature_table_name)
    impacted_cor_1, impacted_cor_2, impacted_cor_3, pristine_cor_1, pristine_cor_2, pristine_cor_3 = find_correction_value(otu_counts)
    
    print(impacted_cor_1, impacted_cor_2, impacted_cor_3) 
    print(pristine_cor_1, pristine_cor_2, pristine_cor_3)
    
    taxa_to_counts = {}
    for taxa, otus in taxa_to_otu_dict.items():
        i1, i2, i3 = 0, 0, 0
        p1, p2, p3 = 0, 0, 0
        
        for otu in otus:
            if otu in otu_counts:
               impacted_1, impacted_2, impacted_3, pristine_1, pristine_2, pristine_3 = otu_counts[otu]
               i1 += impacted_1
               i2 += impacted_2
               i3 += impacted_3
               
               p1 += pristine_1
               p2 += pristine_2
               p3 += pristine_3
           
        taxa_to_counts[taxa] = [i1, i2, i3, p1, p2, p3]
    
    outfile = open(args.output_file+'_counts.tsv', 'w')
    header = ('#taxa\tP1\tP2\tP3\tS1\tS2\tS3\tV1\tV2\n')
    outfile.write(header)
    
    for taxa, cts in taxa_to_counts.items():
        i1, i2, i3, p1, p2, p3 = cts
        outline = ('{taxa}\t{i1}\t{i2}\t{i3}\t{p1}\t{p2}\t{p3}\n').format(taxa=taxa, p1=p1, p2=p2, p3=p3, i1=i1, i2=i2, i3=i3)        
        outfile.write(outline)
    outfile.close()
        
    #store otu data
    taxa_dict = {'total':[0,0], 'observed':0}  
    
    taxa_raw_dict = {}
    
    for taxa, counts in taxa_to_counts.items():
        #
        impacted_1 = counts[0]
        impacted_2 = counts[1]
        impacted_3 = counts[2]
        raw_impacted = [impacted_1, impacted_2, impacted_3]
        #
        pristine_1 = counts[3]
        pristine_2 = counts[4]  
        pristine_3 = counts[5]
        raw_pristine = [pristine_1, pristine_2, pristine_3]
        #      
        if taxa not in taxa_raw_dict:
            taxa_raw_dict[taxa] = [raw_impacted, raw_pristine]
        else:
            for index in range(3):
                taxa_raw_dict[taxa][0][index] += raw_impacted[index]
                taxa_raw_dict[taxa][1][index] += raw_pristine[index]
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
            
        taxa_dict['total'][0] += sum(impacted)
        taxa_dict['total'][1] += sum(pristine)
        
        taxa_dict['observed'] += len([s for s in pristine if s != 1]) + len([s for s in impacted if s != 1])
        
    header = ('taxa\timpacted\tpristine\tvalley\n')
    
    all_outfile_name = ('all_unnormalized_taxa_abundance_{}.tab').format(taxa_cutoff_name)
    all_outfile = open(all_outfile_name, 'w')
    all_outfile.write(header)
    
    for taxa, raw_taxa_array  in taxa_raw_dict.items():
        outline = ('{}\t{}\t{}\n').format(taxa, sum(raw_taxa_array[0]), sum(raw_taxa_array[1]))
        all_outfile.write(outline)
    all_outfile.close()
    
    all_outfile_name = ('all_normalized_taxa_abundance_{}.tab').format(taxa_cutoff_name)
    all_outfile = open(all_outfile_name, 'w')
    all_outfile.write(header)
    
    for taxa, taxa_array in taxa_dict.items():
        if taxa != 'total' and taxa != 'observed':
            outline = ('{}\t{}\t{}\n').format(taxa, sum(taxa_array[0]), sum(taxa_array[1]))
            all_outfile.write(outline)
    all_outfile.close()
    
    #        
    taxa_outfile_name = ('site_specific_{}_enrichment.tab').format(taxa_cutoff_name)
    outfile = open(taxa_outfile_name, 'w')
    #taxa, uid, pval, p_1v2, p_1v3, p_1v4, p_2v3, p_2v4, p_3v4
    header = ('#taxa\tuid\tkruskal-willis\tmedian_log_10_Impacted\tmedian_log_10_Pristine\n')
    outfile.write(header)
    
    plot_top10_taxa(taxa_dict, 'normalized', 0.01)
    plot_top10_taxa(taxa_raw_dict, 'un_normalized', 0.01)
    
    uid = 0
    
    max_obs = 0
    figure_dict = {}
    
    dotplot_dict = {}

    total_i = taxa_dict['total'][0]    
    total_p = taxa_dict['total'][1]
    
    total_observed = taxa_dict['observed']
    
    unique_dict = {}
        
    for taxa in taxa_set:
        bonferroni_corrected_pvalue = 0.05
        if taxa in taxa_dict:
            #
            raw_i_array, raw_i_mean, raw_i_std, raw_i_set, raw_i_obs = return_deets(taxa_raw_dict[taxa][0], taxa_raw_dict[taxa][0])
            raw_p_array, raw_p_mean, raw_p_std, raw_p_set, raw_p_obs = return_deets(taxa_raw_dict[taxa][1], taxa_raw_dict[taxa][1])

            i_array, i_mean, i_std, i_set, i_obs = return_deets(taxa_dict[taxa][0], taxa_dict[taxa][0])
            p_array, p_mean, p_std, p_set, p_obs = return_deets(taxa_dict[taxa][1], taxa_dict[taxa][1])
            
            if p_obs or i_obs:
                if sum(raw_p_array) >= 100 and sum(raw_i_array) == 0:
                    ct = 0
                    for each in raw_p_array:
                        if each >= 10:
                            ct+=1
                            
                    if ct >= 3:
                        unique_dict[taxa]=[np.median(return_log10(raw_i_array)), np.median(return_log10(raw_p_array))] 
                                                
                if sum(raw_i_array) >= 100 and sum(raw_p_array) == 0:
                    ct = 0
                    for each in raw_i_array:
                        if each >= 10:
                            ct+=1
                            
                    if ct >= 3:
                        unique_dict[taxa]=[np.median(return_log10(raw_i_array)), np.median(return_log10(raw_p_array))] 
                                       
            #                                      
            if p_obs or i_obs:
                p_1v2 = 0
                p_sum, i_sum = sum(p_array), sum(i_array)
                p_less, i_less = (total_p - p_sum), (total_i - i_sum)
                
                obs = np.array([[max(p_sum,1), p_less], [max(i_sum,1), i_less]])
                chi2, pval, dof, expected = stats.chi2_contingency(obs, correction=True)
                
                runmode = 'kruskal-willis'
                if runmode == 'kruskal-willis':
                    p_1v2 = run_kruskal(i_set, p_set)
                else:
                    p_1v2 = run_chi2(i_sum, i_less, p_sum, p_less)
                
                pass_criteria = criteria(i_array, p_array, p_1v2,  pval, taxa, 0.2, pval_threshold=bonferroni_corrected_pvalue)
                
                if pass_criteria:
                    log_p_set = return_log10(p_array)
                    log_i_set = return_log10(i_array)
                    
                    max_obs = test_max([log_p_set, log_i_set], max_obs)
                    
                    clog_p_set = log_p_set
                    clog_i_set = log_i_set
                    
                    if (len(clog_p_set)+len(clog_i_set)) >= 3:
                        if clog_i_set:
                            med_log_i = np.median(clog_i_set)
                        else:
                            med_log_i = 0
                            
                        if clog_p_set:
                            med_log_p = np.median(clog_p_set)
                        else:
                            med_log_p = 0
                        
                        if taxa not in figure_dict:
                            pass_dict['figure_dict']+=1
                            figure_dict[taxa] = [clog_i_set, clog_p_set]
                        else:
                            print('error')
                            1/0
                        
                        is_how = simplify_enrichment(taxa, p_1v2, med_log_i, med_log_p)
                        outline = ('{taxa}\t{uid}\t{pval}\t{median_p}\t{median_i}\t{pv12}\t{is_how}\n').format(taxa=taxa, uid=uid, pval=pval, median_p=med_log_p, median_i=med_log_i, pv12=p_1v2, is_how=is_how)
                        print(outline)
                        outfile.write(outline)
                        
                        uid += 1  
                        if taxa not in dotplot_dict:
                            dotplot_dict[taxa] = [med_log_i, med_log_p]
                            
                        else:
                            print('taxa conflict', taxa)
                           
    if unique_dict:
        import plotly.graph_objects as go
        
        taxa_by_site_list = []
        site_abundance_dict = {'I':[], 'P':[]}
        
        #taxon_set = set()
            
        taxa_list = list(unique_dict.keys())
        taxa_list.sort(reverse=True)
        
        for taxa in taxa_list:
            sites = unique_dict[taxa]
            
            taxa_by_site_list.append(taxa)
            site_abundance_dict['I'].append(sites[0])
            site_abundance_dict['P'].append(sites[1])

        outfile_name = ('_unique_raw_10_{}_.pdf').format(taxa_cutoff_name)
            
        fig = go.Figure()
        fig.add_trace(go.Scatter(
            y=taxa_by_site_list,
            x=site_abundance_dict['I'],
            marker=dict(color='rgba(229,43,80, 0.5)', size=10),
            mode="markers",
            name="Impacted",
        ))
                
        fig.add_trace(go.Scatter(
            y=taxa_by_site_list,
            x=site_abundance_dict['P'],
            marker=dict(color='rgba(44, 160, 101, 0.5)', size=10),
            mode="markers",
            name="Pristine",
        ))

        fig.update_layout(title='Site specific taxonomic differentials',
                          xaxis_title="Median Log10 Taxa Abundance",
                          yaxis_title="Taxa",
                          font_size=10,
                          width=1500,
                          height=1500)
        
        fig.write_image(outfile_name)     
        
        
    if figure_dict:
        #build_x_y()

        taxa_list = list(figure_dict.keys())
        taxa_list.sort(reverse=True)
        
        for taxa in taxa_list:
            set_list = figure_dict[taxa]
            global_x_data = []
            global_y_data = []
            global_colors = []
            
            if '/' in taxa:
                taxa = taxa.replace('/','_or_')
            
            p_tag = ('{}_{}_Pristine').format(taxa_cutoff_name, taxa)
            i_tag = ('{}_{}_Impacted').format(taxa_cutoff_name, taxa)
             
            x_data = p_tag, i_tag
            
            p_set = mod_null_set(set_list[1])
            i_set = mod_null_set(set_list[0])
                    
            outfile.close()
                                                           
            y_data = p_set, i_set
            
            colors = 'rgba(44, 160, 101, 0.5)', 'rgba(229,43,80, 0.5)'
            
            global_x_data += x_data
            global_y_data += y_data
            global_colors += colors
            
            fig = go.Figure()
            
            outfile_name = ('site_specific_{}_{}_enrichment.pdf').format(taxa_cutoff_name, taxa)
            print(outfile_name)
            
            for xd, yd, cls in zip(global_x_data, global_y_data, global_colors):
                    fig.add_trace(go.Box(
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
                title=taxa,
                xaxis_title="Sample Site",
                yaxis_title="Log10(Relative Taxa Abundance)",

            )

            fig.write_image(outfile_name)
    
    if dotplot_dict:
        import plotly.graph_objects as go
        
        taxa_by_site_list = []
        site_abundance_dict = {'I':[], 'P':[]}
        
        taxa_list = list(dotplot_dict.keys())
        taxa_list.sort(reverse=True)
            
        for taxa in taxa_list:
            sites = dotplot_dict[taxa]
            taxa_by_site_list.append(taxa)
            site_abundance_dict['I'].append(sites[0])
            site_abundance_dict['P'].append(sites[1])

        outfile_name = ('_site_specific_spread_{}_.pdf').format(taxa_cutoff_name)
            
        fig = go.Figure()
        fig.add_trace(go.Scatter(
            y=taxa_by_site_list,
            x=site_abundance_dict['I'],
            marker=dict(color='rgba(229,43,80, 0.5)', size=10),
            mode="markers",
            name="Impacted",
        ))
                
        fig.add_trace(go.Scatter(
            y=taxa_by_site_list,
            x=site_abundance_dict['P'],
            marker=dict(color='rgba(44, 160, 101, 0.5)', size=10),
            mode="markers",
            name="Pristine",
        ))

        fig.update_layout(title='Site specific taxonomic differentials',
                          xaxis_title="Median Log10 Taxa Abundance",
                          yaxis_title="Taxa",
                          font_size=10,
                          width=1500,
                          height=1500)
        
        fig.write_image(outfile_name)     