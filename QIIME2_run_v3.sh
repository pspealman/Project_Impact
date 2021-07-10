#Loading modules
#module load biom/intel/2.1.5

#starting qiime
singularity shell core_2019.4.sif

#Joining reads
	#scrits_folder contains our custom python scripts
	scripts_folder=/scratch/ps163/Dr_Carolina/Project_Impact/scripts/
	#joined_folder contains the completed joined files
	joined_folder=/scratch/ps163/Dr_Carolina/Project_Impact/combined_fastq
	#qiime_results contains qiime run objects
	qiime_results=/scratch/ps163/Dr_Carolina/Project_Impact/qiime_results/
	#Mapping File contains metadata for samples
	mapping_file=/scratch/ps163/Dr_Carolina/Project_Impact/metadata/MappingFile_Impact.csv

# #Importing 
	qiime tools import \
	 --type 'SampleData[SequencesWithQuality]' \
	 --input-path ${joined_folder}/ \
	 --input-format CasavaOneEightSingleLanePerSampleDirFmt \
	 --output-path ${qiime_results}/demux-joined.qza
#
	qiime dada2 denoise-single \
	 --i-demultiplexed-seqs ${qiime_results}/demux-joined.qza \
	 --p-trim-left 3 \
	 --p-trunc-len 0 \
	 --o-representative-sequences ${qiime_results}/rep-seqs-dada2.qza \
	 --o-table ${qiime_results}/table-dada2.qza \
	 --o-denoising-stats ${qiime_results}/stats-dada2.qza
	 
	 
	###
	# Combine
	###
	m_fastq_dir=/scratch/ps163/Project_Impact/fastq/M/
	p_fastq_dir=/scratch/ps163/Project_Impact/fastq/P/
	combined_dir=/scratch/ps163/Project_Impact/fastq/MP/
	#
	qiime feature-table merge \
	--i-tables ${p_fastq_dir}/pre-otu-table-dada2.qza \
		${m_fastq_dir}/pre-otu-table-dada2.qza \
		${qiime_results}/table-dada2.qza \
	--o-merged-table ${combined_dir}/pre-otu-table-dada2.qza \
	
	qiime feature-table merge-seqs \
	--i-data ${p_fastq_dir}/pre-otu-rep-seqs-dada2.qza \
		${m_fastq_dir}/pre-otu-rep-seqs-dada2.qza \
		${qiime_results}/stats-dada2.qza \
	--o-merged-data ${combined_dir}pre-otu-rep-seqs-dada2.qza
# #	
	# #Phylogeny
	qiime phylogeny align-to-tree-mafft-fasttree \
	  --i-sequences ${qiime_results}/rep-seqs-dada2.qza \
	  --o-alignment  ${qiime_results}/aligned-rep-seqs.qza \
	  --o-masked-alignment  ${qiime_results}/masked-aligned-rep-seqs.qza \
	  --o-tree  ${qiime_results}/unrooted-tree.qza \
	  --o-rooted-tree  ${qiime_results}/rooted-tree.qza
	  
	#to view trees 
	qiime tools export \
	  --input-path ${qiime_results}/unrooted-tree.qza \
	  --output-path ${qiime_results}/exported-unrooted-tree
	  
	#to view trees 
	qiime tools export \
	  --input-path ${qiime_results}/rooted-tree.qza \
	  --output-path ${qiime_results}/exported-rooted-tree
	
	###  
	#Stop here -
	# 1. Check min_depth and max_depth
	# 	from ${qiime_results}/feature_table.qzv
	# 2. edit variables below
	###
	
	min_depth=48600
	max_depth=106700
	
	#Alpha raref.
	qiime diversity alpha-rarefaction \
	  --i-table ${qiime_results}/table-dada2.qza \
	  --i-phylogeny ${qiime_results}/rooted-tree.qza \
	  --p-max-depth ${max_depth} \
	  --m-metadata-file ${mapping_file} \
	  --o-visualization ${qiime_results}/core-metrics-results/alpha-rarefaction.qzv
	
	#double check sampling depth
	qiime diversity core-metrics-phylogenetic \
	  --i-phylogeny ${qiime_results}/rooted-tree.qza \
	  --i-table ${qiime_results}/table-dada2.qza \
	  --p-sampling-depth ${min_depth} \
	  --m-metadata-file ${mapping_file} \
	  --output-dir ${qiime_results}core-metrics-results
	  
	#double check sampling depth
	qiime diversity alpha-group-significance \
	  --i-alpha-diversity ${qiime_results}/core-metrics-results/observed_otus_vector.qza \
	  --m-metadata-file ${mapping_file} \
	  --o-visualization ${qiime_results}/core-metrics-results/observed_otus_vector.qzv
	
	qiime diversity alpha-group-significance \
	  --i-alpha-diversity ${qiime_results}/core-metrics-results/faith_pd_vector.qza \
	  --m-metadata-file ${mapping_file} \
	  --o-visualization ${qiime_results}/core-metrics-results/faith-pd-group-significance.qzv
	
	qiime diversity alpha-group-significance \
	  --i-alpha-diversity ${qiime_results}/core-metrics-results/shannon_vector.qza \
	  --m-metadata-file ${mapping_file} \
	  --o-visualization ${qiime_results}/core-metrics-results/shannon-group-significance.qzv
	 
	 qiime diversity alpha-group-significance \
	  --i-alpha-diversity ${qiime_results}/core-metrics-results/evenness_vector.qza \
	  --m-metadata-file ${mapping_file} \
	  --o-visualization ${qiime_results}/core-metrics-results/evenness-group-significance.qzv
	
	qiime diversity beta-group-significance \
	  --i-distance-matrix  ${qiime_results}/core-metrics-results/weighted_unifrac_distance_matrix.qza \
	  --m-metadata-file ${mapping_file} \
	  --m-metadata-column Site \
	  --o-visualization ${qiime_results}/core-metrics-results/weighted-unifrac-site-significance.qzv \
	  --p-pairwise
	  
	qiime diversity beta-group-significance \
	  --i-distance-matrix  ${qiime_results}/core-metrics-results/unweighted_unifrac_distance_matrix.qza \
	  --m-metadata-file ${mapping_file} \
	  --m-metadata-column Site \
	  --o-visualization ${qiime_results}/core-metrics-results/unweighted-unifrac-site-significance.qzv \
	  --p-pairwise
	  
	qiime diversity beta-group-significance \
	  --i-distance-matrix  ${qiime_results}/core-metrics-results/bray_curtis_distance_matrix.qza \
	  --m-metadata-file ${mapping_file} \
	  --m-metadata-column Site \
	  --o-visualization ${qiime_results}/core-metrics-results/bray_curtis-site-significance.qzv \
	  --p-pairwise
	  
	qiime diversity beta-group-significance \
	  --i-distance-matrix  ${qiime_results}/core-metrics-results/jaccard_distance_matrix.qza \
	  --m-metadata-file ${mapping_file} \
	  --m-metadata-column Site \
	  --o-visualization ${qiime_results}/core-metrics-results/jaccard_site-significance.qzv \
	  --p-pairwise
	  
	qiime diversity beta-group-significance \
	  --i-distance-matrix  ${qiime_results}/core-metrics-results/jaccard_distance_matrix.qza \
	  --m-metadata-file ${mapping_file} \
	  --m-metadata-column Site \
	  --o-visualization ${qiime_results}/core-metrics-results/jaccard_site-significance.qzv \
	  --p-pairwise
	  
	qiime diversity alpha-group-significance \
	  --i-alpha-diversity ${qiime_results}/core-metrics-results/faith_pd_vector.qza \
	  --m-metadata-file ${mapping_file} \
	  --o-visualization ${qiime_results}/core-metrics-results/faith-pd-group-significance.qzv
	  
	qiime diversity alpha-group-significance \
	  --i-alpha-diversity ${qiime_results}/core-metrics-results/evenness_vector.qza \
	  --m-metadata-file ${mapping_file} \
	  --o-visualization ${qiime_results}/core-metrics-results/evenness-group-significance.qzv
	
	qiime diversity beta-group-significance \
	  --i-distance-matrix  ${qiime_results}/core-metrics-results/unweighted_unifrac_distance_matrix.qza \
	  --m-metadata-file ${mapping_file} \
	  --m-metadata-column Site \
	  --o-visualization ${qiime_results}/core-metrics-results/unweighted-unifrac-site-significance.qzv \
	  --p-pairwise

###
	#taxonomy
	qiime feature-classifier classify-sklearn \
	  --i-classifier meta_data/silva-132-99-515-806-nb-scikit0.20-classifier.qza \
	  --i-reads ${qiime_results}/rep-seqs-dada2.qza \
	  --o-classification ${qiime_results}/taxonomy.qza
	 
	qiime metadata tabulate \
	  --m-input-file ${qiime_results}/taxonomy.qza \
	  --o-visualization ${qiime_results}/taxonomy.qzv
	  
	qiime taxa barplot \
	  --i-table ${qiime_results}/table-dada2.qza \
	  --i-taxonomy ${qiime_results}/taxonomy.qza \
	  --m-metadata-file ${mapping_file} \
	  --o-visualization ${qiime_results}/taxa-bar-plots.qzv

###heirarchial
	qiime_results=/scratch/ps163/Dr_Carolina/Project_Impact/qiime_results/
	#Mapping File contains metadata for samples
	mapping_file=/scratch/ps163/Dr_Carolina/Project_Impact/metadata/MappingFile_Impact.csv
	
qiime tools import \
	--type FeatureData[Taxonomy] \
	--input-path ${qiime_results}/edited_taxonomy.tsv \
	--output-path ${qiime_results}/edited_taxonomy.qza
	
qiime taxa collapse \
  --i-table ${qiime_results}/table-dada2.qza \
  --i-taxonomy ${qiime_results}/edited_taxonomy.qza \
  --p-level 5 \
  --o-collapsed-table ${qiime_results}/table-dada2_family.qza
	
qiime gneiss correlation-clustering \
    --i-table ${qiime_results}/table-dada2_family.qza \
	--p-pseudocount 1 \
    --o-clustering ${qiime_results}/hierarchy_family.qza
	
qiime gneiss dendrogram-heatmap \
    --i-table ${qiime_results}/table-dada2_family.qza \
    --i-tree ${qiime_results}/hierarchy_family.qza \
	--m-metadata-file ${mapping_file} \
	--m-metadata-column Site \
	--o-visualization ${qiime_results}/VISUALIZATION_family.pdf

