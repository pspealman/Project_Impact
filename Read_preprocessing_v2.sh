#Loading modules
#module load biom/intel/2.1.5

#starting qiime
singularity shell core_2019.4.sif

#Joining reads
	#scrits_folder contains our custom python scripts
	scripts_folder=/scratch/ps163/Dr_Carolina/Project_Impact/scripts/
	#fastq_folder contains forward and reverse reads:
	#fastq_folder=/scratch/ps163/Dr_Carolina/UFBA_Carolina_Jan_2019-115734619/
	#vsearch_folder contains the results of qiime2 vsearch join
	#vsearch_folder=/scratch/ps163/Dr_Carolina/qiime2/vsearch_join_fastq/
	
	#joined_folder contains the completed joined files
	joined_folder=/scratch/ps163/Dr_Carolina/Project_Impact/joined_fastq/
	#combined_folder contains the completed joined files
	combined_folder=/scratch/ps163/Dr_Carolina/Project_Impact_v2/combined_fastq/
	
	
	
	
	#combine the previously combined reads (see Prev_Read_preprocessing.sh)
	rep_name=P1
		name_1=P1-1
		name_2=P1-2
		name_3=P1-3
		cat ${joined_folder}/${name_1}_0_L001_R1_001.fastq ${joined_folder}/${name_2}_0_L001_R1_001.fastq ${joined_folder}/${name_3}_0_L001_R1_001.fastq > ${combined_folder}/${rep_name}_0_L001_R1_001.fastq 

	rep_name=P2
		name_1=P2-1
		name_2=P2-2
		name_3=P2-3
		cat ${joined_folder}/${name_1}_0_L001_R1_001.fastq ${joined_folder}/${name_2}_0_L001_R1_001.fastq ${joined_folder}/${name_3}_0_L001_R1_001.fastq > ${combined_folder}/${rep_name}_0_L001_R1_001.fastq 
		
	rep_name=P3
		name_1=P3-1
		name_2=P3-2
		name_3=P3-3
		cat ${joined_folder}/${name_1}_0_L001_R1_001.fastq ${joined_folder}/${name_2}_0_L001_R1_001.fastq ${joined_folder}/${name_3}_0_L001_R1_001.fastq > ${combined_folder}/${rep_name}_0_L001_R1_001.fastq 