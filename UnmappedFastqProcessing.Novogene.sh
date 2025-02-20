#!/bin/bash
hg_ref_primary="/czlab/References/GRCh38/Homo_sapiens_assembly38.fasta";
hg_ref_decoy="/czlab/References/GRCh38/GRCh38_plus_decoy/Homo_sapiens_assembly38.fasta";
hg_ref_sponge='/czlab/References/GRCh38/primary_assembly_plus_decoy_sponge/Homo_sapiens_assembly38.plus_sponge.fa'
hg_ref='/czlab/References/GRCh38/primary_assembly_plus_sponge/hg38.primary_assembly_plus_sponge.fa';
hg_ref='/czlab/References/GRCh38/hg38_primary.fa'
picard_dir='/czlab/Software/picard/current/bin'

input_dir=$1;
# for BWA only
THREADS=24
clear_up_script='clearup.sbatch';
	echo "#!/bin/bash" > ${clear_up_script}
while read line
do	
	SM_alias=`echo ${line} | sed "s:^\(\S*\).*:\1:"` ;
	Sample=`echo ${line} | sed "s:${SM_alias}\s*::" | sed "s:\(\S*\).*:\1:"`;	
	RD1=`find -L ${input_dir} -name "${SM_alias}_*1.fq.gz"`;
	RD2=`echo $RD1 | sed "s:1.fq.gz:2.fq.gz:"`;
	Rname=`zless ${RD1} | head -n1`;
	RN=`zless ${RD1} | head -n1 | sed "s:^@::" | sed "s:SRR\S*::" | sed "s:\s*length.*::" | cut -d' ' -f1`
	BC=`zless ${RD1} | head -n1 | sed "s:^@::" | sed "s:SRR\S*::" | sed "s:\s*length.*::" | cut -d' ' -f2 | cut -d':' -f4 | sed "s:+:_:"`
	PU=$(echo ${RN} | cut -d':' -f1)	
	FC=$(echo ${RN} | cut -d':' -f3)
	LN=$(echo ${RN} | cut -d':' -f4)	
	RG=`echo ${PU}.${FC} | sed "s:..[X-Y][X-Y]::"`.${LN}
	PL="NovaSeq"
	runDate="2023-01-10"
	CN="Novogene"
	SM=${Sample}
	LB=${SM}.${BC};
	RGinfo=`echo "ID:${RG}	PL:${PL}	PU:${PU}	LB:${LB}	DT:${runDate}	SM:${SM}	CN:${CN}"`
	sample_name=${LB}.${RG};
	unmapped_bam_file=${sample_name}.unmapped.bam
	
	# CreateUnmappedBam
	tmp_dir=/cluster/singlecellcenter/${sample_name}
	jobScript=${sample_name}.createUnmapped.sbatch
	echo "#!/bin/bash" > ${jobScript}
	echo "#SBATCH -c 2" >> ${jobScript}
	echo "#SBATCH --mem=16G" >> ${jobScript}
	echo "#SBATCH -o ${sample_name}_%j.out" >> ${jobScript}
	echo "#SBATCH -e ${sample_name}_%j.err" >> ${jobScript}

	echo "#Step 0: Clear directory for temporary files" >> ${jobScript}
	echo "mkdir -p ${tmp_dir};" >> ${jobScript}
	echo "rm -rf ${tmp_dir}/*;" >> ${jobScript}
	echo "rm -rf ${tmp_dir};" >> ${clear_up_script}
	echo "#Step 1: Create unmapped bam from fastq files" >> ${jobScript}
	echo "java -Xmx10000m -jar ${picard_dir}/picard.jar FastqToSam TMP_DIR=${tmp_dir} VALIDATION_STRINGENCY=SILENT MAX_RECORDS_IN_RAM=1000000 FASTQ=$RD1 FASTQ2=$RD2 OUTPUT=`pwd`/${sample_name}.unmapped.bam READ_GROUP_NAME=$RG SAMPLE_NAME=$SM LIBRARY_NAME=$LB PLATFORM_UNIT=$PU PLATFORM=$PL SEQUENCING_CENTER=${CN} SORT_ORDER=unsorted &>> ${sample_name}.createUnmapped.log" >> ${jobScript}

	# bwa alignment
	jobScript=${sample_name}.bwa_aln.sbatch
	echo "#!/bin/bash" >> ${jobScript}
	echo "#SBATCH --cpus-per-task=${THREADS}" >> ${jobScript}        # cpu-cores per task (>1 if multi-threaded tasks)
	echo "#SBATCH --mem-per-cpu=4G" >> ${jobScript} 
	echo "#SBATCH -o ${sample_name}_bwa-mem_%j.out" >> ${jobScript} 
	echo "#SBATCH -e ${sample_name}_bwa-mem_%j.err">> ${jobScript} 

	echo "module load bwa/0.7.17-r1188">> ${jobScript} 
	echo "module load samtools/1.9.1">> ${jobScript} 
	BWA="bwa">> ${jobScript} 
	run_Samtools="samtools";>> ${jobScript} 

	bwa_aln="${BWA} mem -t ${THREADS} ${hg_ref} -Y -M">> ${jobScript} 
	
	echo "#Step 1.1: Align ${sample_name}/${sample_name}.1.fastq.gz and ${sample_name}/${sample_name}.2.fastq.gz with bwa mem"	>> ${jobScript} 
	echo "${bwa_aln} ${RD1} ${RD2} 2>> ${sample_name}.bwa_mem.log | $run_Samtools view -b - > ${sample_name}.aligned.bam 2>> ${sample_name}.bwa_mem.log" >> ${jobScript} 

	# Add alignment data to unmapped bam
	jobScript=${sample_name}.HitsPairing.sbatch
	echo "#!/bin/bash" > ${jobScript}
	echo "#SBATCH -c 2" >> ${jobScript}
	echo "#SBATCH --mem=16G" >> ${jobScript}
	echo "#SBATCH -o ${sample_name}_%j.out" >> ${jobScript}
	echo "#SBATCH -e ${sample_name}_%j.err" >> ${jobScript}
	echo "#Step 2: Merge aligned bam with unmapped bam and match paired alignments" >> ${jobScript}
	echo "touch ${sample_name}.paired.bam" >> ${jobScript}		
	sh /czlab/SequenceProcessing/CreatePicardOrGATKJobs.sh "local" "HitsPairing" "TMP_DIR=${tmp_dir} INPUT=${sample_name}.aligned.bam OUTPUT=${sample_name} MATE_PAIR=false DISABLE_INSERT_INFERENCE=false SORT_ORDER=unsorted CHIMERA_KB_MIN=2000 IF_PROPER_FR=true IF_PROPER_RF=false IF_PROPER_FF=false VALIDATION_STRINGENCY=SILENT PREFERRED_ORIENTATION=FR UNMAPPED=${sample_name}.unmapped.bam ADD_RG=true" ${sample_name} "16g" >> ${jobScript}
	done < BarcodeMap.txt
exit

