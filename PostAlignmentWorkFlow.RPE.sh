#!/bin/bash

hg_ref_primary="/czlab/References/GRCh38/Homo_sapiens_assembly38.fasta";
hg_ref_decoy="/czlab/References/GRCh38/GRCh38_plus_decoy/Homo_sapiens_assembly38.fasta";
hg_ref_sponge='/czlab/References/GRCh38/primary_assembly_plus_decoy_sponge/Homo_sapiens_assembly38.plus_sponge.fa'
hg_ref=${hg_ref_sponge};

sample_name=$1;
aligned_bam=${sample_name}.paired.bam
tmp_dir=/cluster/singlecellcenter/${sample_name};


if [[ ! -e ${aligned_bam} ]]; then
return;
fi

jobScript=${sample_name}.create.sbatch

echo "#!/bin/bash"
echo "#SBATCH -c 2"
echo "#SBATCH --mem=16G"
echo "#SBATCH -o ${sample_name}_%j.out"
echo "#SBATCH -e ${sample_name}_%j.err"

echo "#Step 0: Clear directory for temporary files"
echo "mkdir -p ${tmp_dir};"
echo "rm -rf ${tmp_dir}/*;"

#Step 1: Sort aligned bam by read name
echo "#Step 1: Sort aligned bam by read name"
sh /czlab/SequenceProcessing/CreatePicardOrGATKJobs.sh "picard" "SortSam" "TMP_DIR=${tmp_dir} INPUT=${sample_name}.paired.bam OUTPUT=${sample_name}.aligned.bam VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate" ${sample_name}

#Step 2: Mark duplicates
echo "#Step 2: Mark duplicates for primary reads"
if [[ ! -e ${sample_name}.duplicate_metrics ]] ; then
sh /czlab/SequenceProcessing/CreatePicardOrGATKJobs.sh "picard" "MarkDuplicates" "TMP_DIR=${tmp_dir} INPUT=${sample_name}.aligned.bam OUTPUT=${sample_name}.primary-deduped.bam VALIDATION_STRINGENCY=SILENT MAX_RECORDS_IN_RAM=1000000 MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 SORTING_COLLECTION_SIZE_RATIO=0.15 METRICS_FILE=${sample_name}.duplicate_metrics ASSUME_SORTED=true" ${sample_name} "12g" 
fi

#Step 3: Sort primary-deduped bam by read name
echo "#Step 3: Sort primary-deduped bam by read name"
sh /czlab/SequenceProcessing/CreatePicardOrGATKJobs.sh "picard" "SortSam" "TMP_DIR=${tmp_dir} INPUT=${sample_name}.primary-deduped.bam OUTPUT=${sample_name}.primary-deduped.RNsorted.bam VALIDATION_STRINGENCY=SILENT SORT_ORDER=queryname" ${sample_name}

#Step 4: Collect discordant reads	
echo "#Step 5: Get discordant reads"
sh /czlab/SequenceProcessing/CreatePicardOrGATKJobs.sh "local" "GetDiscordantReadsFromQuerynameSortedBam" "TMP_DIR=${tmp_dir} INPUT=${sample_name}.primary-deduped.RNsorted.bam EXPECTED_ORIENTATION=FR MAX_INSERT=2000 DISABLE_INSERT_INFERENCE=false SORT_ORDER=queryname IS_SAMPLE=1000000  MINIMUM_MAPPING_QUALITY=0 MAXIMUM_MISMATCHES=10 OUTPUT_DIR=`pwd`" ${sample_name}

#Step 5: Mark secondary/supplementary duplicates	
echo "#Step 4: Mark secondary/supplementary duplicates"
sh /czlab/SequenceProcessing/CreatePicardOrGATKJobs.sh "local" "MarkDuplicatesForNonPrimary" "TMP_DIR=${tmp_dir} INPUT=${sample_name}.primary-deduped.RNsorted.bam OUTPUT=${sample_name}.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT  CREATE_INDEX=true" ${sample_name} "10g"

#Step 6: Collect read depths
echo "#Step 6: Collect read depths"
run_command=`sh /czlab/SequenceProcessing/CollectReadCounts.sh "${sample_name}.bam" | sed "s:java:java -Djava.io.tmpdir=${tmp_dir}:"`	
sh /czlab/SequenceProcessing/CreatePicardOrGATKJobs.sh "command" "CollectRC" "${run_command}" ${sample_name}  

#Step 7: Collect Allelic Counts
echo "#Step 7: Collect Allelic Counts"
run_command=`sh /czlab/SequenceProcessing/CollectAllelicDepths.RPE.sh "${sample_name}.bam" | sed "s:java:java -Djava.io.tmpdir=${tmp_dir}:"`
sh /czlab/SequenceProcessing/CreatePicardOrGATKJobs.sh "command" "CollectAD" "${run_command}" ${sample_name}  

#Step 8: Clear up
echo "#Step 8: Clear up"
echo "rm -rf ${tmp_dir};"

