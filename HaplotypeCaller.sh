#!/bin/bash

# Binaries
picard_dir='/czlab/chzhang/github/picard-public/dist/';
picard='/czlab/chzhang/github/picard-public/dist/picard.jar';

GATK4_dir='/czlab/chzhang/github/gatk4/build/libs';
GATK4=${GATK4_dir}/gatk.jar;

GATK='/czlab/Software/GATK/current/GenomeAnalysisTK.jar'

# References
ref_dir='/czlab/References/GRCh38/'

# GRCh38 reference
reference=${ref_dir}/Homo_sapiens_assembly38.fasta; 
ref="-R ${reference}";
scatter_intervals_dir=${ref_dir}/hg38_scatter_300

hg_ref_primary="/czlab/References/GRCh38/Homo_sapiens_assembly38.fasta";
hg_ref_decoy="/czlab/References/GRCh38/GRCh38_plus_decoy/Homo_sapiens_assembly38.fasta";
hg_ref_sponge='/czlab/References/GRCh38/primary_assembly_plus_decoy_sponge/Homo_sapiens_assembly38.plus_sponge.fa'

# Standard resources
resource_dir='/czlab/References/GRCh38/Broad';
hapmap_reference=${resource_dir}/resources/hapmap_3.3.hg38.vcf.gz;
hapmapBA_reference=${resource_dir}/resources/hapmap_3.3.hg38.BA.vcf.gz;
omni_reference=${resource_dir}/resources/1000G_omni2.5.hg38.vcf.gz;
kg_snp_reference=${resource_dir}/resources/1000G_phase1.snps.high_confidence.hg38.vcf.gz;
kgBA_snp_reference=${resource_dir}/resources/1000G_phase1.snps.high_confidence.BA.hg38.vcf.gz;
kg_indel_reference=${resource_dir}/resources/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz;
dbsnp_reference=${resource_dir}/Homo_sapiens_assembly38.dbsnp138.vcf;
dbsnp_reference=${resource_dir}/resources/dbsnp_146.hg38.vcf.gz;
kg_snp_reference='/czlab/References/GRCh38/Broad/1000G_phase3_v4_20130502.snp.1percentAF.sites.hg38.vcf.gz';
#RPE1_hets='/singlecellcenter/RPE-1/RPE-1_Genotyping_GRCh38/RPE-1.hets.vcf.gz';

# Input 
NAME=$1
INPUT=`find . -name "${NAME}*.bam" -print0 | sed "s:./: -I ./:g"`

tmp=/cluster/singlecellcenter/${NAME};

# Standard engine parameters
run_GATK="java -Xmx8000m -Djava.io.tmpdir=${tmp} -jar ${GATK}"
GATK_rf="-rf DuplicateRead -rf FailsVendorQualityCheck -rf NotPrimaryAlignment -drf BadMate -rf MappingQualityUnavailable -rf UnmappedRead -rf BadCigar -maxInsert 1000 -rf MateSameStrand";

run_GATK4="java -Xmx8000m -Djava.io.tmpdir=${tmp} -jar ${GATK4}"
GATK4_standard_filter="--read-filter PassesVendorQualityCheckReadFilter --read-filter HasReadGroupReadFilter \
 --read-filter NotDuplicateReadFilter --read-filter MappingQualityAvailableReadFilter"
ProperPair_filter="--read-filter PairedReadFilter \
	--read-filter FragmentLengthReadFilter --max-fragment-length 1000 --read-filter MateOnSameContigOrNoMappedMateReadFilter --read-filter MateDifferentStrandReadFilter" 

# Variant discovery parameters
GATK4_HC_rf="--read-filter MappingQualityReadFilter --minimum-mapping-quality 30 --read-filter OverclippedReadFilter --filter-too-short 25 --read-filter GoodCigarReadFilter --read-filter AmbiguousBaseReadFilter "	
GATK4_HC_Discovery="${run_GATK4} HaplotypeCaller --genotyping-mode DISCOVERY --disable-tool-default-annotations true --annotation DepthPerAlleleBySample --annotation Coverage ${ref} ${GATK4_standard_filter} ${GATK4_HC_rf} \
	--native-pair-hmm-threads 2 --output-mode EMIT_ALL_SITES  --seconds-between-progress-updates 100"

mkdir -p ${NAME}/sbatch
output_dir="`pwd`/${NAME}";
# Commands
for interval in ${scatter_intervals_dir}/*.interval_list
do
	id=`echo ${interval} | sed "s:_of.*interval_list::" | sed "s:.*interval_0*::" `
	OUTPUT=sbatch/${NAME}.gatk4_HC.${id}.sbatch
	echo "${output_dir}/${OUTPUT}" >> ${NAME}.sbatch_jobs.txt;
	echo "#!/bin/bash" > ${output_dir}/${OUTPUT}
	echo "#SBATCH -c 2" >> ${output_dir}/${OUTPUT}
	echo "#SBATCH --mem=10G" >> ${output_dir}/${OUTPUT}
	echo "#SBATCH -e ${output_dir}/${OUTPUT}_%j.out" >> ${output_dir}/${OUTPUT}
	echo "#SBATCH -o ${output_dir}/${OUTPUT}_%j.err" >> ${output_dir}/${OUTPUT}
	echo "mkdir -p ${tmp}" >> ${output_dir}/${OUTPUT}
	echo "${GATK4_HC_Discovery} ${INPUT} -L ${interval} -O ${output_dir}/${NAME}.gatk4_HC.${id}.vcf.gz" >> ${output_dir}/${OUTPUT};
done
vcf_count=`find ${scatter_intervals_dir} -name "*.interval_list" -print | wc -l`;
for id in $(seq 1 ${vcf_count}) ;
do
    echo "${output_dir}/${NAME}.gatk4_HC.${id}.vcf.gz" >> ${NAME}.vcf_list.txt ;
done

OUTPUT=${NAME}.postHC.sbatch
echo "#!/bin/bash" > ${OUTPUT}
echo "#SBATCH -c 2" >> ${OUTPUT}
echo "#SBATCH --mem=10G" >> ${OUTPUT}
echo "#SBATCH -e ${OUTPUT}_%j.out" >> ${OUTPUT}
echo "#SBATCH -o ${OUTPUT}_%j.err" >> ${OUTPUT}
echo "rmdir -rf ${tmp}" >> ${OUTPUT}
echo "/czlab/bin/bcftools concat -f ${NAME}.vcf_list.txt -O z -o ${NAME}.HC_all.vcf.gz" >> ${OUTPUT}
echo "/czlab/bin/tabix ${NAME}.HC_all.vcf.gz" >> ${OUTPUT}
exit

