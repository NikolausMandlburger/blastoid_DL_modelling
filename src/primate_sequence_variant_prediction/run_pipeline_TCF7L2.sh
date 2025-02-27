
#PARAMETERS
#constant paramenter and scripts
HG38=/groups/stark/nikolaus.mandlburger/Projects/blastoid_project/blastoid_DL_modelling/ext_data/reference_genomes/hg38/hg38.chrom.sizes
GTF_ANNOT_FILE=/groups/stark/nikolaus.mandlburger/Projects/blastoid_project/blastoid_DL_modelling/ext_data/reference_genomes/hg38/gencode.v28.primary_assembly.annotation.gtf.gz
EPI=/groups/stark/nikolaus.mandlburger/Projects/blastoid_project/blastoid_DL_modelling/res/pseudobulks/broadct_epi/epiblast_96h_rpm_norm.bw
HYPO=/groups/stark/nikolaus.mandlburger/Projects/blastoid_project/blastoid_DL_modelling/res/pseudobulks/broadct_hypo/hypoblast_96h_rpm_norm.bw
TE=/groups/stark/nikolaus.mandlburger/Projects/blastoid_project/blastoid_DL_modelling/res/pseudobulks/broadct_TE/TE_96h_rpm_norm.bw
pTE=/groups/stark/nikolaus.mandlburger/Projects/blastoid_project/blastoid_DL_modelling/res/pseudobulks/pTE_refined/pTE_96h_rpm_norm.bw #The refined pTE track, slightly different from the track used for training
pTE_FRAGMETSTSV=/groups/stark/nikolaus.mandlburger/Projects/blastoid_project/blastoid_DL_modelling/data/pTE_refined/fragments.tsv.gz #The refined pTE track, slightly different from the track used for training
PEAKCALL_Q=0.4
SELECT_TE_SPEC_PEAKS_SCRIPT=/groups/stark/nikolaus.mandlburger/Projects/blastoid_project/blastoid_DL_modelling/src/primate_sequence_variant_prediction/select_TEspec_peaks.py
FC_THRESHOLD=1.7
TE_MIN_ACC_THRESHOLD=0.5
SPECPEAKS_TO_FASTA_SCRIPT=/groups/stark/nikolaus.mandlburger/Projects/blastoid_project/blastoid_DL_modelling/src/primate_sequence_variant_prediction/filtered_peaks_to_fasta_summitcentered.py
RETRIEVAL_FLANK=500
PREDICTIONS_SCRIPT=/groups/stark/nikolaus.mandlburger/Projects/blastoid_project/blastoid_DL_modelling/src/primate_sequence_variant_prediction/make_candidate_predictions_with_window_shifts.py
MODELDIR=/groups/stark/nikolaus.mandlburger/Projects/blastoid_project/blastoid_DL_modelling/res/model_training/pTE/model_training/trained_models/
OUT_SUPERDIR=/groups/stark/nikolaus.mandlburger/Projects/blastoid_project/blastoid_DL_modelling/res/primate_sequence_variant_prediction


#target gene specific parameters
REGION=chr10:112890945:113301222
TARGETGENE=TCF7L2

CHROM=$(echo $REGION | cut -d':' -f1)
START=$(echo $REGION | cut -d':' -f2)
END=$(echo $REGION | cut -d':' -f3)


#dynamically generated parameters
OUTD=${OUT_SUPERDIR}/${TARGETGENE}/
mkdir -p $OUTD

pTE_PEAKS_IN_REGION=${OUTD}/pTE_${TARGETGENE}_${REGION//:/_}.narrowPeak
pTE_PEAKS_IN_REGION_TE_SPEC=${OUTD}/pTE_${TARGETGENE}_${REGION//:/_}_TEspec.narrowPeak
pTE_NARROW_PEAKS=${OUT_SUPERDIR}/pTE_peaks_Q_${PEAKCALL_Q}/pTE_peaks_Q_${PEAKCALL_Q}_peaks.narrowPeak
pTE_FRAGMET_BED=$OUT_SUPERDIR/pTE_fragments.bed

#UNCOMMENT TO PERFORM PEAK CALLING
#call peaks on whole pTE frags 
#convert fragments to bed
zcat $pTE_FRAGMETSTSV | awk 'BEGIN{OFS="\t"} $1 ~ /^chr/ {print $1,$2,$3}' | sort -k 1,1 -k2,2n -k3,3n > $pTE_FRAGMET_BED

#call peaks on whole pTE
module load build-env/f2022
module load macs3/3.0.1-foss-2022b
macs3 callpeak -f BEDPE -t $pTE_FRAGMET_BED -g hs -n pTE_peaks_Q_${PEAKCALL_Q} -q $PEAKCALL_Q --keep-dup all --outdir ${OUT_SUPERDIR}/pTE_peaks_Q_${PEAKCALL_Q}
module unload macs3/3.0.1-foss-2022b

#select peaks in region
ROI=$(printf "%s\t%s\t%s" ${REGION//:/ })
module load build-env/2020
module load bedtools/2.27.1-foss-2018b
bedtools intersect -a $pTE_NARROW_PEAKS -b <(echo -e "$ROI") > $pTE_PEAKS_IN_REGION

#select peaks with TE specific accessibility
$SELECT_TE_SPEC_PEAKS_SCRIPT -p $pTE_PEAKS_IN_REGION -g $GTF_ANNOT_FILE -t $pTE -e $EPI -y $HYPO -f $FC_THRESHOLD -m $TE_MIN_ACC_THRESHOLD

#retrieve sequence alignments for all TE specific peaks and process tp fasta 
for WINDOW_SHIFT in -30 -20 -10 0 10 20 30;do
    $SPECPEAKS_TO_FASTA_SCRIPT -p $pTE_PEAKS_IN_REGION_TE_SPEC -s $WINDOW_SHIFT -e $RETRIEVAL_FLANK;
done

#run predictions on candidates 
for CAND_FA in $OUTD/*shift_0.fa;do 
    BASENAME=$(basename $CAND_FA);
    CANDIDATE_NAME=${BASENAME%_shift_0.fa};
    PREFIX=$(echo $CANDIDATE_NAME | cut -d'_' -f1-3);
    echo $PREFIX;
    $MYBSUB -m 30 -n preds_${PREFIX} -o ${OUTD}/log -T 00:15:00  "$PREDICTIONS_SCRIPT --chromosome $CHROM --model_dir $MODELDIR --res_dir $OUTD --ROI_name $PREFIX"
done