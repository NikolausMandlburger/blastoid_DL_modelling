#!/groups/stark/nikolaus.mandlburger/.conda/envs/genomics_python_2/bin/python3.12

#Use e.g. env /groups/stark/nikolaus.mandlburger/Projects/blastoid_project/blastoid_DL_modelling/src/envs/genomics_python_env.yml

import argparse 
import pyBigWig
import pandas as pd
import numpy as np
import os
import bioframe as bf


#generate a parser
parser = argparse.ArgumentParser(prog='selectTEspec_peaks', description="Among a set of peaks, select trophectoderm specific ones")

#add arguments
parser.add_argument("-p","--narrow_peaks_file",type=str, help="Path to the narrow peaks file thatshould be filtered (should cover a region of interest)")
parser.add_argument("-g","--gtf_file",type=str, help="Path to genome annotation gtf from which promotor regions are extracted")
parser.add_argument("-t","--TE_bigwig_file",type=str, help="Path to the bigwig file for TE")
parser.add_argument("-e","--epi_bigwig_file",type=str, help="Path to the bigwig file for epi")
parser.add_argument("-y","--hypo_bigwig_file",type=str, help="Path to the bigwig file for hypo")
parser.add_argument("-f","--foldchange_threshold",type=float, help="Threshold for foldchange according to which peaks are selected")
parser.add_argument("-m","--TE_min_accessibility",type=float, help="Threshold for minimum accessibility for TE specific peaks")


                                
#parse and retrieve arguments
args=parser.parse_args()
                                
#access and check individual arguments an assign to variables

if args.narrow_peaks_file!=None:
    narrow_peaks_file = args.narrow_peaks_file
else:
    print("No narrow peaks file given")

if args.gtf_file!=None:
    gtf_file = args.gtf_file
else:
    print("No gtf_file file given, setting to default (/groups/stark/annotations/hg38/gencode.v28.primary_assembly.annotation.gtf.gz)")
    gtf_file = "/groups/stark/annotations/hg38/gencode.v28.primary_assembly.annotation.gtf.gz"

if args.TE_bigwig_file!=None:
    TE_bigwig_file = args.TE_bigwig_file
else:
    print("No TE bigwig file given")

if args.epi_bigwig_file!=None:
    epi_bigwig_file = args.epi_bigwig_file
else:
    print("No epi bigwig file given")

if args.hypo_bigwig_file!=None:
    hypo_bigwig_file = args.hypo_bigwig_file
else:
    print("No hypo bigwig file given")

if args.foldchange_threshold!=None:
    foldchange_threshold = args.foldchange_threshold
else:
    foldchange_threshold = 3
    print("No foldchange threshold given, setting to default (3)")

if args.TE_min_accessibility!=None:
    TE_min_acc = args.TE_min_accessibility
else:
    TE_min_acc = 0.1
    print("No minimum accessibility threshold given, setting to default (0.1)")


#obtain mean score in bigwig 
def fetch_interval_bigwig_score(intervals, bigwig_file,score_type = "mean"):
    bw = pyBigWig.open(bigwig_file)
    interval_scores = []
    for c,s,e in zip(intervals["chrom"],intervals["start"], intervals["end"]):
        interval_score = bw.stats(c, s, e, type=score_type)[0]
        interval_scores.append(interval_score)
    bw.close()
    interval_scores = np.array(interval_scores)
    return interval_scores


#load peaks in region of interest
ROI_narrow_peaks = bf.read_table(
    narrow_peaks_file,
    schema="narrowPeak",
)

#load promotor locations
gtf = pd.read_csv(gtf_file, 
                  sep='\t', 
                  comment='#', 
                  header=None, 
                  names=["seqnames", "source", "type", "start", "end", "score", "strand", "frame", "attribute"])

tss = gtf[gtf['type'] == 'transcript']
tss["TSS"]=tss["start"]
tss.loc[tss["strand"]=="-", "TSS"] = tss["end"]
tss["start_prom"] = tss["TSS"]-100
tss["end_prom"] = tss["TSS"]+100

promotors = tss.loc[:,["seqnames","start_prom","end_prom"]].copy()
promotors.columns = ["chrom","start","end"]

#select peaks that do not overlap with promotors
ROI_narrow_peaks_wo_promotor = bf.setdiff(ROI_narrow_peaks, promotors)

#obtain mean accessibility values in the intervals for TE, epi and hypo

epi_scores =  fetch_interval_bigwig_score(ROI_narrow_peaks_wo_promotor,epi_bigwig_file)
hypo_scores =  fetch_interval_bigwig_score(ROI_narrow_peaks_wo_promotor,hypo_bigwig_file)
TE_scores =  fetch_interval_bigwig_score(ROI_narrow_peaks_wo_promotor,TE_bigwig_file)

ROI_narrow_peaks_wo_promotor.loc[:,"epiblast_ATAC_rpmnorm"] = epi_scores
ROI_narrow_peaks_wo_promotor.loc[:,"hypoblast_ATAC_rpmnorm"] = hypo_scores
ROI_narrow_peaks_wo_promotor.loc[:,"TE_ATAC_rpmnorm"] = TE_scores
ROI_narrow_peaks_wo_promotor.fillna(0, inplace=True)

#calculate mean between epi and hypo scores
ROI_narrow_peaks_wo_promotor.loc[:,"mean_ATAC_rpmnorm_epi_hypo"] = ROI_narrow_peaks_wo_promotor.loc[:,["epiblast_ATAC_rpmnorm", "hypoblast_ATAC_rpmnorm"]].mean(axis=1)
ROI_narrow_peaks_wo_promotor.loc[:,"TE_vs_epihypo_foldchange"] =  ROI_narrow_peaks_wo_promotor["TE_ATAC_rpmnorm"]/ROI_narrow_peaks_wo_promotor["mean_ATAC_rpmnorm_epi_hypo"]

#select peaks that exceed the FC threshold
ROI_narrow_peaks_wo_promotor_highFC = ROI_narrow_peaks_wo_promotor.loc[(ROI_narrow_peaks_wo_promotor["TE_ATAC_rpmnorm"]>TE_min_acc) & (ROI_narrow_peaks_wo_promotor["TE_vs_epihypo_foldchange"]>foldchange_threshold)]

#save highFC peaks
savefile= narrow_peaks_file.replace(".narrowPeak", "_TEspec.narrowPeak")
ROI_narrow_peaks_wo_promotor_highFC.iloc[:,:10].to_csv(savefile, header= False, index = False, sep = "\t")

#usage:
#SCRIPT=/groups/stark/nikolaus.mandlburger/Projects/blastoid_project/src/primate_sequence_variants/enhancer_candidate_selection/select_TEspec_peaks.py
#PEAKS=/groups/stark/nikolaus.mandlburger/Projects/blastoid_project/res/primate_sequence_variants/enhancer_candidate_selection/GCM1/pTE_GCM1_chr6_53065400_53262808.narrowPeak
#TE=/groups/stark/nikolaus.mandlburger/Projects/blastoid_project/res/pseudobulks/Saurabh_definition_2024-09-13/broadct_TE/TE_96h_rpm_norm.bw
#EPI=/groups/stark/nikolaus.mandlburger/Projects/blastoid_project/res/pseudobulks/Saurabh_definition_2024-09-13/broadct_epi/epiblast_96h_rpm_norm.bw
#HYPO=/groups/stark/nikolaus.mandlburger/Projects/blastoid_project/res/pseudobulks/Saurabh_definition_2024-09-13/broadct_hypo/hypoblast_96h_rpm_norm.bw
#THRESHOLD=3
#$SCRIPT -p $PEAKS -t $TE -e $EPI -y $HYPO -f $THRESHOLD