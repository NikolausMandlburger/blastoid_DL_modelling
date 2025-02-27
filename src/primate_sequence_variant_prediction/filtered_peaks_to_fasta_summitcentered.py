#!/groups/stark/nikolaus.mandlburger/.conda/envs/genomics_python_2/bin/python3.12

#Use e.g. env /groups/stark/nikolaus.mandlburger/Projects/blastoid_project/blastoid_DL_modelling/src/envs/genomics_python_env.yml

import pyBigWig
import pandas as pd
import numpy as np
import os
import bioframe as bf
from Bio import AlignIO
from Bio.Seq import Seq
import copy
import requests
import argparse


parser = argparse.ArgumentParser(prog='filtered_peaks_to_fasta', description="Generates sequences for different primate species corresponding to loci from a certain narrowPeaks file")

#argumentslist:
#p, narrowpeaks_file
#s, window_shift
#e, window_flank_extend


parser.add_argument("-p","--filtered_peaks_file",type=str, help="Path to narrowPeaks file", required=True) 
parser.add_argument("-s","--window_shift",type=int, help="Shift of the window center relative to the peak center", default=0)
parser.add_argument("-e","--window_flank_extend",type=int, help="Extension of the window flanks", default=500)

                                
#parse and retrieve arguments
args=parser.parse_args()
                                

    
if args.filtered_peaks_file != None:
    filtered_peaks_file = args.filtered_peaks_file
else:
    print("No narrowpeaks_file given")

if args.window_shift != None:
    window_shift = args.window_shift
    window_shift=int(window_shift)
else:
    print("No window_shift given, using default value 0")

if args.window_flank_extend != None:
    window_flank_extend = args.window_flank_extend
else:
    print("No window_flank_extend given, using default value 500")


def get_cons30_from_UCSC(chrom,start,end, output_file):
    url = "https://genome.ucsc.edu/cgi-bin/hgTables"

    # Define the form data based on your selections
    data = {
        'clade': 'mammal',                      # Clade: Mammal
        'org': 'Human',                         # Genome: Human
        'db': 'hg38',                           # Assembly: GRCh38/hg38
        'hgta_group': 'compara',                # Group: Comparative Genomics
        'hgta_track': 'multiz30way',            # Track: Cons 30 Primates
        'hgta_table': 'multiz30way',            # Table: Multiz Align
        'hgta_regionType': 'range',             # Region: Specify by position
        'position': f'{chrom}:{start}-{end}', # Position (without commas)
        'hgta_outputType': 'maf',               # Output format: MAF - Multiple Alignment Format
        'hgta_doTopSubmit': 'get output',       # Equivalent to clicking "get output"
        'hgta_outFileName': output_file,     # Output filename
        'boolshad.sendToGalaxy': '0',           # Not sending to Galaxy
        'boolshad.sendToGreat': '0',            # Not sending to GREAT
        'hgta_compressType': 'none',            # No compression
        'boolshad.hgta_doSummaryStats': '0'     # Not getting summary statistics
    }

    # Make the POST request
    response = requests.post(url, data=data)

    # Check if the request was successful
    if response.status_code == 200:
        # Save the result to a file
        with open(output_file, "w") as f:
            f.write(response.text)
        print("Data successfully written to file")
    else:
        print(f"Error: Status code {response.status_code}")


species_dict={
    'hg38':["Human","Homo_sapiens"],
    'panTro5':["Chimp","Pan_troglodytes"],
    'panPan2':["Bonobo","Pan_paniscus"],
    'gorGor5':["Gorilla","Gorilla_gorilla_gorilla"],
    'ponAbe2':["Orangutan","Pongo_pygmaeus_abelii"],
    'nomLeu3':["Gibbon","Nomascus_leucogenys"],
    'rheMac8':["Rhesus","Macaca_mulatta"],
    'macFas5':["Crab-eating macaque","Macaca_fascicularis"],
    'macNem1':["Pig-tailed macaque","Macaca_nemestrina"],
    'cerAty1':["Sooty_mangabey", "Cercocebus_atys"],
    'papAnu3':["Baboon", "Papio_anubis"],
    'chlSab2':["Green_monkey","Chlorocebus_sabaeus"],
    'manLeu1':["Drill","Mandrillus_leucophaeus"],
    'nasLar1':["Proboscis_monkey","Nasalis_larvatus"], 
    'colAng1':["Angolan_colobus","Colobus_angolensis_palliatus"],
    'rhiRox1':["Golden_snub-nosed monkey","Rhinopithecus_roxellana"],
    'rhiBie1':["Black_snub-nosed monkey","Rhinopithecus_bieti"],
    "calJac3":["Marmoset","Callithrix_jacchus"],
    "saiBol1":["Squirrel_monkey","Saimiri_boliviensis"],
    "cebCap1":["White-faced_sapajou","Cebus_capucinus_imitator"],
    "aotNan1":["Ma's_night monkey","Aotus_nancymaae"],
    "tarSyr2":["Tarsier","Tarsius_syrichta"],
    "micMur3":["Mouse_lemur","Microcebus_murinus"],
    "proCoq1":["Coquerel's_sifaka","Propithecus_coquereli"],
    "eulMac1":["Black_lemur","Eulemur_macaco"],
    "eulFla1":["Sclater's_lemur","Eulemur_flavifrons"],
    "otoGar3":["Bushbaby","Otolemur_garnettii"],
    "mm10":["Mouse","Mus musculus"],
    "canFam3":["Dog","Canis_lupus_familiaris"],
    "dasNov3":["Armadillo","Dasypus novemcinctus"]
}

def maf_to_fasta(maf_file, fa_file, tsv_file, coordinates):
    #look up which species are in the multiple sequence alignment
    species = []
    with open(maf_file, "r") as maf:
        alignments = AlignIO.parse(maf, "maf")
        for a in alignments:
            for r in a:
                if not r.id in species:
                    species.append(r.id)
    
    #generate contiguous sequences for all the species in the alignment
    species_seqs = {s:Seq('') for s in species}

    with open(maf_file, "r") as maf:
        alignments = AlignIO.parse(maf, "maf")
        for a in alignments:
            for r in a:
                species_seqs[r.id]=species_seqs[r.id]+r.seq.replace("-","")

    #strip species names of chromosome extension
    species_seqs = {k[:k.find(".")]:v for k,v in species_seqs.items()}

    #narrow down to 1001bp central part
    species_seqs_1kb={}

    for k in species_seqs.keys():
        seq=species_seqs[k].upper()
        seq_1kb = copy.copy(seq)
        excess_len = len(seq)-1001
        excess_len_half =int(excess_len/2)
        seq_1kb = seq_1kb[excess_len_half:-excess_len_half]
        if len(seq_1kb)<1001:
            continue
        while len(seq_1kb)>1001:
            seq_1kb=seq_1kb[1:]
        
        if set(seq_1kb)=={'A', 'C', 'G', 'T'}:
            species_seqs_1kb[k]=seq_1kb
        else:
            print(f"Warning: {k} sequence contains non-ACGT characters, skipping")
            continue
    
    #write to fasta
    with open(fa_file,"w") as fa:
        with open(tsv_file, "w") as txt:
            
            txt.write(f"seqname\tpeak_name\tcoordinates\tspecies\tspecies_trivial\n")
            for species,seq in species_seqs_1kb.items():
                species_name=species_dict[species][1]
                species_name_trivial=species_dict[species][0]
                fa.write(f">{peak_name}_shift_{window_shift}_fwd_{species_name}\n")
                fa.write(str(seq)+"\n")
                fa.write(f">{peak_name}_shift_{window_shift}_rv_{species_name}\n")
                fa.write(str(seq.reverse_complement())+"\n")
                
                txt.write(f"{peak_name}_shift_{window_shift}_fwd_{species_name}\t{coordinates}\t{peak_name}_shift_{window_shift}\t{species_name}\t{species_name_trivial}\n")
                txt.write(f"{peak_name}_shift_{window_shift}_rv_{species_name}\t{coordinates}\t{peak_name}_shift_{window_shift}\t{species_name}\t{species_name_trivial}\n")



filtered_peaks_base_path = os.path.dirname(filtered_peaks_file)

filtered_peaks = bf.read_table(filtered_peaks_file, schema="narrowPeak")
input_peaks_name = os.path.basename(filtered_peaks_file).split("_")
target_gene=input_peaks_name[-5]


for chrom,start,rel_summit,name in zip(filtered_peaks.chrom, filtered_peaks.start, filtered_peaks.relSummit, filtered_peaks.name):
    peak_name = "_".join(name.split("_")[-2:])
    peak_name = f"{target_gene}_{peak_name}"

    #1001 bp window around summit (0 based coordinates)
    summit = start+rel_summit
    summit = summit+window_shift
    start_core_window = summit-501
    end_core_window = summit+500
    print(f"Processing {peak_name} at {chrom}:{start_core_window+1}-{end_core_window}")

    peak_maf_file =  os.path.join(filtered_peaks_base_path, f"{peak_name}_{chrom}_{start_core_window+1}_{end_core_window}_shift_{window_shift}.maf")

    #extend flanks 
    start_window_extended = start_core_window - window_flank_extend
    end_window_extended = end_core_window + window_flank_extend


    #retrieve cons30 alignment 
    get_cons30_from_UCSC(chrom, start_window_extended, end_window_extended, peak_maf_file)

    peak_fa_file = peak_maf_file.replace(".maf",".fa")
    peak_tsv_file = peak_maf_file.replace(".maf",".tsv")
    maf_to_fasta(peak_maf_file, peak_fa_file, peak_tsv_file, f"{chrom}:{start_core_window+1}-{end_core_window}") #1-based coordinates


#usage:
#SCRIPT=/groups/stark/nikolaus.mandlburger/Projects/blastoid_project/src/primate_sequence_variants/enhancer_candidate_selection/filtered_peaks_to_fasta.py
#FILTEREDPEAKS=/groups/stark/nikolaus.mandlburger/Projects/blastoid_project/res/primate_sequence_variants/enhancer_candidate_selection/GCM1/pTE_GCM1_chr6_53065400_53262808_TEspec.narrowPeak
#$SCRIPT -p $FILTEREDPEAKS -s 0 -e 500
