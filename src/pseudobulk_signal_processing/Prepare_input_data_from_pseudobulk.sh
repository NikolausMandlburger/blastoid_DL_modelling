#!/bin/bash


#adapted from Monika's script Prepare_input_data_1001_200bp_STARRseq.sh
#which has probably been adapted from Bernardos script /groups/stark/almeida/Projects/Sequence_rules_epigenomic_layers/results/20220203_CNN_DHS_STARRseq/Data/Prepare_input_data.sh

#This script takes as input a genome sizes file and a bigwig file with (rpm) normalised ATAC coverage (not loged)
#and it generates genomic windows of a certain length and stride and calculates the mean0 ATAC coverage at the center regions of the windows.
#It produces a 0-based bed file with the coordinates of the windows and a column with the mean0 ATAC coverage values for each window.

usage() {
  echo -e "Usage: bash Prepare_input_data_from_pseudobulk.sh -i <input.bw> -o <output.bw> -c <genome.chrom.sizes> -h\n\n"
  echo -e "Description: his script takes as input a genome sizes file and a bigwig file with (rpm) normalised ATAC coverage (not loged) and it generates genomic windows of a certain length and stride and calculates the mean0 ATAC coverage at the center regions of the windows.\nIt produces a 0-based bed file with the coordinates of the windows and a column with the mean0 ATAC coverage values for each window.\n"
  echo -e "  -i/--input_bigwig: Path to the bigwig file with (rpm) normalised ATAC coverage. \n"
  echo -e "  -o/--output_dir: Path to output directory\n"
  echo -e "  -c/--chrom_sizes: Path to chrom.sizes file of the genome used, needed for bedtools genomecov (-g). E.g. Two-column tab-separated text file containing assembly
    sequence names and sizes as in https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/\n"
  echo -e "  -h: Display this help message\n"
  exit 1
}

# Parse command-line options
while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--input_bigwig)
      input_pseudobulk_bigwig="$2"
      shift 2
      ;;
    -o|--output_dir)
      output_dir="$2"
      shift 2
      ;;
    -c|--chrom_sizes)
      chrom_sizes_file="$2"
      shift 2
      ;;
    -h|--help)
      usage
      ;;
    *)
      echo "Invalid option: $1" >&2
      usage
      exit 1
      ;;
  esac
done

#checks and setting default values

if [[ -z "$input_pseudobulk_bigwig" ]]; then
  echo "Input pseudobulk bigwig file not provided, exiting";
  exit 1;
fi

if [[ -z "$chrom_sizes_file" ]]; then
  echo "Chrom sizes file file not provided, exiting";
  exit 1;
fi

if [[ -z "$output_dir" ]]; then
  echo "Output dir not provided, setting to current workdir";
  output_dir=.
fi

### Genomic bins
width=1001 # 1001
center=-400 # for 1001 cut -400 on each side
center2=400 # for 1001 --> 400
stride=50 # for 1001 --> 50

sizes_basename=$(basename $chrom_sizes_file)
GENOME=${sizes_basename/.chrom.sizes/}
input_basename=$(basename $input_pseudobulk_bigwig)
INPUT_NAME=${input_basename/.bw/}


echo "Input pseudobulk bigwig file:"
echo $input_pseudobulk_bigwig
echo "Output directory:"
echo $output_dir
echo "Genome:"
echo $GENOME

windows_bed_file=${output_dir}/${GENOME}.${width}bp_${stride}s.windows.bed

#generating genomic windows bed file:
echo "Making windows..."
bedtools makewindows -g $chrom_sizes_file -w $width -s ${stride} | awk -v OFS='\t' '{print $1,$2,$3,$1"_"$2"_"$3,"0","+"}' > $windows_bed_file #REMOVE HEAD FOR PRODUCTION!

# 201bp center coordinates and sort bed file
echo "Restricting windows to center part..."
bedtools slop -i $windows_bed_file -g $chrom_sizes_file -b $center | sort -V -k1,1 -k2,2 > ${output_dir}/Regions_center.bed

#Calculating average over central part of windows...
echo "Calculating average over central part of windows..."
bigWigAverageOverBed $input_pseudobulk_bigwig ${output_dir}/Regions_center.bed ${output_dir}/out.tab
# add column mean0
awk '{print $5}' ${output_dir}/out.tab | paste ${output_dir}/Regions_center.bed - > ${output_dir}/Regions_center_ATAC_mean0.bed 
#resize to 1kb 
echo "Resizing to 1kb ..."
bedtools slop -i ${output_dir}/Regions_center_ATAC_mean0.bed -g $chrom_sizes_file -b $center2 > ${output_dir}/Regions_fulllength_ATAC_mean0.bed
echo "Producing final output file ..."
# add colnames -> create final file
echo -e "chr\tstart\tend\tname\tscore\tstrand\tATAC" | cat - ${output_dir}/Regions_fulllength_ATAC_mean0.bed > Regions_coverage_${INPUT_NAME}.bed
#removing intermediary files
echo "Removing intermediary files ..."
rm ${output_dir}/Regions_center.bed ${output_dir}/out.tab ${output_dir}/Regions_center_ATAC_mean0.bed  ${output_dir}/Regions_fulllength_ATAC_mean0.bed

#HG38=/groups/stark/nikolaus.mandlburger/Projects/blastoid_project/ext_data/reference_genomes/hg38/hg38.chrom.sizes
#BW=/groups/stark/nikolaus.mandlburger/Projects/blastoid_project/res/pseudobulks/R14550_ctrl/pTE/polar_cell_96h_rpm_norm.bw
#SCRIPT=/groups/stark/nikolaus.mandlburger/Projects/blastoid_project/src/scripts/Prepare_input_data_from_pseudobulk.sh
#bash $SCRIPT -i $BW -c $HG38