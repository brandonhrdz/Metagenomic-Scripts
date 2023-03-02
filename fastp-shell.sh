#!/bin/bash
set -e
set -u
set -o pipefail

#Provide:
#1) The file with the name of paired-end file 
#2) The path to the directory where the files are located

sample_names=($(cut -f 1 "$1"))

for sample in ${sample_names[@]}
do
  sample_name=$(basename ${sample} _1_sequence.fastq.gz)
  path=$2
  cd $path

  echo "Limpiando muestra: ${sample_name}"
  
  # Comprueba si los archivos existen antes de procesarlos
  if [ -f "${sample_name}_1_sequence.fastq.gz" ] && [ -f "${sample_name}_2_sequence.fastq.gz" ]; then 
    fastp -i "${sample_name}_1_sequence.fastq.gz" -I "${sample_name}_2_sequence.fastq.gz" \
    -o /hd1/brandon/data/geotraces_descargas/fastp_output/"${sample_name}_1_Q.fastq.gz" -O /hd1/brandon/data/geotraces_descargas/fastp_output/"${sample_name}_2_Q.fastq.gz" \
    -V -w 20 -q 20 -x --detect_adapter_for_pe --poly_x_min_len=40 -g --poly_g_min_len=10 -l 40 -n 15 -P 20 -p -y -Y 30 -f 10 -F 10 --html=/hd1/brandon/data/geotraces_descargas/fastp_output/"${sample_name}.html"
  fi
done
