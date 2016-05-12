
#PBS -S /bin/sh
#PBS -V

#PBS -l nodes=1:ppn=1
#PBS -l pmem=40gb
#PBS -l walltime=8:00:00

#PBS -N mendel-phase

#PBS -m ae
#PBS -M daniel.malmer@colorado.edu

python /projects/Down/Dowellseq/genomes/mendel_phasing/phase_main.py ${INPUT_FILE}

