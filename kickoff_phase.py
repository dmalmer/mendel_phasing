
import os
from glob import glob

wkdir = '/scratch/Users/dama9282/freebayes/'

for filename in glob(wkdir + '*p5qu23.vcf'):
    name = filename.rsplit('/',1)[1].split('.vcf')[0]
    fam_size = '4' if 'familyE' in name else '3'
    if 'familyE' not in name:
        continue
    os.system('qsub -v INPUT_FILE="' + filename + '",FAM_SIZE="' + fam_size + '" -o ${PWD}/oe/ -e ${PWD}/oe/ -N mendel-phase_' + name + ' run_phase.sh')
    break

