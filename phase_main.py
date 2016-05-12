
import sys
from collections import defaultdict

from util import chrom_sort, convert_to_ploidy_alleles, read_ped, read_vcf, output_phased_vcf
from mendel import find_nondisjunct_parent, phase_trio_diploid, phase_trio_trisomy

def phase(diploid_alleles, trisomy_alleles):
    phases_by_chr = defaultdict(list)
    unique_phases = {}

    # diploid phasing
    for child_alleles, dad_alleles, mom_alleles, chrom, pos, child_col_str, dad_col_str, mom_col_str, line_pre_str in diploid_alleles: 
        # inheritance phasing
        child_phased_GT, dad_phased_GT, mom_phased_GT, is_mendelian = phase_trio_diploid(child_alleles, dad_alleles, mom_alleles)
    
        # add phase (pos, child_phase, dad_phase, mom_phase, is_mendelian, line_pre_str)
        phases_by_chr[chrom].append((pos, child_phased_GT + child_col_str, dad_phased_GT + dad_col_str, 
                                     mom_phased_GT + mom_col_str, is_mendelian, line_pre_str))

        k = '_'.join((child_phased_GT, dad_phased_GT, mom_phased_GT))
        if k not in unique_phases:
            unique_phases[k] = (child_phased_GT, dad_phased_GT, mom_phased_GT, child_alleles, dad_alleles, mom_alleles, is_mendelian, '')
    
    # trisomy phasing
    nondisjunct_parent = find_nondisjunct_parent(trisomy_alleles)
    for child_alleles, dad_alleles, mom_alleles, chrom, pos, child_col_str, dad_col_str, mom_col_str, line_pre_str in trisomy_alleles:
        # inheritance phasing
        child_phased_GT, dad_phased_GT, mom_phased_GT, is_mendelian = phase_trio_trisomy(child_alleles, dad_alleles, mom_alleles, nondisjunct_parent)
            
        # add phase (pos, child_phase, dad_phase, mom_phase, is_mendelian, line_pre_str)
        phases_by_chr[chrom].append((pos, child_phased_GT + child_col_str, dad_phased_GT + dad_col_str, 
                                     mom_phased_GT + mom_col_str, is_mendelian, line_pre_str))

        k = '_'.join(('_'.join(child_alleles), '_'.join(dad_alleles), '_'.join(mom_alleles), nondisjunct_parent))
        if k not in unique_phases:
            unique_phases[k] = (child_phased_GT, dad_phased_GT, mom_phased_GT, child_alleles, dad_alleles, mom_alleles, is_mendelian, nondisjunct_parent)
    
    return phases_by_chr, unique_phases


if __name__ == '__main__':
    ped_file = '/scratch/Users/dama9282/phasing/DS_trios.ped'
    output_dir = '/scratch/Users/dama9282/phasing/output/'

    # first input parameter is the vcf filename
    if len(sys.argv) < 2:
        raise Exception('Must pass in input filename as input, eg. "python %s /path/to/file.vcf"' % sys.argv[0]) 
    filename = sys.argv[1]
    print 'file: ' + filename

    # second (optional) input parameter is if we're phasing a quartet rather than a trio
    fam_size = 'trio'
    if len(sys.argv) > 2:
        fam_size = 'quartet' if 'quart' in sys.argv[2].lower() or sys.argv[2] == '4' else 'trio'

    # read ped file
    rels = read_ped(ped_file)
    print rels
    exit(0)

    # phase
    dad_SNPs_by_chr, mom_SNPs_by_chr, child_SNPs_by_chr, vcf_data_by_chr, header = read_vcf(filename, rels)
    diploid_alleles, trisomy_alleles = convert_to_ploidy_alleles(child_SNPs_by_chr, dad_SNPs_by_chr, mom_SNPs_by_chr, vcf_data_by_chr)
    phases_by_chr, unique_phases = phase(diploid_alleles, trisomy_alleles)
            
    # print results
    #for k in unique_phases:
    #    print '\nphase:'
    #    print 'dad_alleles = %i' % unique_phases[k][4]
    #    print 'mom_alleles = %i' % unique_phases[k][5]
    #    print 'child_alleles = %i' % unique_phases[k][3]
    #    print 'dad_GT = %i' % unique_phases[k][1]
    #    print 'mom_GT = %i' % unique_phases[k][2]
    #    print 'child_GT = %i' % unique_phases[k][0]
    #    print 'is_mendelian = %i' % unique_phases[k][6]
    #    print 'nondisjunct_parent = %i' % unique_phases[k][7]

    print '\n\nPHASE COUNTS:'
    dip_het_mend_phased = 0
    dip_het_mend_unphased = 0
    dip_hom_mend_phased = 0
    dip_hom_mend_unphased = 0
    dip_nonmend = 0
    tri_het_mend_phased = 0
    tri_het_mend_unphased = 0
    tri_het_mend_disj_phased = 0
    tri_hom_mend_phased = 0
    tri_hom_mend_unphased = 0
    tri_nonmend = 0
    phased_regions = []
    curr_region = ['', '', '']
    for curr_chr in chrom_sort(phases_by_chr.keys()):
        for pos, child_phased_col, dad_phased_col, mom_phased_col, is_mendelian, line_pre in phases_by_chr[curr_chr]:
            if is_mendelian:
                child_GT = child_phased_col.split(':')[0]
                dad_GT = dad_phased_col.split(':')[0]
                mom_GT = mom_phased_col.split(':')[0]
                alleles = set()
                alleles.update(child_GT.replace('|',',').replace('/',',').split(','))
                alleles.update(dad_GT.replace('|',',').replace('/',',').split(','))
                alleles.update(mom_GT.replace('|',',').replace('/',',').split(','))
                if '/' in child_GT + dad_GT + mom_GT:
                    if len(alleles) > 1:
                        if 'chr21' not in curr_chr:
                            dip_het_mend_unphased += 1
                        else:
                            tri_het_mend_unphased += 1
                            if '|' in child_GT:
                                tri_het_mend_disj_phased += 1
                    else:
                        if 'chr21' not in curr_chr:
                            dip_hom_mend_unphased += 1
                        else:
                            tri_hom_mend_unphased += 1
                    
                    # reset current phased region
                    if curr_region[1] != '':
                        phased_regions.append(curr_region)
                    curr_region = ['', '', '']

                else:
                    if len(alleles) > 1:
                        if 'chr21' not in curr_chr:
                            dip_het_mend_phased += 1
                        else:
                            tri_het_mend_phased += 1
                    else:
                        if 'chr21' not in curr_chr:
                            dip_hom_mend_phased += 1
                        else:
                            tri_hom_mend_phased += 1
                    
                    # add to current phased region
                    curr_region[2] = pos
                    if curr_region[1] == '':
                        curr_region[0] = curr_chr
                        curr_region[1] = pos

            else:
                if 'chr21' not in curr_chr:
                    dip_nonmend += 1
                else:
                    tri_nonmend += 1

    def print_division(print_str, numerator, denominator):
        print print_str,
        try:
            print float(numerator) / denominator
        except ZeroDivisionError:
            print float(numerator), '/ 0 (division by zero exception caught)'

    print ' DIPLOID:'
    print '  dip_het_mend_phased = %i' % dip_het_mend_phased
    print '  dip_het_mend_unphased = %i' % dip_het_mend_unphased
    print '  dip_hom_mend_phased = %i' % dip_hom_mend_phased
    print '  dip_hom_mend_unphased = %i' % dip_hom_mend_unphased
    print '  dip_nonmend = %i' % dip_nonmend
    print_division('  percent non-mend = ', dip_nonmend, dip_het_mend_phased + dip_het_mend_unphased + dip_hom_mend_phased + dip_hom_mend_unphased + dip_nonmend)
    print_division('  percent of HET mend that can be phased = ', dip_het_mend_phased, dip_het_mend_phased + dip_het_mend_unphased)
    print_division('  percent of HOM mend that can be phased = ', dip_hom_mend_phased, dip_hom_mend_phased + dip_hom_mend_unphased)

    print ' TRISOMY:'
    print '  tri_het_mend_phased_21 = %i' % tri_het_mend_phased
    print '  tri_het_mend_unphased_21 = %i' % tri_het_mend_unphased
    print '  tri_het_mend_disj_phased = %i' % tri_het_mend_disj_phased
    print '  tri_hom_mend_phased_21 = %i' % tri_hom_mend_phased
    print '  tri_hom_mend_unphased_21 = %i' % tri_hom_mend_unphased
    print '  tri_nonmend_21 = %i' % tri_nonmend
    print_division('  percent non-mend_21 = ', tri_nonmend, tri_het_mend_phased + tri_het_mend_unphased + tri_hom_mend_phased + tri_hom_mend_unphased + tri_nonmend)
    print_division('  percent of HET mend_21 that can be phased = ', tri_het_mend_phased, tri_het_mend_phased + tri_het_mend_unphased)
    print_division('  percent of HOM mend_21 that can be phased = ', tri_hom_mend_phased, tri_hom_mend_phased + tri_hom_mend_unphased)
    
    #output_phased_vcf(phases_by_file[shortname], headers_by_file[shortname], phased_regions, 
    #                  output_dir + '_'.join(shortname.split('_')[1:4]) + '/', shortname) 

