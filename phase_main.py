
import sys
from collections import defaultdict

from util import read_ped, read_vcf, output_phased_vcf
from mendel import find_nondisjunct_parent, phase_trio_diploid, phase_trio_trisomy


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

    # read vcf file and find nondisjunction parent
    loci_by_chr, header = read_vcf(filename, rels, fam_size)
    trisomy_loci = [locus for loci_list in loci_by_chr.values() for locus in loci_list if locus.ploidy == 3]
    nondisjunct_parent = find_nondisjunct_parent(trisomy_loci) if len(trisomy_loci) > 0 else None

    # phase
    for loci_list in loci_by_chr.values():
        for locus in loci_list:
            if locus.ploidy == 3:
                locus.nondisjunct_parent = nondisjunct_parent
            locus.phase()

    # collect statistics and track phased regions
    phased_regions_by_chr = defaultdict(list)
    dip_het_mend_phased = 0
    dip_het_mend_unphased = 0
    dip_hom_mend_phased = 0
    dip_hom_mend_unphased = 0
    dip_nonmend = 0
    tri_het_mend_phased = 0
    tri_het_mend_sib_phased = 0
    tri_het_mend_unphased = 0
    tri_het_mend_disj_phased = 0
    tri_hom_mend_phased = 0
    tri_hom_mend_unphased = 0
    tri_nonmend = 0
    for chrom, loci_list in loci_by_chr.items():
        curr_region = ['', '', '']
        for locus in loci_list: 
            if locus.is_mendelian:
                alleles = set(locus.child_alleles + locus.dad_alleles + locus.mom_alleles)
                #if '/' in locus.child_phased_GT + locus.dad_phased_GT + locus.mom_phased_GT:
                if locus.is_phased:
                    # phased
                    #if len(alleles) > 1:
                    if locus.is_heterozygous:
                        if locus.ploidy == 2:
                            dip_het_mend_phased += 1
                        else:
                            tri_het_mend_phased += 1
                            if locus.is_phased_by_sibling:
                                tri_het_mend_sib_phased += 1
                    else:
                        if locus.ploidy == 2:
                            dip_hom_mend_phased += 1
                        else:
                            tri_hom_mend_phased += 1
                    
                    # add to current phased region
                    curr_region[2] = locus.pos
                    if curr_region[1] == '':
                        curr_region[0] = locus.chrom
                        curr_region[1] = locus.pos

                else:
                    # unphased
                    #if len(alleles) > 1:
                    if locus.is_heterozygous:
                        if locus.ploidy == 2:
                            dip_het_mend_unphased += 1
                        else:
                            tri_het_mend_unphased += 1
                            if '|' in locus.child_phased_GT:
                                #should always occur (we phase 100% of the disjunction parent alleles), assertion at end
                                tri_het_mend_disj_phased += 1
                    else:
                        #should never occur, assertion at end
                        if locus.ploidy == 2:
                            dip_hom_mend_unphased += 1
                        else:
                            tri_hom_mend_unphased += 1
                    
                    # reset current phased region
                    if curr_region[1] != '':
                        phased_regions_by_chr[locus.chrom].append(curr_region)
                    curr_region = ['', '', '']

            else:
                if locus.ploidy == 2:
                    dip_nonmend += 1
                else:
                    tri_nonmend += 1

                # reset current phased region
                if curr_region[1] != '':
                    phased_regions_by_chr[locus.chrom].append(curr_region)
                curr_region = ['', '', '']

    #TODO: ADD TWO ASSERTIONS

    def print_division(print_str, numerator, denominator):
        print print_str,
        try:
            print float(numerator) / denominator
        except ZeroDivisionError:
            print float(numerator), '/ 0 (division by zero exception caught)'

    print '\n\nPHASE COUNTS:'
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

