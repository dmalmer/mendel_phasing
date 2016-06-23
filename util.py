
import itertools
from collections import defaultdict

from locus import Locus


def chrom_sort(chroms):
    # get characters after 'chr'
    alphanum = [c.split('chr')[1] for c in chroms]
    
    # split up chroms with and without '_', convert numbers to ints for easier sorting
    mapped = []
    unmapped = []
    for an in alphanum:
        if '_' in an:
            tup = (an.split('_')[0], '_'.join(an.split('_')[1:]))
            try:
                unmapped.append((int(tup[0]), tup[1]))
            except ValueError:
                unmapped.append((tup[0], tup[1]))
        else:
            try:
                mapped.append(int(an))
            except ValueError:
                mapped.append(an)

    # sort, and move M's to front
    sorted_chroms = []
    mapped.sort()
    unmapped.sort()
    #  mapped
    if 'M' in mapped:
        sorted_chroms.append('chrM')
    sorted_chroms.extend(['chr' + str(an) for an in mapped if an != 'M'])
    #  unmapped
    sorted_chroms.extend(['chr' + tup[0] + '_' + tup[1] for tup in unmapped if tup[0] == 'M'])
    sorted_chroms.extend(['chr' + str(tup[0]) + '_' + tup[1] for tup in unmapped if tup[0] != 'M'])
    
    return sorted_chroms


def read_ped(filename):
    rels = {}
    with open(filename, 'r') as f:
        for line in f:
            fam, child, father, mother = line.split('\t')[:4]
            if father != '0':
                #if we've already added the proband child, add the sibling
                if fam in rels.keys():
                    rels[fam + '-s'] = (child, father, mother)
                else:
                    rels[fam] = (child, father, mother)
    return rels


def read_vcf(filename, rels, fam_type):
    fam = filename.split('family')[1].split('.')[0]
    header = []
    alleles_by_chr = defaultdict(list)
    with open(filename, 'r') as f_in:
        # header
        line = f_in.readline()
        header.append(line)
        while line[:2] == '##':
            line = f_in.readline()
            header.append(line)

        # get indexes of child/dad/mom/sibling columns
        cols = line.strip().split('\t')
        child_i = cols.index(rels[fam][0])
        dad_i = cols.index(rels[fam][1])
        mom_i = cols.index(rels[fam][2])
        sib_i = cols.index(rels[fam + '-s'][0]) if fam_type == 'quartet' else None
        
        line = f_in.readline()
        while line != '':
            # grab genotypes and vcf info
            cols = line.strip().split('\t')
            nucleotides = [cols[3]] + cols[4].split(',')
            child_GT = cols[child_i].split(':')[0]
            dad_GT = cols[dad_i].split(':')[0]
            mom_GT = cols[mom_i].split(':')[0]
            child_col_str = ':' + ':'.join(cols[child_i].split(':')[1:])
            dad_col_str = ':' + ':'.join(cols[dad_i].split(':')[1:])
            mom_col_str = ':' + ':'.join(cols[mom_i].split(':')[1:])
            sibling_GT = cols[sib_i].split(':')[0] if fam_type == 'quartet' else None
            sibling_col_str = ':' + ':'.join(cols[sib_i].split(':')[1:]) if fam_type == 'quartet' else None
            
            # convert to GT value lists
            dad_vals = dad_GT.replace('/',',').replace('|',',').split(',') if (dad_GT != './.' and dad_GT != '0') else ['0', '0']
            mom_vals = mom_GT.replace('/',',').replace('|',',').split(',') if (mom_GT != './.' and mom_GT != '0') else ['0', '0']
            if 'chr21' not in cols[0]:
                child_vals = child_GT.replace('/',',').replace('|',',').split(',') if (child_GT != './.' and child_GT != '0') else ['0', '0']
            else:
                child_vals = child_GT.replace('/',',').replace('|',',').split(',') if ('./.' not in child_GT and child_GT != '0') else ['0', '0', '0']
            sibling_vals = None
            if fam_type == 'quartet':
                sibling_vals = sibling_GT.replace('/',',').replace('|',',').split(',') if (sibling_GT != './.' and sibling_GT != '0') else ['0', '0']
            
            # convert to nucleotides
            dad_alleles = [nucleotides[int(dad_v)] for dad_v in dad_vals]
            mom_alleles = [nucleotides[int(mom_v)] for mom_v in mom_vals]
            child_alleles = [nucleotides[int(child_v)] for child_v in child_vals]
            sibling_alleles = [nucleotides[int(sib_v)] for sib_v in sibling_vals] if fam_type == 'quartet' else None

            # add to alleles list
            alleles_by_chr[cols[0]].append(Locus(cols[0], int(cols[1]), '\t'.join(cols[:9]), child_alleles, child_col_str,
                                                 dad_alleles, dad_col_str, mom_alleles, mom_col_str, sibling_alleles,
                                                 sibling_col_str))
                        
            line = f_in.readline()

    return alleles_by_chr, header


def output_phased_vcf(filename, loci_by_chr, header, rels, fam_type, phased_regions_by_chr):
    file_prefix = filename.split('.vcf')[0]
    fam = filename.split('family')[1].split('.')[0]

    # write phased vcf file
    with open(file_prefix + '.phased.vcf', 'w') as f:
        # write VCF header up to column header line
        for line in header[:-1]:
            f.write(line)
        
        # the order of the individuals needs to be set by us for the output (it can't just be the same order that we
        #  read it in) because there might be extra individuals in the original VCF that weren't phased)
        # write columns names up to the indv IDs
        f.write('\t'.join(header[-1].split('\t')[:9]))
        # write indv IDs in order: child, dad, mom, sibling (if quartet is used)
        f.write('\t%s\t%s\t%s' % (rels[fam][0], rels[fam][1], rels[fam][2]))
        if fam_type == 'quartet':
            f.write('\t%s' % rels[fam + '-s'][0])
        f.write('\n')

        # write phased positions
        for curr_chr in chrom_sort(loci_by_chr.keys()):
            for locus in loci_by_chr[curr_chr]:
                # write vcf line preceding phased GT info
                f.write(locus.vcf_data)

                # write phased GT columns in order of child, dad, mom, sibling (if quartet is used)
                f.write('\t%s:%s' % (locus.child_phased_GT, locus.child_col_str))
                f.write('\t%s:%s' % (locus.dad_phased_GT, locus.dad_col_str))
                f.write('\t%s:%s' % (locus.mom_phased_GT, locus.mom_col_str))
                if fam_type == 'quartet':
                    f.write('\t%s:%s' % (locus.sibling_phased_GT, locus.sibling_col_str))
                f.write('\n')

    # write bed file with contiguous phased regions
    with open(file_prefix + '.phased-regions.bed', 'w') as f:
        for curr_chr in chrom_sort(phased_regions_by_chr.keys()):
            for chrom, start, end in phased_regions_by_chr[curr_chr]:
                f.write('%s\t%i\t%i\n' % (chrom, start, end))

