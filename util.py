
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


def pairwise(iterable):
    a, b = itertools.tee(iterable)
    next(b, None)
    return itertools.izip(a, b)


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


def read_vcf(filename, rels, fam_size):
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
        sib_i = cols.index(rels[fam + '-s'][0]) if fam_size == 'quartet' else None
        
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
            sibling_GT = cols[sib_i].split(':')[0] if fam_size == 'quartet' else None
            sibling_col_str = ':' + ':'.join(cols[sib_i].split(':')[1:]) if fam_size == 'quartet' else None
            
            # convert to GT value lists
            dad_vals = dad_GT.replace('/',',').replace('|',',').split(',') if dad_GT != './.' else ['0', '0']
            mom_vals = mom_GT.replace('/',',').replace('|',',').split(',') if mom_GT != './.' else ['0', '0']
            if 'chr21' not in cols[0]:
                child_vals = child_GT.replace('/',',').replace('|',',').split(',') if child_GT != './.' else ['0', '0']
            else:
                child_vals = child_GT.replace('/',',').replace('|',',').split(',') if './.' not in child_GT else ['0', '0', '0']
            sibling_vals = None
            if fam_size == 'quartet':
                sibling_vals = sibling_GT.replace('/',',').replace('|',',').split(',') if sibling_GT != './.' else ['0', '0']
            
            # convert to nucleotides
            dad_alleles = [nucleotides[int(dad_v)] for dad_v in dad_vals]
            mom_alleles = [nucleotides[int(mom_v)] for mom_v in mom_vals]
            child_alleles = [nucleotides[int(child_v)] for child_v in child_vals]
            sibling_alleles = [nucleotides[int(sib_v)] for sib_v in sibling_vals] if fam_size == 'quartet' else None

            # add to alleles list
            alleles_by_chr[cols[0]].append(Locus(cols[0], int(cols[1]), '\t'.join(cols[:-3]), child_alleles, child_col_str,
                                                 dad_alleles, dad_col_str, mom_alleles, mom_col_str, sibling_alleles,
                                                 sibling_col_str))
                        
            line = f_in.readline()

    return alleles_by_chr, header


def output_phased_vcf(phased_tups, vcf_header, phased_regions, fdir, shortname):
    with open(fdir + shortname + '.vcf', 'w') as f:
        for line in vcf_header:
            f.write(line)

        #tup = (chr, pos, child_phase, dad_phase, mom_phase, is_mendelian, line_pre)
        for tup in phased_tups:
            f.write(tup[6] + '\t' + '\t'.join(tup[2:5]) + '\n')

    with open(fdir + shortname + '_phased-regions.bed', 'w') as f:
        for tup in phased_regions:
            f.write('\t'.join(tup) + '\n')

