
import itertools
from collections import defaultdict

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


def convert_to_ploidy_alleles(child_SNPs_by_chr, dad_SNPs_by_chr, mom_SNPs_by_chr, vcf_data_by_chr):
    diploid_alleles = []
    trisomy_alleles = []
    for chrom in chrom_sort(vcf_data_by_chr.keys()):
        for (pos, child_alleles, child_col_str), (_, dad_alleles, dad_col_str), (_, mom_alleles, mom_col_str), vcf_data in \
                zip(child_SNPs_by_chr[chrom], dad_SNPs_by_chr[chrom], mom_SNPs_by_chr[chrom], vcf_data_by_chr[chrom]):
                    if len(child_alleles) == 3:
                        trisomy_alleles.append((child_alleles, dad_alleles, mom_alleles, chrom, pos, child_col_str,
                                                dad_col_str, mom_col_str, '\t'.join(vcf_data)))
                    else:
                        diploid_alleles.append((child_alleles, dad_alleles, mom_alleles, chrom, pos, child_col_str,
                                                dad_col_str, mom_col_str, '\t'.join(vcf_data)))
    return diploid_alleles, trisomy_alleles


def pairwise(iterable):
    a, b = itertools.tee(iterable)
    next(b, None)
    return itertools.izip(a, b)


def read_ped(filename):
    rels = {}
    with open(filename, 'r') as f:
        for line in f:
            cols = line.split('\t')
            if cols[2] != '0':
                rels[cols[0]] = (cols[1], cols[2], cols[3])
    return rels


def read_vcf(filename, rels, name=None):
    fam = filename.split('family')[1][0]
    header = []
    dad_SNPs_by_chr = defaultdict(list)
    mom_SNPs_by_chr = defaultdict(list)
    child_SNPs_by_chr = defaultdict(list)
    vcf_data_by_chr = defaultdict(list)
    with open(filename, 'r') as f_in:
        # header
        line = f_in.readline()
        header.append(line)
        while line[:2] == '##':
            line = f_in.readline()
            header.append(line)

        # get indexes of child/dad/mom columns
        cols = line.strip().split('\t')
        child_i = cols.index(rels[fam][0])
        dad_i = cols.index(rels[fam][1])
        mom_i = cols.index(rels[fam][2])
        
        line = f_in.readline()
        while line != '':
            cols = line.strip().split('\t')
            nucleotides = [cols[3]] + cols[4].split(',')
            child_GT = cols[child_i].split(':')[0]
            dad_GT = cols[dad_i].split(':')[0]
            mom_GT = cols[mom_i].split(':')[0]
            child_col_str = ':' + ':'.join(cols[child_i].split(':')[1:])
            dad_col_str = ':' + ':'.join(cols[dad_i].split(':')[1:])
            mom_col_str = ':' + ':'.join(cols[mom_i].split(':')[1:])
            
            dad_alleles = dad_GT.replace('/',',').replace('|',',').split(',') if dad_GT != './.' else ['0', '0']
            mom_alleles = mom_GT.replace('/',',').replace('|',',').split(',') if mom_GT != './.' else ['0', '0']
            if 'chr21' not in cols[0]:
                child_alleles = child_GT.replace('/',',').replace('|',',').split(',') if child_GT != './.' else ['0', '0']
            else:
                child_alleles = child_GT.replace('/',',').replace('|',',').split(',') if './.' not in child_GT else ['0', '0', '0']
            
            if len(set(dad_alleles + mom_alleles + child_alleles)) > 1:
                dad_SNPs_by_chr

                dad_SNPs_by_chr[cols[0]].append((int(cols[1]), [nucleotides[int(dad_a)] for dad_a in dad_alleles], dad_col_str))
                mom_SNPs_by_chr[cols[0]].append((int(cols[1]), [nucleotides[int(mom_a)] for mom_a in mom_alleles], mom_col_str))
                child_SNPs_by_chr[cols[0]].append((int(cols[1]), [nucleotides[int(child_a)] for child_a in child_alleles], child_col_str))
            
                vcf_data_by_chr[cols[0]].append('\t'.join(cols[:-3]))
                        
            line = f_in.readline()

    if name:
        if name == rels[fam][0]:
            return child_SNPs_by_chr
        elif name == rels[fam][1]:
            return dad_SNPs_by_chr
        elif name == rels[fam][2]:
            return mom_SNPs_by_chr
        else:
            raise Exception('Indvidual\'s name not contained in pedigree file')

    return dad_SNPs_by_chr, mom_SNPs_by_chr, child_SNPs_by_chr, vcf_data_by_chr, header


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

