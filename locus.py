
from mendel import phase_trio_diploid, phase_trio_trisomy

# Locus class contains alleles and vcf meta data for each individual
class Locus:
    def __init__(self, chrom, pos, vcf_data, child_alleles, child_col_str, dad_alleles, dad_col_str, mom_alleles, 
                 mom_col_str, sibling_alleles=None, sibling_col_str=None):
        self.chrom = chrom
        self.pos = pos
        self.vcf_data = vcf_data
        self.child_alleles = child_alleles
        self.child_col_str = child_col_str
        self.dad_alleles = dad_alleles
        self.dad_col_str = dad_col_str
        self.mom_alleles = mom_alleles
        self.mom_col_str = mom_col_str
        self.sibling_alleles = sibling_alleles
        self.sibling_col_str = sibling_col_str

        self.ploidy = len(child_alleles)
        self.nondisjunct_parent = None

        self.child_phased_GT = None
        self.dad_phased_GT = None
        self.mom_phased_GT = None
        self.sibling_phased_GT = None

        self.is_mendelian = False
        self.is_phased = False
        self.is_heterozygous = len(set(self.child_alleles)) != 1
        self.is_phased_by_sibling = False


    def phase(self):
        #set phased_GT, is_mendelian, and is_phased_by_sibling variables
        if self.ploidy == 2:
            self.child_phased_GT, self.dad_phased_GT, self.mom_phased_GT, self.is_mendelian = \
                    phase_trio_diploid(self.child_alleles, self.dad_alleles, self.mom_alleles) 
            if self.sibling_alleles is not None:
                self.sibling_phased_GT, _, _, _ = phase_trio_diploid(self.sibling_alleles, self.dad_alleles, self.mom_alleles)
        elif self.ploidy == 3:
            self.child_phased_GT, self.dad_phased_GT, self.mom_phased_GT, self.is_mendelian = \
                    phase_trio_trisomy(self.child_alleles, self.dad_alleles, self.mom_alleles, self.nondisjunct_parent)

            if not self.is_mendelian:
                #if DS trio is not mendelian, set the sibling as unphased
                if self.sibling_alleles is not None:
                    self.sibling_phased_GT = '/'.join(self.sibling_alleles)
            else:
                if self.sibling_alleles is not None:
                    self.sibling_phased_GT, self.dad_phased_GT, self.mom_phased_GT, _ = phase_trio_diploid(self.sibling_alleles, self.dad_alleles, self.mom_alleles)
                if self.sibling_alleles is not None and '/' in self.child_phased_GT and '|' in self.sibling_phased_GT:
                    # quartet sibling phased het position in nondisjunction parent
                    self.is_phased_by_sibling = True
                    if self.nondisjunct_parent == 'dad':
                        self.child_phased_GT = self.dad_phased_GT + self.child_phased_GT[3:]
                    else:
                        self.child_phased_GT = self.child_phased_GT[:2] + self.mom_phased_GT
        else:
            raise Exception('Invalid ploidy in ' + self.__class__.__name__ + ' class: ' + str(self.ploidy))
        
        #set is_phased variable
        if '/' not in self.child_phased_GT:
            self.is_phased = True
                


