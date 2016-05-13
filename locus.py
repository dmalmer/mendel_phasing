
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
        self.nondisjunction_parent = None

        self.child_phased_GT = None
        self.dad_phased_GT = None
        self.mom_phased_GT = None
        self.sibling_phased_GT = None


    def phase(self):
        #self phased_GT variables
        if self.ploidy == 2:
            self.child_phased_GT, self.dad_phased_GT, self.mom_phased_GT = \
                    phase_trio_diploid(self.child_alleles, self.dad_alleles, self.mom_alleles) 
            if sibling_alleles != None:
                self.sibling_phased_GT, _, _ = phase_trio_diploid(self.sibling_alleles, self.dad_alleles, self.mom_alleles)

        elif self.ploidy == 3:
            self.child_phased_GT, self.dad_phased_GT, self.mom_phased_GT = \
                    phase_trio_trisomy(self.child_alleles, self.dad_alleles, self.mom_alleles) 
            #TODO: make sure sibling_alleles == None in null case
            if sibling_alleles != None:
                self.sibling_phased_GT, _, _ = phase_trio_diploid(self.sibling_alleles, self.dad_alleles, self.mom_alleles)
            pass
        else:
            raise Exception('Invalid ploidy in ' + self.__class__.__name__ + ' class: ' + str(self.ploidy))


