
def find_nondisjunct_parent(trisomy_alleles):
    dad_more_shared = 0
    mom_more_shared = 0
    # only care about trio alleles, ignore rest (this can be done much cleaner in python 3.x)
    for child_alleles, dad_alleles, mom_alleles, _, _, _, _, _, _ in trisomy_alleles:
        c_tmp_d = list(child_alleles)
        c_tmp_m = list(child_alleles)
        
        d_curr = 0
        m_curr = 0
        for d in dad_alleles:
            try:
                c_tmp_d.remove(d)
                d_curr += 1
            except ValueError:
                pass
        
        for m in mom_alleles:
            try:
                c_tmp_m.remove(m)
                m_curr += 1
            except ValueError:
                pass

        if d_curr > m_curr:
            dad_more_shared += 1
        elif m_curr > d_curr:
            mom_more_shared += 1
        
    #print 'count more shared:'
    #print 'dad_more_shared = ' + str(dad_more_shared)
    #print 'mom_more_shared = ' + str(mom_more_shared)

    nondisjunct_parent = 'dad' if dad_more_shared > mom_more_shared else 'mom'
    if float(max(dad_more_shared, mom_more_shared)) / (dad_more_shared + mom_more_shared) < .66:
        pass#print 'WARNING: max more_shared is < 300% bigger than other parent'
        #raise Exception('max more_shared is < 300% bigger than other parent')
    
    return nondisjunct_parent


def phase_trio_diploid(child_alleles, dad_alleles, mom_alleles):
    # returns (child_GT, dad_GT, mom_GT, is_mendelian True/False)
    #  each GT can be phased (a|b) or unphased (a/b)
    child_GT = ''
    dad_GT = ''
    mom_GT = ''
    is_mendelian = True

    # get state of parents (0/1/2 homozygous, 0/1/2 shared)
    hom_parents = []
    num_shared_alleles = 0
    if dad_alleles[0] == dad_alleles[1]:
        hom_parents.append('dad')
    if mom_alleles[0] == mom_alleles[1]:
        hom_parents.append('mom')
    if sorted(dad_alleles)[0] == sorted(mom_alleles)[0]:
        num_shared_alleles += 1
    if sorted(dad_alleles)[1] == sorted(mom_alleles)[1]:
        num_shared_alleles += 1

    # two hom, two shared (eg. 0/0, 0/0)
    if len(hom_parents) == 2 and num_shared_alleles == 2:
        # trivially phase mom and dad
        dad_GT = '|'.join(dad_alleles)
        mom_GT = '|'.join(mom_alleles)
        
        # mendelian child must has same alleles as mom and dad
        if child_alleles[0] == child_alleles[1] == dad_alleles[0]:
            child_GT = '|'.join(child_alleles)
        else:
            # non-mendelian
            is_mendelian = False
            child_GT = '/'.join(child_alleles)
   
    # two hom, zero shared (eg. 0/0, 1/1)
    elif len(hom_parents) == 2 and num_shared_alleles == 0:
        # trivially phase mom and dad
        dad_GT = '|'.join(dad_alleles)
        mom_GT = '|'.join(mom_alleles)
        # mendelian child must have one allele from one parent and one from the other
        if child_alleles[0] in dad_alleles and child_alleles[1] in mom_alleles:
            child_GT = child_alleles[0] + '|' + child_alleles[1]
        elif child_alleles[1] in dad_alleles and child_alleles[0] in mom_alleles:
            child_GT = child_alleles[1] + '|' + child_alleles[0]
        else:
            # non-mendelian
            is_mendelian = False
            child_GT = '/'.join(child_alleles)

    # one hom, one shared (eg. 0/0, 0/1)
    elif len(hom_parents) == 1 and num_shared_alleles == 1:
        # if child is homozygous, must have received shared allele from het parent
        if child_alleles[0] == child_alleles[1]:
            # if dad is homozygous, phase mom
            if hom_parents[0] == 'dad':
                # trivially phase dad
                dad_GT = '|'.join(dad_alleles)
                # check that dad passed one allele
                if child_alleles[0] == dad_alleles[0]:
                    # phase mom
                    try:
                        c_i = mom_alleles.index(child_alleles[0])
                        n_i = 0 if c_i == 1 else 1
                        mom_GT = mom_alleles[c_i] + '|' + mom_alleles[n_i]
                        # trivially phase child
                        child_GT = '|'.join(child_alleles)
                    except:
                        # if child allele isn't in mom, non-mendelian
                        is_mendelian = False
                        child_GT = '/'.join(child_alleles)
                        mom_GT = '/'.join(mom_alleles)
                else:
                    # if child didn't receive any alleles from dad, non-mendelian
                    is_mendelian = False
                    child_GT = '/'.join(child_alleles)
                    mom_GT = '/'.join(mom_alleles)
            # if mom is homozygous, phase dad
            else:
                # trivially phase mom
                mom_GT = '|'.join(mom_alleles)
                # check that mom passed one allele
                if child_alleles[0] == mom_alleles[0]:
                    # phase dad
                    try:
                        c_i = dad_alleles.index(child_alleles[0])
                        n_i = 0 if c_i == 1 else 1
                        dad_GT = dad_alleles[c_i] + '|' + dad_alleles[n_i]
                        # trivially phase child
                        child_GT = '|'.join(child_alleles)
                    except:
                        # if child allele isn't in mom, non-mendelian
                        is_mendelian = False
                        child_GT = '/'.join(child_alleles)
                        dad_GT = '/'.join(dad_alleles)
                else:
                    # if child didn't receive any alleles from mom, non-mendelian
                    is_mendelian = False
                    child_GT = '/'.join(child_alleles)
                    dad_GT = '/'.join(dad_alleles)

        #if child is het, must have received non-shared allele from het parent
        else:
            # if dad is homozygous, phase child with mom
            if hom_parents[0] == 'dad':
                # trivially phase dad
                dad_GT = '|'.join(dad_alleles)
                try:
                    # phase child
                    d_i = child_alleles.index(dad_alleles[0])
                    m_i = 0 if d_i == 1 else 1
                    child_GT = child_alleles[d_i] + '|' + child_alleles[m_i]
                    # phase mom
                    c_i = mom_alleles.index(child_alleles[m_i])
                    n_i = 0 if c_i == 1 else 1
                    mom_GT = mom_alleles[c_i] + '|' + mom_alleles[n_i]
                except:
                    # either:
                    #  1. neither of child's alleles came from dad, non-mendelian
                    #  2. neither of mom's alleles are in child, non-mendelian
                    is_mendelian = False
                    child_GT = '/'.join(child_alleles)
                    mom_GT = '/'.join(mom_alleles)
            # if mom is homozygous, phase child with dad
            else:
                # trivially phase mom
                mom_GT = '|'.join(mom_alleles)
                try:
                    # phase child
                    m_i = child_alleles.index(mom_alleles[0])
                    d_i = 0 if m_i == 1 else 1
                    child_GT = child_alleles[d_i] + '|' + child_alleles[m_i]
                    # phase dad
                    c_i = dad_alleles.index(child_alleles[d_i])
                    n_i = 0 if c_i == 1 else 1
                    dad_GT = dad_alleles[c_i] + '|' + dad_alleles[n_i]
                except:
                    # either:
                    #  1. neither of child's alleles came from mom, non-mendelian
                    #  2. neither of dad's alleles are in child, non-mendelian
                    is_mendelian = False
                    child_GT = '/'.join(child_alleles)
                    dad_GT = '/'.join(dad_alleles)
        
    # one hom, zero shared (eg. 0/0, 1/2)
    elif len(hom_parents) == 1 and num_shared_alleles == 0:
        # if dad is homozygous, phase child with mom
        if hom_parents[0] == 'dad':
            # trivially phase dad
            dad_GT = '|'.join(dad_alleles)
            try:
                # phase child
                d_i = child_alleles.index(dad_alleles[0])
                m_i = 0 if d_i == 1 else 1
                child_GT = child_alleles[d_i] + '|' + child_alleles[m_i]
                # phase mom
                c_i = mom_alleles.index(child_alleles[m_i])
                n_i = 0 if c_i == 1 else 1
                mom_GT = mom_alleles[c_i] + '|' + mom_alleles[n_i]
            except:
                # either:
                #  1. neither of child's alleles came from dad, non-mendelian
                #  2. neither of mom's alleles are in child, non-mendelian
                is_mendelian = False
                child_GT = '/'.join(child_alleles)
                mom_GT = '/'.join(mom_alleles)
        # if mom is homozygous, phase child with dad
        else:
            # trivially phase mom
            mom_GT = '|'.join(mom_alleles)
            try:
                # phase child
                m_i = child_alleles.index(mom_alleles[0])
                d_i = 0 if m_i == 1 else 1
                child_GT = child_alleles[d_i] + '|' + child_alleles[m_i]
                # phase dad
                c_i = dad_alleles.index(child_alleles[d_i])
                n_i = 0 if c_i == 1 else 1
                dad_GT = dad_alleles[c_i] + '|' + dad_alleles[n_i]
            except:
                # either:
                #  1. neither of child's alleles came from dad, non-mendelian
                #  2. neither of mom's alleles are in child, non-mendelian
                is_mendelian = False
                child_GT = '/'.join(child_alleles)
                dad_GT = '/'.join(dad_alleles)
    
    # zero hom, two shared (eg. 0/1, 0/1)
    elif len(hom_parents) == 0 and num_shared_alleles == 2:
        # if child is homozygous, check mendelian and phase trivially
        if child_alleles[0] == child_alleles[1]:
            try:
                # phase dad
                dc_i = dad_alleles.index(child_alleles[0])
                dn_i = 0 if dc_i == 1 else 1
                dad_GT = dad_alleles[dc_i] + '|' + dad_alleles[dn_i]
                # phase mom
                mc_i = mom_alleles.index(child_alleles[0])
                mn_i = 0 if mc_i == 1 else 1
                mom_GT = mom_alleles[mc_i] + '|' + mom_alleles[mn_i]
                # trivially phase child
                child_GT = '|'.join(child_alleles)
            except:
                # child allele not in either mom or dad, non-mendelian
                is_mendelian = False
                child_GT = '/'.join(child_alleles)
                dad_GT = '/'.join(dad_alleles)
                mom_GT = '/'.join(mom_alleles)
        # if child is heterozygous, check mendelian then give up (can't phase)
        else:
            if not ((child_alleles[0] in dad_alleles and child_alleles[1] in mom_alleles) or \
                    (child_alleles[1] in dad_alleles and child_alleles[0] in mom_alleles)):
                # non-mendelian
                is_mendelian = False
            child_GT = '/'.join(child_alleles)
            dad_GT = '/'.join(dad_alleles)
            mom_GT = '/'.join(mom_alleles)

    # zero hom, one shared (eg. 0/1, 0/2)
    elif len(hom_parents) == 0 and num_shared_alleles == 1:
        # if child is homozygous, must have inherited shared allele from both parents
        if child_alleles[0] == child_alleles[1]:
            try:
                # trivially phase child
                child_GT = '|'.join(child_alleles)
                # phase dad
                dc_i = dad_alleles.index(child_alleles[0])
                dn_i = 0 if dc_i == 1 else 1
                dad_GT = dad_alleles[dc_i] + '|' + dad_alleles[dn_i]
                # phase mom
                mc_i = mom_alleles.index(child_alleles[0])
                mn_i = 0 if mc_i == 1 else 1
                mom_GT = mom_alleles[mc_i] + '|' + mom_alleles[mn_i]
            except ValueError:
                # if child allele isn't in mom or dad,  non-mendelian
                is_mendelian = False
                child_GT = '/'.join(child_alleles)
                dad_GT = '/'.join(dad_alleles)
                mom_GT = '/'.join(mom_alleles)
        # if child is heterozygous
        else:
            # get parent indexes of shared/unshared alleles
            #  (no need for exception catch, guaranteed to have one shared allele)
            shared_allele = dad_alleles[0] if dad_alleles[0] in mom_alleles else dad_alleles[1]
            ds_i = dad_alleles.index(shared_allele)
            dn_i = 0 if ds_i == 1 else 1
            ms_i = mom_alleles.index(shared_allele)
            mn_i = 0 if ms_i == 1 else 1
            # if one of the parents passed the shared allele, phase using who passed the unshared allele
            if shared_allele in child_alleles:
                # get child indexes of shared/unshared alleles
                cs_i = child_alleles.index(shared_allele)
                cn_i = 0 if cs_i == 1 else 1
                # check if mom passed unshared allele
                if child_alleles[cn_i] == mom_alleles[mn_i]:
                    child_GT = child_alleles[cs_i] + '|' + child_alleles[cn_i]
                    dad_GT = dad_alleles[ds_i] + '|' + dad_alleles[dn_i]
                    mom_GT = mom_alleles[mn_i] + '|' + mom_alleles[ms_i]
                # check if dad passed unshared allele
                elif child_alleles[cn_i] == dad_alleles[dn_i]:
                    child_GT = child_alleles[cn_i] + '|' + child_alleles[cs_i]
                    dad_GT = dad_alleles[dn_i] + '|' + dad_alleles[ds_i]
                    mom_GT = mom_alleles[ms_i] + '|' + mom_alleles[mn_i]
                else:
                    # neither parent passed unshared allele,  non-mendelian
                    is_mendelian = False
                    child_GT = '/'.join(child_alleles)
                    dad_GT = '/'.join(dad_alleles)
                    mom_GT = '/'.join(mom_alleles)

            # if neither parent passed shared allele, check mendelian and phase
            else:
                if child_alleles[0] in (dad_alleles[dn_i], mom_alleles[mn_i]) and \
                        child_alleles[1] in (dad_alleles[dn_i], mom_alleles[mn_i]):
                    cd_i = child_alleles.index(dad_alleles[dn_i])
                    cm_i = child_alleles.index(mom_alleles[mn_i])
                    child_GT = child_alleles[cd_i] + '|' + child_alleles[cm_i]
                    dad_GT = dad_alleles[dn_i] + '|' + dad_alleles[ds_i]
                    mom_GT = mom_alleles[mn_i] + '|' + mom_alleles[ms_i]
                else:
                    # non-mendelian
                    is_mendelian = False
                    child_GT = '/'.join(child_alleles)
                    dad_GT = '/'.join(dad_alleles)
                    mom_GT = '/'.join(mom_alleles)

    # zero hom, zero shared (eg. 0/1, 2/3)
    elif len(hom_parents) == 0 and num_shared_alleles == 0:
        # check mendelian and phase trivially
        if child_alleles[0] in dad_alleles and child_alleles[1] in mom_alleles:
            # phase child
            child_GT = child_alleles[0] + '|' + child_alleles[1]
            # phase dad
            dc_i = dad_alleles.index(child_alleles[0])
            dn_i = 0 if dc_i == 1 else 1
            dad_GT = dad_alleles[dc_i] + '|' + dad_alleles[dn_i]
            # phase mom
            mc_i = mom_alleles.index(child_alleles[1])
            mn_i = 0 if mc_i == 1 else 1
            mom_GT = mom_alleles[mc_i] + '|' + mom_alleles[mn_i]
        elif child_alleles[1] in dad_alleles and child_alleles[0] in mom_alleles:
            # phase child
            child_GT = child_alleles[1] + '|' + child_alleles[0]
            # phase dad
            dc_i = dad_alleles.index(child_alleles[1])
            dn_i = 0 if dc_i == 1 else 1
            dad_GT = dad_alleles[dc_i] + '|' + dad_alleles[dn_i]
            # phase mom
            mc_i = mom_alleles.index(child_alleles[0])
            mn_i = 0 if mc_i == 1 else 1
            mom_GT = mom_alleles[mc_i] + '|' + mom_alleles[mn_i]
        else:
            # non-mendelian
            is_mendelian = False
            child_GT = '/'.join(child_alleles)
            dad_GT = '/'.join(dad_alleles)
            mom_GT = '/'.join(mom_alleles)
    
    return (child_GT, dad_GT, mom_GT, is_mendelian)


def phase_trio_trisomy(child_alleles, dad_alleles, mom_alleles, nondisjunct_parent):
    # returns (child_GT, nondisj_GT, disj_GT, is_mendelian True/False)
    #  GT of parent that passed extra chrom can't be phased
    #  GT of parent that didn't can be phased (a|b) or unphased (a/b)
    #  child GT can be completely phased (a|b|c), partially phased (a|b/c), or unphased (a/b/c)
    #   the left-most child allele will is from the dad and the right-most child allele is from the mom,
    #   the middle child allele is from whichever parent passed the extra chromosome
    child_GT = ''
    dad_GT = ''
    mom_GT = ''
    is_mendelian = True

    # use parent alleles in terms of non-disjunction alleles and disjunction alleles
    nondisjunct_alleles = dad_alleles if nondisjunct_parent == 'dad' else mom_alleles
    disjunct_alleles = mom_alleles if nondisjunct_parent == 'dad' else dad_alleles

    # check for mendelian inheritance
    c_tmp = list(child_alleles)
    for a in nondisjunct_alleles:
        try:
            c_tmp.remove(a)
        except ValueError:
            is_mendelian = False
            child_GT = '/'.join(child_alleles)
            dad_GT = '|'.join(dad_alleles) if dad_alleles[0] == dad_alleles[1] else '/'.join(dad_alleles) 
            mom_GT = '|'.join(mom_alleles) if mom_alleles[0] == mom_alleles[1] else '/'.join(mom_alleles) 
            
            return (child_GT, dad_GT, mom_GT, is_mendelian)
    else:
        if c_tmp[0] in disjunct_alleles:
            # must be mendelian, phase
            #  passed disjunction allele is whatever is left in c_tmp
            disjunct_GT = c_tmp[0] + '|'
            disjunct_GT += disjunct_alleles[0] if c_tmp[0] == disjunct_alleles[1] else disjunct_alleles[1]

            #  non-disjunction alleles can only be phased if they're homozygous
            if nondisjunct_alleles[0] == nondisjunct_alleles[1]:
                nondisjunct_GT = '|'.join(nondisjunct_alleles)
            else:
                nondisjunct_GT = '/'.join(nondisjunct_alleles)
            
        else:
            is_mendelian = False
            child_GT = '/'.join(child_alleles)
            dad_GT = '|'.join(dad_alleles) if dad_alleles[0] == dad_alleles[1] else '/'.join(dad_alleles) 
            mom_GT = '|'.join(mom_alleles) if mom_alleles[0] == mom_alleles[1] else '/'.join(mom_alleles) 
            
            return (child_GT, dad_GT, mom_GT, is_mendelian)

    # assign GTs
    if nondisjunct_parent == 'dad':
        child_GT = nondisjunct_GT + '|' + disjunct_GT[0]
        dad_GT = nondisjunct_GT
        mom_GT = disjunct_GT
    else:
        child_GT = disjunct_GT[0] + '|' + nondisjunct_GT
        dad_GT = disjunct_GT
        mom_GT = nondisjunct_GT
    
    return (child_GT, dad_GT, mom_GT, is_mendelian)
     
