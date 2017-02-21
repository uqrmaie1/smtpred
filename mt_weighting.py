#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# This script combines a number of single-trait SNP effects (or individual predictors) into a series of multi-trait SNP effects (or summary statistics)
#
# Robert Maier
# 30/01/2017
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------

import re
import sys
import code
import numpy
import pandas
import logging
import os.path
import argparse
import scipy.stats
import ldsc_wrapper


logging.getLogger().setLevel(logging.INFO)
logger = logging.getLogger()
ch = logging.StreamHandler()
#ch.setLevel(logging.INFO)
formatter = logging.Formatter('%(message)s')
ch.setFormatter(formatter)
logger.handlers = [ch]

logger.info('-------------------------------------------------------------------\n' + \
            ' Multi-trait weighting                            February 2017    \n' + \
            ' Robert Maier, Peter Visscher, Matt Robinson                       \n' + \
            ' Contact: r.maier@uq.edu.au                                        \n' + \
            '-------------------------------------------------------------------\n')


class Values:
    # global variables are bad
    ids = None
    ntraits = None
    snpindex = None
    numindex = None
    a2missing = None
    basenames = None
    plinkscore = None


#-----------------------------------------------------------------
# parse parameters
#-----------------------------------------------------------------

def parse_arguments():

    parser = argparse.ArgumentParser()

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--h2', nargs='+', type=float, help='SNP-heritability estimates for each trait')
    group.add_argument('--h2file', type=str, help='file with SNP-heritability estimates for each trait')

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--n', nargs='+', type=float, help='sample size for each trait')
    group.add_argument('--nfile', type=str, help='file with sample size for each trait')

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--rg', nargs='+', type=float, help='genetic correlation estimates for each pair of traits')
    group.add_argument('--rgfile', type=str, help='file with genetic correlation estimates for each pair of traits')

    parser.add_argument('--meff', default='90000', type=float, help='effective number of markers')
    parser.add_argument('--mtot', default='1e6', type=float, help='total number of markers')
    parser.add_argument('--blup', action='store_true', help='flag to indicate that SNP effects (or individual scores) are BLUP effects, not GWAS (OLS) effects')
    parser.add_argument('--alltraits', action='store_true', help='flag to indicate if multi-trait effects should be returned for all traits (otherwise only first trait)')
    parser.add_argument('--skipidcheck', action='store_true', help='flag to indicate that IDs among input file should not be merged. Will increase speed if all input files have same IDs in same order.')
    parser.add_argument('--out', type=str, default='./multi_trait', required=False, help='output file prefix')

    filegroup = parser.add_mutually_exclusive_group()
    filegroup.add_argument('--scorefiles', nargs='+', type=str, help='single-trait profile score files which should be weighted')
    filegroup.add_argument('--betafiles', nargs='+', type=str, help='single-trait beta files which should be weighted')
    filegroup.add_argument('--scorepath', type=str, help='path to single-trait profile score files which should be weighted; assumes files end in .profile')
    filegroup.add_argument('--betapath', type=str, help='path to single-trait beta files which should be weighted; assumes files end in .profile')
    

    try:
        args = parser.parse_args()
    except IOError, msg:
        parser.error(str(msg))

    args.out = args.out + ('multi_trait' if args.out.endswith('/') else '')
    fh = logging.FileHandler(args.out + '.log', mode='w')
    logger.addHandler(fh)

    return args


#-----------------------------------------------------------------
# check input
#-----------------------------------------------------------------

def check_h2_n_rg(h2s, ns, rgs):

    """
    Checks if correct number of input values have been provided
    Args: h2s: SNP heritability for each trait
          ns: sample size for each trait
          rgs: genetic correlation for each pair of traits
    """
    
    Values.ntraits = len(h2s)

    if len(ns) != Values.ntraits:
        raise ValueError("You have to provide a heritability estimate and sample size for each trait!")

    if len(rgs) != Values.ntraits * (Values.ntraits - 1) / 2:
        raise ValueError("You have to provide a genetic correlation estimate for each trait!")


#-----------------------------------------------------------------
# read h2, rg and N files
#-----------------------------------------------------------------

def read_h2(file, nam):
    """
    Args: file: file with h2 estimates. First column trait name, second column h2 estimate
          nam: ordered trait names
    Returns: list of h2 in order of nam
    Ignores everything after 2nd column (SE)
    """
    h2dict = {}
    with open(file, 'r') as f:
        for l in f:
            h2dict[l.split()[0]] = float(l.split()[1])

    if set(h2dict.iterkeys()) != set(nam):
        raise ValueError("h2 file has to have one value for each trait and has to have matching trait names!")
    h2 = []
    for name in nam:
        h2.append(h2dict[name])
    return h2

def read_n(file, nam):
    """
    Args: file: file with sample sizes. First column trait name, second column sample size
          nam: ordered trait names
    Returns: list of n in order of nam
    Ignores everything after 2nd column
    """
    ndict = {}
    with open(file, 'r') as f:
        for l in f:
            ndict[l.split()[0]] = float(l.split()[1])
    if set(ndict.iterkeys()) != set(nam):
        raise ValueError("N file has to have one value for each trait and has to have matching trait names!")
    n = []
    for name in nam:
        n.append(ndict[name])
    return n

def read_rg(file, nam):
    """
    Reads rg file. Ignores everything after 3rd column (SE). Fills in zero for missing pairs.
    Args: file: name of file containing trait1, trait2, rg in each trait
          nam: ordered trait names
    Returns: vector of genetic correlations in correct order
    """
    
    rgdict = {}
    with open(file, 'r') as f:
        for l in f:
            rgdict[(l.split()[0], l.split()[1])] = float(l.split()[2])

    for k in rgdict.keys():
        rgdict[(k[1], k[0])] = rgdict[k]

    rg = []
    for i, n1 in enumerate(nam):
        for j, n2 in enumerate(nam):
            if j > i:
                if (n1, n2) not in rgdict.iterkeys():
                    rg.append(0)
                else:
                    rg.append(rgdict[(n1, n2)])

    return rg

def get_names_from_nfile(file):
    """
    reads nfile, returns list of names
    """
    nam = []
    with open(file, 'r') as f:
        for l in f:
            nam.append(l.split()[0])
    return nam


#-----------------------------------------------------------------
# read single-trait sumstats files
#-----------------------------------------------------------------


def read_files(files, args):
    """
    Args: files: list of file names
          args: command line arguments
    Returns: list of pandas DataFrames with SNP or score data
    """

    with open(files[0], 'r') as f:
        first_line = f.readline()
    nam = first_line.split()

    Values.plinkscore = False
    Values.a2missing = False
    if args.scorefiles or args.scorepath:
        Values.plinkscore = True
        if not re.match(r'\s+FID\s+IID\s+PHENO\s+CNT\s+CNT2\s+SCORE\s+', first_line):
            raise ValueError('Expected PLINK score header in file ' + f)
    

    st_data = []

    if Values.plinkscore:
        logger.info('Files in PLINK score format found. Weighting individual effects...')
        Values.numindex = 5 # SCORE column
        for fil in files:
            logger.info('reading ' + fil + '...')
            newdat = pandas.read_csv(fil, sep='\s+')
            st_data.append(newdat)

        id0 = [str(a) + ' ' + str(b) for a, b in zip(st_data[0].ix[:,0], st_data[0].ix[:,1])]

        if args.skipidcheck:
            logger.info('skipping ID check...')
            Values.ids = id0

        else:
            logger.info('merging IDs...')
            
            Values.ids = set(id0)
            idlists = [id0]
            for i, newdat in enumerate(st_data[1:]):
                print '\rmerging file ' + str(i+2) + ' of ' + str(len(st_data)) + '...',
                sys.stdout.flush()

                idlists.append([str(a) + ' ' + str(b) for a, b in zip(newdat.ix[:,0], newdat.ix[:,1])])
                newids = set(idlists[i+1])
                Values.ids = set([x for x in Values.ids if x in newids])
            print

            # subset to common IDs
            for i in range(len(st_data)):

                newids = idlists[i]
                indices = [j for j, item in enumerate(newids) if item in Values.ids]
                st_data[i] = st_data[i].ix[indices,:]
            # bring in order of first file
            Values.ids = [x for x in id0 if x in Values.ids]

            logger.info(str(len(Values.ids)) + ' unique IDs found.')


    else:
        logger.info('Assuming SNP effects...')
        
        snpnam = ['snp', 'snpid', 'rsid', 'rs']
        for i in range(len(nam)):
            na = nam[i].lower()
            if na in snpnam:
                Values.snpindex = i
            if na == 'a1':
                a1index = i
            if na == 'a2':
                a2index = i
            if na in ['b', 'beta']:
                Values.numindex = i

        if Values.snpindex == None:
            raise ValueError('SNP column missing')

        if a1index == None:
            raise ValueError('Column A1 missing')

        if a2index == None:
            logger.info('Column A2 missing')
            Values.a2missing = True
            a2index = a1index

        
        for fil in files:
            logger.info('reading ' + fil + '...')
            newdat = pandas.read_csv(fil, sep='\s+')
            st_data.append(newdat)



        id0 = [str(snp) + ' ' + str(a1) + ' ' + str(a2) for snp, a1, a2 in zip(st_data[0].ix[:,Values.snpindex], st_data[0].ix[:,a1index], st_data[0].ix[:,a2index])]

        if args.skipidcheck:
            logger.info('skipping ID check...')
            Values.ids = id0

        else:
            logger.info('merging IDs...')

            Values.ids = set(id0)
            idlists = [id0]

            for i, newdat in enumerate(st_data[1:]):
                print '\rmerging file ' + str(i+2) + ' of ' + str(len(st_data)) + '...',
                sys.stdout.flush()

                idlists.append([str(snp) + ' ' + str(a1) + ' ' + str(a2) for snp, a1, a2 in zip(newdat.ix[:,Values.snpindex], newdat.ix[:,a1index], newdat.ix[:,a2index])])
                newids = set(idlists[i+1])
                Values.ids = set([x for x in Values.ids if x in newids])
            print

            logger.info(str(len(Values.ids)) + ' unique and matching SNPs found.')

            # Subset to set of SNPs which occur in all sumstats files and have matching allels (no flipping performed)
            for i in range(len(st_data)):
                newids = idlists[i]
                indices = [j for j, item in enumerate(newids) if item in Values.ids]
                st_data[i] = st_data[i].ix[indices,:]
            # bring in order of first file
            Values.ids = [x for x in id0 if x in Values.ids]
            # A2 not used after here, so it is removed.
            Values.ids = [re.sub(' [a-zA-Z]+$', '', x) for x in Values.ids]

            
    return st_data


#-----------------------------------------------------------------
# write output
#-----------------------------------------------------------------

def write_weights(weights, expected_variances, nam, file_prefix):
    """
    Args: weights: ntraits * ntraits numpy array of weights for each trait
          expected_variances: ntraits vector of variances which each trait should have for the weights to be correct
    """

    with open(file_prefix + '.variances', 'w') as f:
        for i in range(len(nam)):
            f.write(nam[i] + '\t')
            f.write(str(expected_variances[i]) + '\n')
    logger.info('Variances written to "' + file_prefix + '.variances''"')
    
    with open(file_prefix + '.weights', 'w') as f:
        f.write('\t' + '\t'.join(nam) + '\n')
        for i in range(len(nam)):
            f.write(nam[i] + '\t')
            f.write('\t'.join([ str(x) for x in weights[i,:]]) + '\n')
    logger.info('Weights written to "' + file_prefix + '.weights''"')
    

def write_output(mt_effects, file_prefix, args):

    if not args.alltraits:
        Values.ntraits = 1
        mt_effects = mt_effects[:, 0, None]
    if Values.plinkscore:
        ending = '.score'
        with open(file_prefix + ending, 'w') as f:
            f.write('FID\tIID\t')
            for i in range(Values.ntraits):
                f.write(Values.basenames[i] + '\t')
            f.write('\n')

            for i in range(len(Values.ids)):
                f.write(Values.ids[i].replace(' ', '\t') + '\t' + '\t'.join( map(str, mt_effects[i,]) ) + '\n')
        
        logger.info('Multi-trait scores written to "' + file_prefix + ending + '"')

    else:
        ending = '.beta'
        with open(file_prefix + ending, 'w') as f:
            f.write('snpid\tA1\t')
            for i in range(Values.ntraits):
                f.write(Values.basenames[i] + '\t')
            f.write('\n')

            for i in range(len(Values.ids)):
                f.write(Values.ids[i].replace(' ', '\t') + '\t' + '\t'.join( map(str, mt_effects[i,]) ) + '\n')

        logger.info('Multi-trait SNP effects written to "' + file_prefix + ending + '"')

#-----------------------------------------------------------------
# calculate weights
#-----------------------------------------------------------------

def r2_blup(n, h2, meff):
    """
    Args: n: sample size
          h2: SNP heritability
          meff: effective number of markers
    Returns: expected squqred correlation between BLUP predictor and phenotype; is at the same time expected variance of BLUP predictor; see Daetwyler 2008
    """

    k = meff/n
    return ( (k + h2) - numpy.sqrt( (k+h2)**2 - 4*k*h2**2) ) / (2*k)

def get_gcovmat(h2, rg):
    """
    Args: h2: vector with SNP heritabilities
          rg: vector with genetic correlations
    Returns: numpy trait by trait array with h2 on diagonal and genetic covariance on offdiagnoals
    """
    mat = numpy.zeros((len(h2), len(h2)))
    mat[numpy.triu_indices(len(h2), 1)] = rg
    mat = mat + mat.T
    mat = mat * numpy.sqrt(numpy.outer(h2, h2))
    numpy.fill_diagonal(mat, h2)
    return numpy.array(mat)


# When input files are score files, not beta files, mtot may be unknown.
# Here mtot=1e6 is assumed. The absolute value of the expected variances for each trait is irrelevant for the multi-trait weighting, so it doesn't matter too much what this value is, expecially if M > N.

def ols_variances(n, h2, mtot=1e6):
    """
    returns expected OLS beta variances for standardised genotypes and phenotypes (h2/mtot + 1/n)
    Args: n: number of individuals
          h2: SNP heritability
          mtot: total number of markers
    Returns: list of expected variances
    """
    variances = [x[1]/mtot + 1/x[0] for x in zip(n, h2)]
    return variances


def blup_variances(n, h2, meff):
    """
    returns expected BLUP beta variances for standardised genotypes and phenotypes (R2/meff)
    Args: n: number of individuals
          h2: SNP heritability
          meff: effective number of markers
    Returns: list of expected variances
    """
    blup_variances = numpy.array([r2_blup(n[i], h2[i], meff) for i in range(len(n))]) / meff
    return blup_variances

# the weights are very insensitive to changes in mtot
def ols_weights(n, h2, rg, mtot=1e6):
    """
    Args: n: vector with sample size for each trait
          h2: vector with SNP heritabilities
          rg: vector with rg for each pair of traits (3 traits: 1,2; 1,3; 2,3)
          mtot: total number of markers (doesn't change result much)
    Returns: ntraits * ntraits array with ols weights. weights in each row are for are for a multi-trait predictor of the trait in this row
    """
    ntraits = len(n)
    gcovmat = get_gcovmat(h2, rg)
    V = gcovmat / mtot
    numpy.fill_diagonal(V, ols_variances(n, h2, mtot))
    C = gcovmat / mtot

    weights = numpy.zeros([ntraits, ntraits])
    for i in range(ntraits):
        nonzero = V[i,] != 0
        Vi = V[numpy.array(numpy.where(nonzero)[0])[:, None], nonzero]
        Vinv = numpy.linalg.inv(Vi)
        weights[i, nonzero] = numpy.dot(Vinv, C[i, nonzero])

    return weights

def blup_weights(n, h2, rg, meff):
    """
    Args: n: vector with sample size for each trait
          h2: vector with SNP heritabilities
          rg: vector with rg for each pair of traits (3 traits: 1,2; 1,3; 2,3)
          meff: effective number of markers (doesn't change result much)
    Returns: ntraits * ntraits array with blup weights. weights in each row are for are for a multi-trait predictor of the trait in this row
    """

    ntraits = len(n)
    gcovmat = get_gcovmat(h2, rg)

    rsqs = numpy.array([r2_blup(n[i], h2[i], meff) for i in range(ntraits)])
    V = gcovmat * numpy.outer(rsqs, rsqs) / (numpy.outer(h2, h2) * meff)
    numpy.fill_diagonal(V, rsqs/meff)
    C = gcovmat * rsqs / (numpy.array(h2) * meff)

    weights = numpy.zeros([ntraits, ntraits])
    
    for i in range(ntraits):
        nonzero = V[i,] != 0
        Vi = V[numpy.array(numpy.where(nonzero)[0])[:, None], nonzero]
        Vinv = numpy.linalg.inv(Vi)
        weights[i, nonzero] = numpy.dot(Vinv, C[i, nonzero])
   
    return weights


def effective_weights(weights, expected_variances):
    """
    When operating on z-scores, effective weights are normal weights times expected standard deviations
    Args: weights: vector of weights
          expected_variances: vector of expected_variances
    Returns: vector of effective variances
    """

    effw = (weights * numpy.sqrt(expected_variances))
    return effw / sum(effw)



#-----------------------------------------------------------------
# apply weights
#-----------------------------------------------------------------

def get_st_matrix(st_data, column):
    mat = numpy.zeros([ len(st_data), numpy.shape(st_data[0])[0] ])
    for i in range(len(st_data)):
        mat[i,:] = st_data[i].ix[:, column]
    return mat


def fix_variance(mat, expected_variances):
    """
    Single trait effects may not have variances which follow expectations from theory.
    Here, they are standardised and multiplied with their expectations.
    Variances of scores (individual effects) are usually around mtot * var(beta). This
    is ignored here. What matters here is that the variances of different traits are
    correct relative to each other.
    Args: mat: matrix of SNP effects or individual scores
          expected_variances: vector with expected variances for each trait
    Returns: matrix of SNP effects or individual scores with fixed variances
    """
    mat = scipy.stats.zscore(mat, axis=1)
    for i in range(len(expected_variances)):
        mat[i, :] *= numpy.sqrt(expected_variances[i])
    return mat


def apply_weights(mat, weights):
    """
    Args: mat: numpy array with single trait SNP effects or individual scores
          weights: ntraits * ntraits numpy array with weights for each trait
    Returns: numpy array of multi trait SNP effects or individual scores
    """
    mt_effects = numpy.zeros([numpy.shape(mat)[1], Values.ntraits])
    for i in range(Values.ntraits):
        for j in range(Values.ntraits):
            mt_effects[:, i] += mat[j, ] * weights[i, j]
    return mt_effects


def get_mt_effects(st_data, weights, expected_variances):
    """
    Args: st_data: list of pandas DataFrames with single-trait SNP effects or individual scores
          weights: ntraits * ntraits numpy array with weights for each trait
          expected_variances: theoretical variance which each trait should have given n, h2 and rg
    Returns: numpy matrix of multi-trait SNP-effects or multi-trait individual scores
    """

    mat = get_st_matrix(st_data, Values.numindex)
    
    mat = fix_variance(mat, expected_variances)
    return apply_weights(mat, weights)



#-----------------------------------------------------------------
# main
#-----------------------------------------------------------------


def main():
    
    args = parse_arguments()

    # extract trait names
    if args.scorefiles != None:
        Values.basenames = map(lambda x: os.path.splitext(os.path.basename(x))[0], args.scorefiles)
    elif args.betafiles != None:
        Values.basenames = map(lambda x: os.path.splitext(os.path.basename(x))[0], args.betafiles)
    elif args.scorepath != None:
        Values.basenames = [re.sub('.profile', '', x) for x in sorted(os.listdir(args.scorepath))]
    elif args.betapath != None:
        Values.basenames = [re.sub('.txt', '', x) for x in sorted(os.listdir(args.betapath))]
    elif args.nfile != None:
        Values.basenames = get_names_from_nfile(args.nfile)
    else:
        Values.basenames = ['trait_' + str(i) for i in range(len(args.n))]

    
    # get parameters
    h2s = args.h2 or read_h2(args.h2file, Values.basenames)
    ns  = args.n  or read_n(args.nfile, Values.basenames)
    rgs = args.rg or read_rg(args.rgfile, Values.basenames)
    
    check_h2_n_rg(h2s, ns, rgs)

    # calculate weights
    if args.blup:
        logger.info('calculating BLUP variances and weights...')
        weights = blup_weights(ns, h2s, rgs, args.meff)
        expected_variances = blup_variances(ns, h2s, args.meff)
    else:
        logger.info('calculating OLS variances and weights...')
        weights = ols_weights(ns, h2s, rgs, mtot=1e6)
        expected_variances = ols_variances(ns, h2s, mtot=1e6)

    logger.debug(effective_weights(weights[0,:], expected_variances))

    write_weights(weights, expected_variances, Values.basenames, args.out)

    # if files with beta values or scores are provided, calculate multi-trait beta values or scores
    if args.scorefiles or args.betafiles or args.scorepath or args.betapath:

        path = args.scorepath or args.betapath
        files = args.scorefiles or args.betafiles or [path + '/' + x for x in sorted(os.listdir(path))]
        st_data = read_files(files, args)

        logger.info('calculate multi-trait ' + ('scores' if args.scorefiles or args.scorepath else 'SNP effects' ) + '...')
        mt_effects = get_mt_effects(st_data, weights, expected_variances)
        #code.interact(local=dict(globals(), **locals()))
        write_output(mt_effects, args.out, args)


if __name__ == "__main__":
    main()

