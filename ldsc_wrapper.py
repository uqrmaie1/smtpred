#!/usr/bin/env python

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# This script is a simple wrapper around LDSC to simplify the application of multi-trait weighting, which requires h2 and rG estimates
#
# Robert Maier
# 31/01/2017
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------

from __future__ import print_function
import subprocess
import argparse
import logging
import pandas
import numpy
import code
import os
import re


logger = logging.getLogger()
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(message)s')
ch.setFormatter(formatter)
logger.handlers = [ch]

logger.info('-------------------------------------------------------------------\n' + \
            ' LDSC wrapper                                     February 2017    \n' + \
            ' Robert Maier, Peter Visscher, Naomi Wray, Matt Robinson           \n' + \
            ' Contact: r.maier@uq.edu.au                                        \n' + \
            '-------------------------------------------------------------------\n')


def parse_arguments():

    parser = argparse.ArgumentParser()
    
    parser.add_argument('--out', type=str, required=False, default='./', help='output path')

    parser.add_argument('--n', nargs='+', type=float, required=False, help='sample size for each trait, if not in sumstats files')
    parser.add_argument('--usealln', action='store_true', help='flag to that all SNPs should be used to infer N for each trait (not just the first 100)')
    parser.add_argument('--files', nargs='+', type=str, required=False, help='file names of summary statistics for which to calculate h2 and rg')
    parser.add_argument('--snplist', type=str, required=False, help='path of the ldsc snplist file')
    parser.add_argument('--ref_ld', type=str, required=False, help='path of the ldsc ref_ld files')
    parser.add_argument('--w_ld', type=str, required=False, help='path of the ldsc w_ld files')

    group_extract = parser.add_mutually_exclusive_group(required=True)
    group_extract.add_argument('--extract', type=str, help='Path of log files of ldsc output. With this option ldsc is not run, but the h2 and rg estimates are extracted from already existing ldsc log files.')
    group_extract.add_argument('--ldscpath', type=str, help='path of the ldsc program (munge_sumstats.py and ldsc.py)')

    try:
        args = parser.parse_args()
    except IOError as msg:
        parser.error(str(msg))

    fh = logging.FileHandler(args.out + '/ldsc_wrapper.log')
    logger.addHandler(fh)

    return args


#-----------------------------------------------------------------
# read sumstats
#-----------------------------------------------------------------


def get_n(args):
    if args.usealln:
        nrows = None
    else:
        nrows = 100

    n = {}
    for i, f in enumerate(args.files):
        nam = os.path.splitext(os.path.basename(f))[0]
        if args.n and len(args.n) == len(args.files): # get n from command line input
            n[nam] = args.n[i]
        else: # get n from sumstats files
            print(f)
            newdat = pandas.read_csv(f, sep='\s+', nrows=nrows)
            if not 'N' in list(newdat):
                raise ValueError('No sample sizes provided and file ' + f + ' is missing column "N"!')
            n[nam] = numpy.median(newdat['N'])
    return n


#-----------------------------------------------------------------
# run LDSC
#-----------------------------------------------------------------

def munge_one(py, munge_sumstats, sumstats, n, snplist, out):
    if not sumstats:
        raise ValueError('Please provide sumstats file!')
    if not snplist:
        raise ValueError('Please provide SNP list file!')
    command = py + ' ' + munge_sumstats + ' --sumstats ' + sumstats + ' --N ' + str(n) + ' --merge-alleles ' + snplist + ' --out ' + out
    errorcode = subprocess.call(command, shell=True) 
    if errorcode:
        raise Exception('Error running "munge_sumstats.py"!')

def ldsc_all(py, ldsc, basenames, ref_ld, w_ld, out):
    if not ref_ld:
        raise ValueError('Please ref_ld directory!')
    if not w_ld:
        raise ValueError('Please w_ld directory!')
    sumstats = [out + '/' + s + '.sumstats.gz' for s in basenames]
    for i, nam in enumerate(basenames[:]):
        sums = ','.join([sumstats[i]] + sumstats[i:])       
        command = py + ' ' + ldsc + ' --rg ' + sums + ' --ref-ld-chr ' + ref_ld + ' --w-ld-chr ' + w_ld + ' --out ' + out + '/' + nam
        errorcode = subprocess.call(command, shell=True)
        if errorcode:
            raise Exception('Error running "ldsc.py"!')


#-----------------------------------------------------------------
# read results
#-----------------------------------------------------------------

def extract_h2_rg(ldsc_logfiles):
    """
    Args: ldsc_logfiles: directory of ldsc log files with h2 and rg estimates
    Returns: h2out: dictionary with tuples of h2, se estimates for each trait. If multiple
                  estimates are present for a trait, it returns the median of the h2 estimates
             df: Pandas DataFrame with ldsc rg results
    """
    logger.info('extract h2 and rg...')
    logger.debug('reading ldsc logfiles:' + '\n'.join(ldsc_logfiles))
    df = None
    h2s = {}
    ses = {}
    for i, fil in enumerate(ldsc_logfiles):
        nam = None
        f = open(fil, 'r')
        for l in f.readlines():
            spl = l.split()
            # h2 header
            if re.match(r'^Reading summary statistics from ', l) != None and nam == None:
                nam = re.sub('.+/', '', re.sub('.sumstats.gz ...\n', '', l))
            # h2 line
            if re.match(r'^Total Observed scale h2', l) != None:
                h2 = float(re.sub(' .+', '', re.sub('.+: ', '', l)))
                se = float(re.sub('\).*', '', re.sub('.+\(', '', l)))
                if nam not in h2s:
                    h2s[nam] = [h2]
                    ses[nam] = [se]
                else:
                    h2s[nam] += [h2]
                    ses[nam] += [se]
            # rg header
            if re.match(r'\s*p1', l) != None:
                if df is None:
                    df = pandas.DataFrame(columns=spl)
            # rg line
            if len(spl) == 12:
                rgline = True
                for stri in spl[2:]:
                    try:
                        float(stri)
                    except ValueError:
                        rgline = False
                if rgline:
                    df.loc[df.shape[0] + 1] = l.split()
        f.close()
    medianpos = [numpy.argsort(x)[len(x)//2] for x in h2s.itervalues() ]
    h2out = {}
    # order of dict is guaranteed if not modified
    for i, (k, v) in enumerate(h2s.iteritems()):
        h2out[k] = (h2s[k][medianpos[i]], ses[k][medianpos[i]])
    return h2out, df

def extract_h2(ldsc_logfiles):
    logger.info('extract h2...')
    logger.debug('reading ldsc logfiles:' + '\n'.join(ldsc_logfiles))
    h2s = {}
    for i, fil in enumerate(ldsc_logfiles):
        nam = None
        with open(fil, 'r') as f:
            l = f.readline()
            while l != '':
                if re.match(r'^Reading summary statistics from ', l) != None and nam == None:
                    nam = re.sub('.+/', '', re.sub('.sumstats.gz ...\n', '', l))
                if re.match(r'^Total Observed scale h2', l) != None:
                    h2 = float(re.sub(' .+', '', re.sub('.+: ', '', l)))
                    se = float(re.sub('\).*', '', re.sub('.+\(', '', l)))
                    h2s[nam] = (h2, se)
                    break
                l = f.readline()
    return h2s


#-----------------------------------------------------------------
# write rg, h2 and N files
#-----------------------------------------------------------------

def write_rg(rg_df, filename):
    with open(filename, 'w') as f:
        for index, row in rg_df.iterrows():
            n1 = os.path.basename(row[0]).split('.')[0]
            n2 = os.path.basename(row[1]).split('.')[0]
            if n1 != n2:
                f.write('\t'.join([n1, n2, row['rg'], row['se'], row['p']]) + '\n')
    logger.info('rG written to "' + filename + '"')

def write_h2(h2s, filename):
    with open(filename, 'w') as f:
        for nam in sorted(h2s.iterkeys()):
            h2 = h2s[nam]
            f.write(nam + '\t' + str(h2[0]) + '\t' + str(h2[1]) + '\n')
    logger.info('h2 written to "' + filename + '"')

def write_n(ns, filename):
    with open(filename, 'w') as f:
        for nam in sorted(ns.iterkeys()):
            n = ns[nam]
            f.write(nam + '\t' + str(n) + '\n')
    logger.info('n written to "' + filename + '"')





#-----------------------------------------------------------------
# main
#-----------------------------------------------------------------

def main():

    args = parse_arguments()
        
    

    #ldscpath = '/ibscratch/wrayvisscher/robert/programs/ldsc/'

    #snplist = '/ibscratch/wrayvisscher/robert/data/ldscores/w_hm3.snplist'
    #ref_ld = '/ibscratch/wrayvisscher/robert/data/ldscores/eur_w_ld_chr/'
    #w_ld = '/ibscratch/wrayvisscher/robert/data/ldscores/eur_w_ld_chr/'

    if not args.extract:
        py = 'python'

        munge_sumstats = args.ldscpath + 'munge_sumstats.py'
        ldsc = args.ldscpath + 'ldsc.py'

        basenames = map(lambda x: os.path.splitext(os.path.basename(x))[0], args.files)
        ns = get_n(args)
        
        for i, f in enumerate(args.files):
            munge_one(py, munge_sumstats, f, ns[basenames[i]], args.snplist, args.out + '/' + basenames[i])

        ldsc_all(py, ldsc, basenames, args.ref_ld, args.w_ld, args.out)

        h2s, rg_df = extract_h2_rg([ args.out + '/' + s + '.log' for s in basenames ])

        write_n(ns, args.out + '/ldsc_ns.txt')

    else:
        logfiles = [args.extract + '/' + x for x in os.listdir(args.extract) if x.endswith('.log')]
        h2s, rg_df = extract_h2_rg(logfiles)
        #h2s = extract_h2(logfiles)

    
    write_rg(rg_df, args.out + '/ldsc_rgs.txt')
    write_h2(h2s, args.out + '/ldsc_h2s.txt')
    

if __name__ == "__main__":
    main()













