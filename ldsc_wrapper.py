#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# This script is a simple wrapper around LDSC to simplify the application of multi-trait weighting, which requires h2 and rG estimates
#
# Robert Maier
# 31/01/2017
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------

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
            ' Robert Maier, Peter Visscher, Matt Robinson                       \n' + \
            ' Contact: r.maier@uq.edu.au                                        \n' + \
            '-------------------------------------------------------------------\n')


def parse_arguments():

    parser = argparse.ArgumentParser()
    
    parser.add_argument('--out', type=str, required=False, default='./', help='output path')

    group = parser.add_mutually_exclusive_group(required=False)
    
    group_run = group.add_argument_group()
    group_run.add_argument('--n', nargs='+', type=float, required=False, help='sample size for each trait, if not in sumstats files')
    group_run.add_argument('--usealln', action='store_true', help='flag to that all SNPs should be used to infer N for each trait (not just the first 100)')
    group_run.add_argument('--files', nargs='+', type=str, required=False, help='file names of summary statistics for which to calculate h2 and rg')
    group_run.add_argument('--ldscpath', type=str, required=False, help='path of the ldsc program (munge_sumstats.py and ldsc.py)')
    group_run.add_argument('--snplist', type=str, required=False, help='path of the ldsc snplist file')
    group_run.add_argument('--ref_ld', type=str, required=False, help='path of the ldsc ref_ld files')
    group_run.add_argument('--w_ld', type=str, required=False, help='path of the ldsc w_ld files')

    group_extract = group.add_argument_group()
    group_run.add_argument('--extract', type=str, help='Path of log files of ldsc output. With this option ldsc is not run, but the h2 and rg estimates are extracted from already existing ldsc log files.')

    try:
        args = parser.parse_args()
    except IOError, msg:
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
            print f
            newdat = pandas.read_csv(f, sep='\s+', nrows=nrows)
            if not 'N' in list(newdat):
                raise ValueError('No sample sizes provided and file ' + f + ' is missing column "N"!')
            n[nam] = numpy.median(newdat['N'])
    return n


#-----------------------------------------------------------------
# run LDSC
#-----------------------------------------------------------------

def munge_one(py, munge_sumstats, sumstats, n, snplist, out):
    command = py + ' ' + munge_sumstats + ' --sumstats ' + sumstats + ' --N ' + str(n) + ' --merge-alleles ' + snplist + ' --out ' + out
    subprocess.call(command, shell=True)


def ldsc_all(py, ldsc, basenames, ref_ld, w_ld, out):

    sumstats = [out + '/' + s + '.sumstats.gz' for s in basenames]
    for i, nam in enumerate(basenames[:]):
        sums = ','.join([sumstats[i]] + sumstats[i:])       
        command = py + ' ' + ldsc + ' --rg ' + sums + ' --ref-ld-chr ' + ref_ld + ' --w-ld-chr ' + w_ld + ' --out ' + out + '/' + nam
        subprocess.call(command, shell=True)


#-----------------------------------------------------------------
# read results
#-----------------------------------------------------------------

def extract_rg(ldsc_logfiles):
    logger.info('extract rg...')
    logger.debug('reading ldsc logfiles:' + '\n'.join(ldsc_logfiles))
    df = None
    for i, fil in enumerate(ldsc_logfiles):
        with open(fil, 'r') as f:
            l = f.readline()
            while l != '':
                if re.match(r'\s*p1', l) != None:
                    if df is None:
                        df = pandas.DataFrame(columns=l.split())
                    l = f.readline()
                    while not re.match(r'^\s+$', l) and len(l.split()) == 12:
                        df.loc[df.shape[0] + 1] = l.split()
                        l = f.readline()
                    break
                l = f.readline()
    return df

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

        rg_df = extract_rg([ args.out + '/' + s + '.log' for s in basenames ])
        h2s = extract_h2([ args.out + '/' + s + '.log' for s in basenames ])

        write_n(ns, args.out + '/ldsc_ns.txt')

    else:
        logfiles = [args.extract + '/' + x for x in os.listdir(args.extract) if x.endswith('.log')]
        rg_df = extract_rg(logfiles)
        h2s = extract_h2(logfiles)

    #code.interact(local=dict(globals(), **locals()))
    write_rg(rg_df, args.out + '/ldsc_rgs.txt')
    write_h2(h2s, args.out + '/ldsc_h2s.txt')
    

if __name__ == "__main__":
    main()













