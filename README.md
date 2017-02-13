# MTweighting
A program to combine SNP effects or individual scores from multiple-traits according to their sample size, h2 and rg.

Table of Contents
=================

* [Simple Example](#simple-example)
* [General process](#general-process)
* [Input formats](#input-formats)
    * [SNP effect files](#snp-effect-files)
    * [Score files](#score-files)
    * [Sample size file](#sample-size-file)
    * [h2 file](#h2-file)
    * [rg file](#rg-file)
* [Output formats](#output-formats)
    * [SNP effect file](#snp-effect-file)
    * [Score file](#score-file)
    * [Weights](#weights)
    * [Variances](#variances)
* [Additional options](#additional-options)
* [Converting OLS effects to SBLUP effects](#converting-ols-effects-to-sblup-effects)
* [Further examples](#further-examples)
    * [Using ldsc wrapper to get h2 and rg](#using-ldsc-wrapper-to-get-h2-and-rg)
    * [Using ldsc wrapper to extract h2 and rg](#using-ldsc-wrapper-to-extract-h2-and-rg)
    * [Weighting OLS SNP effects](#weighting-ols-snp-effects)
    * [Weighting SBLUP individual scores](#weighting-sblup-individual-scores)
    * [Weighting OLS individual score (profile scores) files](#weighting-ols-individual-score-profile-scores-files)

Simple example
==============

Let's say we want to combine traitA, traitB and traitC to create a more accurate predictor for traitA. It is assumed that single-trait predictors for traitA, traitB and traitC already exist, and that N, h2 and rg are known and are 1e5, 0.5 and 0.5, respectively.

```
python mt_weighting.py \
  --h2 0.5 0.5 0.5 \
  --rg 0.5 0.5 0.5 \
  --n 1e5 1e5 1e5 \
  --scorefiles individual_scores/OLS/traitA.profile \
               individual_scores/OLS/traitB.profile \
               individual_scores/OLS/traitC.profile \
  --out individual_scores/wMT-OLS/
```

This will create a file "multi_trait.score" with columns FID, IID and the multi-trait profile score.

General process
===============

When combining multiple traits to create a more powerful predictor, it is possible to either generate a weighted sum of SNP effects and use those for weighting, or to first create individual scores for each trait (using the ```PLINK --score``` function) and then create a weighted sum of those. If the set of SNPs for each trait is similar, the resulting predictor will also be very similar.

By default it is assumed that single trait SNP effects are OLS estimates (GWAS profile scores). If instead of OLS estimates they are BLUP or SBLUP estimates, the ```--blup``` option can be used to calculate the appropriate weights. Even though the weights will be different under this option, the resulting weighted sum will be very similar, because of changes in the expected variance of the SNP effects or individual scores.

The examples below can be recreted using the files in the data directory. However, since this data is based on traits with low rg, it will not necessarily increase prediction accuracy.

Input formats
=============

## SNP effect files

Weighting is performed on SNP effects, if the option ```--betafiles``` or ```--betapath``` is specified. SNP effect files for each trait all have to be in the same format, and have to have a header line with three required fields: SNP ID (called "snp", "snpid", "rs", "rsid"; case insensitive), effect allele (called "a1"; case insensitive) and SNP effect (called "beta" or "b"; case insensitive). SNP IDs will be matched on their ID and effect allele "a1", and optionally on "a2" if it exists. "a1" (and "a2") have to match exactly among traits, otherwise the SNP will not be used.

## Score files

Weighting is performed on individual scores, if the option ```--scorefiles``` or ```--scorepath``` is specified. Score files have to be in the format of the output of ```PLINK --score``` (```.profile``` files) ```--scorepath``` assumes that all files in this directory are PLINK score files.


## Sample size file

File that contains sample size of each trait (option ```--nfile```). This file has no header and two columns: Trait and sample size. Alternatively sample size input can be provided directly using the option ```--n```.

## h2 file

File that contains SNP heritability estimates of each trait (option ```--h2file```). This file has no header and two columns: Trait and SNP heritability. Alternatively SNP heritability input can be provided directly using the option ```--h2```.

## rg file

File that contains genetic correlation (rg) estimates of each trait (option ```--rgfile```). This file has no header and three columns: Trait 1, Trait 2 and SNP heritability. Alternatively genetic correlation input can be provided directly using the option ```--rg```.


Output formats
==============

## SNP effect file

If SNP effects have been provided as input, the file ```multi_trait.beta``` contains the multi-trait SNP effects. It has columns for SNP ID, effect allele and multi-trait beta for the trait of interest, which is assumed to be the first trait provided. If multi-trait SNP effects for all traits are of interest, the option ```--alltraits``` will result in one column for each trait in the input files. 

## Score file

If individual scores have been provided as input, the file ```multi_trait.score``` contains the multi-trait individual scores. It has columns for FID, IID and multi-trait scores for the trait of interest, which is assumed to be the first trait provided. If multi-trait individual scores for all traits are of interest, the option ```--alltraits``` will result in one column for each trait in the input files.

## Weights

```multi_trait.weights``` will contain the weights that are used to combine traits. The header line of the file contains the traits that are used to create a multi-trait predictor. Each line contains the weights for creating a multi-trait predictor for one trait, with the first column containing the trait name and the other columns the weights for eaach trait. If ```--alltraits``` is specified, the file will have one line for each trait.

## Variances

```multi_trait.variances``` will contain expected variances for each trait. This is necessary becasue the weights assume that the variances of the SNP effects are exactly identical to their expectations. Since that is not always the case, each trait is scaled to its expected variance before weighting. For OLS effects the expected variance for each trait is h2/mtot + 1/n. For BLUP effects the expected variance for each trait is R2/meff, where R2 = h2/(1+meff*(1-R2)/(n*h2)). Despite the differences in weights and expected variances between OLS and BLUP effects, the combined effect of both will mostly cancel out and the specification of the ```--blup``` option will not change the weighted output substantially.

Additional options
==================

## ```--alltraits```

This option specifies that multi-trait weighting should be performed for all traits, rather than just for the first trait.


## ```--blup```

This option specifies that the input SNP effect or individual scores are estimated using BLUP, rather than OLS (GWAS estimates). This will affect both weights and expected variances, and have a small effect on the resulting multi-trait SNP effects and individual scores.

## ```--skipidcheck```

This option skips the ID check, which ensures that IDs will be correctly matched across SNP effect or individual scores. It results in a speedup, but requires that all input files have the same IDs (and reference alleles, in case of SNP effects).

## ```--mtot```

This option specifies the total number of markers, which is needed for the calculation of expected variances and weights if the ```--blup``` option is not specified. If ```--mtot``` is not specified, it is set to 1e6 by default.

## ```--meff```

This option specifies the effective number of markers, which is needed for the calculation of expected variances and weights if the ```--blup``` option is specified. If ```--meff``` is not specified, it is set to 90000 by default.

## ```--out```

This option specfies the location of the output files. If a path is given, the output files will have the prefix "multi_trait". If a path and file name prefix is given, that file name prefix will be used instead.


LDSC wrapper
===============

Multi trait weighting requires SNP heritability estimates for each trait and rg estimates for each pair of traits. If only summary statistics are available, LD score regression can be used to estimate these parameters. ```ldsc_wrapper.py``` is a wrapper around LD score regression, and in addition it extracts the parameters of interest from the LD score regression output files and saves them in the format used by ```mt_weighting.py```.

If LD score regression has already been run, the ```--extract``` option can be used to only process LD score regression output files.

```ldsc_wrapper.py``` is a helper script to more conveniently obtain the input parameters for ```mt_weighting.py``` and has not been tested extensively.

## Further ```ldsc_wrapper.py``` options

### ```--ldscpath```

Location of ```munge_sumstats.py``` and ```ldsc.py```.

### ```--snplist```, ```--ref_ld``` and ```--w_ld```

These are input parameters for LD score regression. Please refer to the LD score regression documentation for more information.

### ```--usealln```

By default, only the median of the first 100 lines is used to determine the sample size for each trait. This option will make sure the sample size is determined based on the median of all values.

Converting OLS effects to SBLUP effects
=======================================

Multi-trait weighting can be applied to both OLS (GWAS) effects, as well as BLUP effects, which often result in higher prediction accuracy. Typically BLUP effects require individual level genotype data, but ```GCTA --cojo-sblub``` allows to transform OLS effects into BLUP-like (SBLUP) effects, requiring only summary statistics and an LD reference panel:

```bash
lambda=5000000
for sumstats in `ls snp_effects/OLS/`; do
      gcta --bfile testset/test \
           --cojo-file snp_effects/OLS/${sumstats} \
           --cojo-sblup ${lambda} \
           --cojo-wind 2000 \
           --thread-num 20 \
           --out snp_effects/SBLUP/`basename ${sumstats} .txt`
      awk '{print $1, $2, $4}' snp_effects/SBLUP/`basename ${sumstats} .txt`.sblup.cojo > snp_effects/SBLUP/`basename ${sumstats}`
done
```



The shrinkage parameter lambda should be M * (1-h2)/h2, where M is the total number of (overlapping) markers and h2 is the SNP heritability of that trait.

Summary statistics input files have to these columns: SNP, A1, A2, freq, b, se, p, N; and should include a header line. For more information see http://cnsgenomics.com/software/gcta/cojo.html.

The reference panel genotype file should be in PLINK binary format. For more information see https://www.cog-genomics.org/plink2/formats#bed.

The analysis can be sped up by parallelizing over chromosomes.


Further examples
================


## Using ldsc wrapper to get h2 and rg

```
python ldsc_wrapper.py \
    --out ldsc/ \
    --files snp_effects/OLS/traitA.txt \
            snp_effects/OLS/traitB.txt \
            snp_effects/OLS/traitC.txt \
    --ldscpath /path/to/ldsc/ \
    --snplist /path/to/w_hm3.snplist \
    --ref_ld /path/to/eur_w_ld_chr/ \
    --w_ld /path/to/eur_w_ld_chr/ \
```

Will result in these files:
```
traitA.sumstats.gz
traitB.sumstats.gz
traitC.sumstats.gz

traitA.log
traitB.log
traitC.log

ldsc_ns.txt
ldsc_rgs.txt
ldsc_h2s.txt
```

## Using ldsc wrapper to extract h2 and rg 

If LDSC output files already exist, the following can be used to extract h2 and rg

```
python ldsc_wrapper.py \
    --extract ldsc/ \
    --out ldsc/
```

This will result in these files:

```
ldsc_rgs.txt
ldsc_h2s.txt
```


## Weighting OLS SNP effects

```
python mt_weighting.py \
  --h2file ldsc/ldsc_h2s.txt \
  --rgfile ldsc/ldsc_rgs.txt \
  --nfile ldsc/ldsc_ns.txt \
  --betapath snp_effects/OLS/ \
  --out snp_effects/wMT-OLS/ \
  --alltraits
```

Will result in these files:

```
multi_trait.variances
multi_trait.weights
multi_trait.beta
multi_trait.log
```

Equivalently, parameters can be specified on the command line rather than in files

```bash
h2s=`awk '{printf $2 " "}' ldsc/ldsc_h2s.txt`
rgs=`awk '{printf $3 " "}' ldsc/ldsc_rgs.txt`
ns=`awk '{printf $2 " "}' ldsc/ldsc_ns.txt`

python mt_weighting.py \
  --h2 $h2s \
  --rg $rgs \
  --n $ns \
  --betapath snp_effects/OLS/ \
  --out snp_effects/wMT-OLS/direct_input \
  --alltraits
 ```

The order of h2, rg, n has to be in same order as files, which is alphabetical, if a path is provided
rg order for 3 traits: 1-2, 1-3, 2-3
for n traits:
```R
apply(combn(1:n, 2), 2, function(x) paste(x, collapse='-'))
```

This produces the same output as before:
```
diff snp_effects/wMT-OLS/direct_input.beta snp_effects/wMT-OLS/multi_trait.beta
```

To get only weights, don't provide any SNP effect or score files:
```
python mt_weighting.py \
  --h2 0.2 0.5 0.8 0.8 \
  --rg 0.8 0.8 0.8 0.5 0.5 0.5 \
  --n 1e4 1e5 1e5 1e6 \
  --out snp_effects/wMT-OLS/weights_only \
  --alltraits
```

This results in these files:

```
weights_only.variances
weights_only.weights
weights_only.log
```

## Weighting SBLUP individual scores

```
python mt_weighting.py \
  --h2file ldsc/ldsc_h2s.txt \
  --rgfile ldsc/ldsc_rgs.txt \
  --nfile ldsc/ldsc_ns.txt \
  --scorefiles individual_scores/SBLUP/traitA.profile \
               individual_scores/SBLUP/traitB.profile \
               individual_scores/SBLUP/traitC.profile \
  --out individual_scores/wMT-SBLUP/ \
  --alltraits \
  --blup
```

This results in these files:
```
multi_trait.variances
multi_trait.weights
multi_trait.score
multi_trait.log
```

Not specifying the BLUP option gives different weights and variances, but a very similar combined weighting
If meff is set to larger values than the default of 90000, the results will become indistinguishable. 

```
python mt_weighting.py \
  --h2file ldsc/ldsc_h2s.txt \
  --rgfile ldsc/ldsc_rgs.txt \
  --nfile ldsc/ldsc_ns.txt \
  --scorefiles individual_scores/SBLUP/traitA.profile \
               individual_scores/SBLUP/traitB.profile \
               individual_scores/SBLUP/traitC.profile \
  --out individual_scores/wMT-SBLUP/multi_trait_olsweighting \
  --alltraits
```

```R
blupweights = read.table('individual_scores/wMT-SBLUP/multi_trait.score', h=T)
olsweights = read.table('individual_scores/wMT-SBLUP/multi_trait_olsweighting.score', h=T)
diag(cor(blupweights[,-1:-2], olsweights[,-1:-2]))
#                     traitA                      traitB                    traitC
#                  0.9976447                   0.9996590                 0.9998719
```

Create individual scores using PLINK --score after multi-trait weighting:

```
plink2 --bfile testset \
       --score <( tail -n +2 snp_effects/wMT-OLS/multi_trait.beta )
       --out individual_scores/wMT-OLS/traitA
```


## Weighting OLS individual score (profile scores) files

```
python mt_weighting.py \
  --h2file ldsc/ldsc_h2s.txt \
  --rgfile ldsc/ldsc_rgs.txt \
  --nfile ldsc/ldsc_ns.txt \
  --scorefiles individual_scores/OLS/traitA.profile \
               individual_scores/OLS/traitB.profile \
               individual_scores/OLS/traitC.profile \
  --out sample_data/individual_scores/wMT-OLS/traitA
```

This demonstrates that first creating individual score files for each trait and then combining individual scores leads to similar results as first combining SNP effects for all traits and then creating individual scores from SNP effects.

```R
combine_score = read.table('individual_scores/wMT-OLS/traitA.score', h=T)
combine_beta = read.table('individual_scores/wMT-OLS/traitA.profile', h=T)
cor(combine_score[,3], combine_beta$SCORE)
# [1] 0.9573187
```











