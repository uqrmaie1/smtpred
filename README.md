# mt_weighting
Program to combine SNP effects or individual scores from multiple-traits according to their sample size, h2 and rg

## Simple example

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

## General process

When combining multiple traits to create a more powerful predictor, it is possible to either generate a weighted sum of SNP effects and use those for weighting, or to first create individual scores for each trait (using the ```PLINK --score``` function) and then create a weighted sum of those. If the set of SNPs for each trait is similar, the resulting predictor will also be very similar.

By default it is assumed that single trait SNP effects are OLS estimates (GWAS profile scores). If instead of OLS estimates they are BLUP or SBLUP estimates, the ```--blup``` option can be used to calculate the appropriate weights. Even though the weights will be different under this option, the resulting weighted sum will be very similar, because of changes in the expected variance of the SNP effects or individual scores.


## Input formats

### Score files

Weighting is performed on individual scores, if the option ```--scorefiles``` or ```--scorepath``` is specified. Score files have to be in the format of the output of ```PLINK --score``` (```.profile``` files) ```--scorepath``` assumes that all files in this directory are PLINK score files.

### SNP effect files

Weighting is performed on SNP effects, if the option ```--betafiles``` or ```--betapath``` is specified. SNP effect files for each trait all have to be in the same format, and have to have a header line with three required fields: SNP ID (called "snp", "snpid", "rs", "rsid"; case insensitive), effect allele (called "a1"; case insensitive) and SNP effect (called "beta" or "b"; case insensitive). SNP IDs will be matched on their ID and effect allele "a1", and optionally on "a2" if it exists. "a1" (and "a2") have to match exactly among traits, otherwise the SNP will not be used.

### Sample size file

File that contains sample size of each trait (option ```--nfile```). This file has no header and two columns: Trait and sample size. Alternatively sample size input can be provided directly using the option ```--n```.

### h2 file

File that contains SNP heritability estimates of each trait (option ```--h2file```). This file has no header and two columns: Trait and SNP heritability. Alternatively sample size input can be provided directly using the option ```--h2```.

### rg file

File that contains genetic correlation (rg) estimates of each trait (option ```--rgfile```). This file has no header and three columns: Trait 1, Trait 2 and SNP heritability. Alternatively sample size input can be provided directly using the option ```--rg```.


## Output formats

### SNP effect file

If SNP effects have been provided as input, the file ```multi_trait.beta``` contains the multi-trait SNP effects. It has columns for SNP ID, effect allele and multi-trait beta for the trait of interest. By default, the trait of interest is the first trait. 

### Score file

### variances

### weights


## Options

### ```--alltraits```

### ```--blup```

## LDSC wrapper


## GCTA --cojo-sblup




















