# mt_weighting
Program to combine SNP effects or individual scores from multiple-traits according to their sample size, h2 and rg

## Simple example

The goal is to combine traitA, traitB and traitC to create a more accurate predictor for traitA. It is assumed that single-trait predictors for traitA, traitB and traitC already exist, and that N, h2 and rg are known and are 1e5, 0.5 and 0.5, respectively.

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


## Input formats

### Score files

### SNP effect files

### sample size file

### h2 file

### rg file


## Output formats

### SNP effect file

### Score file

### variances

### weights


## Options

### ```--alltraits```

### ```--blup```

## LDSC wrapper


## GCTA --cojo-sblup




















