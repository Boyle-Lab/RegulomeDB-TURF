# RegulomeDB-TURF

* Precalculated TURF scores for GWAS variants from GWAS Catalog after LD expansion (with R<sup>2</sup> threshold of 0.6) are under `GWAS_variants/`.

* All ASB SNVs and MPRA variants used for training and evaluation are under `Training/`.

<<<<<<< HEAD
* The script used to train the final TURF models `train_TURF_ASB.py` is under `Final_models/`.  We included the following features on each column as described in Supplemental Table S2:

| Column name       | Description |
| :----------- | :----------- |
| **Generic features** ||
| CHIP      | TF binding sites from ChIP-seq|
| DNASE   | DNase I hypersensitive sites from DNase-seq|
| PWM   | TF motifs from PWM matching|
| FOOTPRINT   | DNase footprints|
| EQTL_2   |eQTLs|
| PWM_matched   | DNase footprints with matched TF ChIP-seq peaks|
| FOOTPRINT_matched | TF motifs from PWM matching with matched TF ChIP-seq peaks|
| IC_change | Information content change of two alleles in PWM matching|
| IC_matched_change | Information content change of two alleles in PWM matching with matched TF ChIP-seq peaks        |
| funsig |Functional significance score from DeepSEA      |
| ChIP_quantile1,ChIP_quantile2,ChIP_quantile3,ChIP_max,ChIP_var   |Quantiles (25%,50%,75% and 100%) and variance of ChIP-seq signals across all available ChIP-seq experiments from ENCODE       |
| **Tissue-specific features** ||
|H3K4me1_tissueSp|H3K4me1 peaks from ChIP-seq|
|H3K4me3_tissueSp|H3K4me3 peaks from ChIP-seq|
|H3K27ac_tissueSp|H3K27ac peaks from ChIP-seq|
|H3K36me3_tissueSp|H3K36me3 peaks from ChIP-seq|
|H3K27me3_tissueSp|H3K27me3 peaks from ChIP-seq|
|DNASE_tissueSp|DNase I hypersensitive sites from DNase-seq|
|FOOTPRINT_tissueSp|DNase footprints|

The features were retrieved from the RegulomeDB web server as json objects (for example: https://regulomedb.org/regulome-search/?regions=chr1:39492461-39492462&genome=GRCh37&format=json), and processed with custom script `Training/Features_generation/get_RegDB_features.py`.

* An example input feature file `example_input_features.txt` and a python script for calculating TURF scores `predict_TURF.py` are under `Final_models/`.
=======
* An example of input feature file `example_input_features.txt` and a python script for calculating TURF scores `predict_TURF.py` are under `Final_models/`.
