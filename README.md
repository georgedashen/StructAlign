# StructAlign-evaluator
Benchmarking protein structure alignment algorithms on several databases
* SCOP200: 200 non-homologous proteins from the SCOP database used in the TM-align paper (TM-score)
* Malisam: 130 difficult non-homologous proteins in different families from SCOP (Accuracy, TM-score)
* Malidup: 241 difficult proteins with internal duplicated structures (Accuracy, TM-score)
* HOMSTRAD: 590 single-domain non-homologous protein structure alignments (Accuracy, TM-score)
* CAFA3-MF: molecular function prediction using homology search (Fmax, AUPR)
* ?, phylogeny dataset, reconstruct evolutionary distances

## Environment
```
* python=3.7.16
* numpy=1.21.5
```

## Database download
Before running any codes, please download all databases for benchmarking following the instruction in the `database` folder.

## Accuracy evaluation
First retrieve the alignment pattern from the result and arrange it as the following:
```
ppakRPEQGLLRLRKGLD--lYANLRPAQIF--DVDILVVREltGNMFGDILSDEASQLTgs----igMLPSASLGe-----------graMYEPIHGS
-ftyEEVLAFEERLEREAeapSLYTVEHKVDfpVEHCYEKAL--GAEGVEEVYRRGLAQRhalpfeadGVVLKLDDltlwgelgytaraprFALAYKFP
```

Then use run the following script to calculate the accuracy given a ground truth alignment:
```
python accuracy.py <groundtruth.ali> <predict.ali>
```

The outputs contain accuracy, recall, and precision score. When there are lowercase letters in the `<predict.ali>` file, which stand for unaligned or low-confident positions, three additional scores: accuracy, recall, and precision are calculated by taking the confidence into account. Noted that letters in the `<groundtruth.ali>` file are all uppercase, and `<predict.ali>` file with only uppercase output results in identical accuracy, recall and precision scores.

## TM-score evaluation
First download **TM-align** and compile it according to the intruction on [Zhang's lab](https://zhanggroup.org/TM-align/) website.

