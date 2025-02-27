## Malidup and Malisam

The Malisam and Malidup databases are downloaded from http://prodata.swmed.edu/malisam/ and http://prodata.swmed.edu/malidup/, respectively. (see the [DeepBLAST](https://github.com/flatironinstitute/deepblast) github)

Extract the compressed files and rename them into folder `Malisam` and `Malidup`, respectively, and **place them in the project folder** but NOT the `database` folder.

The complete PDB structure files are provided as the PDB entries. The aligned regions of each proteins in a pair is stored as separated files. The ground truth alignment pattern of the aligned regions are provided in `*_*.aln` and `*_*.manual.ali`, with the latter simplified from the former. Results from **TM-align**, **DALI**, and **fast** are also provided, with `*_*.<algm>.ali` as the aligned pattern and `*_*.<algm>.pdb` as the superimposed positions.

For DaliLite to work appropriately on createing data, edit the `import.pl` script on **line 102** by changing the single quote to double quote on the `system` command:
```
# original: system('mkdir $dir')
system("mkdir $dir")
```

## SCOP140 (homology detection)
The SCOP140 is downloaded according to the instructions on the [Dali](http://ekhidna2.biocenter.helsinki.fi/dali/README.benchmark) website.

## SwissTree (phylogeny reconstruction)
Full data and the filtered low-identity data used in this study is provided in **SwissTree.tar.gz** on [Zenodo](https://zenodo.org/records/14938229) and **SwissTree_cluster.tar.gz** in this folder, respectively. You can also find the data on the [foldtree](https://github.com/DessimozLab/fold_tree) github page.

## CAFA3-MF (function inference)
Data is provided in **CAFA3_MF.tar.gz** on [Zenodo](https://zenodo.org/records/14938229). You can also find the data as described in the [TEMPROT](https://doi.org/10.1186/s12859-023-05375-0) paper or on [zenodo](https://zenodo.org/records/7409660).
You can download the protein structures used in CAFA3 using the download API provided on the [AlphaFoldDB](https://alphafold.com/#/default/get_predictions_api_prediction__qualifier__get) website. UniProt IDs for target proteins and query proteins are provided in **mf_train_uniprot_list.txt** and **mf_test_uniprot_list.txt**, respectively. Noted that some predictions are not available in AlphaFoldDB v4 and it results in 31740 structures instead of 32421 for the target dataset. Structures of 15 query proteins are also missing. We predicted the missing query proteins in the test set by AlphaFold3 (truncated within 5000 aa.), converted the CIF-formatted models to PDB format, and provided them in the compressed file.
