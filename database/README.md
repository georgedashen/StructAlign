The Malisam and Malidup databases are downloaded from http://prodata.swmed.edu/malisam/ and http://prodata.swmed.edu/malidup/, respectively. (see the [DeepBLAST](https://github.com/flatironinstitute/deepblast) github)

Extract the compressed files and rename them into folder `Malisam` and `Malidup`, respectively.

The complete PDB structure files are provided as the PDB entries. The aligned regions of each proteins in a pair is stored as separated files. The ground truth alignment pattern of the aligned regions are provided in `*_*.aln` and `*_*.manual.ali`, with the latter simplified from the former. Results from **TM-align**, **DALI**, and **fast** are also provided, with `*_*.<algm>.ali` as the aligned pattern and `*_*.<algm>.pdb` as the superimposed positions.
