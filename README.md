## epitopes package

Processing and Feature Extraction of epitope data from the Immune Epitope Database ([IEDB](http://iedb.org))

## Authors: 
- Felipe Campelo - [f.campelo@aston.ac.uk](mailto:f.campelo@aston.ac.uk)
- Jodie Ashford - [ashfojsm@aston.ac.uk](ashfojsm@aston.ac.uk)
- Francisco Lobo - [franciscolobo@ufmg.br](franciscolobo@ufmg.br)

## Current version:
- Consolidation of linear B-cell epitopes from the full IEDB database
- Extraction of observations for feature calculation using a sliding window representation.
- Calculation of summary peptide-based features (using package [Peptides](https://cran.r-project.org/package=Peptides)).

## Upcoming:
- Calculation of more features: cojoint triads and dipeptide frequencies (code tested, to be released to the package in late May/2020).

## Usage  
1. Download the Complete Database Export (XML format) from [https://www.iedb.org/database_export_v3.php](https://www.iedb.org/database_export_v3.php)  
2. Run routine `get_LBCE()` to extract the linear B-cell epitopes.
3. Run `get_proteins()` to retrieve protein information from [GenBank](https://www.ncbi.nlm.nih.gov/genbank/)
4. Run `assemble_windowed_dataframe()` to assemble a tidy dataframe based on a sliding window representation.
5. Run `calc_features()` to calculate a variety of statistical and 
physiochemical features for each window.


***

Please report any bugs or suggestions directly in the repository [Issues](https://github.com/fcampelo/epitopes/issues) page.
