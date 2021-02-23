# epitopes package
<!-- badges: start -->
  [![R-CMD-check](https://github.com/fcampelo/epitopes/workflows/R-CMD-check/badge.svg)](https://github.com/fcampelo/epitopes/actions)
  <!-- badges: end -->

Processing and Feature Extraction of epitope data from the Immune Epitope Database ([IEDB](http://iedb.org))

## Authors: 
- Felipe Campelo - [f.campelo@aston.ac.uk](mailto:f.campelo@aston.ac.uk), [fcampelo@gmail.com](mailto:fcampelo@gmail.com)
- Jodie Ashford - [ashfojsm@aston.ac.uk](mailto:ashfojsm@aston.ac.uk)

## Contributors:
- Francisco Lobo - [franciscolobo@ufmg.br](mailto:franciscolobo@ufmg.br)
- Joao Reis-Cunha - [jaumlrc@gmail.com](jaumlrc@gmail.com)

## Current version:
- Consolidation of linear B-cell epitopes from the full IEDB database
- Data splitting for machine learning based on protein dissimilarity
- Extraction of observations for feature calculation using a sliding window representation.
- Calculation of features.

## Upcoming:
- Modelling and plotting

## Usage  
1. Download the Complete Database Export (XML format) from [https://www.iedb.org/database_export_v3.php](https://www.iedb.org/database_export_v3.php)  
2. Run routine `get_LBCE()` to extract the linear B-cell epitopes.
3. Run `get_proteins()` to retrieve protein information from [GenBank](https://www.ncbi.nlm.nih.gov/genbank/)
4. Run `assemble_windowed_dataframe()` to assemble a tidy dataframe based on a sliding window representation.
5. Run `calc_features()` to calculate a variety of statistical and 
physiochemical features for each window.


***

Please report any bugs or suggestions directly in the repository [Issues](https://github.com/fcampelo/epitopes/issues) page.
