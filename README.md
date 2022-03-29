# Identifying Engaging Bird Species and Traits with Community Science Observations

This repository contains all code needed to replicate the analyses performed in "Identifying Engaging Bird Species and Traits with Community Science Observations" by Stoudt, Goldstein, and de Valpine.

### Components

The main code to be executed is contained in 8 code files, numbered 00-07.

The following data inputs are required but not published with this directory. They may be obtained from their original creators. These are:

- *The eBird Basic Dataset (EBD)*. Needed by script #00. Can be downloaded at [the eBird download page](https://ebird.org/data/download/)
- *The iNaturalist bird dataset*. Needed by script #01. Can be downloaded on [GBIF](https://doi.org/10.15468/ab3s5x)
- *The EltonTraits dataset*. Needed by script #03. [See the publication](http://doi.wiley.com/10.1890/13-1917.1).
- *The supplementary trait dataset produced by Schuetz and Johnston (2019)*. Needed by script #03. [See the publication](https://doi.org/10.1073/pnas.1820670116).
- *IUCN Red List species endangerment status*. Needed by script #03. Available on [GBIF](https://www.gbif.org/dataset/19491596-35ae-4a91-9a98-85cf505f1bd3); place
the entire folder in the `raw_data` directory.

They must be placed in the `raw_data` directory and modify the relevant code to match the file name provided.


### Other notes

This code takes a long time to run, especially files #00 (processing raw eBird) and #04-#05 (which do the actual data fitting). The former takes several days and the latter took the authors more than 3 weeks.

Throughout the data processing pipeline, the R package `taxalight` is used to associate species' scientific names between eBird, iNaturalist, and the trait datasets. Some manual adjustments were made to correct species that weren't automatically merged across platforms. To view or modify these, see the function `manual_get_ids` in `helper_code/contHexDetectionRate_fn.R`.