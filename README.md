# supercentenarian
Code accompanying the paper "*Human mortality at extreme age*"

The supplementary material is split in three folders:

- The folder `scripts` contains the main script, `Semi_supercentenarian.R`, in which a description of the other source files can be found. The file `Semi_supercentenarian_fn.R` contains auxiliary functions (likelihoods, optimisation routine, etc.) called by the user. The **R** code assumes the presence of three databases, `italcent.rda`, `IDL.rda` and `francent.rda` in the repository.

- The `data` folder contains descriptions of the database and scripts for preprocessing the raw data and recreating the `.rda` data files given the original database. 
    - The Istat data can be bought for a small fee from the
National Institute of Statistics by registering at the Contact Center (https://contact.istat.it) and mentioning the semi-supercentenarian Survey and Marco Marsili as contact person. A description of the spreadsheet, entitled `ISTAT_Italian_SSC_1896-1910.xlsx`, is provided in `ISTAT_description.md`.
    - The other two datasets can be obtained freely by registering on http://www.supercentenarians.org/; the supplementary material provides detailed instruction for obtaining the `idl_complete.csv` file. This spreadsheet can be converted into `IDL.rda` and `francent.rda` using the script `data/IDL2021_preprocess.R`.

- The folder `auxiliary` contains metadata to enhance reproducibility. The data are updated on an ongoing basis and may fail to match the version used in the paper. To palliate this problem, and assuming no records are deleted in the future, we provide two files with the dummy identifiers of IDL individuals whose records were analysed in the paper: `IDL2021_id.txt` and `France_id.txt` give the identifiers of respectively supercentenarian (excluding French and Japanese records) and French records (both semi and supercentenarians) used. The `IDL2021_truncation_bounds.txt` provides the calendar dates for the data collection (these were extracted from the IDL country-specific metadata). Please email Leo Belzile if you need further help with extraction or formatting of the database from the IDL.

The code and supplementary materials are made available under the CC-BY-4.0 License. The list of **R** packages (with version numbers) used at the time of writing is given in the `renv.lock` file.
