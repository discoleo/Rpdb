
# Rpdb 2.4.0

> started work;

## Updates to Existing Functions

- Function read.pdb:
  - Bug fixed: process only FieldName CRYST1;
  - New arg: verbose = TRUE;
  - Warn if pdb is marked as obsolete;
  - Multi-line titles: import only the text;
  - Extract resolution from pdb;
  - Renaming arg CRYST1 to CRYSTAL;
- Function is.crystal: as replacement to is.cryst1;
- Function visualize.coords: starting some refactoring / optimization;
- Function toSymbols: started refactoring;


# Rpdb 2.3.4

- new maintainer;

## Updates of Existing Functions

- code updated to use rgl >= 1.1.3;

## Documentation

- corrected various bugs with Roxygen2;
- improved compliance with CRAN-checks;
