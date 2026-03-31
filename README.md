# Rpdb: Read, Write, Visualize and Manipulate PDB Files

- Reanimated the Rpdb package from the CRAN archive;
- Original author: Julien Idé;

## Changes

- see [News.md] for all changes;

## Usage

```
library(Rpdb)
# or: devtools::load_all(_path_)

# Download some legacy PDB files:
x = read.pdb(file.choose())
visualize(x, type="l")
```

Note:
- Visualize is incomplete;
- Started redesigning the Bounding box vs PBC box;

<!-- badges: start -->
[![R-CMD-check](https://github.com/discoleo/Rpdb/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/discoleo/Rpdb/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->
