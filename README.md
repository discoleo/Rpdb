# Rpdb: Read, Write, Visualize and Manipulate PDB Files

- Reanimated the Rpdb package from the CRAN archive;
- Original author: Julien Idé;

## Chamges

- [fixed] deprecated rgl functions;
- [fixed] some Roxygen bugs in centre.R;

## Usage

	library(Rpdb)
	# or: devtools::load_all(_path_)
	
	x = read.pdb(file.choose())
	visualize(x, type="p")

Note:
- the default visualize is broken;
- the protein is NOT centered in the bounding box;
