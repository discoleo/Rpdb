
### Molecular Structure


is.structure = function(x) {
	if(inherits(x, "structure")) return(TRUE);
	return(FALSE);
}

as.pdb.structure = function(helix = NULL, sheet = NULL, ...) {
	x = list(Helix = helix, Sheet = sheet);
	class(x) = c("structure", "list");
	return(x);
}

# Note:
# - Positions are coded using 4 digits;
# - Duplicate positions are disambiguated using the Insrtion Code;
# - Insertion codes use letters: potential for 26*2 more positions;
# - Sense of Strand:
#   0 = First strand;
#   1 = Parallel to previous, -1 = Anti-parallel to previous;
extract.pdb.structure = function(x, recname) {
	### Helix:
	idHelix = which(recname == "HELIX ");
	isHelix = length(idHelix) > 0;
	helix   = NULL;
	if(isHelix) {
		sH = x[idHelix];
		hID   = trim(substr(sH,  7, 10));
		hName = trim(substr(sH, 11, 14));
		hAAS  = trim(substr(sH, 15, 18));
		hChS  = trim(substr(sH, 19, 20));
		hAA1  = trim(substr(sH, 21, 25));
		hAAE  = trim(substr(sH, 27, 30)); # 1 Space before?
		hChE  = trim(substr(sH, 31, 32));
		hAA2  = trim(substr(sH, 33, 37)); # 1 Space after?
		helix = data.frame(ID = as.integer(hID), Name  = hName,
			chS = hChS, chE = hChE,
			AAS = hAAS, idAAS = as.integer(hAA1),
			AAE = hAAE, idAAE = as.integer(hAA2));
	}
	### Sheet:
	# https://www.wwpdb.org/documentation/file-format-content/format33/sect5.html#SHEET
	idSheet = which(recname == "SHEET ");
	isSheet = length(idSheet) > 0;
	sheet   = NULL;
	if(isSheet) {
		sSh = x[idSheet];
		sID   = trim(substr(sSh,  8, 10)); # Strand ID
		sName = trim(substr(sSh, 12, 14)); # Sheet ID
		sNStr = trim(substr(sSh, 15, 16)); # No. Strands
		# Start:
		sAAS  = trim(substr(sSh, 18, 20)); # AA Code
		sChS  = trim(substr(sSh, 22, 22)); # Chain ID
		sAASn = trim(substr(sSh, 23, 26)); # Seq No.
		sAASI = trim(substr(sSh, 27, 27)); # Insertion Code
		# End:
		sAAE  = trim(substr(sSh, 29, 31)); # AA Code
		sChE  = trim(substr(sSh, 32, 33)); # Chain ID
		sAAEn = trim(substr(sSh, 34, 37)); # Seq No
		sAAEI = trim(substr(sSh, 38, 38)); # Insertion Code
		sType = trim(substr(sSh, 39, 40)); # Sense of Strand
		### Registration
		# Current Strand:
		atomRc   = trim(substr(sSh, 42, 45)); # Atom
		insertRc = trim(substr(sSh, 55, 55)); # Insertion
		# Previous Strand:
		atomRp   = trim(substr(sSh, 42, 45)); # Atom
		insertRp = trim(substr(sSh, 70, 70)); # Insertion
		# TODO: remaining fields;
		sheet = data.frame(ID = as.integer(sID), Name  = sName,
			Type = as.integer(sType), # Sense of Strand
			chS = sChS, chE = sChE,
			AAS = sAAS, idAAS = as.integer(sAASn), InsertS = sAASI,
			AAE = sAAE, idAAE = as.integer(sAAEn), InsertE = sAAEI,
			# Registration
			atomRc   = atomRc,
			insertRc = insertRc,
			atomRp   = atomRp,
			insertRp = insertRp
			);
	}
	#
	structMol = as.pdb.structure(helix, sheet);
	return(structMol);
}
