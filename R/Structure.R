
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
	idSheet = which(recname == "SHEET ");
	isSheet = length(idSheet) > 0;
	sheet   = NULL;
	if(isSheet) {
		sSh = x[idSheet];
		sID   = trim(substr(sSh,  7, 10));
		sName = trim(substr(sSh, 11, 16));
		sAAS  = trim(substr(sSh, 17, 20));
		sChS  = trim(substr(sSh, 21, 22));
		sAA1  = trim(substr(sSh, 23, 26));
		sAAE  = trim(substr(sSh, 28, 31)); # 1 Space before?
		sChE  = trim(substr(sSh, 32, 33));
		sAA2  = trim(substr(sSh, 34, 37)); # 1 Space after?
		sType = trim(substr(sSh, 38, 40));
		# TODO: remaining fields;
		sheet = data.frame(ID = as.integer(sID), Name  = sName,
			Type = as.integer(sType),
			chS = sChS, chE = sChE,
			AAS = sAAS, idAAS = as.integer(sAA1),
			AAE = sAAE, idAAE = as.integer(sAA2));
	}
	#
	structMol = as.pdb.structure(helix, sheet);
	return(structMol);
}
