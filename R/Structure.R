
### Molecular Structure


extract.pdb.structure = function(x, recname) {
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
	structMol = list(Helix = helix);
	return(structMol);
}
