


# TODO: better strategy;
isProtein = function(x) {
	isPr = "CA" %in% x$atoms$elename;
	return(isPr);
}
isProteinA = function(x) {
	isPr = "CA" %in% x$elename;
	return(isPr);
}

### Backbone:
asBackbone = function(x) {
	# TODO: check continuity
	FUN = function(tmp) {
		id = match(c("N", "CA", "C", "O"), tmp$elename);
		id = tmp$eleid[id];
		id = c(id, NA, id[3]);
	}
	if(inherits(x, "atoms")) {
		atoms = x;
	} else if(inherits(x, "pdb")) atoms = x$atoms;
	# Process each chain:
	ch = chains(x);
	idBB = lapply(ch, function(ch) {
		tmp = atoms[atoms$chainid == ch, c("resid", "eleid", "elename", "Hetero")];
		tmp = tmp[tmp$Hetero == FALSE, c("resid", "eleid", "elename")];
		if(nrow(tmp) == 0) return(NULL);
		#
		idBB = tapply(tmp, tmp$resid, FUN);
		idBB = unlist(idBB);
		idBB = data.frame(idBB[- length(idBB)], idBB[-1]);
		isNA = is.na(idBB[,1]) | is.na(idBB[,2]);
		idBB = idBB[! isNA, , drop = FALSE];
		names(idBB) = c("eleid.1", "eleid.2");
		return(idBB);
	});
	idBB = do.call(rbind, idBB);
	connect.default(idBB);
}
