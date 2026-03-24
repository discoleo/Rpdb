


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
		id = c(id, id[3]);
	}
	atoms = x$atoms;
	ch = chains(x);
	idBB = lapply(ch, function(ch) {
		tmp = atoms[atoms$chainid == ch, c("resid", "eleid", "elename")];
		if(! isProteinA(tmp)) return(NULL);
		idBB = tapply(tmp, tmp$resid, FUN);
		idBB = unlist(idBB);
		idBB = data.frame(idBB[- length(idBB)], idBB[-1]);
		names(idBB) = c("eleid.1", "eleid.2");
		return(idBB);
	});
	idBB = do.call(rbind, idBB);
	connect.default(idBB);
}
