


# TODO: better strategy;
isProtein = function(x) {
	isPr = "CA" %in% x$atoms$elename;
	return(isPr);
}
isProteinA = function(x) {
	isPr = "CA" %in% x$elename;
	return(isPr);
}
# TODO
isNucleicAcid.pdb = function(x) {
	nc = unique(x$atoms$resname);
	isNc = any(nc %in% c("A","C","G", "T", "U"));
	return(isNc);
}
isNucleicAcid.atoms = function(x) {
	nc = unique(x$resname);
	isNc = any(nc %in% c("A","C","G", "T", "U"));
	return(isNc);
}
isNucleicAcid.character = function(x) {
	nc = unique(x);
	isNc = any(nc %in% c("A","C","G", "T", "U"));
	return(isNc);
}
# Special AA:
isAASpecial = function(x) {
	x %in% c("HZP", "PTR", "TYS");
}

### Backbone:
asBackbone = function(x) {
	# TODO: check continuity
	FUN = function(tmp) {
		id = match(c("N", "CA", "C", "O"), tmp$elename);
		id = tmp$eleid[id];
		id = c(id, NA, id[3]);
	}
	FUN.Nc = backbone.nucleic();
	if(inherits(x, "atoms")) {
		atoms = x;
	} else if(inherits(x, "pdb")) atoms = x$atoms;
	# Modified AA: may require better checking if true AA;
	isAA = atoms$Hetero == FALSE;
	isAS = isAASpecial(atoms$resname[atoms$Hetero]);
	isAA[atoms$Hetero] = isAS;
	# Duplicate resid:
	lvl = unique(atoms$insert);
	lvl = sort(lvl[nchar(lvl) > 0]);
	if(length(lvl) == 0) {
		hasDuplicates = FALSE;
	} else {
		hasDuplicates = TRUE;
		lvl = c("", lvl);
		atoms$insert = ordered(atoms$insert, levels = lvl);
	}
	# Process each chain:
	ch = chains(x);
	idBB = lapply(ch, function(ch) {
		isCh = atoms$chainid == ch & isAA;
		tmp  = atoms[isCh, c("resid", "insert", "eleid", "elename")];
		if(nrow(tmp) == 0) return(NULL);
		#
		if(isNucleicAcid.character(atoms$resname[isCh])) {
			FUN = FUN.Nc;
		}
		#
		idGroup = if(hasDuplicates) list(tmp$insert, tmp$resid) else tmp$resid;
		idBB = tapply(tmp, idGroup, FUN, simplify = FALSE);
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

residues = function(x, chains = NULL, hetero = NULL) {
	if(inherits(x, "pdb")) x = x$atoms;
	ch = chains(x);
	if(! is.null(chains)) {
		idCh = match(chains, ch);
		isNA = is.na(idCh);
		if(any(isNA)) {
			warning("The chains: ",
				paste(chains[isNA], collapse = ", "), " do NOT exist!");
			idCh = idCh[! isNA];
		}
		ch = ch[idCh];
		if(length(ch) == 0) return();
	}
	doHetero = ! is.null(hetero);
	tmp = lapply(ch, function(ch) {
		isCh = x$chainid == ch;
		if(doHetero) isCh = isCh & (x$Hetero == hetero);
		tmp = x[isCh, c("resname", "resid"), drop = FALSE];
		tmp = unique(tmp);
		if(nrow(tmp) > 0) tmp$chainid = ch;
		return(tmp);
	})
	tmp = do.call(rbind, tmp);
	return(tmp);
}

# 5' -PO2-O-CH2-C4-C3-O-> 3'
backbone.nucleic = function(x) {
	FUN = function(x) {
		if(nrow(x) == 0) return(NULL);
		atom = x$elename;
		atom = sub("\\*", "'", atom);
		id = match(
			c("O5'", "C5'", "C4'", "C3'", "O3'",
			"P", "O1P", "O2P",  "C2'", "C1'", "O4'"), atom);
		id = x$eleid[id];
		id = c(id[6:7], NA, id[6], id[8], NA, # PO4
			id[6], id[-(6:11)], NA,
			id[4], id[9:11], id[3], NA, # Pentose-Ring
			id[5]);
		return(id);
	}
}

connect.nucleic.pdb = function(x, ...) {
	connect.nucleic.atoms(atoms(x), ...);
}

connect.nucleic.atoms = function(x, ...) {
	atoms = x;
	# Older variants:
	atoms$elename = sub("\\*", "'", atoms$elename);
	# Duplicate resid:
	lvl = unique(atoms$insert);
	lvl = sort(lvl[nchar(lvl) > 0]);
	if(length(lvl) == 0) {
		hasDuplicates = FALSE;
	} else {
		hasDuplicates = TRUE;
		lvl = c("", lvl);
		atoms$insert = ordered(atoms$insert, levels = lvl);
	}
	isNc = atoms$resname %in% c("A","C","G", "T", "U");
	# Process each chain:
	ch = chains(atoms);
	idBB = lapply(ch, function(ch) {
		isCh = atoms$chainid == ch & isNc;
		tmp  = atoms[isCh, c("resid", "insert", "eleid", "elename", "resname")];
		if(nrow(tmp) == 0) return(NULL);
		#
		idGroup = if(hasDuplicates) list(tmp$insert, tmp$resid) else tmp$resid;
		idAll   = tapply(tmp, idGroup, connectNc, simplify = FALSE);
		if(is.null(idAll)) return(NULL);
		# TODO:
		idNc = lapply(idAll, function(x) rbind(x$BB, x$N, x$Jn));
		idNc = do.call(rbind, idNc);
		# BB-Joins:
		idJb = lapply(idAll, function(x) x$Jb);
		idJb = unlist(idJb);
		if(! is.null(idJb) && length(idJb) > 3) {
			idJb = data.frame(
				E1 = idJb[seq(2, length(idJb) - 1, by = 2)],
				E2 = idJb[seq(3, length(idJb), by = 2)]);
			idNc = rbind(idNc, idJb);
		}
		names(idNc) = c("eleid.1", "eleid.2");
		return(idNc);
	});
	idBB = do.call(rbind, idBB);
}

# TODO: RNA vs DNA
connectNc = function(x) {
		if(nrow(x) == 0) return(NULL);
		# BackBone:
		id = match(
			c("C1'", "C2'", "C3'", "C4'", "C5'",
			"O5'", "O4'", "O3'",
			"P", "O1P", "O2P"), x$elename);
		id = x$eleid[id];
		idJ = id[1]; Jb = c(id[9], id[8]);
		idB = data.frame(
			E1 = id[c(9, 9, 9, 1,2,3,4,4,7, 5, 3)],
			E2 = id[c(6,10,11, 2,3,4,5,7,1, 6, 8)]);
		# N-Base:
		nc  = x$resname[1];
		nmN = switch(nc,
			"A" = c("N1", "C2", "N3", "C4", "C5", "C6", # Cycle 1
				"N7", "C8", "N9",    "C6", "N6"), # Cycle 2, Substitutions
			"G" = c("N1", "C2", "N3", "C4", "C5", "C6",
				"N7", "C8", "N9",    "C2", "N2",  "C6", "O6"),
			"T" = c("N1", "C2", "N3", "C4", "C5", "C6",
				"C2", "O2",  "C4", "O4",  "C5", "C7"), # Subst
			"U" = c("N1", "C2", "N3", "C4", "C5", "C6",
				"C2", "O2",  "C4", "O4"),
			"C" = c("N1", "C2", "N3", "C4", "C5", "C6",
				"C2", "O2",  "C4", "N4")
		)
		idN = match(nmN, x$elename);
		idN = x$eleid[idN];
		# Cycle: always C6 - N1;
		cycNc = data.frame(E1 = idN[1:6], E2 = idN[c(2:6, 1)]);
		isAG  = nc == "A" || nc == "G";
		if(isAG) {
			id1 = c(5,7,8,9,  10);
			id2 = c(7,8,9,4,  11);
			if(nc == "G") {
				id1 = c(id1, 12); id2 = c(id2, 13);
			}
		} else {
			id1 = c(7, 9);
			id2 = c(8,10);
			if(nc == "T") {
				id1 = c(id1, 11); id2 = c(id2, 12);
			}
		}
		cyc2  = data.frame(E1 = idN[id1], E2 = idN[id2]);
		cycNc = rbind(cycNc, cyc2);
		idJnc = if(isAG) 9 else 1;
		lst = list(BB = idB, N = cycNc,
			Jb = Jb,
			Jn = data.frame(E1 = idJ, E2 = idN[idJnc]));
		return(lst);
}
