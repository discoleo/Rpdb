

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

