#' Centres-of-Geometry: Rolling Window
#' 
#' Computes centres-of-geometry using a rolling window on a polypeptide chain
#' 
#' \sQuote{centres.pp.roll} is a generic function to compute the centres-of-geometry
#' using a rolling window on an object containing atomic coordinates for
#' a polypeptide chain.
#'
#'
#' @usage centres.pp.roll(x, window, ...)
#'   
#' @param x an R object containing atomic coordinates.
#' @param window size of the rolling window, specified as number of amino-acids;
#' @param chain apply only to the respective chains, by default all chains;
#' @param na.rm a logical value indicating whether NA values should be stripped
#'   before the computation proceeds.
#' @param \dots further arguments passed to or from other methods.
#' 
#' @seealso \code{\link{coords}}, \code{\link{atoms}}, \code{\link{pdb}}.
#'
#' @return Return an object of class \sQuote{coords} containing the coordinates
#'   of the centres.
#' 

#' @name centres.pp.roll
#' @export centres.pp.roll
centres.pp.roll <- function(x, window, ...) {
	UseMethod("centres.pp.roll");
}

#' @rdname centres.pp.roll
#' @exportS3Method centres.pp.roll default
centres.pp.roll.default = function(x, window = 34, na.rm = TRUE, ...) {
	ids = unique(x$resid[x$elename == "CA"]);
	len = length(ids);
	if(len == 0) {
		xyz = data.frame(x = numeric(0), y = numeric(0), z = numeric(0));
		return(as.coords(xyz));
	}
	if(len <= window) {
		xyz = data.frame(
			x = mean(x$x1, na.rm=na.rm),
			y = mean(x$x2, na.rm=na.rm),
			z = mean(x$x3, na.rm=na.rm));
		return(as.coords(xyz));
	}
	ids = sort(ids);
	maxID = ids[len] - window;
	xyz = lapply(seq(ids[1], maxID, by = 1), function(id) {
		idw = which(ids >= id & ids <= id + window);
		idw = ids[idw];
		isW = x$resid %in% idw;
		xyz = data.frame(
			x = mean(x$x1[isW], na.rm=na.rm),
			y = mean(x$x2[isW], na.rm=na.rm),
			z = mean(x$x3[isW], na.rm=na.rm));
	});
	xyz = do.call(rbind, xyz);
	return(invisible(as.coords(xyz)));
}

#' @rdname centres.pp.roll
#' @method centres.pp.roll atoms
#' @exportS3Method centres.pp.roll atoms
centres.pp.roll.atoms = function(x, window = 34, chain = NULL, na.rm = TRUE, ...) {
	if(is.null(chain)) {
		chain = unique(x$chainid);
	}
	FUN = function(chain) {
		xyz = x[x$chainid == chain, c("x1","x2","x3", "elename", "resid"),
			drop = FALSE];
		xyz = centres.pp.roll.default(xyz, window=window, na.rm=na.rm, ...);
		if(nrow(xyz) > 0) xyz$chainid = chain;
		return(xyz);	
	}
	xyz = lapply(chain, FUN);
	xyz = do.call(rbind, xyz);
	invisible(xyz);
}
