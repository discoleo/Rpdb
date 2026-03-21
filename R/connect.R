#' Create \sQuote{connect} Object
#' 
#' Creates an object of class \sQuote{connect} containing the IDs of bonded atoms
#' defining the connectivity of a molecular system.
#' 
#' \code{connect} is a generic function to create objects of class 
#' \sQuote{connect}. The purpose of this class is to store CONNECT records from 
#' PDB files, indicating the connectivity of a molecular system.\cr
#' The default method creates a \code{connect} object from its different components, i.e.: 
#'   \code{eleid.1} and \code{eleid.2}. Both arguments have to be specified.\cr 
#' The S3 method for object of class \sQuote{coords} determines the connectivity 
#'   from atomic coordinates. A distance matrix is computed, then, for each pair 
#'   of atoms the distance is compared to a bounding distance computed from atomic 
#'   radii. If this distance is lower than the bounding distance then the atoms 
#'   are assumed to be connected.\cr
#' The S3 method for object of class \sQuote{pdb} first uses the element names
#'   to search for atomic radii in the \code{elements} data set.
#'   Then atomic coordinates and radii are passed to \code{connect.coords}.\cr
#'   If \code{by.block == TRUE}, a grid is defined to determined the connectivity by block.
#'   The method is slow but allows to deal with very large systems.\cr
#' The S3 method for object of class \sQuote{character} converts the raw string
#'   from a pdb file into a \code{connect} object.
#' \code{is.connect} tests if an object is of class \sQuote{connect},
#'   i.e. if it has a \dQuote{class} attribute equal to \code{connect}.
#' 
#' @return \code{connect} returns a two-column data.frame of class
#' \sQuote{connect} whose rows contain the IDs of bonded atoms. The columns of
#' this data.frame are described below: \item{eleid.1}{a integer vector
#' containing the elements IDs defining the connectivity of the system.} 
#' \item{eleid.2}{a integer vector containing the elements IDs defining the
#' connectivity of the system.} \cr\cr \code{is.connect} returns TRUE if \code{x}
#' is an object of class \sQuote{coords} and FALSE otherwise.
#' 
#' @param \dots arguments passed to methods.
#' @param eleid.1 a integer vector containing the IDs of bonded atoms.
#' @param eleid.2 a integer vector containing the IDs of bonded atoms.
#' @param x,atoms an R object containing atomic coordinates.
#' @param radii a numeric vector containing atomic radii used to find neigbours.
#' @param safety a numeric value used to extend the atomic radii.
#' @param by.block a logical value indicating whether the connectivity has to be
#'   determine by block (see details).
#' @param pdbRec the raw text lines from the pdb file.
#'   
#' @seealso \code{\link{pdb}}
#' 
#' @examples 
#' # If atom 1 is connected to atom 2, 3, 4 and 5
#' # then we can prepare the following 'connect' object:
#' x <- connect(rep(1,4), 2:5)
#' print(x)
#' is.connect(x)
#' 
#' # Compute connectivity from coordinates
#' x <- read.pdb(system.file("examples/PCBM_ODCB.pdb", package="Rpdb"), CONNECT = FALSE)
#' x$connect
#' x$connect <- connect(x)
#' x$connect
#' 
#' @keywords classes
#' 

#' @name connect
# TODO: remove
conect.default = function(eleid.1, eleid.2, ...) {
	if(missing(eleid.2)) {
		if(inherits(eleid.1, c("connect", "conect"))) {
			return(eleid.1);
		} else if(inherits(eleid.1, "data.frame")) {
			class(eleid.1) = c("connect", "data.frame");
			return(eleid.1);
		} else
			stop("Missing eleid2!");
	}
	connect.default(eleid.1, eleid.2, ...);
}

#' @name connect
#' @export
connect <- function(...)
  UseMethod("connect")

#' @rdname connect
#' @export
connect.default <- function(eleid.1, eleid.2, ...)
{
  if(is.null(eleid.1) & is.null(eleid.2))
    return(NULL)
  eleid.1 <- as.integer(eleid.1)
  eleid.2 <- as.integer(eleid.2)
  con <- data.frame(eleid.1, eleid.2)
#   con[con$eleid.1 > con$eleid.2,] <- rev(con[con$eleid.1 > con$eleid.2,])
  con <- con[order(con$eleid.2),]
  con <- con[order(con$eleid.1),]
  rownames(con) <- NULL
  class(con) = c("connect", "data.frame");
  if(nrow(con) == 0)
    con <- NULL
  return(con)
}

#' @rdname connect
#' @export
connect.coords <- function(x, radii = 0.75, safety = 1.2, by.block = FALSE, ...) {
  if(!is.coords(x))
    stop("'x' must be an object of class 'coords'")

  radii <- radii*safety
  data <- cbind(x, radii)

  findCon <- function(data) {
    nat <- nrow(data)
    if(nat==0) return(NULL)
    r <- sqrt(
        outer(data$x1, data$x1, "-")^2 +
        outer(data$x2, data$x2, "-")^2 +
        outer(data$x3, data$x3, "-")^2
    )
    bond.dist <- outer(data$radii, data$radii, "+")
    M <- lower.tri(r) & (r < bond.dist)
    if(all(!M)) return(NULL)
    eleid <- matrix(rownames(data), nrow = nat, ncol = nat)   
    eleid.1 <- as.integer(t(eleid)[M])
    eleid.2 <- as.integer(  eleid [M]) 
    return(connect.default(eleid.1, eleid.2))
  }

  if(!by.block) {
    con <- findCon(data)
  }
  else {
    get.con <- function(shift = c(0,0,0), x, radii, step) {
      coords.range <- range(x)
      coords.range <- t(t(coords.range) - shift)

      x.cuts <- seq(coords.range["min","x"], coords.range["max","x"] + step, step)
      y.cuts <- seq(coords.range["min","y"], coords.range["max","y"] + step, step)
      z.cuts <- seq(coords.range["min","z"], coords.range["max","z"] + step, step)
    
      x.index <- cut(x$x1, x.cuts, include.lowest = TRUE)
      y.index <- cut(x$x2, y.cuts, include.lowest = TRUE)
      z.index <- cut(x$x3, z.cuts, include.lowest = TRUE)
      
      f <- interaction(x.index, y.index, z.index)
      
      data <- cbind(x, radii)
      
      con <- split(data, f)
      con <- lapply(con, findCon)
      con <- do.call(rbind, con)

      return(con)
    }

    width <- 10
    shift <- expand.grid(0:1,0:1,0:1) * width/2;
	
	# TODO: check of correct
    con = apply(shift, 1, get.con, x, radii, width)
    con = unique(do.call(rbind, con))
    con = connect(con$eleid.1, con$eleid.2);
    rownames(con) = NULL;
  }
  
  return(con)
}

#' @rdname connect
#' @export
connect.pdb <- function(x, safety = 1.2, by.block = FALSE, ...) {
  symb <- toSymbols(x$atom$elename)
  symb[is.na(symb)] <- "Xx"
  rcov <- Rpdb::elements[match(symb, Rpdb::elements[,"symb"]), "rcov"]
  con <- connect(coords(x), rcov, safety, by.block)
  con <- connect(x$atoms$eleid[con$eleid.1], x$atoms$eleid[con$eleid.2])
  return(con)
}


#' @rdname connect
#' @export
connect.character = function(pdbRec, atoms, ...) {
	C0 = as.integer(substr(pdbRec,  7, 11));
	C1 = as.integer(substr(pdbRec, 12, 16));
	C2 = as.integer(substr(pdbRec, 17, 21));
	C3 = as.integer(substr(pdbRec, 22, 26));
	C4 = as.integer(substr(pdbRec, 27, 31));
	C0 = rep(C0, 4);
	C1 = c(C1,C2,C3,C4);
	C0 = subset(C0, !is.na(C1));
	C1 = subset(C1, !is.na(C1));
	connect = connect.default(eleid.1 = C0, eleid.2 = C1);
	tokeep  = connect$eleid.1 %in% atoms$eleid & connect$eleid.2 %in% atoms$eleid;
	connect = subset(connect, tokeep);
	return(connect);
}

#' @rdname connect
#' @export
is.connect <- function(x)
{
  to.return = inherits(x, c("connect", "conect"));
  return(to.return)
}
