stopQuietly <- function() {
  opt <- options(show.error.messages = FALSE)
  on.exit(options(opt))
  stop()
}

install_missing <- function(packages, repos='http://cran.us.r-project.org') {
	new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
	if(length(new.packages)) install.packages(new.packages, repos=repos) #install missing using selected mirror (default cran)
	lapply(packages, library, character.only = TRUE) # import installed
}

gDims <- function(a) {
	return(c(NROW(a),NCOL(a)));
}

compareNA <- function(cond) {
	# This function returns TRUE wherever elements are the same, including NA's,
	# and false everywhere else.
	cond[is.na(cond)] <- FALSE
	return(cond)
}

# recursiveBinary <- function(n, l) { #Finds binary numbers in string form given length l and n on bits (n<=l)
# 	# print(paste(n, " ", l))
# 	if (n==0){
# 		return(0);
# 	}else if (n == l){
# 		return(bitwShiftL(1,l)-1);
# 	}else{
# 		return(c(recursiveBinary(n,l-1), bitwShiftL(1,l-1)+recursiveBinary(n-1,l-1)));
# 	}
# }

recursiveBinary <- function(l) {
	return(sapply(1:2^l, function(x) substrRight(dec2bin(x-1), l)))
}

dec2bin <- function(x) paste(as.integer(rev(intToBits(x))), collapse = "") # decimal to binary string

substrRight <- function(x, n){ #Substring from right
  substr(x, nchar(x)-n+1, nchar(x))
}

readInput <- function(st) { #Reads input
	inp <- "";
	if(interactive()){
		inp <- readline(st);
	}else{
		fin <- file("stdin")
		open(fin)
		cat(st)
		inp <- readLines(fin,1);
		close(fin)

	}
	return(inp);
}


fileExists <- function(dname, errorMessage="file not found.") {
	if (!dir.exists(dname)){
		print(errorMessage)
		stopQuietly()
	}
	return(dname)
}

dist_no_na <- function(mat) {
    edist <- dist(mat)
    edist[which(is.na(edist)|!is.finite(edist))] <- max(edist, na.rm=TRUE) * 1.1 
    return(edist)
}

cv <- function(dat, na.rm=FALSE) {
	return(sd(dat, na.rm=na.rm)/mean(dat, na.rm=na.rm))
}

mergeList <- function(original, new) {
	# Will add fields of new to original if they dont exist in original
	for (i in names(original)) {
		new[[i]] = original[[i]]
	}
	return(new)
}