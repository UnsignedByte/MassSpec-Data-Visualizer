stopQuietly <- function() {
  opt <- options(show.error.messages = FALSE)
  on.exit(options(opt))
  stop()
}

install_missing <- function(packages, repos='http://cran.us.r-project.org', github=c()) {
	new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
	github.packages <- github[names(github) %in% new.packages]
	if(length(github.packages)){
		if (!('devtools' %in% packages)){
			install_missing('devtools')
			library(devtools)
		}
		sapply(github.packages, install_github)
	}
	if(length(new.packages)) install.packages(setdiff(new.packages, names(github)), repos=repos) #install missing using selected mirror (default cran)
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

getGroups <- function(data, fids, groupids) {
	groups <- matrix(0L, nrow=NROW(data),ncol=length(groupids)); # get final test groups averaged
	for(x in 1:length(groupids)){
		fnames <- fids$ID[fids$Test_Group==groupids[x]] # get file ids in this group
		colnames <- paste("x_OfSpectra", fnames, sep="_");
		selected <- data[,colnames];
		if (NCOL(selected) > 1){
			selected <- rowMeans(selected, na.rm=TRUE);
		}
		groups[,x] = selected;
	}
	return(groups)
}

null_na <- function(x) ifelse(is.null(x), NA, x)

dynamicCeil <- function(x) { # Rounds all numbers up, keeping 1 significant figure
	if(x==0) return(0)
	factor = 10^floor(log10(abs(x))) #smallest power of 10 larger than x
	return(ceiling(x/factor)*factor)
}
dynamicFloor <- function(x) { # Rounds all numbers down, keeping 1 significant figure
	if(x==0) return(0)
	factor = 10^floor(log10(abs(x))) #smallest power of 10 larger than x
	return(floor(x/factor)*factor)
}