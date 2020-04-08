stopQuietly <- function() {
  opt <- options(show.error.messages = FALSE)
  on.exit(options(opt))
  stop()
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

recursiveBinary <- function(n, l) { #Finds binary numbers in string form given length l and n on bits (n<=l)
	# print(paste(n, " ", l))
	if (n==0){
		return(0);
	}else if (n == l){
		return(bitwShiftL(1,l)-1);
	}else{
		return(c(recursiveBinary(n,l-1), bitwShiftL(1,l-1)+recursiveBinary(n-1,l-1)));
	}
}

dec2bin <- function(x) paste(as.integer(rev(intToBits(x))), collapse = "") # decimal to binary string

substrRight <- function(x, n){ #Substring from right
  substr(x, nchar(x)-n+1, nchar(x))
}

library(VennDiagram)
library(RColorBrewer)

futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger");

fin <- file("stdin")
open(fin)

if(interactive()){
	dataset.name <- readline("Dataset name: ")
}else{
	cat("Dataset name:")
	dataset.name <- readLines(fin,1)
}

close(fin)

homedir <- file.path("Results",dataset.name);

if (!dir.exists(homedir)){
	print(paste("Dataset", dataset.name, "cannot be found. Please run combinedHM to generate HeatMap data before this program is run."))
	stopQuietly()
}

setwd(homedir);

cutoff <- 0.8; # the proportion of the max in the test groups at which it cuts off including the other rows

fids <- read.csv("fileIDs.csv");

dataset.groupids <- unique(fids$Test_Group) #unique test groups
palette <- brewer.pal(length(dataset.groupids), "Pastel2")
# if (length(dataset.groupids) == 1){
# 	print("Not enough test groups to compare!")
# 	stopQuietly();
# }

hms <- list.files(path="HeatMap"); #All heatmaps generated

for(hm in hms){
	f <- read.csv(file.path("HeatMap", hm));
	f <- f[f$Row_Type==1,];
	groups <- matrix(0L, nrow=NROW(f),ncol=length(dataset.groupids)); # get final test groups averaged
	for(x in 1:length(dataset.groupids)){
		fnames <- fids$ID[fids$Test_Group==dataset.groupids[x]] # get file ids in this group
		colnames <- paste(rep("x_OfSpectra", length(fnames)), fnames, sep="_");
		selected <- f[,colnames];
		if (NCOL(selected) > 1){
			selected <- rowMeans(selected, na.rm=TRUE);
		}
		groups[,x] = selected;
	}

	if (length(dataset.groupids) > 1){
		maxgroup <- apply(groups, 1, function(x) max(x, na.rm=TRUE));
	}else{
		maxgroup <- groups;
	}
	# print(dim(groups))
	# sapply(1:length(dataset.groupids), 1, function(x) groups[>cutoff*maxgroup)
	sets <- lapply(seq_len(length(dataset.groupids)), function(x) which(compareNA(groups[,x]>cutoff*maxgroup)));
	# sets <- lapply(sets, function(x) find(x));
	# sets = sets==TRUE;
	dir.create("VennDiagram", showWarnings = FALSE)

	v <- venn.diagram(
		x = sets,
		category.names = paste("Testgroup_", dataset.groupids, sep=""),
		filename = file.path("VennDiagram", paste(unlist(strsplit(hm, ".", fixed=TRUE))[1], "png", sep=".")),
		output=TRUE,

		# Output features
		imagetype="png",
		height = 1800,
		width = 1800,
		resolution = 300,
		compression = "lzw",
		units = 'px',

		# Circles
        lwd = 2,
        lty = 'solid',
        fill = palette,
        
        # Numbers
        cex = 0.6,
        fontface = "bold",
        fontfamily = "sans",
        
        # Set names
        cat.cex = 0.6,
        cat.fontface = "bold",
        cat.fontfamily = "sans",
        cat.default.pos = "outer"
	)

	overlap <- calculate.overlap(sets);
	cnames <- c();

	for(i in 1:length(dataset.groupids)){
		cnames <- c(cnames, recursiveBinary(i,length(dataset.groupids)));
	}
	cnames <- sapply(cnames, function(x){substrRight(dec2bin(x),length(dataset.groupids))})
	names(overlap) <- cnames;
	max_l <- max(lengths(overlap))
	overlapdf <- rapply(overlap, function(x) 'length<-'(x, max_l), how="list")
	write.csv(overlapdf, file=file.path("VennDiagram", paste(unlist(strsplit(hm, ".", fixed=TRUE))[1], "csv", sep=".")))
}
