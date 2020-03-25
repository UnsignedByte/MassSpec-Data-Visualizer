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

library(VennDiagram)

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
# if (length(dataset.groupids) == 1){
# 	print("Not enough test groups to compare!")
# 	stopQuietly();
# }

hms <- list.files(path="HeatMap"); #All heatmaps generated

for(hm in hms){
	f <- read.csv(file.path("HeatMap", hm));
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
	venn.diagram(
	  x = sets,
	  category.names = paste("Testgroup_", dataset.groupids, sep=""),
	  filename = file.path("VennDiagram", paste(unlist(strsplit(hm, ".", fixed=TRUE))[1], "png", sep=".")),
	  output=TRUE
	)
}
