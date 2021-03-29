source("utils/rtools.r");

list.packages = c("Rcpp", "stringr", "jsonlite")
install_missing(list.packages)

sourceCpp('utils/parseParams.cpp')

params <- list()

params <- mergeList(parseParams('volcano.r'), params);

if (!('name' %in% names(params))) {
	params$name <- readInput("dataset name:");
}

homedir <- fileExists(file.path("Results",params$name), paste("Dataset", params$name, "cannot be found. Please run combinedHM to generate HeatMap data before this program is run."));

setwd(homedir);


if(dir.exists("Volcano")) {
	unlink("Volcano", recursive=TRUE)
}
dir.create("Volcano")

fids <- read.csv("fileIDs.csv");

dataset.groupids <- unique(fids$Test_Group) #unique test groups

hms <- list.files(path="HeatMap/TestGroups"); #All heatmaps generated

jsonData <- list(VennDiagram = list(), ClusterHeatMap = list())

parsedHM <- str_match(hms, "(.+)\\.csv")[,2]

if (length(dataset.groupids) > 50){
	message("Number of test groups exceeds 50. 2d statistics tests will be calculated unless a list of pairs is provided.")
	stopQuietly();
}

# json list
jsonData <- list()
# all pairs of files
cnames <- combn(dataset.groupids, 2)
message(paste("Generating ", choose(NCOL(cnames), 2), " volcano plots."))

for(hmid in 1:length(hms)){
	hm <- hms[hmid];
	hmname <- unlist(strsplit(hm, ".", fixed=TRUE))[1];
	f <- read.csv(file.path("HeatMap", "TestGroups", hm));
	pvals <- read.csv(file.path("StatTests", paste(hmname, "_", params$statTest, ".csv", sep="")));
	dir.create(file.path("Volcano", hmname))
	for(pairI in 1:NCOL(cnames)){
		pairname <- paste(cnames[1,pairI], "_", cnames[2,pairI], sep="");
		dir.create(file.path("Volcano", hmname, pairname))
		group <- data.frame(matrix(NA, nrow=NROW(f), ncol=4))
		group[,1:2] <- f[,paste("x_OfSpectra", cnames[,pairI], sep="_")]
		group[,3] <- sapply(1:NROW(f), function(i) log2(group[i,2]/group[i,1]))
		group[,4] <- pvals[,paste("X", pairname, sep="")]
		names(group) <- c(paste("x_OfSpectra", cnames[,pairI], sep="_"), "log2foldchange", params$statTest)
		write.csv(group, file=file.path("Volcano", hmname, pairname, "raw.csv"))
	}
}