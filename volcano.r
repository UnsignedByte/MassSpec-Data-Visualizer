source("utils/rtools.r");

list.packages = c("Rcpp", "stringr", "jsonlite", "ggplot2")
install_missing(list.packages)

sourceCpp('utils/parseParams.cpp')

params <- list(
	foldchangethreshold = 1,
	pthreshold = 0.05
)

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
message(paste("Generating", NCOL(cnames), "volcano plots."))

classify <- function(foldchange,pval){ # filter out up/downregulated
	if (is.na(foldchange) || is.na(pval)) return("invalid")
	if (10^(-pval) < params$pthreshold){
		if(foldchange < -params$foldchangethreshold){
			return("upregulated")
		}else if (foldchange > params$foldchangethreshold) {
			return("downregulated")
		}
	}
	return("none")
}

for(hmid in 1:length(hms)){
	hm <- hms[hmid];
	hmname <- unlist(strsplit(hm, ".", fixed=TRUE))[1];
	f <- read.csv(file.path("HeatMap", "TestGroups", hm));
	pvals <- read.csv(file.path("StatTests", paste(hmname, "_", params$statTest, ".csv", sep="")));
	dir.create(file.path("Volcano", hmname))
	for(pairI in 1:NCOL(cnames)){
		pairname <- paste(cnames[1,pairI], cnames[2,pairI], sep="_");
		message(paste("Generating raws for groups", cnames[1,pairI], "and", cnames[2,pairI]));
		dir.create(file.path("Volcano", hmname, pairname))
		group <- data.frame(matrix(NA, nrow=NROW(f), ncol=6))
		group[,1:2] <- f[,paste("x_OfSpectra", cnames[,pairI], sep="_")]
		group[,3] <- sapply(1:NROW(f), function(i) log2(group[i,2]/group[i,1]))
		group[which(!is.finite(group[,3])),3] <- NA;
		group[,4] <- pvals[,paste("X", pairname, sep="")]
		group[,5] <- sapply(group[,4], function(i) -log10(i));
		group[,6] <- sapply(1:NROW(f), function(i) classify(group[i,3],group[i,5]))
		names(group) <- c(paste("x_OfSpectra", cnames[,pairI], sep="_"), "log2foldchange", params$statTest, paste("-log10(", params$statTest, ")", sep=""), "significance")
		message("Creating plot")
		pdf(file.path("Volcano", hmname, pairname, "plot.pdf"))
		print(
			ggplot(data=group, aes_string(x=names(group)[3], y=paste("-log10(", params$statTest, ")", sep=""), color="significance"))
			+geom_point()
			+scale_y_continuous(limits=c(0,dynamicCeil(max(group[,5], na.rm=TRUE))),expand=c(0,0))
			+scale_x_continuous(limits=c(dynamicFloor(min(group[,3], na.rm=TRUE)), dynamicCeil(max(group[,3], na.rm=TRUE))),expand=c(0,0))
			+scale_colour_manual(values = c("black","red","grey","black"))
		)
		dev.off()
		write.csv(group, file=file.path("Volcano", hmname, pairname, "raw.csv"))
		message("Finished")
	}
}