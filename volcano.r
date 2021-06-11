source("utils/rtools.r");

list.packages = c("Rcpp", "stringr", "jsonlite", "ggplot2")
install_missing(list.packages)

sourceCpp('utils/parseParams.cpp')

params <- list(
	wantedCol="x_OfSpectra",
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

jsonData <- list();

parsedHM <- str_match(hms, "(.+)\\.csv")[,2]

if (length(dataset.groupids) > 50){
	message("Number of test groups exceeds 50. 2d statistics tests will be calculated unless a list of pairs is provided.")
	stopQuietly();
}

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
	f <- f[f$Row_Type==1,]; # take only class values
	rownames(f) <- NULL; # Reset row names
	pvals <- read.csv(file.path("StatTests", paste(hmname, "_", params$statTest, ".csv", sep="")));
	pvals <- pvals[pvals$Row_Type==1,]; # take only class values
	rownames(pvals) <- NULL; # Reset row names
	significance <- read.csv(file.path("Significance", hmname, "raw.csv"));
	significance$volcano = vector(mode="character", length=NROW(significance));
	# print(significance)
	dir.create(file.path("Volcano", hmname))
	jsonData$Volcano[[hmid]] <- list(name = hmname, graph=list(), raw=list());
	for(pairI in 1:NCOL(cnames)){
		pairname <- paste(cnames[1,pairI], cnames[2,pairI], sep="_");
		message(paste("Generating raws for groups", cnames[1,pairI], "and", cnames[2,pairI]));
		dir.create(file.path("Volcano", hmname, pairname))
		group <- data.frame(matrix(NA, nrow=NROW(f), ncol=6))
		group[,1:2] <- f[,paste(params$wantedCol, cnames[,pairI], sep="_")]
		group[,3] <- sapply(1:NROW(f), function(i) log2(group[i,2]/group[i,1]))
		group[which(!is.finite(group[,3])),3] <- NA;
		group[,4] <- pvals[,paste("X", pairname, sep="")]
		group[,5] <- sapply(group[,4], function(i) -log10(i));
		group[,6] <- sapply(1:NROW(f), function(i) classify(group[i,3],group[i,5])) # Label upregulated and downregulated

		for(i in 1:NROW(f)) {
			significance[significance$Rank_Number==i,"volcano"] = paste(significance[significance$Rank_Number==i,"volcano"], switch(group[i,6],
				"upregulated" = paste("Classified upregulated between group pairs ", cnames[1,pairI], "&", cnames[2,pairI], "\n", sep=""),
				"downregulated" = paste("Classified downregulated between groups pairs ", cnames[1,pairI], "&", cnames[2,pairI], "\n", sep=""),
				"none" = ""
				), sep="");
		}

		names(group) <- c(paste(params$wantedCol, cnames[,pairI], sep="_"), "log2foldchange", params$statTest, paste("-log10(", params$statTest, ")", sep=""), "significance")
		message("Creating plot")
		outsvg <- file.path("Volcano", hmname, pairname, "plot.svg");
		ggsave(
			file=outsvg,
			plot=ggplot(data=group, aes_string(x=names(group)[3], y=paste("-log10(", params$statTest, ")", sep=""), color="significance"))
					+geom_point()
					+scale_y_continuous(limits=c(0,dynamicCeil(max(group[,5], na.rm=TRUE))),expand=c(0,0))
					+scale_x_continuous(limits=c(dynamicFloor(min(group[,3], na.rm=TRUE)), dynamicCeil(max(group[,3], na.rm=TRUE))),expand=c(0,0))
					+scale_colour_manual(values = c("#69a048","red","#999998","#69a048"))
					+ggtitle(paste("Volcano plot comparing", params$wantedCol, "for test groups", cnames[1,pairI], "and", cnames[2,pairI]))
		)
		jsonData$Volcano[[hmid]]$graph[[pairname]] <- readChar(outsvg, file.info(outsvg)$size)
		write.csv(group, file=file.path("Volcano", hmname, pairname, "raw.csv"), row.names=FALSE)
		jsonData$Volcano[[hmid]]$raw[[pairname]] <- group;
		message("Finished")
	}
	write.csv(significance, file=file.path("Significance", hmname, "raw.csv"), row.names=FALSE);
}
write(toJSON(jsonData, auto_unbox=TRUE), file=file.path("Raws", "volcano.json"));