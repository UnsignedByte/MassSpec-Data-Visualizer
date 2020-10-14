source("utils/rtools.r");

list.packages = c("stats", "utils", "Rcpp", "stringr", "jsonlite")
install_missing(list.packages)

sourceCpp('utils/parseParams.cpp')

params <- list()

params$stats <- list(
	"wilcoxon"=function(x,y) tryCatch(wilcox.test(x,y)$p.value, error=function(cond) return(NaN))
)

params <- mergeList(parseParams('statTests.r'), params);

if (!('name' %in% names(params))) {
	params$name <- readInput("dataset name:");
}

homedir <- fileExists(file.path("Results",params$name), paste("Dataset", params$name, "cannot be found. Please run combinedHM to generate HeatMap data before this program is run."));

setwd(homedir);

dir.create("statTests", showWarnings = FALSE)

fids <- read.csv("fileIDs.csv");

dataset.groupids <- unique(fids$Test_Group) #unique test groups

hms <- list.files(path="HeatMap"); #All heatmaps generated

jsonData <- list(VennDiagram = list(), ClusterHeatMap = list())

parsedHM <- str_match(hms, "(.+)\\.csv")[,2]

if (length(dataset.groupids) > 10){
	message("Number of test groups exceeds 10. 2d statistics tests cannot be calculated unless a list of pairs is provided.")
	stopQuietly();
}

# json list
jsonData <- list(statTests=list())
# all pairs of files
cnames <- combn(dataset.groupids, 2)



for(hmid in 1:length(hms)){
	hm <- hms[hmid]
	f <- read.csv(file.path("HeatMap", hm));
	groups <- list();

	for(x in 1:length(dataset.groupids)){
		fnames <- fids$ID[fids$Test_Group==dataset.groupids[x]] # get file ids in this group
		colnames <- paste("x_OfSpectra", fnames, sep="_");
		selected <- f[,colnames];
		groups[[x]] <- as.data.frame(selected);
	}

	statTables <- list()

	for (stat in names(params$stats)) {
		statTables[[stat]] <- data.frame(matrix(NA, nrow = NROW(f), ncol = NCOL(cnames)));
		for(col in 1:NCOL(cnames)) {
			pair <- cnames[,col];
			colnames(statTables[[stat]])[[col]] <- paste(pair, collapse="_");
			for(row in 1:NROW(f)) {
				p <- lapply(pair, function(i) as.numeric(groups[[i]][row,]))
				statTables[[stat]][row,col] <- params$stats[[stat]](p[[1]], p[[2]]);
			}
		}
		write.csv(statTables[[stat]], file=file.path("StatTests", paste(unlist(strsplit(hm, ".", fixed=TRUE))[1], '_', stat, '.csv', sep='')))
	}
	jsonData$statTests[[hm]] <- statTables;
	# print(statTables);
}
write(toJSON(jsonData, auto_unbox=TRUE), file=file.path("Raws", "statTests.json"));