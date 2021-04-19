source("utils/rtools.r");

list.packages = c("stats", "utils", "Rcpp", "stringr", "jsonlite")
install_missing(list.packages)

sourceCpp('utils/parseParams.cpp')

params <- list()

params$twoStats <- list( # stats comparing 2 test groups
	"wilcoxon"=function(x,y) tryCatch(wilcox.test(x,y)$p.value, error=function(cond) return(NaN))
)

params$multiStats <- list( # stats comparing value to test group
	"anova"=function(x,y) null_na(summary(aov(x~y))[[1]][1,5])
)

params <- mergeList(parseParams('statTests.r'), params);

if (!('name' %in% names(params))) {
	params$name <- readInput("dataset name:");
}

homedir <- fileExists(file.path("Results",params$name), paste("Dataset", params$name, "cannot be found. Please run combinedHM to generate HeatMap data before this program is run."));

setwd(homedir);

dir.create("StatTests", showWarnings = FALSE)
dir.create("StatsHeatMap", showWarnings = FALSE)

fids <- read.csv("fileIDs.csv");

dataset.groupids <- unique(fids$Test_Group) #unique test groups

hms <- list.files(path="HeatMap/Files"); #All heatmaps generated

parsedHM <- str_match(hms, "(.+)\\.csv")[,2]

if (length(dataset.groupids) > 50){
	message("Number of test groups exceeds 50. 2d statistics tests will be calculated unless a list of pairs is provided.")
	stopQuietly();
}

# json list
jsonData <- list()
# all pairs of files
cnames <- combn(dataset.groupids, 2)

for(hmid in 1:length(hms)){
	hm <- hms[hmid]
	hmname <- unlist(strsplit(hm, ".", fixed=TRUE))[1];
	f <- read.csv(file.path("HeatMap", "Files", hm));
	groups <- list();

	for(x in 1:length(dataset.groupids)){
		fnames <- fids$ID[fids$Test_Group==dataset.groupids[x]] # get file ids in this group
		colnames <- paste("x_OfSpectra", fnames, sep="_");
		selected <- f[,colnames];
		groups[[x]] <- as.data.frame(selected);
	}

	colnames <- paste("x_OfSpectra", 1:length(fids$ID), sep="_")
	ungrouped <- f[,colnames];
	colnames(ungrouped) <- fids$Test_Group;

	statTables <- list()

	for (stat in names(params$twoStats)) {
		statTables[[stat]] <- data.frame(matrix(NA, nrow = NROW(f), ncol = NCOL(cnames)));
		for(col in 1:NCOL(cnames)) {
			pair <- cnames[,col];
			colnames(statTables[[stat]])[[col]] <- paste(pair, collapse="_");
			for(row in 1:NROW(f)) {
				p <- lapply(pair, function(i) as.numeric(groups[[i]][row,]))
				statTables[[stat]][row,col] <- params$twoStats[[stat]](p[[1]], p[[2]]);
			}
		}
		statTables[[stat]] <- cbind(f[c('Rank_Number','Protein_Name','Gene_Name')], statTables[[stat]], f['Row_Type'])
		write.csv(statTables[[stat]], file=file.path("StatTests", paste(hmname, '_', stat, '.csv', sep='')))
	}

	multiName <- "MultiDim"

	statTables[[multiName]] <- data.frame(matrix(NA, nrow = NROW(f), ncol = length(params$multiStats)));
	for (col in 1:length(params$multiStats)) {
		stat <- names(params$multiStats)[[col]];
		colnames(statTables[[multiName]])[[col]] <- stat;
		for(row in 1:NROW(f)) {
			statTables[[multiName]][row,col] <- params$multiStats[[stat]](as.numeric(colnames(ungrouped)), as.numeric(ungrouped[row,]));
		}
	}
	write.csv(cbind(f[,names(mtcars)!="Row_Type"], statTables[[multiName]], f['Row_Type']), file=file.path("StatsHeatMap", paste(hmname, '.csv', sep='')))
	statTables[[multiName]] <- cbind(f[c('Rank_Number','Protein_Name','Gene_Name')], statTables[[multiName]], f['Row_Type'])
	write.csv(statTables[[multiName]], file=file.path("StatTests", paste(hmname, '_', multiName, '.csv', sep='')))
	jsonData$StatTests[[hmid]] <- list(name=hmname, data=statTables);
	# print(statTables);
}
write(toJSON(jsonData, auto_unbox=TRUE), file=file.path("Raws", "statTests.json"));
