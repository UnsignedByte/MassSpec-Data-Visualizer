source("utils/rtools.r");

list.packages = c("Rcpp", "stringr", "jsonlite", "ggplot2")
install_missing(list.packages)

sourceCpp('utils/parseParams.cpp')

params <- list(
	wantedCol="x_OfSpectra",
	foldchangethreshold = 1,
	pthreshold = 0.05
)

params <- mergeList(parseParams('violin.r'), params);

if (!('name' %in% names(params))) {
	params$name <- readInput("dataset name:");
}

homedir <- fileExists(file.path("Results",params$name), paste("Dataset", params$name, "cannot be found. Please run combinedHM to generate HeatMap data before this program is run."));

setwd(homedir);

if(dir.exists("Violin")) {
	unlink("Violin", recursive=TRUE)
}
dir.create("Violin")

fids <- read.csv("fileIDs.csv");

dataset.groupids <- unique(fids$Test_Group) #unique test groups

hms <- list.files(path="HeatMap/TestGroups"); #All heatmaps generated

jsonData <- list();

parsedHM <- str_match(hms, "(.+)\\.csv")[,2]

if (length(dataset.groupids) > 50){
	message("Number of test groups exceeds 50. 2d statistics tests will be calculated unless a list of pairs is provided.")
	stopQuietly();
}

for(hmid in 1:length(hms)){
	hm <- hms[hmid];
	hmname <- unlist(strsplit(hm, ".", fixed=TRUE))[1];
	f <- read.csv(file.path("HeatMap", "Files", hm));
	significance <- read.csv(file.path("Significance", hmname, "raw.csv"), stringsAsFactors=FALSE);
	significance[is.na(significance)] <- "";
	dir.create(file.path("Violin", hmname))
	jsonData$Violin[[hmid]] <- list(name = hmname, graph=list(), raw=list());

	for(idx in 1:NROW(f))
	{
		if (f$Row_Type[[idx]] == 1)
		{
			rank <- f$Rank[[idx]];

			nname = gsub("/", ",", significance$Gene_Name[rank]);

			dataset <- data.frame(count=NA, group=NA);
			for(fidx in 1:NROW(fids))
			{
				dataset[fidx,] <- c(f[idx, paste(params$wantedCol, fidx, sep="_")], fids$Test_Group[[fidx]]);
			}

			if (any(nchar(as.character(significance[rank,-(1:3)]))))
			{
				dataset$group = as.factor(dataset$group)
				dir.create(file.path("Violin", hmname, nname), recursive=TRUE);
				outsvg <- file.path("Violin", hmname, nname, "plot.svg");
				# print(dataset);
				ymax <- dynamicCeil(max(dataset$count, na.rm=TRUE));
				ggsave(
					file=outsvg,
					plot=ggplot(data=dataset, aes(x=group, y=count))
							+geom_violin(trim=FALSE)
							+geom_dotplot(binaxis='y', stackdir='center', binwidth = 1, dotsize=0.3)
							+scale_y_continuous(limits=c(0,ymax),expand=c(0,0))
							# +scale_x_discrete(limits=c(0,length(dataset.groupids)))
				)
				jsonData$Violin[[hmid]]$graph[[significance$Gene_Name[rank]]] <- readChar(outsvg, file.info(outsvg)$size);
			}
		}
	}
}
write(toJSON(jsonData, auto_unbox=TRUE), file=file.path("Raws", "violin.json"));