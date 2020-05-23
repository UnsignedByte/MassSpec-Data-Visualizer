library(VennDiagram)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(measurements)

source("utils/rtools.r");

futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger");

dataset.name <- readInput("Dataset name:");

homedir <- fileExists(file.path("Results",dataset.name), paste("Dataset", dataset.name, "cannot be found. Please run combinedHM to generate HeatMap data before this program is run."));

setwd(homedir);

cutoff <- 0.8; # the proportion of the max in the test groups at which it cuts off including the other rows
# threshold between counted/not in venn

fids <- read.csv("fileIDs.csv");

dataset.groupids <- unique(fids$Test_Group) #unique test groups
palette <- brewer.pal(length(dataset.groupids), "Pastel2")
# if (length(dataset.groupids) == 1){
# 	print("Not enough test groups to compare!")
# 	stopQuietly();
# }

hms <- list.files(path="HeatMap"); #All heatmaps generated

imsize <- 2800
hmcount <- 64 #number of rows to take for heatmap
colors <- colorRamp2(c(0, hmcount/2, hmcount), c("blue", "white", "red"))(seq(0, hmcount));
linespacing <- 1.5;
mm2in <- 0.0393701;

for(hm in hms){
	f <- read.csv(file.path("HeatMap", hm));
	f <- f[f$Row_Type==1,]; #remove contaminated
	groups <- matrix(0L, nrow=NROW(f),ncol=length(dataset.groupids)); # get final test groups averaged
	for(x in 1:length(dataset.groupids)){
		fnames <- fids$ID[fids$Test_Group==dataset.groupids[x]] # get file ids in this group
		colnames <- paste("x_OfSpectra", fnames, sep="_");
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
	dir.create("ClusterHeatMap", showWarnings = FALSE)

	f <- f[,paste("x_OfSpectra", fids$ID,sep="_")]
	colnames(f) <- fids$ID;

	colnames(groups) <- dataset.groupids;
	# print(names(groups))

	combinedf <- data.matrix(cbind(f, groups))
	# print(percentile(c(0,0,0,1,2,3)))
	combinedf <- combinedf[1:hmcount,]
	combinedf <- apply(combinedf,2,rank, na.last = "keep", ties.method = "max");
	# combinedf <- log2(combinedf)
	# print(combinedf)

	outfname <- paste(unlist(strsplit(hm, ".", fixed=TRUE))[1], "pdf", sep=".");

	rowgp <- gpar(fontsize = 7);
	colgp <-gpar(fontsize = 7);

	minwidth <- unit(6, "cm");
	minheight <- unit(10, "cm");

	hmheight <- max(max_text_height(
	        rownames(combinedf), 
	        gp = colgp
	    )*hmcount*linespacing, minheight)

	hm2 <- Heatmap(
		combinedf[,1:NCOL(f)],
		name="File Id",
		column_title="File Id",
		row_title="Peptide Rank",
		col=colors,
		row_names_gp = rowgp,
		column_names_gp = colgp,
		width = max(max_text_height(
	        rownames(combinedf), 
	        gp = rowgp
	    )*NCOL(f)*linespacing, minwidth),
		height=hmheight
	)+Heatmap(
		combinedf[,NCOL(f)+1:NCOL(groups)],
		name="Test Group",
		column_title="Test Group",
		col=colors,
		row_names_gp = rowgp,
		column_names_gp = colgp,
		width = max(max_text_height(
	        rownames(combinedf), 
	        gp = rowgp
	    )*NCOL(groups)*linespacing, minwidth),
		height=hmheight
	)
	# print(attributes(hm2))
	plot <- draw(hm2)
	# attr(plot, "layout")$page_size[1] <- attr(plot, "ht_list_param")$width

	hmw <- as.double(ComplexHeatmap:::width(plot))*mm2in
	hmh <- as.double(ComplexHeatmap:::height(plot))*mm2in
	# print(hmw)
	# print(conv_unit(303.1907, "mm", "inch"))
	pdf(file= file.path("ClusterHeatMap", outfname), 
		width=hmw, 
		height=hmh)
	print(plot)
	dev.off()

	outfname <- paste(unlist(strsplit(hm, ".", fixed=TRUE))[1], "png", sep=".");

	v <- venn.diagram(
		x = sets,
		category.names = paste("Testgroup_", dataset.groupids, sep=""),
		filename = file.path("VennDiagram", outfname),
		output=TRUE,

		# Output features
		imagetype="png",
		height = imsize,
		width = imsize,
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
	cnames <- sapply(cnames, function(x){substrRight(dec2bin(x),length(dataset.groupids))});
	names(overlap) <- cnames;
	max_l <- max(lengths(overlap));
	overlapdf <- rapply(overlap, function(x) 'length<-'(x, max_l), how="list");
	write.csv(overlapdf, file=file.path("VennDiagram", paste(unlist(strsplit(hm, ".", fixed=TRUE))[1], "csv", sep=".")));
}
