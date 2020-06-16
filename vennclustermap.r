source("utils/rtools.r");

list.packages = c("VennDiagram", "ComplexHeatmap", "circlize", "RColorBrewer", "measurements", "svglite", "stringr", "jsonlite", "base64enc")
install_missing(list.packages)

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

# Disable Rplots.pdf from generating;
pdf(NULL)

jsonData <- list(VennDiagram = list(), ClusterHeatMap = list())

parsedHM <- str_match(hms, "(.+)\\.csv")[,2]

# Function to autogen heatmap

genHM.rowgp <- gpar(fontsize = 7);
genHM.colgp <-gpar(fontsize = 7);

genHM.minwidth <- unit(6, "cm");
genHM.minheight <- unit(10, "cm");


# Helper function to generate a heatmap given name
genHM <- function(loc, name, data, fnum){

	print(paste('Generating Heatmap', name))
	# use percentilism
	data <- apply(data,2,rank, na.last = "keep", ties.method = "max")

	hmheight <- max(max_text_height(
	    rownames(data), 
	    gp = genHM.colgp
	  )*hmcount*linespacing, genHM.minheight)
	rowhmheight <- max_text_height(
	    rownames(data), 
	    gp = genHM.rowgp
	  )

	hm2 <- Heatmap(
		data[,1:fnum],
		name="File Id",
		column_title="File Id",
		row_title="Peptide Rank",
		col=colors,
		row_names_gp = genHM.rowgp,
		column_names_gp = genHM.colgp,
		width = max(rowhmheight*fnum*linespacing, genHM.minwidth),
		height=hmheight
	)+Heatmap(
		data[,(fnum+1):NCOL(data)],
		name="Test Group",
		column_title="Test Group",
		col=colors,
		row_names_gp = genHM.rowgp,
		column_names_gp = genHM.colgp,
		width = max(rowhmheight*(NCOL(data)-fnum)*linespacing, genHM.minwidth),
		height=hmheight
	)
	hm2plot <- draw(hm2)
	# attr(plot, "layout")$page_size[1] <- attr(plot, "ht_list_param")$width

	hmw <- as.double(ComplexHeatmap:::width(hm2plot))*mm2in
	hmh <- as.double(ComplexHeatmap:::height(hm2plot))*mm2in
	# print(hmw)
	# print(conv_unit(303.1907, "mm", "inch"))
	outfname <- file.path("ClusterHeatMap", paste(loc, name, sep='_'));
	outsvg <- stringSVG(print(hm2plot), 
		width=hmw, 
		height=hmh)
	write(outsvg, file= paste(outfname, "svg", sep="."))

	return(list(name=name, data=as.character(outsvg)))
}

for(hmid in 1:length(hms)){
	hm <- hms[hmid]
	f <- read.csv(file.path("HeatMap", hm));
	
	f <- f[Reduce("&", list(f$Row_Type==1, f$Contaminant == 0)), ]; #take only uncontamiated classes
	# f <- f[f$contaminant==0,]
	
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


	rownames(f) <- f$Rank_Number #save ranks of each
	f.geneNames <- f[,c("Rank_Number", "Gene_Name")] #save gene names
	f <- f[,paste("x_OfSpectra", fids$ID, sep="_")]
	colnames(f) <- fids$ID;

	colnames(groups) <- dataset.groupids;

	combinedf <- data.matrix(cbind(f, groups))


	# AAAAA unfinsihed  A A A A must set hmid

	hm.functions <- list(
		"ranks"=function(i, data) i,
		"stdev"=function(i, data) -sd(data[i,], na.rm=TRUE) # negative is used because rank() puts low number with the highest rank (lowest # rank)
	)

	hm.list <- list()

	for(hm.name in names(hm.functions)){
		data <- combinedf[rank(sapply(1:NROW(combinedf), hm.functions[[hm.name]], combinedf), ties.method = "first"),]
		data <- data[1:hmcount,]

		hm.list[[length(hm.list)+1]] <- genHM(parsedHM[hmid], paste(hm.name, 'ids', sep='_'), data, NCOL(f))
		rownames(data) <- sapply(rownames(data), function(v) f.geneNames[f.geneNames$Rank_Number==as.numeric(v),"Gene_Name"])
		hm.list[[length(hm.list)+1]] <- genHM(parsedHM[hmid], paste(hm.name, 'names', sep='_'), data, NCOL(f))
	}

	jsonData[['ClusterHeatMap']][[hmid]] <- list(name=parsedHM[hmid], sheets=hm.list);
	
	outfname <- paste(unlist(strsplit(hm, ".", fixed=TRUE))[1], "png", sep=".");

	# Try to generate PWA, if fail too bad
	if (length(sets) > 5){
		message("Number of test groups exceeds 5, the maximum number for the venn diagram. As a result, no venn diagram will be generated.")
	}else{
		venn.diagram(
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

		max_l <- max(lengths(overlap));
		overlap <- as.data.frame(matrix(unlist(lapply(overlap, 'length<-', max_l)), nrow=max_l))
		colnames(overlap) <- cnames;
		# print(overlap)
		overlap <- data.frame(lapply(overlap, function(x) ifelse(is.na(x), "", x)));
		# print(overlap)

		# ifelse(v=="", "", f.geneNames[f.geneNames$Rank_Number==v,"Gene_Name"])

		# convert list of IDs to names
		overlapnames <- data.frame(apply(overlap, c(1,2), function(v) ifelse(v=="", "", as.character(f.geneNames[f.geneNames$Rank_Number==as.numeric(v),"Gene_Name"]))))


		# add venn diagram data to json
		jsonData[['VennDiagram']][[hmid]] <- list(
			Name=parsedHM[hmid],
			img = paste("<img src=\"", 
									dataURI(mime = "image/png", encoding = "base64", file = file.path("VennDiagram", outfname)), 
									"\" class=\"vennImg\" alt=\"Venn Diagram\"/>", sep=""),
			raw = list(
				Id=overlap,
				Name=overlapnames
			)
		)

		#maximum length of any of the venn diagram areas (used to equalize for saving to csv)

		# save id and name csv files
		write.csv(overlap, file=file.path("VennDiagram", paste(parsedHM[hmid], "_ids.csv", sep="")), na = "");
		write.csv(overlapnames, file=file.path("VennDiagram", paste(parsedHM[hmid], "_names.csv", sep="")), na = "");
	}
}

# save json raw file
write(toJSON(jsonData, auto_unbox=TRUE), file=file.path("Raws", "vennclustermap.json"));