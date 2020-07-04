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
if (length(dataset.groupids) <= 8) {
	palette <- brewer.pal(length(dataset.groupids), "Pastel2")
}
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

	message(paste('Generating Heatmap', name))
	# use percentilism

	outfname <- file.path("ClusterHeatMap", paste(loc, name, sep='_'));
	write.csv(data, file=paste(outfname, 'csv', sep='.'), na="")

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

	dir.create("VennDiagram", showWarnings = FALSE)
	dir.create("ClusterHeatMap", showWarnings = FALSE)


	rownames(f) <- f$Rank_Number #save ranks of each
	f.geneNames <- f[,c("Rank_Number", "Gene_Name")] #save gene names
	f <- f[,paste("x_OfSpectra", fids$ID, sep="_")]
	colnames(f) <- fids$ID;

	colnames(groups) <- dataset.groupids;

	combinedf <- data.matrix(cbind(f, groups))

	# Return TRUE if include, FALSE if not (or any other truthy/falsey value)
	hm.prefunctions <- list(
		"ranks"=function(i, data) TRUE, # just include all
		"stdev"=function(i, data) i<NROW(data)/10, # top 10% by rank
		"coeffvar"=function(i, data) i<NROW(data)/10 # same as above
	)

	# Return a value by which to sort the list (low->high)
	hm.functions <- list(
		"ranks"=function(i, data) i, # sort by rank
		"stdev"=function(i, data) -sd(data[i,], na.rm=TRUE), # negative is used because rank() puts low number with the highest rank (lowest # rank)
		"coeffvar"=function(i, data) -cv(data[i,], na.rm=TRUE) # same as above but with coeff of variation
	)

	hm.list <- list()

	for(hm.name in names(hm.functions)){
		hm.func <- function(i, data){
			if(hm.prefunctions[[hm.name]](i, data)){ # if prefunction says to include, return parsed value
				return(hm.functions[[hm.name]](i, data))
			}
			return(NA) # return NA if it shouldnt be included (will be put last)
		}
		data <- sapply(1:NROW(combinedf), hm.func, combinedf) # get sort order first
		data <- data[!is.na(data)] # remove NAs (things that we dont want to include)
		data <- combinedf[rank(data, ties.method = "first"),] # convert list of sorted to actual data
		data <- data[1:min(NROW(data), hmcount),] # take the first <hmcount> so the table isnt too big

		hm.list[[length(hm.list)+1]] <- genHM(parsedHM[hmid], paste(hm.name, 'ids', sep='_'), data, NCOL(f))
		rownames(data) <- sapply(rownames(data), function(v) f.geneNames[f.geneNames$Rank_Number==as.numeric(v),"Gene_Name"])
		hm.list[[length(hm.list)+1]] <- genHM(parsedHM[hmid], paste(hm.name, 'names', sep='_'), data, NCOL(f))
	}

	jsonData[['ClusterHeatMap']][[hmid]] <- list(name=parsedHM[hmid], sheets=hm.list);

	# for sets of >10, generating overlap data will not only take too long but also be virtually useless due to 2^n possible overlap configs
	if (length(dataset.groupids) > 10){
		message("Number of test groups exceeds 10. Overlap data and Venn Diagram images will not be calculated.")
		next
	}

	# create sets for venn diagram later usage
	sets <- sapply(1:length(dataset.groupids), function(x) which(compareNA(groups[,x]>cutoff*maxgroup)));

	# calculate overlap sections of the venn diagram
	cnames <- recursiveBinary(length(dataset.groupids));
	overlap <- replicate(length(cnames), c())
	names(overlap) <- cnames;

	# loop through each protein and assign it to the right section
	for(prot in 1:NROW(combinedf)){
		p.x <- paste(sapply(sets, function(set) as.numeric(is.element(prot, set))), collapse="")
		overlap[[p.x]][[length(overlap[[p.x]])+1]] <- prot
	}

	# convert list of vectors with different lengthsto a dataframe
	max_l <- max(lengths(overlap));
	# overlap <- as.data.frame(matrix(unlist(lapply(overlap, 'length<-', max_l)), nrow=max_l))
	overlap <- do.call(cbind, lapply(overlap, 'length<-', max_l))
	
	# Replace NA values with empty
	overlap <- data.frame(apply(overlap, c(1, 2), function(x) ifelse(is.na(x), "", x)));

	# ifelse(v=="", "", f.geneNames[f.geneNames$Rank_Number==v,"Gene_Name"])

	# convert list of IDs to names
	overlapnames <- data.frame(apply(overlap, c(1,2), function(v) ifelse(v=="", "", as.character(f.geneNames[f.geneNames$Rank_Number==as.numeric(v),"Gene_Name"]))))

	# add venn diagram data to json
	jsonData[['VennDiagram']][[hmid]] <- list(
		Name=parsedHM[hmid],
		raw = list(
			Id=overlap,
			Name=overlapnames
		)
	)

	#maximum length of any of the venn diagram areas (used to equalize for saving to csv)

	# save id and name csv files
	message("Saving Venn Diagram Raws")
	write.csv(overlap, file=file.path("VennDiagram", paste(parsedHM[hmid], "_ids.csv", sep="")), na = "");
	write.csv(overlapnames, file=file.path("VennDiagram", paste(parsedHM[hmid], "_names.csv", sep="")), na = "");

	# too many datasets to generate venn diagram
	if (length(sets) > 5){
		message("Number of test groups exceeds 5, the maximum number for the venn diagram. As a result, no venn diagram image will be generated.")
		next
	}

	outfname <- paste(unlist(strsplit(hm, ".", fixed=TRUE))[1], "png", sep=".");

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
	message("Saving Venn Diagram Image")
	jsonData[['VennDiagram']][[hmid]][['img']] <- paste("<img src=\"", 
							dataURI(mime = "image/png", encoding = "base64", file = file.path("VennDiagram", outfname)), 
							"\" class=\"vennImg\" alt=\"Venn Diagram\"/>", sep="")
}

# save json raw file
message("Saving Final JSON Raw")
write(toJSON(jsonData, auto_unbox=TRUE), file=file.path("Raws", "vennclustermap.json"));
