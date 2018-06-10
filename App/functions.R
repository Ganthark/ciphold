##This file contains all the function used in the shiny application, grouped in sections.
##Many functions remain untouched from the original code found on this link:
##https://github.com/nolanlab/scaffold/tree/multiFilesClustering
##They are not commented.
##The new functions are on the top of their respective sections.


##GLOBAL FUNCTIONS

chooseDir <-  function() {
  OS <- Sys.info()["sysname"]
    if (OS=="Windows") {
    Dir <-
      choose.dir(default = "", caption = "Select a Folder for saving:")
    }
    else if (OS=="Linux") {
      Dir <- tk_choose.dir(default = "", caption = "Select a Folder for saving:")
    }
    else {
      Dir <- choose.mac.dir()
    }
  return(Dir)
}
#Function used to select a folder via interface on a mac OS system.
choose.mac.dir <- function() {
	system("osascript -e 'tell app \"R\" to POSIX path of (choose folder with prompt \"Select a Folder for saving:\")' > /tmp/R_folder",
			intern = FALSE, ignore.stderr = TRUE)
	p <- system("cat /tmp/R_folder && rm -f /tmp/R_folder", intern = TRUE)
	return(ifelse(length(p), p, NA))
}

##CLUSTERING
options(stringsAsFactors = F)

#Function used to organize data. Groups are processed depending whether they contain multiple flow frames or not.
get_flow_frames_from_groups <- function(flow.frames, files.list)
{
  v <- list()
  v <- lapply(files.list, function(x) {
    if (length(x) > 1) {
      w <- list()
      w <- lapply(x, function(y) {
        j <- flow.frames[names(flow.frames) == y]
        names(j) <- y
        return(j)
      })
      names(w) <- x
      return (w)
    }
    else {
      u <- flow.frames[names(flow.frames) == x]
      names(u) <- x
      return(u)
    }
  })
  return(v)
}

#This function is used to process the flow frames for each group, data are compensed, they are transformed according to user's choice and then clustered.
process_files_groups <- function(flow.frames, col.names, num_clusters, num_samples, cofactor, downsample.to, output.dir, transform.data, output_type ="legacy")
{

	cluster_data <- function(tab, col.names, k, algorithm = "", transform.data, ...)
	{

		m <- as.matrix(tab[, col.names])
		if(algorithm == "clara")
		{
			print("Performing clara clustering")
			groups <- clara(m, k, ...)$clustering
		}

		else if(algorithm == "hierarchical")
		{
			print("Performing hierarchical clustering")
			dend <- hclust(dist(m), ...)
			groups <- cutree(dend, k)
		}

		print("Clustering done")
		print(Sys.time())
		tab <- cbind(tab, groups)
		return(tab)
	}
	tab <- NULL
	orig.data <- NULL
	temp <- NULL
	if (length(names(flow.frames)) > 1) {
		f <- names(flow.frames[[1]])
	}
	else {
		f <- names(flow.frames)
	}
	for(i in c(1:length(flow.frames)))
	{
		if (class(flow.frames[[1]]) == "flowFrame")
			elt <- flow.frames[[1]]
		else
			elt <- flow.frames[[i]][[1]]
		temp.orig.data <- exprs(elt)
		temp.tab <- convert_fcs(elt, cofactor, transform.data = transform.data)
		colnames(temp.tab) <- pData(parameters(elt))$desc

		if(any(is.na(colnames(temp.tab))))
		{
			w <- is.na(colnames(temp.tab))
			colnames(temp.tab)[w] <- pData(parameters(elt))$name[w]
		}

		temp.tab <- as.matrix(temp.tab)
		temp.tab[temp.tab < 0] <- 0

		if(downsample.to > 0)
		{
			x <- NULL
			if(nrow(temp.tab) <= downsample.to)
			{
				print("Number of events smaller than downsampling target, taking all events")
				x <- 1:nrow(temp.tab)
			}
			else
			{
				print(sprintf("Predownsampling to %d events", downsample.to))
				x <- sample(1:nrow(temp.tab), size = downsample.to)
			}
			temp.tab <- temp.tab[x,]
			temp.orig.data <- temp.orig.data[x,]
		}

		temp.tab <- as.data.frame(temp.tab, check.names = F, stringsAsFactors = F)

		temp.tab <- data.frame(temp.tab, sample = f, check.names = F, stringsAsFactors = F)
		temp.orig.data <- data.frame(temp.orig.data, sample = f, check.names = F, stringsAsFactors = F)
		tab <- rbind(tab, temp.tab)
		orig.data <- rbind(orig.data, temp.orig.data)
	}

	m <- cluster_data(tab, col.names, k = num_clusters, algorithm = "clara", sampsize = min(nrow(tab), 1000), samples = num_samples)
	colnames(m) <- gsub("groups", "cellType", colnames(m))
	orig.data <- cbind(orig.data, cellType = m[, "cellType"])

	temp <- get_stats_by_sample(m)
	m <- data.frame(m, check.names = F)
	orig.data <- data.frame(orig.data, stringsAsFactors = FALSE, check.names = FALSE)
	write_clustering_output(base.name = f, tab.medians = temp, clustered.data =  m, output.dir =  output.dir)
}

#Function that actually compensates data, then transforms them.
convert_fcs <- function(f, cofactor, transform.data = "")
{
  comp <- grep("SPILL", names(description(f)), value = T)
  if(length(comp) > 0)
  {
    print("Found compensation matrix, applying...")
    comp <- description(f)[comp][[1]]
    if(is.character(comp))
    {
      comp <- strsplit(comp, ",")[[1]]
      num.channels <- as.numeric(comp[1])
      m <- matrix(nrow = num.channels, byrow = T, data = as.numeric(comp[(num.channels + 2):length(comp)]))
      colnames(m) <- comp[2:(1 + num.channels)]
      comp <- m
    }
    f <- compensate(f, spillover = comp)
  }

  #Any transformation method can be added here, but you must also modify server.R and ui.R accordingly.
  if(transform.data == "A"){
    tab <- exprs(f)
    m <- as.matrix(tab)
    m <- asinh(m / cofactor)
  }
  else if (transform.data == "L") {
    fTrans <- logiclTransformCiphe(f) # transform value
    tab <- exprs(fTrans)
    m <- as.matrix(tab)
  }
  else {
    tab <- exprs(f)
    m <- as.matrix(tab)
  }
  col.names <- colnames(m)
  tab <- data.frame(m)
  names(tab) <- col.names
  return(tab)
}

#Function used to transform data via Logicle transform, determining automaticly the better factor.
logiclTransformCiphe <- function(flow.frame)
{
  # no.transform <- c("FSC-A","FSC-H","FSC-W","SSC-A","SSC-H","SSC-W","Time","Flag")
  # markers.transform <- colnames(flow.frame)[colnames(flow.frame)%in%no.transform == FALSE]
  if(is.null(flow.frame@description[["SPILL"]])){
    markers.transform <- colnames(flow.frame)
  } else {
    markers.transform <- colnames(flow.frame@description[["SPILL"]])
  }


  list.index <- names(unlist(lapply(markers.transform, function(x) return(which(flow.frame@description==x)))))
  list.index <- gsub("N","", list.index)
  list.index <- gsub("\\$P","", list.index)

  if(!is.null(flow.frame@description[[paste0("$P",list.index[1],"MS")]]))
  {
    r.values <- unlist(lapply(list.index, function(x)
      as.integer(flow.frame@description[[paste0("$P",x,"MS")]]))
    )
  } else if(!is.null(flow.frame@description[[paste0("P",list.index[1],"MS")]]))
  {
    r.values <- unlist(lapply(list.index, function(x)
      as.integer(flow.frame@description[[paste0("P",x,"MS")]]))
    )
  } else
  {
    r.values <- rep(90, length(list.index))
  }

  w.values <- (4.5-log10(262143/abs(r.values)))/2
  w.values[which(w.values<0)] <- 0.5
  w.values[which(is.infinite(w.values))] <- 0.5

  for(t in 1:length(markers.transform)){
    lgcl <- logicleTransform(w=w.values[t])
    flow.frame <- transform(flow.frame, transformList(markers.transform[t],lgcl))
  }

  return(flow.frame)
}


cluster_fcs_files_groups <- function(flow.frames, files.list,  col.names, num_clusters, num_samples,
                                     cofactor, downsample.to, output_dir = NULL, transform.data)
{

  lapply(flow.frames,
         process_files_groups, col.names = col.names, num_clusters = num_clusters, num_samples = num_samples,
         cofactor = cofactor, downsample.to = downsample.to, output_type = output_type, output.dir = output_dir, transform.data = transform.data)

  return(files.list)
}


get_cluster_groups_table <- function(v, key) {
  tags$table(class = "table table-hover table-striped",
             tags$tr(tags$td(
               v[1],
               tags$button(class = "btn btn-xs btn-warning pull-right", onClick = sprintf("Shiny.onInputChange('clusteringui_remove_clustering_group', {'key':'%s', 'x':Math.random()})", key),
                           tags$span(class = "glyphicon glyphicon-trash")
               )
             )),
             ifelse(length(v > 1),
                    tagList(lapply(tail(v, n = -1), function(x) {tags$tr(tags$td(x))})),
                    tagList()
             )
  )
}

write_clustering_output <- function(base.name, tab.medians, clustered.data, output.dir, output.type = "legacy")
{
	if(output.type == "legacy")
	{
	  OS <- Sys.info()["sysname"]
    if (OS=="Windows") {
  		write.table(tab.medians, file = paste(output.dir, "\\", base.name, ".clustered.txt", sep = ""), row.names = F, sep = "\t", quote = F)
  		my_save(clustered.data, paste(output.dir, "\\", base.name, ".clustered.all_events.RData", sep = ""))
    }
    else if (OS=="Linux") {
      write.table(tab.medians, file = paste(output.dir, "/", base.name, ".clustered.txt", sep = ""), row.names = F, sep = "\t", quote = F)
	  	my_save(clustered.data, paste(output.dir, "/", base.name, ".clustered.all_events.RData", sep = ""))
    }
    else {
      write.table(tab.medians, file = paste(output.dir, "/", base.name, ".clustered.txt", sep = ""), row.names = F, sep = "\t", quote = F)
	  	my_save(clustered.data, paste(output.dir, "/", base.name, ".clustered.all_events.RData", sep = ""))
    }

	}
	else if(output.type == "directory")
	{
		clustered.data.dir <- "clustered.data"
		txt.file.name <- paste(base.name, ".clustered.txt", sep = "")
		full.path <- file.path(output.dir, clustered.data.dir, txt.file.name)
		dir.create(full.path, recursive = T)

		write.table(tab.medians, file.path(output.dir, txt.file.name, sep = ""),
					row.names = F, sep = "\t", quote = F)
		ddply(clustered.data, ~cellType, function(x) {
			saveRDS(x, file = file.path(full.path, sprintf("cluster_%d.RData", x$cellType[1])))
		})
	}
}

get_stats_by_sample <- function(tab)
{
	tab.medians <- ddply(tab, ~cellType, colwise(median, is.numeric))
	tab.medians.by.sample <- ddply(tab, ~cellType * sample, colwise(median, is.numeric))
	pop.size <- ddply(tab, ~cellType, nrow)
	names(pop.size) <- gsub("V1", "popsize", names(pop.size))
	pop.size.by.sample <- ddply(tab, ~cellType * sample, nrow)
	names(pop.size.by.sample) <- gsub("V1", "popsize", names(pop.size.by.sample))
	tab.medians <- merge(tab.medians, pop.size, by = "cellType")
	tab.medians.by.sample <- merge(tab.medians.by.sample, pop.size.by.sample, by = c("cellType", "sample"), all.x = T)


	#Rotate the by.sample table
	temp <- melt(tab.medians.by.sample, id = c("cellType", "sample"))
	temp$variable <- paste(temp$variable, temp$sample, sep = "@")
	temp$sample <- NULL
	temp <- cast(temp, cellType~variable)


	ret <- merge(tab.medians, temp, by = "cellType", all.x = T)

	return(ret)

}

my_save <- function(obj, f_name)
{
	con <- file(f_name, "wb")
	serialize(obj, con, ascii = F)
	close(con)
}

get_fcs_col_names <- function(f)
{
	fcs.file <- f
	ret <- as.vector(pData(parameters(fcs.file))$desc)

	if(any(is.na(ret)))
	{
		w <- is.na(ret)
		ret[w] <- as.vector(pData(parameters(fcs.file))$name[w])
	}

	return(ret)
}



##ANALYSIS

#Runs the gated analysis, mainly like the original code but with a few tricks due to how data are processed in this version.
run_analysis_gated <- function(flow.frames, clusteredFiles, all_events, outputDir, col.names, cofactor, transform.data, ...)
{
	gated.flow.frames <- lapply(c(1:length(flow.frames)), function(i) {flow.frames[[i]]@exprs})
	print(sprintf("Markers used for SCAFFoLD: %s", paste(col.names, collapse = ", ")))
	print(paste("Using as reference", names(flow.frames)[[1]], sep = " "))
	gated_data <- load_attractors_from_gated_data(flow.frames, all_events, cofactor, transform.data, col.names)
	tab.attractors <- gated_data$tab.attractors
	att.labels <- gated_data$cellType_key$population
	G.attractors <- NULL
	ret <- process_files(clusteredFiles,  G.attractors, tab.attractors, att.labels, col.names, scaffold.mode = "gated", all_events, ...)
	ret <- c(list(scaffold.col.names = col.names, landmarks.data = gated_data$downsampled.data), ret)
	print(sprintf("%s.scaffold", names(flow.frames)[[1]], sep = "/"))
	my_save(ret, paste(outputDir, sprintf("%s.scaffold", names(clusteredFiles)[[1]]), sep = "/"))
	return(names(clusteredFiles))
}

#Runs the existing analysis, mainly like the original code but with a few tricks due to how data are processed in this version.
run_analysis_existing <- function(.scaffold, refClusteredFiles, clusteredTables, rdata, outputDir, col.names, names.mapping = NULL, ...)
{
	files.list <- clusteredTables
	print(sprintf("Markers used for SCAFFoLD: %s", paste(col.names, collapse = ", ")))
	print(paste("Using as reference", refClusteredFiles[1], sep = " "))
	ref.scaffold.file <- refClusteredFiles[2]
	ref.scaffold.data <- .scaffold
	ref.scaffold.markers <- ref.scaffold.data$scaffold.col.names
	View(ref.scaffold.markers)
	l <- load_existing_layout(ref.scaffold.data)
	tab.attractors <- l$tab.attractors
	G.attractors <- l$G.attractors
	att.labels <- V(G.attractors)$Label
	ret <- process_files(files.list, G.attractors, tab.attractors, att.labels, col.names,
						 scaffold.mode = "existing", ref.scaffold.markers = ref.scaffold.markers, names.mapping = names.mapping, all_events = rdata, ...)
	ret <- c(list(scaffold.col.names = col.names, landmarks.data = ref.scaffold.data$landmarks.data), ret)
	View(ret)
	my_save(ret, paste(outputDir, sprintf("%s.scaffold", basename(refClusteredFiles[1])), sep = "/"))
	return(files.list)
}

#Function used to process files, only used for the existing analysis from the "map dataset" tab.
process_files <- function(clusteredFiles, G.attractors, tab.attractors, att.labels, col.names, scaffold.mode, all_events = NULL,
                          ref.scaffold.markers = NULL, names.mapping = NULL, ew_influence = NULL,
                          col.names.inter_cluster = NULL, ...)
{
  ret <- list(graphs = list(), clustered.data = list())
  map_names <- names_map_factory(names.mapping)
  for(i in c(1:length(clusteredFiles)))
  {
    print(paste("Processing", names(clusteredFiles)[[i]], sep = " "))
    tab <- clusteredFiles[[i]]
    names(tab) <- map_names(names(tab))
    col.names.inter_cluster <- map_names(col.names.inter_cluster)

    if(scaffold.mode == "existing")
    {
      #Some markers in the reference scaffold file have been designated
      #for mapping, but they are missing from the sample files
      if(any(is.na(names(names.mapping))))
        tab <- add_missing_columns(tab, col.names, fill.data = 0)
      if(is.null(ew_influence))
        ew_influence <- ceiling(sum(!is.na(names(names.mapping))) / 3)
    }
    else
    {
      if(is.null(ew_influence))
        ew_influence <- ceiling(length(col.names) / 3)
    }

    tab <- tab[!apply(tab[, col.names], 1, function(x) {all(x == 0)}),]
    names(tab) <- gsub("cellType", "groups", names(tab))
    names(tab) <- gsub("^X", "", names(tab))
    print(sprintf("Running with Edge weight: %f", ew_influence))
    res <- process_data(tab, G.attractors, tab.attractors,
                        col.names = col.names, att.labels = att.labels, already.clustered = T, ew_influence = ew_influence,
                        col.names.inter_cluster = col.names.inter_cluster, ...)
    G.complete <- get_highest_scoring_edges(res$G.complete)
    clustered.data <- all_events[[i]]
    names(clustered.data) <- map_names(names(clustered.data))
    clustered.data <- downsample_by(clustered.data, "cellType", 1000)

    ret$graphs[names(clusteredFiles)[[i]]] <- list(G.complete)
    ret$clustered.data[names(clusteredFiles)[[i]]] <- list(clustered.data)

    G.attractors <- res$G.attractors
  }
  ret <- c(ret, list(dataset.statistics = get_dataset_statistics(ret)))
  return(ret)
}

#The code for the unsupervised analysis is here as is, in case of further implementation of this functionnality.
# run_analysis_unsupervised <- function(clusteredFiles, clusteredFilesNames, outputDir, ref.file, col.names, ...)
# {
# 	files.list <- list.files(path = working.dir, pattern = "*.clustered.txt$")
# 	files.list <- files.list[files.list != ref.file]
# 	print(sprintf("Markers used for SCAFFoLD: %s", paste(col.names, collapse = ", ")))
# 	files.list <- c(ref.file, files.list)
# 	print(paste("Using as reference", files.list[1], sep = " "))
# 	files.list <- paste(working.dir, files.list, sep = "/")
#
# 	temp <- get_attractors_from_graph_clustering(files.list[[1]], col.names)
# 	tab.attractors <- temp$tab.attractors
# 	att.labels <- temp$att.labels
# 	G.attractors <- NULL
# 	ret <- process_files(files.list, G.attractors, tab.attractors, att.labels, col.names, scaffold.mode = "unsupervised", ...)
# 	ret <- c(list(scaffold.col.names = col.names), ret)
# 	my_save(ret, paste(working.dir, sprintf("%s.scaffold", ref.file), sep = "/"))
# 	return(files.list)
# }

get_file_couples <- function(v, key, tag) {
  tags$table( class = "table table-hover table-striped",
              tags$tr(tags$td(
                v[1],
                tags$button(class = "btn btn-xs btn-warning pull-right", onClick = sprintf("Shiny.onInputChange('%s', {'key':'%s', 'x':Math.random()})", tag, key),
                            tags$span(class = "glyphicon glyphicon-trash")
                )
              )),
              ifelse(length(v > 1),
                     tagList(lapply(tail(v, n = -1), function(x) {tags$tr(tags$td(x))})),
                     tagList()
              )
  )
}

load_attractors_from_gated_data <- function(flow.frames, all_events, cofactor, transform.data = "", col.names, ...)
{
	res <- NULL
	for(i in c(1:length(flow.frames)))
	{

		population <- tail(strsplit(names(flow.frames)[i], "_")[[1]], n = 1)

		fcs <- flow.frames[[i]]

		tab <- convert_fcs(fcs, cofactor, transform.data)
		if(!all(pData(parameters(fcs))$desc == " "))
			colnames(tab) <- pData(parameters(fcs))$desc
		else
			colnames(tab) <- pData(parameters(fcs))$name

		if(any(is.na(colnames(tab))))
		{
			w <- is.na(colnames(tab))
			colnames(tab)[w] <- pData(parameters(fcs))$name[w]
		}

		tab <- as.matrix(tab)
		tab[tab < 0] <- 0
		tab <- as.data.frame(tab)

		tab <- cbind(tab, population, stringsAsFactors = F)
		res <- rbind(res, tab)
	}

	downsampled.data <- downsample_by(res, "population", 1000)
	names(downsampled.data) <- gsub("population", "cellType", names(downsampled.data))

	#Change cellType to be numbers
	k <- unique(res$population)
	k <- data.frame(population = k, cellType = seq_along(k), stringsAsFactors = F)
	res <- merge(res, k)
	res <- res[, grep("population", names(res), invert = T)]
	res <- ddply(res, ~cellType, colwise(median))
	return(list(downsampled.data = downsampled.data, tab.attractors = res, cellType_key = k))
}



downsample_by <- function(tab, col.name, size)
{
	print(sprintf("Downsampling to %d events", size))
	return(ddply(tab, col.name, function(x, size)
	{
		if(nrow(x) <= size){
			return(x)
		}
		else {
			return(x[sample(1:nrow(x), size),])
		}
	}, size = size))
}

names_map_factory <- function(names.map)
{
	function(v)
	{
		sel <- v %in% names(names.map)
		if(any(sel))
			v[sel] <- names.map[v[sel]]
		return(v)

	}
}

process_data <- function(tab, G.attractors = NULL, tab.attractors = NULL, col.names = NULL, att.labels = NULL, dist.thresh = 0.7,
						 already.clustered = FALSE, inter.cluster.connections = FALSE, col.names.inter_cluster = NULL, inter_cluster.weight_factor = 0.7, ew_influence,
						 overlap_method = NULL)
{

	if(!already.clustered)
	{
		tab <- cluster_data(tab, col.names)
		tab.clustered <- ddply(tab, ~groups, colwise(median))
	}
	else
		tab.clustered <- tab

	if(is.null(col.names.inter_cluster) || col.names.inter_cluster == "")
		col.names.inter_cluster = col.names
	if(is.null(G.attractors))
	{
		G.attractors <- build_graph(tab.attractors, col.names)

		G.complete <- add_vertices_to_attractors_graph(G.attractors, tab.clustered, tab.attractors, col.names, dist.thresh)
		G.complete <- complete.forceatlas2(G.complete, first.iter = 50000,
										   overlap.iter = 20000, ew_influence = ew_influence, overlap_method = "repel")
		if(inter.cluster.connections)
		{
			print("Adding inter-cluster connections with markers:")
			print(col.names.inter_cluster)
			print(sprintf("Weight factor:%f", inter_cluster.weight_factor))
			G.complete <- add_inter_clusters_connections(G.complete, col.names.inter_cluster, weight.factor = inter_cluster.weight_factor)
			G.complete <- complete.forceatlas2(G.complete, first.iter = 50000, overlap.iter = 20000,
											   ew_influence = ew_influence, overlap_method = overlap_method)
		}
		V(G.attractors)$x <- V(G.complete)$x[1:vcount(G.attractors)]
		V(G.attractors)$y <- V(G.complete)$y[1:vcount(G.attractors)]
	}
	else
	{
		G.complete <- add_vertices_to_attractors_graph(G.attractors, tab.clustered, tab.attractors, col.names, dist.thresh)

		fixed <- rep(FALSE, vcount(G.complete))
		fixed[1:vcount(G.attractors)] <- TRUE

		G.complete <- complete.forceatlas2(G.complete, first.iter = 50000, overlap.iter = 20000,
										   overlap_method = "repel", ew_influence = ew_influence, fixed = fixed)
		if(inter.cluster.connections)
		{
			print("Adding inter-cluster connections with markers:")
			print(col.names.inter_cluster)
			print(sprintf("Weight factor:%f", inter_cluster.weight_factor))
			G.complete <- add_inter_clusters_connections(G.complete, col.names.inter_cluster, weight.factor = inter_cluster.weight_factor)
			G.complete <- complete.forceatlas2(G.complete, first.iter = 50000, overlap.iter = 20000,
											   overlap_method = overlap_method, ew_influence = ew_influence, fixed = fixed)
		}

	}
	G.complete <- add_attractors_labels(G.complete, att.labels)
	V(G.complete)$name <- gsub(".fcs", "", V(G.complete)$name)
	return(list(G.attractors = G.attractors, G.complete = G.complete, tab.attractors = tab.attractors, tab = tab, col.names = col.names))
}

build_graph <- function(tab, col.names, filtering_T = 0.8)
{
	m <- as.matrix(tab[, col.names])
	row.names(m) <- tab$cellType
	dd <- cosine_similarity_matrix(m)
	diag(dd) <- 0
	dd[is.na(dd)] <- 0 #This can happen if one of the attractors has all 0's for the markers of interest

	if(filtering_T >= 1)
		dd <- filter_similarity_matrix_by_rank(dd, filtering_T)
	else
		dd <- filter_similarity_matrix(dd, filtering_T)
	G <- graph.adjacency(dd, mode = "undirected", weighted = T)
	n.vertices <- length(V(G))
	lay <- layout.kamada.kawai(G)
	colnames(lay) <- c("x", "y")
	G <- set.vertex.attribute(G, name = "x", value = lay[, "x"])
	G <- set.vertex.attribute(G, name = "y", value = lay[, "y"])
	for(i in names(tab))
		G <- set.vertex.attribute(G, name = i, value = tab[, i])

	return(G)
}

cosine_similarity_matrix <- function(m)
{
	ret <- t(apply(m, 1, function(x, m) {cosine_similarity_from_matrix(x, m)}, m = m))
	return(ret)
}

cosine_similarity_from_matrix <- function(v, m)
{
	m <- as.matrix(m[, names(v), drop = F])
	ret <- apply(m, 1, function(x, v) {return(crossprod(x, v)/sqrt(crossprod(x) * crossprod(v)))}, v)
	return(ret)
}

filter_similarity_matrix <- function(m, T)
{
	ret <- t(apply(m, 1, function(x)
	{
		if(max(x) <= T)
			x[x < max(x)] <- 0
		else
			x[x < T] <- 0
		return(x)
	}))
	return(ret)
}

filter_similarity_matrix_by_rank <- function(m, T)
{
	ret <- t(apply(m, 1, function(x)
	{
		r <- rank(x, ties.method = "first")
		r <- max(r) - r + 1
		x[r > T] <- 0
		return(x)
	}))
	return(ret)
}

add_vertices_to_attractors_graph <- function(G, tab.clustered, tab.median, col.names, dist.thresh = 0.7)
{
	dd <- get_distances_from_attractors(tab.clustered, tab.median, col.names, dist.thresh)
	n <- nrow(dd)
	num.vertices <- length(V(G))
	G <- add.vertices(G, n)
	v.seq <- (num.vertices + 1):length(V(G))
	V(G)[v.seq]$name <- as.character(v.seq)
	V(G)[v.seq]$Label <- paste("c", tab.clustered$groups, sep = "")
	row.names(dd) <- as.character(v.seq)
	for(i in 1:nrow(dd))
	{
		v <- dd[i,]
		v <- v[v > 0]
		if(length(v) > 0)
		{
			e.list <- c(rbind(as.character(num.vertices + i), names(v)))
			G <- G + edges(e.list, weight = v)
		}
	}

	#	weight <- E(G)$weight
	#	E(G)$weight <- weight ^ 10
	maxx <- maxy <- rep(Inf, vcount(G))
	minx <- miny <- rep(-Inf, vcount(G))

	maxx[1:num.vertices] <- minx[1:num.vertices] <- V(G)$x[1:num.vertices]
	maxy[1:num.vertices] <- miny[1:num.vertices] <- V(G)$y[1:num.vertices]
	lay <- layout.kamada.kawai(G, minx = minx, maxx = maxx, miny = miny, maxy = maxy)
	colnames(lay) <- c("x", "y")
	G <- set.vertex.attribute(G, name = "x", value = lay[, "x"])
	G <- set.vertex.attribute(G, name = "y", value = lay[, "y"])

	V(G)[1:num.vertices]$type <- 1 #attractor
	V(G)[(num.vertices + 1):vcount(G)]$type <- 2 #cell

	for(i in names(tab.clustered))
		G <- set.vertex.attribute(G, name = i, index = (num.vertices + 1):vcount(G), value = tab.clustered[, i])

	G <- set_visual_attributes(G)
	return(G)
}

distance_from_attractor_hard_filter <- function(dd, tab, col.names, thresh = 0.5)
{
	tab <- tab[, col.names]
	w <- apply(tab[,col.names], 1, function(x, thresh) {all(x < thresh)}, thresh = thresh)
	if(any(w))
		print("Hard removing some connections to unstained landmarks")
	dd[, w] <- 0
	return(dd)
}


get_distances_from_attractors <- function(m, tab, col.names, dist.thresh)
{
	att <- as.matrix(tab[, col.names])
	row.names(att) <- as.character(1:nrow(tab))
	m <- as.matrix(m[, col.names])
	dd <- t(apply(m, 1, function(x, att) {cosine_similarity_from_matrix(x, att)}, att))
	dd <- distance_from_attractor_hard_filter(dd, tab, col.names, thresh = 1)

	dist.thresh <- quantile(dd, probs = 0.85, na.rm = T)
	dist.thresh <- max(c(dist.thresh, 0.5))



	dd[is.na(dd)] <- 0 #This can happen if one of the attractors has all 0's for the markers of interest
	dd <- filter_similarity_matrix(dd, dist.thresh)
	return(dd)
}

set_visual_attributes <- function(G)
{
	att <- V(G)$type == 1
	V(G)$r <- 79
	V(G)$g <- 147
	V(G)$b <- 222
	V(G)$size <- 10

	V(G)[att]$r <- 255
	V(G)[att]$g <- 117
	V(G)[att]$b <- 128
	V(G)[att]$size <- 20

	E(G)$r <- 180
	E(G)$g <- 180
	E(G)$b <- 180

	return(G)
}

layout.forceatlas2 <- function(G, ew_influence = 1, kgrav = 1, iter = 1000, prevent.overlap = FALSE, fixed = rep(FALSE, vcount(G)), stopping_tolerance = 0.001, barnes_hut = FALSE)
{
	if(vcount(G) >= 2000)
		barnes_hut <- TRUE
	if(vcount(G) > 2000)
		stopping_tolerance <- 0.01
	else if(vcount(G) > 800)
		stopping_tolerance <- 0.005
	else
		stopping_tolerance <- 0.001

	if(is.null(get.vertex.attribute(G, "x")))
	{
		lay <- cbind(x = rnorm(vcount(G)), y = rnorm(vcount(G)))
	}
	else
	{
		lay <- cbind(x = V(G)$x, y = V(G)$y)
	}


	#This is only used with prevent.overlap
	if(is.null(get.vertex.attribute(G, "size")))
		V(G)$size <- rep(10, vcount(G))
	mass <- 1 + degree(G)
	F_att <- (E(G)$weight ^ ew_influence)
	edge_list <- get.edgelist(G, names = F) - 1 #This is gonna be used in the C code where the indexing is 0-based

	avg_displ <- numeric(iter)
	max_displ <- numeric(iter)

	print(system.time(layout_forceatlas2Cpp(lay, F_att, mass, V(G)$size, edge_list, avg_displ,
											kgrav,  iter, prevent.overlap, fixed, max_displ, stopping_tolerance, barnes_hut)))

	return(list(lay = lay, avg_displ = avg_displ, max_displ = max_displ))
}

complete.forceatlas2 <- function(G, first.iter = 1000, overlap.iter, overlap_method = NULL, ...)
{

	print("First iteration")
	ret <- layout.forceatlas2(G, prevent.overlap = FALSE, iter = first.iter, ...)
	lay <- ret$lay
	#plot(ret$avg_displ, type = "l")
	#lines(ret$max_displ, col = "red")
	G <- set.vertex.attribute(G, name = "x", value = lay[, 1])
	G <- set.vertex.attribute(G, name = "y", value = lay[, 2])
	if(!is.null(overlap_method))
	{
		if(overlap_method == "repel")
		{
			print("Second iteration with prevent overalp")
			ret <- layout.forceatlas2(G, prevent.overlap = TRUE, iter = overlap.iter, ...)
			lay <- ret$lay
			if(any(is.na(lay)))
			{
				print("Prevent overlap iteration failed")
			}
			#plot(ret$avg_displ, type = "l")
			#lines(ret$max_displ, col = "red")
			else
			{
				G <- set.vertex.attribute(G, name = "x", value = lay[, 1])
				G <- set.vertex.attribute(G, name = "y", value = lay[, 2])
			}
		}
		else if(overlap_method == "expand")
			G <- adaptive_expand(G, overlap.iter)
	}
	return(G)
}

# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

layout_forceatlas2Cpp <- function(lay, F_att_orig, mass, nodes_size, edge_list, avg_displ, kgrav, iter, prevent_overlap, fixed, max_displ, stopping_tolerance, barnes_hut) {
	layout_forceatlas2Cpp(lay, F_att_orig, mass, nodes_size, edge_list, avg_displ, kgrav, iter, prevent_overlap, fixed, max_displ, stopping_tolerance, barnes_hut)
}

add_inter_clusters_connections <- function(G, col.names, weight.factor)
{
	tab <- get_vertex_table(G)
	tab <- tab[tab$type == 2,]
	m <- as.matrix(tab[, col.names])
	row.names(m) <- tab$name
	dd <- cosine_similarity_matrix(m)
	diag(dd) <- 0
	dd[is.na(dd)] <- 0

	dist.thresh <- quantile(dd, probs = 0.85, na.rm = T)
	dist.thresh <- max(c(dist.thresh, 0.5))
	dd <- filter_similarity_matrix(dd, dist.thresh)
	dd <- filter_similarity_matrix_by_rank(dd, 3)

	e.list <- NULL

	for(i in 1:nrow(dd))
	{
		v <- dd[i,]
		v <- v[v > 0]
		if(length(v) > 0)
		{
			e.list <- rbind(e.list, data.frame(a = tab[i, "name"], b = names(v), weight = v, stringsAsFactors = FALSE))

		}
	}
	temp <- as.matrix(e.list[, c("a", "b")])
	temp <- t(apply(temp, 1, sort))
	e.list <- data.frame(temp, weight = e.list$weight, stringsAsFactors = FALSE)
	names(e.list)[1:2] <- c("a", "b")
	e.list <- e.list[!duplicated(e.list[, c("a", "b")]),]
	e.list.igraph <- c(t(as.matrix(e.list[, c("a", "b")])))

	#G <- G + edges(e.list, weight = (v ^ 30))
	G <- G + edges(e.list.igraph, weight = e.list$weight * weight.factor)
	return(G)
}

adaptive_expand <- function(G, max.iter)
{
    print("Starting adaptive expansion")
    x <- V(G)$x
    y <- V(G)$y
    m <- cbind(x, y)
    ss <- outer(V(G)$size, V(G)$size, "+")

    for(i in 1:max.iter)
    {
        dd <- as.matrix(dist(m), method = "euclidean")
        dd <- dd - ss
        dd <- dd[upper.tri(dd)]
        if(all(dd >= 0))
            break
        else
            m <- m * 1.2
    }

    print(sprintf("Expansion stopped at iteration: %d", i))
    V(G)$x <- m[, "x"]
    V(G)$y <- m[, "y"]

    return(G)
}

get_vertex_table <- function(G)
{
	att <- list.vertex.attributes(G)
	ret <- NULL

	for(a in att)
	{
		d <- data.frame(get.vertex.attribute(G, a), stringsAsFactors = FALSE)
		if(is.null(ret))
			ret <- d
		else
			ret <- cbind(ret, d, stringsAsFactors = FALSE)
	}
	names(ret) <- att
	return(ret)
}

add_attractors_labels <- function(G, v)
{
	V(G)$name[1:length(v)] <- V(G)$Label[1:length(v)] <- v
	return(G)
}

get_highest_scoring_edges <- function(G)
{
	#Remove inter-cluster edges for this calculation
	e <- get.edges(G, E(G))
	E(G)$edge_type <- "cluster_to_landmark"
	e <- cbind(V(G)$type[e[,1]], V(G)$type[e[,2]])
	to.remove <- (e[,1] == 2) & (e[,2] == 2)
	E(G)$edge_type[(e[,1] == 2) & (e[,2] == 2)] <- "inter_cluster"
	g.temp <- delete.edges(G, E(G)[to.remove])

	V(g.temp)$highest_scoring_edge <- 0
	for(i in 1:vcount(g.temp))
	{
		if(V(g.temp)$type[i] == 2)
		{
			sel.edges <- incident(g.temp, i)
			max.edge <- sel.edges[which.max(E(G)[sel.edges]$weight)]
			V(g.temp)$highest_scoring_edge[i] <- max.edge
			E(G)$edge_type[max.edge] <- "highest_scoring"
		}
	}
	V(G)$highest_scoring_edge <- V(g.temp)$highest_scoring_edge
	return(G)
}

my_load <- function(f_name)
{
	con <- file(f_name, "rb")
	retval <- unserialize(con)
	close(con)
	return(retval)
}

get_dataset_statistics <- function(dataset)
{
	graphs <- dataset$graphs
	ret <- NULL
	for(G in graphs)
	{
		V(G)$popsize.relative <- V(G)$popsize / sum(V(G)$popsize, na.rm = T)
		tab <- get.data.frame(G, what = "vertices")
		max.vals <- sapply(tab, function(x) {if(is.numeric(x)) return(max(x, na.rm  =T))})
		max.vals <- max.vals[!sapply(max.vals, is.null)]
		for(i in 1:length(max.vals))
		{
			var.name <- names(max.vals)[i]
			if(!(var.name %in% names(ret)) || max.vals[[i]] > ret[[var.name]])
				ret[var.name] <- max.vals[i]
		}
	}
	return(list(max.marker.vals = ret))

}

get_attractors_from_graph_clustering <- function(f_name, col.names)
{
	tab <- read.table(f_name, header = T, sep = "\t", check.names = F, stringsAsFactors = F)
	print("Running graph based clustering")
	tab <- tab[, grep("cellType|popsize|sample", colnames(tab), invert = T)]
	tab <- tab[!(apply(tab[, col.names], 1, function(x) {all(x == 0)} )),]
	G <- build_graph(tab, col.names, filtering_T = 10)
	temp.G <- G
	E(temp.G)$weight <- E(temp.G)$weight * 100
	cc <- multilevel.community(temp.G)
	V(G)$Modularity.Class <- cc$membership


	ret <- get_vertex_table(G)
	ret$x <- ret$y <- NULL
	ret <- ret[, c("Modularity.Class", colnames(tab))]
	colnames(ret) <- gsub("Modularity.Class", "cellType", colnames(ret))

	ret <- ddply(ret, ~cellType, colwise(median))
	ret <- ret[order(ret[, "cellType"]),]

	return(list(tab.attractors = ret, att.labels = c(paste("community", as.character(1:nrow(ret)), sep = "_"))))

}

load_existing_layout <- function(scaffold.data)
{
	G <- scaffold.data$graphs[[1]]
	G <- induced.subgraph(G, V(G)$type == 1, impl = "copy_and_delete")
	tab <- get_vertex_table(G)
	V(G)$name <- 1:vcount(G)
	return(list(G.attractors = G, tab.attractors = tab))
}


# Map ---------------------------------------------------------------------

reactiveNetwork <- function (outputId)
{
	HTML(paste("<div id=\"", outputId, "\" class=\"shiny-network-output\"><svg /></div>", sep=""))
}

get_numeric_vertex_attributes <- function(sc.data, sel.graph)
{
	G <- sc.data$graphs[[sel.graph]]
	d <- get.data.frame(G, what = "vertices")
	#Don't consider attributes which are only present in the landmarks
	d <- d[d$type == 2,]
	num <- sapply(d, function(x) {is.numeric(x) && !any(is.na(x))})
	v <- list.vertex.attributes(G)[num]
	v <- v[grep("@", v, invert = T)]
	exclude <- c("x", "y", "cellType", "type", "groups", "r", "g", "b", "size", "DNA1", "DNA2", "BC1", "BC2", "BC3", "BC4", "BC5", "BC6", "Time", "Cell_length", "Cisplatin", "beadDist", "highest_scoring_edge")
	return(v[!(v %in% exclude)])
}

combine_marker_sample_name <- function(sel.marker, active.sample)
{
	if(active.sample == "All" || active.sample == "Absolute" || sel.marker == "Default")
		return(sel.marker)
	else
		return(paste(sel.marker, active.sample, sep = "@"))

}

get_sample_names <- function(sc.data, sel.graph)
{
	G <- sc.data$graphs[[sel.graph]]
	s <- list.vertex.attributes(G)
	s <- grep("@", s, value = T)
	ret <- sapply(strsplit(s, "@"), function (x) {x[[2]]})
	return(unique(ret))
}

get_graph <- function(sc.data, sel.graph, node.size.attr, min.node.size, max.node.size, landmark.node.size)
{
	G <- sc.data$graphs[[sel.graph]]
	edges <- data.frame(get.edgelist(G, names = F) - 1)
	colnames(edges) <- c("source", "target")
	svg.width <- 1200
	svg.height <- 800
	svg.center <- c(svg.width / 2, svg.height / 2)

	x <- V(G)$x
	y <- V(G)$y

	y <- -1 * y
	x <- x + abs(min(x))
	y <- y + abs(min(y))
	num.landmarks <- sum(V(G)$type == 1)
	trans <- get_graph_centering_transform(x[V(G)$type == 1], y[V(G)$type == 1], svg.width, svg.height)

	x <- (x / trans$scaling) - trans$offset.x
	y <- (y / trans$scaling) - trans$offset.y

	vertex.size <- get_vertex_size(sc.data, sel.graph, svg.width, node.size.attr, min.node.size, max.node.size, landmark.node.size)
	edges <- cbind(edges, x1 = x[edges[, "source"] + 1], x2 = x[edges[, "target"] + 1])
	edges <- cbind(edges, y1 = y[edges[, "source"] + 1], y2 = y[edges[, "target"] + 1])
	edges <- cbind(edges, id = 1:nrow(edges))
	edges <- cbind(edges, is_highest_scoring = 0)
	edges <- cbind(edges, edge_type = "")
	#Set as true for the highest scoring edges of type 2 vertices
	edges[, "is_highest_scoring"][V(G)$highest_scoring_edge[V(G)$type == 2]] <- 1
	if("edge_type" %in% list.edge.attributes(G)) #Old graphs did not have this
		edges[, "edge_type"] <- E(G)$edge_type
	ret <- list(names = V(G)$Label, size = vertex.size / trans$scaling, type = V(G)$type, highest_scoring_edge = V(G)$highest_scoring_edge, X = x, Y = y)
	ret <- c(ret, edges = list(edges))
	return(ret)
}

get_color_for_marker <- function(sc.data, sel.marker, rel.to.sample, sel.graph, active.sample, color.scaling,
								 stats.type, colors.to.interpolate, color.under, color.over, color.scale.limits = NULL, color.scale.mid = NULL)
{
	G <- sc.data$graphs[[sel.graph]]
	if(sel.marker == "Default")
	{
		ret <- rep("#4F93DE", vcount(G))
		ret[V(G)$type == 1] <- "#FF7580"
		return(list(color.vector = ret, color.scale.lim = NULL))
	}
	else
	{
		v <- get.vertex.attribute(G, combine_marker_sample_name(sel.marker, active.sample))

		f <- colorRamp(colors.to.interpolate, interpolate = "linear")

		if(rel.to.sample != "Absolute")
		{
			rel.to.marker <- combine_marker_sample_name(sel.marker, rel.to.sample)
			if(stats.type == "Difference")
				v <- v - (get.vertex.attribute(G, rel.to.marker))
			else if(stats.type == "Ratio")
				v <- v / (get.vertex.attribute(G, rel.to.marker))
			v[is.infinite(v)] <- NA
		}
		color.scale.lim <- NULL
		if(color.scaling == "local")
			color.scale.lim <- list(min = min(v, na.rm = T), max = max(v, na.rm = T))
		if(!is.null(color.scale.limits))
		{
			under <- v < color.scale.limits[1]
			over <- v > color.scale.limits[2]
			v[under] <- color.scale.limits[1]
			v[over] <- color.scale.limits[2]
			if(is.null(color.scale.mid))
				v <-  scales::rescale(v)
			else
				v <- scales::rescale_mid(v, mid = color.scale.mid)
			v <- f(v)
			v <- apply(v, 1, function(x) {sprintf("rgb(%s)", paste(round(x), collapse = ","))})
			v[under] <- sprintf("rgb(%s)", paste(col2rgb(color.under), collapse = ","))
			v[over] <- sprintf("rgb(%s)", paste(col2rgb(color.over), collapse = ","))

		}
		else
		{
			v <- f(scales::rescale(v)) #colorRamp needs an argument in the range [0, 1]
			v <- apply(v, 1, function(x) {sprintf("rgb(%s)", paste(round(x), collapse = ","))})
		}
		return(list(color.vector = v, color.scale.lim = color.scale.lim))
	}
}

get_number_of_cells_per_landmark <- function(sc.data, sel.graph)
{
	G <- sc.data$graphs[[sel.graph]]
	land <- V(G)[V(G)$type == 1]$Label
	ee <- get.edgelist(G)
	ee <- ee[V(G)[V(G)$type == 2]$highest_scoring_edge,]
	vv <- V(G)[as.numeric(ee[,2])]
	popsize <- V(G)[vv]$popsize
	dd <- data.frame(Landmark = ee[,1], popsize)
	dd <- ddply(dd, ~Landmark, function(x) {sum(x["popsize"])})
	dd <- cbind(dd, Percentage = dd$V1 / sum(dd$V1))
	names(dd) <- c("Landmark", "Cells", "Percentage")
	dd$Percentage <- signif(dd$Percentage * 100, digits = 4)
	return(dd)
}

get_summary_table <- function(sc.data, sel.graph, sel.nodes)
{
	G <- sc.data$graphs[[sel.graph]]
	col.names <- get_numeric_vertex_attributes(sc.data, sel.graph)
	tab <- get.data.frame(G, what = "vertices")
	temp <-tab[tab$Label %in% sel.nodes,]
	ret <- temp[, col.names]
	ret <- rbind(ret, apply(ret, 2, median, na.rm = T))
	popsize <- data.frame(Cells = temp$popsize, Percentage = temp$popsize / sum(tab$popsize[tab$type == 2]))
	popsize <- rbind(popsize, colSums(popsize))
	ret <- cbind(popsize, ret)
	ret <- data.frame(Label = c(temp$Label, "Summary"), ret)
	ret$Percentage <- signif(ret$Percentage * 100, digits = 4)
	return(ret)
}

plot_cluster <- function(data, clusters, graph.name, col.names, pool.cluster.data, plot.type)
{
	G <- data$graphs[[graph.name]]
	gated_data <- data$landmarks.data
	clustered_data <- data$clustered.data[[graph.name]]

	names(clustered_data) <- gsub("^X", "", names(clustered_data))
	names(gated_data) <- gsub("^X", "", names(gated_data))

	#This only works if the col.names are actually present in the clustered.data
	#TODO: figure out a consistent way to deal with panel mismatches

	common.names <- col.names[(col.names %in% names(clustered_data)) & (col.names %in% names(gated_data))]
	clustered_data <- clustered_data[, c(col.names, "cellType")]
	gated_data <- gated_data[, c(common.names, "cellType")]
	gated_data <- add_missing_columns(gated_data, col.names, fill.data = NA)
	#Select only the landmark nodes that are connected to these clusters
	land <- V(G)[nei(V(G)$Label %in% clusters)]$Label
	land <- V(G)[(V(G)$Label %in% land) & V(G)$type == 1]$Label
	temp <- gated_data[gated_data$cellType %in% land,]
	clus.num <- as.numeric(gsub("c", "", clusters))
	temp.clustered <- clustered_data[clustered_data$cellType %in% clus.num, ]
	if(pool.cluster.data)
		temp.clustered$cellType <- "Clusters"
	temp <- rbind(temp, temp.clustered)
	p <- NULL
	if(plot.type == "Scatterplot")
	{
		p <- density_scatterplot(temp, x_name = col.names[1], y_name = col.names[2], grouping = "cellType")
	}
	else
	{
		temp <- melt(temp, id.vars = "cellType")
		temp$variable <- as.factor(temp$variable)
		if(plot.type == "Density")
		{
			p <- ggplot(aes(x = value, color = cellType), data = temp) + geom_density() + facet_wrap(~variable, scales = "free")
		}
		else if(plot.type == "Boxplot")
		{
			p <- ggplot(aes(x = variable, fill = cellType, y = value), data = temp) + geom_boxplot()
		}
	}
	plot(p)
	return(p)
}

density_scatterplot  <- function(tab, x_name, y_name, grouping)
{
	m <- ddply(tab, grouping, function(m, x_name, y_name)
	{
		colramp <- grDevices::colorRampPalette(c("black", "red", "yellow"))
		dens.col <- grDevices::densCols(m[, x_name], m[, y_name], colramp = colramp)
		return(data.frame(m, dens.col = dens.col))
	}, x_name = x_name, y_name = y_name)

	maxx <- max(m[, x_name], na.rm = T) + 0.5
	maxy <- max(m[, y_name], na.rm = T) + 0.5

	(p <- ggplot(aes_string(x = x_name, y = y_name, color = "dens.col", size = 1), data = m)
		+ facet_wrap(grouping)
		+ geom_point()
		+ scale_colour_identity()
		+ scale_size_identity()
		+ xlim(0, maxx)
		+ ylim(0, maxy)
	)

	return(p)
}

export_clusters <- function(working.dir, sel.graph, sel.nodes)
{
	d <- gsub(".txt$", ".all_events.RData", sel.graph)
	d <- file.path(working.dir, d)
	d <- my_load(d)
	clus <- as.numeric(gsub("c", "", sel.nodes))
	d <- d[d$cellType %in% clus,]
	f <- flowFrame(as.matrix(d))
	p <- sprintf("scaffold_export_%s_", gsub(".fcs.clustered.txt", "", sel.graph))
	outname <- tempfile(pattern = p, tmpdir = working.dir, fileext = ".fcs")
	write.FCS(f, outname)
}

get_graph_centering_transform <- function(x, y, svg.width, svg.height)
{
	padding <- 50
	G.width <- max(x) - min(x)
	G.height <- max(y) - min(y)
	scaling <- max(c(G.width / (svg.width - (padding * 2)), G.height / (svg.height - (padding * 2))))

	x <- x / scaling
	y <- y / scaling

	offset.y <- min(y) - padding
	graph.x.center <- (min(x) + max(x)) / 2
	offset.x <- graph.x.center - (svg.width / 2)

	return(list(offset.x = offset.x, offset.y = offset.y, scaling = scaling))


}

get_vertex_size <- function(sc.data, sel.graph, figure.width, node.size.attr, min.node.size, max.node.size, landmark.node.size)
{
	G <- sc.data$graphs[[sel.graph]]
	size.attr <- get.vertex.attribute(G, node.size.attr)
	ret <- size.attr / sum(size.attr, na.rm = T)
	ret <- rescale_size(max.node.size, min.node.size, sc.data$dataset.statistics$max.marker.vals[["popsize.relative"]], ret)
	ret[V(G)$type == 1] <- landmark.node.size
	return(ret)
}

rescale_size <- function(max.size, min.size, max.val, x)
{
	return(((max.size - min.size) * x) / max.val + min.size);
}

add_missing_columns <- function(m, col.names, fill.data)
{
	v <- col.names[!(col.names %in% colnames(m))]
	print(sprintf("Adding missing columns: %s", paste(v, collapse = ", ")))
	ret <- matrix(nrow = nrow(m), ncol = length(v), data = fill.data)
	colnames(ret) <- v
	ret <- data.frame(m, ret, check.names = F)
	return(ret)
}

returnOrder <- function(inputId, vars) {
	tagList(
		singleton(tags$head(tags$script(src = 'sort.js'))),
		singleton(tags$head(tags$link(rel = 'stylesheet', type = 'text/css', href = 'sort.css'))),
		HTML(html_list(vars, inputId)),
		tags$script(paste0("$(function() {$( '#",inputId,"' ).sortable({placeholder: 'ui-state-highlight'}); $( '#",inputId,"' ).disableSelection(); });"))
	)
}

updateReturnOrder <- function(session, inputId, vars)
{
	session$sendInputMessage(inputId, list(value = vars))
}

html_list <- function(vars, id) {
	hl <- paste0("<ul id=\'",id,"\' class='stab'>")
	for(i in vars) hl <- paste0(hl, "<li class='ui-state-default stab'><span class='label'>",i,"</span></li>")
	paste0(hl, "</ul>")
}
