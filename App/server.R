##This file contains all the server-side process of the shiny application, grouped by tabs.

# Source, library, options ------------------------------------------------

library(Rcpp)
library(cluster)
library(flowCore)
library(ggplot2)
library(igraph)
library(plyr)
library(reshape)
library(shiny)
library(scales)
library(grDevices)
library(parallel)
library(jsonlite)
library(colourpicker)
library(tcltk)

options(expressions = 5e5,
        shiny.maxRequestSize = 20 * 1024 ^ 3) #increase max file upload size to 20 Gb

source("functions.R")
sourceCpp("forceatlas2.cpp")


# ShinyServer -------------------------------------------------------------


shinyServer(function(input, output, session) {
  # Reactives ---------------------------------------------------------------

  #Definition of the variables needed for each tab of the shiny app
  r <- reactiveValues(
    ##GLOBAL
    del_button_output = NULL,
    outputDirectory = NULL,
    inputDirectory = getwd(),

    ##CLUSTERING
    flow.frames = NULL,
    files.id = NULL,
    clustering.groups = NULL,
    datapaths = NULL,
    transform.data = NULL,
    cofactor = NULL,
    clustering_start = "NO",
    clustering_check = "NO",

    ##ANALYSIS
    scaffoldCoupleFiles = NULL,
    analysisFiles = list(),
    couplesid = NULL,
    inputFiles = NULL,

    a.cofactor = NULL,
    a.transform.data = NULL,

    loaded.rdata = list(),
    clustered.tables = NULL,
    clustered.txt = NULL,
    gated.flow.frames = NULL,
    gated.datapaths = NULL,
    gated.files.id = NULL,

    analysis_start = "NO",
    analysis_check = "NO",

    ##MAPPING
    file.scaffold = NULL,

    ##MAP DATASET
    d.file.scaffold.dataset = NULL,
    d.file.clustered.ref = NULL,
    d.files.clustered.input = NULL,
    d.files.rdata = list(),
    d.files.clustered.tables = NULL,
    d.clustered.txt = NULL,
    d.files.clustered.dataset = list(),
    d.files.clustered.id = NULL
  )

  # Clustering --------------------------------------------------------------
  #Server-side code running for the clustering tab

  #Function triggering when a file is input. It reads the file and appends it to the list of flow frames.
  observeEvent(input$fcsFiles, {
    progress <- Progress$new()
    progress$set(message = "Read Upload file", value = 1)
    filesname <-
      c(names(r$flow.frames),
        as.vector(input$fcsFiles$name))
    r$datapaths <- input$fcsFiles$datapath
    new.flow.frames <-
      lapply(as.vector(r$datapaths), function(x) {
        return(read.FCS(x, emptyValue = FALSE))
      })
    r$flow.frames  <- c(r$flow.frames, new.flow.frames)
    names(r$flow.frames) <- filesname
    progress$close()
  })

  #Function triggering when a clustering group needs to be append to the list of clustering groups.
  observeEvent(input$clusteringui_add_clustering_group, {
    files_list <- isolate({
      input$clusteringui_files_list
    })
    r$clustering.groups <-
      c(r$clustering.groups, setNames(list(files_list), files_list[1]))
  })

  #Usually users want to put all different files onto different groups, this function allows that.
  observeEvent(input$clusteringui_add_all_groups, {
    for (name in names(r$flow.frames)) {
      r$clustering.groups <-
        c(r$clustering.groups, setNames(list(name), name[1]))
    }
  })

  #Usually users want to add all markers for their clustering, this function allows it.
  observeEvent(input$add_all_markers_clustering, {
    if (!is.null(input$clusteringui_file_for_markers) &&
        grepl("*.fcs$", input$clusteringui_file_for_markers))
    {
      a <- (names(r$flow.frames) == input$clusteringui_file_for_markers)
      v <- as.vector(r$flow.frames[a][[1]]@parameters$name)
      updateSelectInput(
        session,
        "clusteringui_markers",
        selected = v,
        choices = c("", v)
      )
    }
  })

  #Function running when the user wants to select where to save files. Depending on the OS different functions are used because no one is actually portable on every OS.
  observeEvent(input$saveClustering, {
   r$outputDirectory <- chooseDir()
   if (is.na(r$outputDirectory)) {
     r$outputDirectory <- NULL
   }
  })

  #Function triggering whenever the user changes the transformation or the cofactor. The function stocks informations about the current selection.
  observe({
    if (input$clustering_transform == "Logicle") {
      r$transform.data <- "L"
      r$cofactor <- NULL
    }
    else if (input$clustering_transform == "Asinh") {
      r$transform.data <- "A"
      r$cofactor <- input$clusteringui_asinh_cofactor
    }
    else {
      r$transform.data <- F
      r$cofactor <- NULL
    }
  })

  #Function ensuring every change in the asinh cofactor will be saved.
  observeEvent(input$clustering_asinh_cofactor, {
    r$cofactor <- input$clusteringui_asinh_cofactor
  })

  #Function creating ui buttons to remove files from the list of FCS files.
  observe({
    if (length(r$flow.frames) < 1) {
      return(NULL)
    }
    r$files.id <-
      as.vector(unlist(lapply(c(
        1:length(r$flow.frames)
      ), function(x) {
        return(paste0(sample(letters, x + 1, replace = TRUE), collapse = ""))
      })))
    del_button_output <- lapply(c(1:length(r$files.id)), function(x) {
      del_button_name <- paste0("delButton_", r$files.id[x])
      del_button_object <-
        actionButton(del_button_name,
                     "",
                     icon = icon(name = "trash", lib = "glyphicon"))
      return(del_button_object)
    })
    do.call(tagList, del_button_output)
    output$buttonDel <- renderUI({
      del_button_output
    })
  })

  #Function triggering when a delete button is used, deleting the corresponding file from the list.
  observe({
    if (length(r$flow.frames) < 1) {
      return(NULL)
    }
    lapply(r$files.id, function(i) {
      observeEvent(input[[paste0("delButton_", i)]], {
        if (length(r$files.id) <= 1) {
          return(NULL)
        }
        index <- which(i == r$files.id)
        r$flow.frames <- r$flow.frames[-index]
      })
    })
  })

  #Updates the choice of markers to select when a reference file is choosed.
  observe({
    if (!is.null(input$clusteringui_file_for_markers) &&
        grepl("*.fcs$", input$clusteringui_file_for_markers))
    {
      a <- (names(r$flow.frames) == input$clusteringui_file_for_markers)
      v <- get_fcs_col_names(r$flow.frames[a][[1]])
      updateSelectInput(session, "clusteringui_markers", choices = v)
    }
  })

  #Removes groups from the list when the user press the corresponding button.
  observe({
    key <- input$clusteringui_remove_clustering_group$key
    if (!is.null(key) && key != "")
      isolate({
        r$clustering.groups[key] <- NULL
      })
  })

  #Shows the actual FCS files read by the app.
  output$flow_frames <- renderTable(expr = {
    if (!is.null(r$flow.frames)) {
      table1 <- as.matrix(names(r$flow.frames))
      colnames(table1) <- c("Selected Files:")
      return(table1)
    }
  },
  colnames = T,
  width = "100%")

  # output$listFiles <- renderPrint(names(r$flow.frames))

  #The two next functions are used to check wheter all minimum conditions are met before launching the clustering.
  observeEvent(input$clusteringui_start, {
    r$clustering_start <- "GO"
  })
  observe({
    print(!is.null(r$outputDirectory))
    print(r$outputDirectory != "")
    print(!is.null(r$flow.frames))
    print(!is.null(r$clustering.groups))
    print(!is.null(input$clusteringui_markers))
    print("-------------")
    if (!is.null(r$outputDirectory) &&
        r$outputDirectory != "" &&
        !is.null(r$flow.frames) &&
        !is.null(r$clustering.groups) &&
        !is.null(input$clusteringui_markers)) {
      r$clustering_check <- "GO"
    }
    else
    r$clustering_check <- "NO"
  })

  #When all conditions are met, executes the clustering and prints various informations about it.
  output$clusteringui_dialog <- renderText({
    if (r$clustering_start == "GO" && r$clustering_check =="GO") {
        isolate({
          col.names <- input$clusteringui_markers
          print(Sys.time())
          flow.frames <-
            get_flow_frames_from_groups(flow.frames = r$flow.frames,
                                        files.list = r$clustering.groups)
          files.analyzed <-
            cluster_fcs_files_groups(
              flow.frames = flow.frames,
              files.list = r$clustering.groups,
              col.names = col.names,
              num_clusters = input$clusteringui_num_clusters,
              num_samples = input$clusteringui_num_samples,
              cofactor = r$cofactor,
              downsample.to = input$clusteringui_downsample_to,
              output_dir = r$outputDirectory,
              transform.data = r$transform.data
            )
          ret <-
            sprintf(
              "Clustering completed with markers %s\n",
              paste(input$clusteringui_markers, collapse = " ")
            )
          ret <-
            paste(ret, sprintf(
              "Files analyzed:\n%s",
              paste(files.analyzed, collapse = "\n")
            ), sep = "")
          updateSelectInput(session,
                            "analysisui_reference",
                            choices = c(
                              "",
                              list.files(path = r$outputDirectory,
                                         pattern = "*.clustered.txt$")
                            ))
          return(ret)})}
    else if (r$clustering_check == "NO") {
      r$clustering_start <- "NO"
      return("Some fields are empty. Please check yout inputs.")
    }
    else
      return("Clustering ready.")
  })

  #Generates the selection widget to select the reference file for markers from input files.
  output$clusteringui1 <- renderUI({
    selectInput(
      "clusteringui_file_for_markers",
      "Load marker names from file",
      choices = c("", names(r$flow.frames)),
      width = "100%"
    )
  })

  #Generates the list of input files for the user to see and updates it whenever a new file is added.
  output$clusteringui2 <- renderUI({
    selectInput(
      "clusteringui_files_list",
      label = "File list",
      choices = names(r$flow.frames),
      selectize = F,
      multiple = T,
      width = "100%"
    )
  })

  #Generates and updates the list of groups for the clustering.
  output$clusteringui3 <- renderUI({
    dd <- r$clustering.groups
    return(tagList(mapply(
      get_cluster_groups_table, dd, names(dd), SIMPLIFY = F
    )))
  })

  #Shows the user the actual saving folder.
  output$saveClusteringFolder <- renderText({
    r$outputDirectory
  })



  # Analysis ----------------------------------------------------------------

  #Process the input of clustered files. Checks whether they are by pair or not, stocks the names and data into the corresponding variables.
  observeEvent(input$scaffoldCoupleFiles, {
    analysis.temp <- NULL
    couple.temp <- NULL
    liste.temp <- as.list(input$scaffoldCoupleFiles$name)

    clusteredList <-
      lapply(input$scaffoldCoupleFiles$name, function(x) {
        grep(pattern = ".clustered.txt", x, value = TRUE)
      })
    clusteredList2 <-
      lapply(clusteredList, function(x) {
        gsub(pattern = ".clustered.txt", replacement = "", x)
      })
    dataList <-
      lapply(input$scaffoldCoupleFiles$name, function(x) {
        grep(pattern = ".clustered.all_events.RData", x, value = TRUE)
      })
    dataList2 <-
      lapply(dataList, function(x) {
        gsub(pattern = ".clustered.all_events.RData", replacement = "", x)
      })
    for (i in clusteredList2) {
      for (j in dataList2) {
        if (identical(i, character(0)) == FALSE &&
            identical(j, character(0)) == FALSE && i == j) {
          analysis.temp[[length(analysis.temp) + 1]] <- i
        }
      }
    }

    couple.temp <- lapply(analysis.temp, function (x) {
      a <- NULL
      for (i in liste.temp) {
        if (identical(grep(pattern = x, i, value = TRUE), character(0)) == FALSE) {
          a <- rbind(a, grep(pattern = x, i, value = TRUE))
        }
      }
      return(as.vector(a))
    })

    r$analysisFiles <- c(r$analysisFiles, analysis.temp)
    r$scaffoldCoupleFiles <- c(r$scaffoldCoupleFiles, couple.temp)
    r$couplesid <-
      as.list(unlist(lapply(c(
        1:length(r$analysisFiles)
      ), function(x) {
        return(paste0(sample(letters, x + 1, replace = TRUE), collapse = ""))
      })))
    names(r$analysisFiles) <- r$couplesid
    names(r$couplesid) <- r$couplesid
    names(r$scaffoldCoupleFiles) <- r$couplesid

    r$clustered.txt <-
      (lapply(r$scaffoldCoupleFiles, function(x) {
        return(x[[2]])
      }))
    names(r$clustered.txt) <-
      (lapply(r$scaffoldCoupleFiles, function(x) {
        return(x[[2]])
      }))

    inputi <- NULL
    inputi <- input$scaffoldCoupleFiles$datapath
    names(inputi) <- input$scaffoldCoupleFiles$name

    clustered <- NULL
    rdata <- NULL
    for (elt in r$scaffoldCoupleFiles) {
      clustered <- c(clustered, elt[2])
      rdata <- c(rdata, elt[1])
    }

    for (elt in clustered) {
      if (!is.null(inputi)) {
        r$clustered.tables[[length(r$clustered.tables) + 1]] <-
          read.table(
            inputi[names(inputi) == elt],
            header = T,
            sep = "\t",
            check.names = F
          )
      }
    }

    names(r$clustered.tables) <- clustered

    for (elt in rdata) {
      if (!is.null(inputi)) {
        r$loaded.rdata[[length(r$loaded.rdata) + 1]] <-
          my_load(inputi[names(inputi) == elt])
      }
    }
    names(r$loaded.rdata) <- rdata
  })

  #When a deleting button for the file pairs is triggered, deletes the files from the different variables.
  observeEvent(input$analTrash$key, {
    key <- input$analTrash$key
    if (!is.null(key) && key != "") {
      r$scaffoldCoupleFiles[key] <- NULL
      r$analysisFiles[key] <- NULL
      r$couplesid <-
        as.list(unlist(lapply(c(
          1:length(r$analysisFiles)
        ), function(x) {
          return(paste0(sample(letters, x + 1, replace = TRUE), collapse = ""))
        })))
      if (length(r$scaffoldCoupleFiles) > 0) {
        names(r$analysisFiles) <- r$couplesid
        names(r$couplesid) <- r$couplesid
        names(r$scaffoldCoupleFiles) <- r$couplesid
      }
      r$clustered.txt <-
        (lapply(r$scaffoldCoupleFiles, function(x) {
          return(x[[2]])
        }))
      names(r$clustered.txt) <-
        (lapply(r$scaffoldCoupleFiles, function(x) {
          return(x[[2]])
        }))

    }
  })

  #Process the input of gated files.
  observeEvent(input$addFCSFiles, {
    progress <- Progress$new()
    progress$set(message = "Read Upload file", value = 1)
    filesname <-
      c(names(r$gated.flow.frames),
        as.vector(input$addFCSFiles$name))
    r$gated.datapaths <- input$addFCSFiles$datapath
    new.flow.frames <-
      lapply(as.vector(r$gated.datapaths), function(x) {
        return(read.FCS(x, emptyValue = FALSE))
      })
    r$gated.flow.frames  <- c(r$gated.flow.frames, new.flow.frames)
    names(r$gated.flow.frames) <- filesname
    progress$close()
  })

  #Adds all markers into the corresponding field when the user press the button.
  observeEvent(input$add_all_markers_analysis, {
    if (!is.null(input$reference_gated))
    {
      tab <-
        (r$clustered.tables[names(r$clustered.tables) == input$reference_gated])
      updateSelectInput(
        session,
        "analysisui_markers",
        selected = names(tab[[1]]),
        choices = c("", names(tab[[1]]))
      )
    }
  })

  #Adds all markers into the corresponding field when the user press the button.
  observeEvent(input$add_all_markers_inter_analysis, {
    if (!is.null(input$reference_gated))
    {
      tab <-
        (r$clustered.tables[names(r$clustered.tables) == input$reference_gated])
      updateSelectInput(
        session,
        "analysisui_markers_inter_cluster",
        selected = names(tab[[1]]),
        choices = c("", names(tab[[1]]))
      )
    }
  })

  #Shows the user the saving folder selected.
  output$saveAnalysisFolder <- renderText({
    r$outputDirectory
  })

  #Allows the user to choose a saving folder.
  observeEvent(input$saveAnalysis, {
    r$outputDirectory <-
      choose.dir(caption = "Select a Folder for saving")
  })

  #Checks any changes to the transformation field and updates the coresponding variables.
  observeEvent(input$analysis_transform, {
    if (input$analysis_transform == "Logicle") {
      r$a.transform.data <- "L"
      r$a.cofactor <- NULL
    }
    else if (input$analysis_transform == "Asinh") {
      r$a.transform.data <- "A"
      r$a.cofactor <- input$analysis_asinh_cofactor
    }
    else {
      r$a.transform.data <- F
      r$a.cofactor <- NULL
    }
  })

  #Checks any update of the cofactor for the transformation.
  observeEvent(input$analysis_asinh_cofactor, {
    r$a.cofactor <- input$analysis_asinh_cofactor
  })

  #Updates the choice of markers when a reference file is selected.
  observe({
    tab <-
      (r$clustered.tables[names(r$clustered.tables) == input$reference_gated])
    if (input$reference_gated != "") {
      updateSelectInput(session, "analysisui_markers", choices = names(tab[[1]]))
      updateSelectInput(session,
                        "analysisui_markers_inter_cluster",
                        choices = names(tab[[1]]))
    }
    else {
      updateSelectInput(session, "analysisui_markers", choices = c(""))
      updateSelectInput(session,
                        "analysisui_markers_inter_cluster",
                        choices = c(""))
    }
  })

  #Updates the choice of reference files when new files are added.
  observe({
    v <- r$clustered.txt
    updateSelectInput(session, "reference_gated", choices = c("", v))

  })

  #Function creating ui buttons to remove files from the list of gated FCS files.
  observe({
    if (length(r$gated.flow.frames) < 1) {
      return(NULL)
    }
    r$gated.files.id <-
      as.vector(unlist(lapply(c(
        1:length(r$gated.flow.frames)
      ), function(x) {
        return(paste0(sample(letters, x + 1, replace = TRUE), collapse = ""))
      })))
    del_button_output_gated <-
      lapply(c(1:length(r$gated.files.id)), function(x) {
        del_button_name <- paste0("gatedDelButton_", r$gated.files.id[x])
        del_button_object <-
          actionButton(del_button_name,
                       "",
                       icon = icon(name = "trash", lib = "glyphicon"))
        return(del_button_object)
      })
    do.call(tagList, del_button_output_gated)
    output$buttonDelGated <- renderUI({
      del_button_output_gated
    })
  })

  #Function triggering when a delete button is used, deleting the corresponding file from the list.
  observe({
    if (length(r$gated.flow.frames) < 1) {
      return(NULL)
    }
    lapply(r$gated.files.id, function(i) {
      observeEvent(input[[paste0("gatedDelButton_", i)]], {
        if (length(r$gated.files.id) <= 1) {
          return(NULL)
        }
        index <- which(i == r$gated.files.id)
        r$gated.flow.frames <- r$gated.flow.frames[-index]
      })
    })
  })

  #Function to initiate the graphic output of the couple files, processing them when they exist.
  output$scaffoldCoupleFiles <- renderUI({
    if (is.null(r$scaffoldCoupleFiles) ||
        !length(r$scaffoldCoupleFiles)) {
      return("No files selected...")
    }
    else {
      dd <- r$scaffoldCoupleFiles
      id <- r$couplesid
      return(tagList(
        mapply(get_file_couples, dd, id, "analTrash", SIMPLIFY = F)
      ))
    }
  })

  #Function used to visualize selected gated FCS files.
  output$gatedFiles <- renderTable({
    if (!is.null(r$gated.flow.frames)) {
      return(names(r$gated.flow.frames))
    }
    else {
      return("No files selected...")
    }
  }, width = "100%")

  #The two next functions are used to check wheter all minimum conditions are met before launching the SCAFFoLD analysis.
  observeEvent(input$analysisui_start, {
    r$analysis_start = "GO"
  })
  observe({
    if (!is.null(r$outputDirectory) &&
        (r$outputDirectory != "NA") &&
        !is.null(r$gated.flow.frames) &&
        !is.null(r$clustered.tables) &&
        !is.null(input$analysisui_markers)) {
      r$analysis_check <- "GO"
    }
    else
      r$analysis_check <- "NO"
  })

  #When all conditions are met, executes the analysis and prints various informations about it.
  output$analysisui_empty <- renderText({

    if (!is.null(input$analysisui_start) &&
        input$analysisui_start != 0) {
      print(input$analysisui_start)

        isolate({
          if (!is.null(input$reference_gated) &&
              input$reference_gated != "" &&
              !is.null(input$analysisui_markers) &&
              length(input$analysisui_markers) > 0)
          {
            files.analyzed <- NULL
            ew_influence <- NULL
            if (!is.null(input$analysisui_ew_influence_type)
                && input$analysisui_ew_influence_type == 'Fixed')
            {
              if (!is.null(input$analysisui_ew_influence))
                ew_influence <- input$analysisui_ew_influence
            }
            {
              files.analyzed <-
                run_analysis_gated(
                  r$gated.flow.frames,
                  r$clustered.tables,
                  r$loaded.rdata,
                  r$outputDirectory,
                  input$analysisui_markers,
                  inter.cluster.connections = input$analysisui_inter_cluster_connections,
                  col.names.inter_cluster = input$analysisui_markers_inter_cluster,
                  cofactor = r$a.cofactor,
                  ew_influence = ew_influence,
                  inter_cluster.weight_factor = input$analysisui_inter_cluster_weight,
                  overlap_method = "repel",
                  transform.data = r$a.transform.data
                )
            }
          }
          # updateSelectInput(session, "graphui_dataset", choices = c("", list.files(path = working.directory, pattern = "*.scaffold$")))
          ret <-
            sprintf(
              "Analysis completed with markers %s\n",
              paste(input$analysisui_markers, collapse = " ")
            )
          ret <-
            paste(ret, sprintf(
              "Files analyzed:\n%s",
              paste(files.analyzed, collapse = "\n")
            ), sep = "")
          return(ret)
        })}
    else if (r$analysis_check == "NO") {
      r$analysis_start <- "NO"
      return("Some fields are empty. Please check yout inputs.")
    }
    else
      return("")
  })

  #Shows the user the actual saving folder.
  output$saveAnalysisFolder <- renderText({
    r$outputDirectory
  })



  # Map ---------------------------------------------------------------------

  #When a .scaffold file is browsed, processes it.
  scaffold_data <- reactive({
    file_name <- input$graphui_dataset
    if (!is.null(file_name) && file_name != "")
    {
      print("Loading data...")
      data <- my_load(file_name$datapath)
      updateSelectInput(session,
                        "graphui_selected_graph",
                        choices = c("", names(data$graphs)))
      return(data)
    }
    else
      return(NULL)
  })

  #Every other function from this "Map" section remains as the original found on https://github.com/nolanlab/scaffold/tree/multiFilesClustering.

  output$graphui_mainnet <- reactive({
    ret <- get_main_graph()
    if (!is.null(ret))
    {
      ret$color <- get_color()
      ret$trans_to_apply <- isolate({
        input$graphui_cur_transform
      })
    }
    return(ret)
  })

  get_main_graph <- reactive({
    sc.data <- scaffold_data()
    if (!is.null(sc.data) &&
        !is.null(input$graphui_selected_graph) &&
        input$graphui_selected_graph != "")
    {
      attrs <-
        get_numeric_vertex_attributes(sc.data, input$graphui_selected_graph)
      node.size.attr <-
        combine_marker_sample_name("popsize", input$graphui_active_sample)

      isolate({
        sel.marker <- NULL
        if (input$graphui_marker %in% attrs)
          sel.marker <- input$graphui_marker
        else
          sel.marker <- "Default"
        updateSelectInput(
          session,
          "graphui_marker",
          choices = c("Default", attrs),
          selected = sel.marker
        )
        updateSelectInput(
          session,
          "graphui_markers_to_plot",
          choices = attrs,
          selected = attrs
        )
        sample.names <-
          get_sample_names(sc.data, input$graphui_selected_graph)
        updateSelectInput(
          session,
          "graphui_active_sample",
          choices = c("All", sample.names),
          selected = input$graphui_active_sample
        )
        updateSelectInput(
          session,
          "graphui_stats_relative_to",
          choices = c("Absolute", sample.names),
          selected = input$graphui_stats_relative_to
        )
      })
      return(
        get_graph(
          sc.data,
          input$graphui_selected_graph,
          node.size.attr,
          input$graphui_min_node_size,
          input$graphui_max_node_size,
          input$graphui_landmark_node_size
        )
      )
    }
    else
      return(NULL)
  })

  read_color_scale_info <- reactive({
    return(
      list(
        sel.marker = input$graphui_marker,
        color.scale.lim = input$graphui_color_scale_lim,
        color.scale.mid = input$graphui_color_scale_mid
      )
    )
  })

  get_color_scale <- reactive({
    #This code only updates the color scales
    sc.data <- scaffold_data()
    if (is.null(sc.data) || is.null(get_main_graph()))
      return(NULL)
    sel.marker <- input$graphui_marker
    rel.to <- input$graphui_stats_relative_to
    color.scaling <- input$graphui_color_scaling
    stats.type <- input$graphui_stats_type
    isolate({
      color <- NULL
      if (sel.marker != "")
      {
        #Colors are not really important here, only included because they need to be passed to the function
        min.color <- input$graphui_color_min
        mid.color <- input$graphui_color_mid
        max.color <- input$graphui_color_max
        under.color <- input$graphui_color_under
        over.color <- input$graphui_color_over
        color <-
          get_color_for_marker(
            sc.data,
            sel.marker,
            rel.to,
            input$graphui_selected_graph,
            input$graphui_active_sample,
            color.scaling,
            stats.type,
            colors.to.interpolate = c(min.color, mid.color, max.color),
            under.color,
            over.color
          )
        if (!is.null(color$color.scale.lim)
            && !(is.null(color.scaling)) && color.scaling == "local")
        {
          updateSliderInput(
            session,
            "graphui_color_scale_lim",
            min = color$color.scale.lim$min,
            max = color$color.scale.lim$max,
            step = 0.1,
            value = c(
              color$color.scale.lim$min,
              color$color.scale.lim$max
            )
          )
          updateSliderInput(
            session,
            "graphui_color_scale_mid",
            min = color$color.scale.lim$min,
            max = color$color.scale.lim$max,
            step = 0.1,
            value = mean(
              c(
                color$color.scale.lim$min,
                color$color.scale.lim$max
              )
            )
          )
        }
      }
    })
  })

  get_color <- reactive({
    #This code does the actual coloring
    get_color_scale()
    color.scale.info <- read_color_scale_info()
    min.color <- input$graphui_color_min
    mid.color <- input$graphui_color_mid
    max.color <- input$graphui_color_max
    under.color <- input$graphui_color_under
    over.color <- input$graphui_color_over
    color.scale.lim <- color.scale.info$color.scale.lim
    colors.to.interpolate <- NULL
    color.scale.mid <- NULL
    if (input$graphui_color_number == 3)
    {
      colors.to.interpolate <- c(min.color, mid.color, max.color)
      color.scale.mid <- color.scale.info$color.scale.mid
    }
    else
      colors.to.interpolate <- c(min.color, max.color)
    return(isolate({
      sel.marker <- color.scale.info$sel.marker

      color.vector <- NULL
      active.sample <- input$graphui_active_sample
      rel.to <- input$graphui_stats_relative_to
      color.scaling <- input$graphui_color_scaling
      stats.type <- input$graphui_stats_type

      if (sel.marker != "")
      {
        sc.data <- scaffold_data()
        if (!is.null(sc.data))
        {
          color <-
            get_color_for_marker(
              sc.data,
              sel.marker,
              rel.to,
              input$graphui_selected_graph,
              active.sample,
              color.scaling,
              stats.type,
              colors.to.interpolate = colors.to.interpolate,
              under.color,
              over.color,
              color.scale.limits = color.scale.lim,
              color.scale.mid = color.scale.mid
            )
          color.vector <- color$color.vector
        }
      }
      return(color.vector)
    }))
  })



  output$graphui_mainnet <- reactive({
    ret <- get_main_graph()
    if (!is.null(ret))
    {
      ret$color <- get_color()
      ret$trans_to_apply <- isolate({
        input$graphui_cur_transform
      })
    }
    return(ret)
  })

  output$graphui_table <- renderDataTable({
    sc.data <- scaffold_data()
    if (!is.null(sc.data) &&
        !is.null(input$graphui_selected_graph) &&
        input$graphui_selected_graph != "")
    {
      if (is.null(input$graphui_selected_nodes) ||
          length(input$graphui_selected_nodes) == 0)
      {
        get_number_of_cells_per_landmark(scaffold_data(),
                                         input$graphui_selected_graph)
      }
      else
      {
        get_summary_table(
          scaffold_data(),
          input$graphui_selected_graph,
          input$graphui_selected_nodes
        )
      }
    }
  }, options = list(
    scrollX = TRUE,
    searching = FALSE,
    scrollY = "800px",
    paging = FALSE,
    info = FALSE,
    processing = FALSE
  ))

  output$graphui_dialog1 <- reactive({
    sc.data <- scaffold_data()
    ret <- ""
    if (!is.null(sc.data))
      ret <-
      sprintf("Markers used for SCAFFoLD: %s",
              paste(sc.data$scaffold.col.names, collapse = ", "))
    return(ret)
  })


  output$graphui_plot = renderPlot({
    p <- NULL
    if (!is.null(input$graphui_plot_clusters) &&
        input$graphui_plot_clusters != 0)
    {
      isolate({
        col.names <- input$graphui_markers_to_plot
        if ((length(col.names) >= 1) &&
            (length(input$graphui_selected_nodes) >= 1))
          p <-
            plot_cluster(
              scaffold_data(),
              input$graphui_selected_nodes,
              input$graphui_selected_graph,
              input$graphui_markers_to_plot,
              input$graphui_pool_cluster_data,
              input$graphui_plot_type
            )
      })
    }
  })

  observe({
    if (!is.null(input$graphui_reset_colors) &&
        input$graphui_reset_colors != 0)
    {
      session$sendCustomMessage(type = "reset_colors", "none")
    }
  })

  observe({
    if (!is.null(input$graphui_reset_graph_position) &&
        input$graphui_reset_graph_position != 0)
    {
      session$sendCustomMessage(type = "reset_graph_position", "none")
    }
  })

  observe({
    if (!is.null(input$graphui_toggle_landmark_labels) &&
        input$graphui_toggle_landmark_labels != 0)
    {
      display <-
        ifelse(input$graphui_toggle_landmark_labels %% 2 == 0,
               "",
               "none")
      session$sendCustomMessage(type = "toggle_label", list(target = "landmark", display = display))
    }
  })

  observe({
    display_edges <- input$graphui_display_edges
    session$sendCustomMessage(type = "toggle_display_edges", display_edges)
  })

  observe({
    if (!is.null(input$graphui_toggle_cluster_labels) &&
        input$graphui_toggle_cluster_labels != 0)
    {
      display <-
        ifelse(input$graphui_toggle_cluster_labels %% 2 == 0,
               "none",
               "")
      session$sendCustomMessage(type = "toggle_label", list(target = "cluster", display = display))
    }
  })

  observe({
    display <- tolower(input$graphui_node_size)
    session$sendCustomMessage(type = "toggle_node_size", list(display = display))
  })


  observe({
    if (!is.null(input$graphui_toggle_node_size) &&
        input$graphui_toggle_node_size != 0)
    {
      display <-
        ifelse(input$graphui_toggle_node_size %% 2 == 0,
               "proportional",
               "default")
      session$sendCustomMessage(type = "toggle_node_size", list(display = display))
    }
  })
  # Map dataset -------------------------------------------------------------

  # observe({
  # 	if(!is.null(input$graphui_color_scaling) && input$graphui_color_scaling == "global")
  # 	{
  # 		updateSliderInput(session, "graphui_color_scale_lim", min = input$graphui_color_scale_min,
  # 						  max = input$graphui_color_scale_max)
  # 	}
  # })


  #Allows the user to select a folder where to save output files.
  observeEvent(input$saveDataset, {
    r$outputDirectory <-
      choose.dir(caption = "Select a Folder for saving")
  })

  #Updates marker choice when a reference clustered file is selected.
  observe({
    if (!is.null(input$mappingui_reference) &&
        input$mappingui_reference != "") {
      tab <-
        names(r$d.files.clustered.tables[names(r$d.files.clustered.tables) == input$mappingui_reference][[1]])
      updateSelectInput(session,
                        "mappingui_sample_clustered_file_markers",
                        choices = tab)
      updateSelectInput(session, "mappingui_markers_inter_cluster", choices = tab)
    }
    else {
      tab <- ""
      updateSelectInput(session,
                        "mappingui_sample_clustered_file_markers",
                        choices = tab)
      updateSelectInput(session, "mappingui_markers_inter_cluster", choices = tab)
    }
  })

  #Checks input pair clustered files and adds them to the different related variables as well as the list visible for the user.
  observeEvent(input$mappingui_added_files, {
    analysis.temp <- NULL
    couple.temp <- NULL
    liste.temp <- as.list(input$mappingui_added_files$name)

    clusteredList <-
      lapply(input$mappingui_added_files$name, function(x) {
        grep(pattern = ".clustered.txt", x, value = TRUE)
      })
    clusteredList2 <-
      lapply(clusteredList, function(x) {
        gsub(pattern = ".clustered.txt", replacement = "", x)
      })
    dataList <-
      lapply(input$mappingui_added_files$name, function(x) {
        grep(pattern = ".clustered.all_events.RData", x, value = TRUE)
      })
    dataList2 <-
      lapply(dataList, function(x) {
        gsub(pattern = ".clustered.all_events.RData", replacement = "", x)
      })
    for (i in clusteredList2) {
      for (j in dataList2) {
        if (identical(i, character(0)) == FALSE &&
            identical(j, character(0)) == FALSE && i == j) {
          analysis.temp[[length(analysis.temp) + 1]] <- i
        }
      }
    }
    couple.temp <- lapply(analysis.temp, function (x) {
      a <- NULL
      for (i in liste.temp) {
        if (identical(grep(pattern = x, i, value = TRUE), character(0)) == FALSE) {
          a <- rbind(a, grep(pattern = x, i, value = TRUE))
        }
      }
      return(as.vector(a))
    })
    r$d.files.clustered.dataset <-
      c(r$d.files.clustered.dataset, analysis.temp)
    r$d.files.clustered.input <-
      c(r$d.files.clustered.input, couple.temp)
    r$d.files.clustered.id <-
      as.list(unlist(lapply(c(
        1:length(r$d.files.clustered.dataset)
      ), function(x) {
        return(paste0(sample(letters, x + 1, replace = TRUE), collapse = ""))
      })))
    if (length(r$d.files.clustered.input) > 0) {
      names(r$d.files.clustered.id) <- r$d.files.clustered.id
      names(r$d.files.clustered.dataset) <- r$d.files.clustered.id
      names(r$d.files.clustered.input) <- r$d.files.clustered.id


      r$d.clustered.txt <-
        (lapply(r$d.files.clustered.input, function(x) {
          return(x[[2]])
        }))
      names(r$d.clustered.txt) <-
        (lapply(r$d.files.clustered.input, function(x) {
          return(x[[2]])
        }))

      inputi <- NULL
      inputi <- input$mappingui_added_files$datapath
      names(inputi) <- input$mappingui_added_files$name

      clustered <- NULL
      rdata <- NULL
      for (elt in couple.temp) {
        clustered <- c(clustered, elt[2])
        rdata <- c(rdata, elt[1])
      }

      for (elt in clustered) {
        if (!is.null(inputi)) {
          r$d.files.clustered.tables[[length(r$d.files.clustered.tables) + 1]] <-
            read.table(
              inputi[names(inputi) == elt],
              header = T,
              sep = "\t",
              check.names = F
            )
        }
      }
      names(r$d.files.clustered.tables) <- clustered

      for (elt in rdata) {
        if (!is.null(inputi)) {
          r$d.files.rdata[[length(r$d.files.rdata) + 1]] <-
            my_load(inputi[names(inputi) == elt])
        }
      }
      names(r$d.files.rdata) <- rdata
    }
  })

  #Checks for delete buttons triggering to remove related pairs.
  observeEvent(input$mappingTrash$key, {
    key <- input$mappingTrash$key
    if (!is.null(key) && key != "") {
      r$d.files.clustered.input[key] <- NULL
      r$d.files.clustered.dataset[key] <- NULL
      r$d.files.clustered.id <-
        as.list(unlist(lapply(c(
          1:length(r$d.files.clustered.dataset)
        ), function(x) {
          return(paste0(sample(letters, x + 1, replace = TRUE), collapse = ""))
        })))
      if (length(r$d.files.clustered.input) > 0) {
        names(r$d.files.clustered.dataset) <- r$d.files.clustered.id
        names(r$d.files.clustered.id) <- r$d.files.clustered.id
        names(r$d.files.clustered.input) <- r$d.files.clustered.id
      }
      r$d.clustered.txt <-
        (lapply(r$d.files.clustered.input, function(x) {
          return(x[[2]])
        }))
      names(r$d.clustered.txt) <-
        (lapply(r$d.files.clustered.input, function(x) {
          return(x[[2]])
        }))

    }
  })

  #Check for the input of a .scaffold file and updates the related variables and fields.
  observeEvent(input$mappingui_ref_scaffold_file, {
    if (!is.null(input$mappingui_ref_scaffold_file) &&
        input$mappingui_ref_scaffold_file != "")
    {
      file_name <- input$mappingui_ref_scaffold_file$datapath
      sc.data <- my_load(file_name)
      r$d.file.scaffold.dataset <- sc.data
      updateSelectInput(
        session,
        "mappingui_ref_scaffold_file_markers",
        choices = sc.data$scaffold.col.names
      )
    }
  })

  #Updates list of markers for clustered files in the corresponding box.
  observe({
    if (!is.null(input$mappingui_sample_clustered_file_markers) &&
        length(input$mappingui_sample_clustered_file_markers > 0))
    {
      updateReturnOrder(
        session,
        "mappingui_clustered_markers_list",
        input$mappingui_sample_clustered_file_markers
      )
    }
  })

  #Updates the list of markers for SCAFFoLD file in the corresponding box.
  observe({
    if (!is.null(input$mappingui_ref_scaffold_file_markers) &&
        length(input$mappingui_ref_scaffold_file_markers > 0))
    {
      updateReturnOrder(
        session,
        "mappingui_ref_markers_list",
        input$mappingui_ref_scaffold_file_markers
      )
    }
  })

  #Updates the list of maps selectable to load markers.
  observe({
    v <- r$d.clustered.txt
    updateSelectInput(session, "mappingui_reference", choices = c("", v))
  })


#Updtaes the list of valid selected pairs.
  output$mappingui_valid_couples <- renderUI({
    if (is.null(r$d.files.clustered.input) ||
        !length(r$d.files.clustered.input)) {
      return("No files selected...")
    }
    else {
      dd <- r$d.files.clustered.input
      id <- r$d.files.clustered.id
      return(tagList(
        mapply(get_file_couples, dd, id, "mappingTrash", SIMPLIFY = F)
      ))
    }
  })

  #Shows the user the current saving folder.
  output$saveDatasetFolder <- renderText({
    r$outputDirectory
  })

  #Runs the analysis and processes the files.
  output$mappingui_dialog <- renderText({
    if (!is.null(input$mappingui_start) && input$mappingui_start != 0)
      isolate({
        col.names <- input$mappingui_clustered_markers_list
        ref.col.names <- input$mappingui_ref_markers_list
        names.map <- ref.col.names
        #Missing values (i.e. non-mapped markers) are filled with NA
        names(names.map) <- col.names
        ew_influence <- NULL
        if (!is.null(input$mappingui_ew_influence_type) &&
            input$mappingui_ew_influence_type == 'Fixed')
        {
          if (!is.null(input$mappingui_ew_influence))
            ew_influence <- input$mappingui_ew_influence
        }

        run_analysis_existing(
          r$d.file.scaffold.dataset,
          input$mappingui_reference,
          r$d.files.clustered.tables,
          r$d.files.rdata,
          r$outputDirectory,
          input$mappingui_ref_markers_list,
          inter.cluster.connections = input$mappingui_inter_cluster_connections,
          names.map = names.map,
          col.names.inter_cluster = input$mappingui_markers_inter_cluster,
          inter_cluster.weight_factor = input$mappingui_inter_cluster_weight,
          overlap_method = input$mappingui_overlap_method,
          ew_influence = ew_influence
        )


        # updateSelectInput(session, "graphui_dataset", choices = c("", list.files(path = working.directory, pattern = "*.scaffold$")))
        ret <-
          sprintf(
            "Analysis completed with markers %s\n",
            paste(input$mappingui_ref_scaffold_fil, collapse = " ")
          )
        if (!is.null(names.map))
          ret <-
          paste(ret, sprintf(
            "Mapping: %s -> %s\n",
            paste(names.map, collapse = " "),
            paste(names(names.map), collapse = " ")
          ), sep = "")
        return(ret)
      })
  })

})
