##This file contains the ui of the shiny application.
##Each section represents a tab of the interface.

source("functions.R")
library(shinythemes)

shinyUI(navbarPage( theme = shinytheme("darkly"), #Theme can be changed here.
	"SCAFFOLD - CIPHOLD v0.2",

	# Run clustering ----------------------------------------------------------

	tabPanel(
		"Run clustering",
		fluidPage(
		  fluidRow(
			column(
				6,
				uiOutput("clusteringui1")
			),
			column(
				6,
				fileInput(
					"fcsFiles",
					"Choose FCS files",
					accept = ".fcs",
					multiple = TRUE
				)
			)
		),
		fluidRow(
			column(6,
				   selectInput(
				   	"clusteringui_markers",
				   	"Choose the markers for clustering",
				   	choices = c(""),
				   	multiple = T,
				   	width = "100%"
				   ),
				   actionButton("add_all_markers_clustering", "Add All")
			),
			column(6, id = "align_button",
				   tableOutput("flow_frames"),
				   uiOutput("buttonDel")
				   )
		),
		fluidRow(column(
			6,
			uiOutput("clusteringui2"),
			column(6,
				actionButton("clusteringui_add_clustering_group", "Add clustering group")
				),
			column(6,
				actionButton("clusteringui_add_all_groups", "Add all files into groups")
				)
		),
		column(6,
			   p(strong("Groups selected:")),
			   uiOutput("clusteringui3"))),
		fluidRow(
			column(
				8,
				numericInput(
					"clusteringui_downsample_to",
					"Pre-downsample data to (0 means no pre-downsampling, only valid for multiple files clustering)",
					value = 0
				),
				numericInput(
					"clusteringui_num_clusters",
					"Number of clusters",
					value = 200,
					min = 1,
					max = 2000
				),
				numericInput(
					"clusteringui_num_samples",
					"Number of samples",
					value = 50,
					min = 1
				),
				selectInput(
					"clustering_transform",
					"Choose a transformation",
					choices = c("None", "Asinh", "Logicle")
				),
				conditionalPanel(
					condition <- "input.clustering_transform == 'Asinh'",
					numericInput("clusteringui_asinh_cofactor", "asinh cofactor", value = 5)
				),
				br(),
				br(),
				actionButton("saveClustering", "Saving Folder"),
				verbatimTextOutput("saveClusteringFolder"),
				br(),
				br(),
				actionButton("clusteringui_start", "Start clustering"),
				br(),
				br(),
				conditionalPanel(
					condition <- "$('html').hasClass('shiny-busy')",
					br(),
					p(strong("Processing data...please wait."))
				),
				conditionalPanel(
					condition <-
						"!$('html').hasClass('shiny-busy') && input.clusteringui_start > 0",
					br(),
					p(strong("Data processing is complete!"))
				)
			)
		),
		fluidRow(column(
			12,
			verbatimTextOutput("clusteringui_dialog"),
			br(),
			br(),
			br(),
			br(),
			br(),
			br()
		))

	)),

	# Scaffold analysis -------------------------------------------------------

	tabPanel(
		"Run SCAFFoLD analysis",
		fluidPage(
			fluidRow(column(
				6,
				fileInput(
					"scaffoldCoupleFiles",
					"Add your .txt/.RData pairs",
					multiple = TRUE,
					accept = c(".clustered.txt", ".RData")
				)
			),
			column(
				6,
				p(strong("Valid selected pairs:")),
				uiOutput("scaffoldCoupleFiles")
			)),
			fluidRow(column(
				6,
				fileInput(
					"addFCSFiles",
					"Add your Gated Files",
					multiple = TRUE,
					accept = ".fcs"
				)
			),
			column(
				6, id = "align_button",
				p(strong("Selected gated files:")),
				uiOutput("gatedFiles"),
				uiOutput("buttonDelGated")
			)),
			fluidRow(
				column(
					6,
					selectInput(
						"reference_gated",
						"Choose a reference file (to load markers)",
						choices = c(""),
						multiple = F,
						width = "100%"
					),
					selectInput(
						"analysisui_markers",
						"Choose the markers for SCAFFoLD",
						choices = c(""),
						multiple = T,
						width = "100%"
					),
					actionButton("add_all_markers_analysis", "Add All"),
					br(),
					br(),
					selectInput(
						"analysisui_ew_influence_type",
						"Edge weight influence",
						choices = c("Proportional", "Fixed"),
						width = "100%"
					),
					conditionalPanel(
						condition = "input.analysisui_ew_influence_type == 'Fixed'",
						numericInput("analysisui_ew_influence", "Specifiy Edge weight value", 12),
						br()
					),
					checkboxInput(
						"analysisui_inter_cluster_connections",
						"Add inter-cluster connections",
						value = FALSE
					),
					conditionalPanel(
						condition = "input.analysisui_inter_cluster_connections == true",
						selectInput(
							"analysisui_markers_inter_cluster",
							"Markers for inter-cluster connections (if different)",
							choices = c(""),
							multiple = T,
							width = "100%"
						),
						actionButton("add_all_markers_inter_analysis", "Add All"),
						br(),
						br(),
						numericInput(
							"analysisui_inter_cluster_weight",
							"Weight factor for inter-cluster connections",
							0.7,
							min = 0,
							max = 10,
							step = 0.1
						),
						br()
					),
					selectInput(
						"analysis_transform",
						"Choose a transformation",
						choices = c("None", "Asinh", "Logicle")
					),
					conditionalPanel(
						condition <- "input.analysis_transform == 'Asinh'",
						numericInput("analysis_asinh_cofactor", "asinh cofactor", value = 5)
					),
					br(),
					br(),
					actionButton("saveAnalysis", "Saving Folder"),
					verbatimTextOutput("saveAnalysisFolder"),
					br(),
					br(),
					actionButton("analysisui_start", "Start analysis"),
					br(),
					br(),
					conditionalPanel(
						condition <- "$('html').hasClass('shiny-busy')",
						br(),
						p(strong("Processing data...please wait."))
					),
					conditionalPanel(
						condition <-
							"!$('html').hasClass('shiny-busy') && input.analysisui_start > 0",
						br(),
						p(strong("Data processing is complete!"))
					)
				)
			),
			fluidRow(column(
				12,
				verbatimTextOutput("analysisui_empty"),
				br(),
				br(),
				br(),
				br(),
				br(),
				br()
			))
		)
	),

	# Map exploration ---------------------------------------------------------


	tabPanel("Map exploration",
			 fluidPage(
			 	fluidRow(
			 		column(
			 			9,
			 			tags$head(tags$script(src = "d3.min.js")),
			 			tags$head(tags$script(src = "graph.js")),
			 			tags$head(tags$script(src = "rect_select.js")),
			 			singleton(tags$head(
			 				tags$link(rel = 'stylesheet', type = 'text/css', href = 'rect_select.css')
			 			)),
			 			singleton(tags$head(
			 				tags$link(rel = 'stylesheet', type = 'text/css', href = 'graph.css')
			 			)),
			 			conditionalPanel(
			 				condition = !is.null("output.graphui_mainnet"),
			 				reactiveNetwork(outputId = "graphui_mainnet")
			 			),

			 			dataTableOutput("graphui_table")
			 		),
			 		column(
			 			3,

			 			fileInput(
			 				"graphui_dataset",
			 				"Choose a dataset",
			 				multiple = F,
			 				accept = ".scaffold"
			 			),
			 			selectizeInput(
			 				"graphui_selected_graph",
			 				"Choose a graph:",
			 				choices = c(""),
			 				width = "100%"
			 			),
			 			selectizeInput(
			 				"graphui_active_sample",
			 				"Active sample",
			 				choices = c("All"),
			 				width = "100%"
			 			),
			 			selectInput(
			 				"graphui_marker",
			 				"Nodes color:",
			 				choices = c("Default"),
			 				width = "100%"
			 			),
			 			fluidRow(column(
			 				6,
			 				selectInput(
			 					"graphui_stats_type",
			 					"Stats type",
			 					choices = c("Ratio", "Difference")
			 				)
			 			),
			 			column(
			 				6,
			 				selectInput(
			 					"graphui_stats_relative_to",
			 					"Stats relative to:",
			 					choices = c("Absolute"),
			 					width = "100%"
			 				)
			 			)),
			 			selectInput(
			 				"graphui_color_scaling",
			 				"Color scaling:",
			 				choices = c("global", "local"),
			 				width = "100%"
			 			),
			 			h4("Colors for scale"),
			 			selectInput("graphui_color_number", "Number of colors", choices = c(2, 3)),
			 			fluidRow(
			 				column(
			 					6,
			 					colourpicker::colourInput("graphui_color_under", "Under:", value = "#FFFF00")
			 				),
			 				column(
			 					6,
			 					colourpicker::colourInput("graphui_color_over", "Over:", value = "#0000FF")
			 				)
			 			),
			 			fluidRow(
			 				column(
			 					4,
			 					colourpicker::colourInput("graphui_color_min", "Min:", value = "#E7E7E7")
			 				),
			 				column(
			 					4,
			 					conditionalPanel(
			 						condition = "input.graphui_color_number == 3",
			 						colourpicker::colourInput("graphui_color_mid", "Mid:", value = "#E7E7E7")
			 					)
			 				),
			 				column(
			 					4,
			 					colourpicker::colourInput("graphui_color_max", "Max:", value = "#E71601")
			 				)
			 			),
			 			conditionalPanel(
			 				condition = "input.graphui_color_number == 3",
			 				sliderInput(
			 					"graphui_color_scale_mid",
			 					"Color scale midpoint",
			 					min = 0.0,
			 					max = 5.0,
			 					value = 2.5,
			 					round = -2,
			 					step = 0.1,
			 					sep = ""
			 				)
			 			),
			 			sliderInput(
			 				"graphui_color_scale_lim",
			 				"Color scale limits",
			 				min = 0.0,
			 				max = 5.0,
			 				value = c(0.0, 5.0),
			 				round = -2,
			 				step = 0.1,
			 				sep = ""
			 			),
			 			fluidRow(column(
			 				6,
			 				numericInput("graphui_color_scale_min", "Color scale min:", 0)
			 			),
			 			column(
			 				6,
			 				numericInput("graphui_color_scale_max", "Color scale max:", 5)
			 			)),
			 			fluidRow(column(
			 				6,
			 				selectInput(
			 					"graphui_node_size",
			 					"Nodes size:",
			 					choices = c("Proportional", "Default"),
			 					width = "100%"
			 				)
			 			),
			 			column(
			 				6,
			 				numericInput(
			 					"graphui_min_node_size",
			 					"Minimum node size",
			 					2,
			 					min = 0,
			 					max = 1000
			 				)
			 			)),
			 			fluidRow(column(
			 				6,
			 				numericInput(
			 					"graphui_max_node_size",
			 					"Maximum node size",
			 					60,
			 					min = 0,
			 					max = 1000
			 				)
			 			),
			 			column(
			 				6,
			 				numericInput(
			 					"graphui_landmark_node_size",
			 					"Landmark node size",
			 					8,
			 					min = 0,
			 					max = 1000
			 				)
			 			)),
			 			selectInput(
			 				"graphui_display_edges",
			 				"Display edges:",
			 				choices = c("All", "Highest scoring", "Inter cluster", "To landmark"),
			 				width = "100%"
			 			),
			 			br(),
			 			actionButton("graphui_reset_graph_position", "Reset graph position"),
			 			br(),
			 			actionButton("graphui_toggle_landmark_labels", "Toggle landmark labels"),
			 			br(),
			 			actionButton("graphui_toggle_cluster_labels", "Toggle cluster labels"),
			 			br(),
			 			actionButton("graphui_plot_clusters", "Plot selected clusters"),
			 			checkboxInput("graphui_pool_cluster_data", "Pool cluster data", value = FALSE),
			 			br(),
			 			selectInput(
			 				"graphui_plot_type",
			 				"Plot type:",
			 				choices = c("Density", "Boxplot", "Scatterplot"),
			 				width = "100%"
			 			),
			 			selectInput(
			 				"graphui_markers_to_plot",
			 				"Markers to plot in cluster view:",
			 				choices = c(""),
			 				multiple = T,
			 				width = "100%"
			 			),
			 			verbatimTextOutput("graphui_dialog1")
			 		)
			 	),
			 	fluidRow(column(12,
			 					plotOutput("graphui_plot")))
			 )),


	# Map dataset -------------------------------------------------------------

	tabPanel("Map dataset",
			 fluidPage(
			 	tags$head(tags$script(src = "jquery-ui.min.js")),
			 	singleton(tags$head(
			 		tags$link(rel = 'stylesheet', type = 'text/css', href = 'custom.css')
			 	)),
			 	fluidRow(
			 		column(6,
			 			   fileInput("mappingui_ref_scaffold_file", "Select reference SCAFFoLD file", accept = ".scaffold")
			 		),
			 		column(
			 			6,
			 			fileInput(
			 				"mappingui_added_files",
			 				"Add your .txt/.RData couples",
			 				multiple = TRUE,
			 				accept = c(".clustered.txt", ".RData")
			 			),
			 			p(strong("Valid selected couples:")),
			 			uiOutput("mappingui_valid_couples"),
			 			selectInput(
			 				"mappingui_reference",
			 				"Choose a reference file (to load markers)",
			 				choices = c(""),
			 				multiple = F,
			 				width = "100%"
			 			)
			 		)),
			 	fluidRow(
			 		column(6,
			 					selectInput("mappingui_ref_scaffold_file_markers", "Select the markers to include in the mapping", choices = c(""), multiple = T, width = "100%"),
			 					actionButton("add_mappingui_ref_scaffold_file_markers", "Add All"),
			 					br(), br(),
			 					wellPanel(returnOrder("mappingui_ref_markers_list", c(""))),
			 					br(), br(), br()
			 	),
			 	column(6,
			 			   selectInput("mappingui_sample_clustered_file_markers", "Select the markers to include in the mapping", choices = c(""), multiple = T, width = "100%"),
			 			   actionButton("add_mappingui_sample_clustered_file_markers", "Add All"),
			 			   br(), br(),
			 			   wellPanel(returnOrder("mappingui_clustered_markers_list", c(""))),
			 			   br(), br(), br()
			 		)
			 	),
			 	fluidRow(
			 		column(12,
			 			   selectInput("mappingui_ew_influence_type", "Edge weight influence", choices = c("Proportional", "Fixed")),
			 			   conditionalPanel(
			 			   	condition = "input.mappingui_ew_influence_type == 'Fixed'",
			 			   	numericInput("mappingui_ew_influence", "Specifiy Edge weight value", 12), br()
			 			   ),
			 			   selectInput("mappingui_overlap_method", "Overlap resolution method", choices = c("repel", "expand")),
			 			   checkboxInput("mappingui_inter_cluster_connections", "Add inter-cluster connections", value = FALSE),
			 			   conditionalPanel(
			 			   	condition = "input.mappingui_inter_cluster_connections == true",
			 			   	selectInput("mappingui_markers_inter_cluster", "Markers for inter-cluster connections (if different)", choices = c(""), multiple = T, width = "100%"),
			 			   	numericInput("mappingui_inter_cluster_weight", "Weight factor for inter-cluster connections", 0.7, min = 0, max = 10, step = 0.1), br()
			 			   )
			 		)
			 	),
			 	fluidRow(
			 		column(12,
			 			   actionButton("saveDataset", "Saving Folder"),
			 			   verbatimTextOutput("saveDatasetFolder"), br(), br()
			 		)),
			 		fluidRow(
			 		  column(12,
			 			   # actionButton("mappingui_start", "Start analysis"), br(), br(),
			 			   # conditionalPanel(
			 			   #   condition <- "!is.null(output.test1)",
			 			   #   br(),
			 			   #   p(strong("Processing data...please wait."))
			 			   # ),
			 			   # conditionalPanel(
			 			   #   condition <- "!is.null(output.mapping_state)",
			 			   #   br(),
			 			   #   strong(textOutput("mapping_state")), br()
			 			   # )),
			 			   # verbatimTextOutput("mappingui_dialog"), br(), br(), br(), br(), br(), br()

			 			   actionButton("mappingui_start", "Start analysis"), br(), br(),
			 			   conditionalPanel(
			 			     condition <- "$('html').hasClass('shiny-busy')",
			 			     br(),
			 			     p(strong("Processing data...please wait."))
			 			   ),
			 			   conditionalPanel(
			 			     condition <- "!$('html').hasClass('shiny-busy') && input.mappingui_start > 0",
			 			     br(),
			 			     p(strong("Data processing is complete!"))
			 			   ),
			 			   verbatimTextOutput("mappingui_dialog"), br(), br(), br(), br(), br(), br() )
			 		)

			 ))


	# tabPanel("Unsupervised map"),
	# tabPanel("Edit SCAFFoLD file")
))
