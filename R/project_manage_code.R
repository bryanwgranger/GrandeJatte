#create project directory
#library(Seurat)
#library(dplyr)
#library(yaml)
#library(ggplot2)

#' Create a project directory
#'
#' @param project_title A string
#' @param project_title A file path
#' @return NULL
#' @examples
#' CreateProjectDirectory("smith", "/Users/asmith01/projects")
#' CreateProjectDirectory("miller", "~/projects")
#' @export
CreateProjectDirectory <- function(project_title, home_directory){

  main_paths <- list(input_path = paste0(home_directory, "/", project_title, "/input_data/"), 
                     input_path_counts = paste0(home_directory, "/", project_title, "/input_data/count"),
                     main_analysis = paste0(home_directory, "/", project_title, "/main_analysis"), 
                     seurat_objects = paste0(home_directory, "/", project_title, "/seurat_objects"), 
                     scripts = paste0(home_directory, "/", project_title, "/scripts"),
                     config = paste0(home_directory, "/", project_title, "/config"))
  
  for (directory in main_paths) {
    if(!dir.exists(directory)){dir.create(directory,recursive = T)} 
  }

  analysis_paths <- list(cell_type = paste0(main_paths$main_analysis, "/cell_type"), 
                              diff_exp = paste0(main_paths$main_analysis, "/differential_expression"), 
                              viz_qc = paste0(main_paths$main_analysis, "/viz_qc"), 
                              viz_clustering = paste0(main_paths$main_analysis, "/viz_clustering"), 
                              trajectory = paste0(main_paths$main_analysis, "/trajectory"),
                              markers = paste0(main_paths$main_analysis, "/markers"),
                              rna_velocity = paste0(main_paths$main_analysis, "/rna_velocity"))
  
  for (directory in analysis_paths) {
    if(!dir.exists(directory)){dir.create(directory,recursive = T)} 
  }
  
  #CREATE YAML FILE
  yml <- list(project_title = project_title,
              home_directory = home_directory,
              project_directory = paste0(home_directory, "/", project_title),
              main_paths = main_paths,
              analysis_paths = analysis_paths)
  write_yaml(yml, file = paste0(main_paths$config, '/project_config.yml'))
  cat('Project YAML file created at', main_paths$config)
}

#LoadSavePaths <- function() {
#  y_file_path <- "config/project_config.yml"
#  data <- yaml.load_file(y_file_path)
#  return(data)
#}

#' Saves a Seurat object to the predefined directory
#'
#' @param seurat_obj A Seurat object
#' @param seurat_filename A filename ending in .rds
#' @param most_recet A boolean
#' @return NULL
#' @examples
#' SaveSeurat(seurat, "05_smith_diff_exp.rds", most_recent = TRUE)
#' SaveSeurat(DATA, "DATA_unfiltered.rds", most_recent = FALSE)
SaveSeurat <- function(seurat_obj, seurat_filename, most_recent = TRUE) {
  y_file_path <- "config/project_config.yml"
  config_file <- yaml.load_file(y_file_path)
  seurat_dir <- config_file[['save_paths']][['seurat_objects']]
  save_path <- paste0(seurat_dir, "/", seurat_filename)
  if (most.recent == TRUE){
    config_file$most_recent_seurat_save = save_path
  }
  saveRDS(seurat_obj, paste0(paths[['save_paths']][['seurat_objects']], "/", seurat_filename))
}

#' Opens the most recently saved Seurat using the SaveSeurat function
#'
#' @return A Seurat object
#' @examples
#' OpenLastSeurat()
#' @export
OpenLastSeurat <- function() {
  y_file_path <- "config/project_config.yml"
  config_file <- yaml.load_file(y_file_path)
  if ("most_recent_seurat_save" %in% names(config_file)) {
    return(readRDS(config_file$most_recent_seurat_save))
  } else {
    stop("Cannot find last saved Seurat Object in config file.\nMake sure to SaveSeurat with condition most.recent = TRUE")
  }
  
}

###########
##VIZ
###########

#save viz general

#' Saves a visualization
#'
#' @param category A string
#' @param filename A file path with an extension
#' @param plot The plot to save, default is last_plot()
#' @param width A number
#' @param height A number
#' @return NULL
#' @importFrom ggplot2 ggsave last_plot
#' @examples
#' SaveViz("qc", "percent_mito.pdf", width = 10, height = 10)
#' SaveViz("other", "bar_plot.png", plot = p1)
#' @export
SaveViz <- function(category, filename, plot = last_plot(), width = 7, height = 7) {
  y_file_path <- "config/project_config.yml"
  config_file <- yaml.load_file(y_file_path)
  
  if (category == "qc") {
    save_dir = config_file$analysis_paths$viz_qc
    
  } else if (category == "cluster"){
    save_dir = config_file$analysis_paths$viz_clustering
    
  } else if (category == "de"){
    save_dir = config_file$analysis_paths$diff_exp
    
  } else if (category == "other"){
    save_dir = config_file$analysis_paths
    
  } else {
    stop('Save category not recognized. Please use one of [qc, cluster, de, or other]')
  }
  ggsave(paste0(save_dir, "/", filename), plot = plot, width = width, height = height)
}

#' Saves a quality control visualization
#'
#' @param filename A file path with an extension
#' @param plot The plot to save, default is last_plot()
#' @param width A number
#' @param height A number
#' @return NULL
#' @importFrom ggplot2 ggsave last_plot
#' @examples
#' SaveVizQC("percent_mito.pdf", width = 10, height = 10)
#' SaveVizQC("bar_plot.png", plot = p1)
#' @export
SaveVizQC <- function(filename, plot = last_plot(), width = 7, height = 7) {
  y_file_path <- "config/project_config.yml"
  config_file <- yaml.load_file(y_file_path)
  save_dir = config_file$analysis_paths$viz_qc
  ggsave(paste0(save_dir, "/", filename, width = width, height = height))
  
}

#' Saves a clustering visualization
#'
#' @param filename A file path with an extension
#' @param plot The plot to save, default is last_plot()
#' @param width A number
#' @param height A number
#' @return NULL
#' @importFrom ggplot2 ggsave last_plot
#' @examples
#' SaveVizCluster("umap_.pdf", width = 10, height = 10)
#' SaveVizCluster("tsne_plot.png", plot = p1)
#' @export
SaveVizCluster <- function(filename, plot = last_plot(), width = 7, height = 7) {
  y_file_path <- "config/project_config.yml"
  config_file <- yaml.load_file(y_file_path)
  save_dir = config_file$analysis_paths$viz_clustering
  ggsave(paste0(save_dir, "/", filename, width = width, height = height))
  
}

#' Saves a differential expression visualization
#'
#' @param filename A file path with an extension
#' @param plot The plot to save, default is last_plot()
#' @param width A number
#' @param height A number
#' @return NULL
#' @importFrom ggplot2 ggsave last_plot
#' @examples
#' SaveVizDE("top_genes.pdf", width = 10, height = 10)
#' SaveVizDE("downregulated_genes.png", plot = p1)
#' @export
SaveVizDE <- function(filename, plot = last_plot(), width = 7, height = 7) {
  y_file_path <- "config/project_config.yml"
  config_file <- yaml.load_file(y_file_path)
  save_dir = config_file$analysis_paths$diff_exp
  ggsave(paste0(save_dir, "/", filename, width = width, height = height))
}





