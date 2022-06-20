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
CreateProjectDirectory <- function(project_title, home_directory, template_file_path = "none") {
  if (template_file_path == 'none') {
    template_file_path <- 'template.yml'
  }
  
  template_file <- yaml.load_file(template_file_path)
  
  new_config = list()
  
  # Start with base_fields
  new_config[['project_title']] = project_title
  new_config[['home_directory']] = home_directory
  
  if ("project_directory" %in% template_file$base_fields) {
    new_config[['project_directory']] = paste0(home_directory, "/", project_title)
  }
  
  # temporarily set other base fields to undefined. These will be updated later.
  for (field in template_file$base_fields) {
    if (!(field %in% c('project_title', 'home_directory', 'project_directory'))){
      new_config[[field]] = "undefined"
    }
  }
  
  # Parse and create main_paths
  main_paths_list <- list()
  project_dir <- new_config[["project_directory"]]
  for (path in template_file$main_paths){
    main_paths_list[[path]] = paste0(project_dir, "/", path)
  }
  new_config[['main_paths']] <- main_paths_list
  
  CreateDirectories(new_config[['main_paths']])
  
  # Parse and create analysis paths
  analysis_paths_list <- list()
  for (path in template_file$analysis_paths){
    analysis_paths_list[[path]] = paste0(project_dir, "/main_analysis/", path)
  }
  new_config[['analysis_paths']] <- analysis_paths_list
  CreateDirectories(new_config[['analysis_paths']])
  
  #save as a new YAML file
  config_dir <- new_config[['main_paths']][['config']]
  write_yaml(new_config, file = paste0(config_dir, "/", project_title, '_config.yml'))
  
  cat("Project directory created.")
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





