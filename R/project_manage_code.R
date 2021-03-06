#create project directory
#library(Seurat)
#library(dplyr)
#library(yaml)
#library(ggplot2)

#' Create a project directory
#'
#' @param project_title A string
#' @param home_directory A file path
#' @param template_file_path A file path to yml file
#' @return NULL
#' @importFrom yaml write_yaml yaml.load_file
#' @examples
#' CreateProjectDirectory("smith", "/Users/asmith01/projects")
#' CreateProjectDirectory("miller", "~/projects")
#' @export
CreateProjectDirectory <- function(project_title, home_directory, template_file_path) {

  template_file <- yaml.load_file(template_file_path)
  #check template file for project - and add if not there
  
  if (project_title %in% names(template_file$projects)) {
    stop('Project with that name exists.')
  } else {
    template_file$projects[[project_title]] <- paste0(home_directory, "/", project_title)
    write_yaml(template_file, template_file_path)
  }
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
  
  for (directory in main_paths_list) {
    if(!dir.exists(directory)){dir.create(directory,recursive = T)} 
  }
  
  # Parse and create analysis paths
  analysis_paths_list <- list()
  for (path in template_file$analysis_paths){
    analysis_paths_list[[path]] = paste0(project_dir, "/main_analysis/", path)
  }
  new_config[['analysis_paths']] <- analysis_paths_list
  for (directory in analysis_paths_list) {
    if(!dir.exists(directory)){dir.create(directory,recursive = T)} 
  }
  
  #add path to global template file (for colors, settings, etc.)
  new_config[['config_path']] <- template_file_path

  #save as a new YAML file
  config_dir <- new_config[['main_paths']][['config']]
  write_yaml(new_config, file = paste0(config_dir, "/project_config.yml"))
  
  cat("Project directory created.")
}

#' Create a default global template for a home directory
#' This should be run once when GrandeJatte is installed
#' 
#' @param home_directory A file path
#' @importFrom yaml write_yaml
#' @examples
#' CreateGlobalTemplate("/Users/asmith01/projects")
#' CreateGlobalTemplate("~/projects")
#' @export
CreateGlobalTemplate <- function(home_directory){
  if (file.exists(paste0(home_directory, "/grandejatte_template.yml"))) {
    stop('A global template file already exists at:', home_directory, "/grandejatte_template.yml")
  }
  default_template = list(
    base_fields = list('project_title', 'home_directory', 'project_directory', 'samples', 'most_recent_seurat'),
    main_paths = list('main_analysis', 'config', 'input_path', 'scripts', 'seurat_objects'),
    analysis_paths = list('cell_type', 'differential_expression', 'viz_qc', 'viz_clustering',
                          'trajectory', 'markers', 'other'),
    projects = list())
  write_yaml(default_template, file = paste0(home_directory, "/grandejatte_template.yml"))
  cat("Template YML file written to ", home_directory, "/grandejatte_template.yml")
}

#LoadSavePaths <- function() {
#  y_file_path <- "config/project_config.yml"
#  data <- yaml.load_file(y_file_path)
#  return(data)
#}

#' Sets the current project
#'
#' @param project_title The project title - this should be found in the global template file
#' @importFrom yaml yaml.load_file
#' @return NULL
#' @examples
#' SetProject('granger')
SetProject <- function(project_title, home_directory) {
  template_file_path <- paste0(home_directory, "/grandejatte_template.yml")
  if(file.exists(template_file_path)){
    template_file <- yaml.load_file(template_file_path)
  } else {
    stop('Cannot find global template file.')
  }
  if(project_title %in% names(template_file$projects)) {
    setwd(template_file$projects[[project_title]])
  } else {
    stop('This project does not exist in the global template file.')
  }
}

#' Saves a Seurat object to the predefined directory
#'
#' @param seurat_obj A Seurat object
#' @param seurat_filename A filename ending in .rds
#' @param most_recet A boolean
#' @importFrom yaml yaml.load_file
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
#' @importFrom yaml yaml.load_file
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
#' @importFrom yaml yaml.load_file
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
#' @importFrom yaml yaml.load_file
#' @examples
#' SaveVizQC("percent_mito.pdf", width = 10, height = 10)
#' SaveVizQC("bar_plot.png", plot = p1)
#' @export
SaveVizQC <- function(filename, plot = last_plot(), width = 7, height = 7) {
  y_file_path <- "config/project_config.yml"
  config_file <- yaml.load_file(y_file_path)
  save_dir = config_file$analysis_paths$viz_qc
  ggsave(paste0(save_dir, "/", filename), width = width, height = height)
  
}

#' Saves a clustering visualization
#'
#' @param filename A file path with an extension
#' @param plot The plot to save, default is last_plot()
#' @param width A number
#' @param height A number
#' @return NULL
#' @importFrom ggplot2 ggsave last_plot
#' @importFrom yaml yaml.load_file
#' @examples
#' SaveVizCluster("umap_.pdf", width = 10, height = 10)
#' SaveVizCluster("tsne_plot.png", plot = p1)
#' @export
SaveVizCluster <- function(filename, plot = last_plot(), width = 7, height = 7) {
  y_file_path <- "config/project_config.yml"
  config_file <- yaml.load_file(y_file_path)
  save_dir = config_file$analysis_paths$viz_clustering
  ggsave(paste0(save_dir, "/", filename), width = width, height = height)
  
}

#' Saves a differential expression visualization
#'
#' @param filename A file path with an extension
#' @param plot The plot to save, default is last_plot()
#' @param width A number
#' @param height A number
#' @return NULL
#' @importFrom ggplot2 ggsave last_plot
#' @importFrom yaml yaml.load_file
#' @examples
#' SaveVizDE("top_genes.pdf", width = 10, height = 10)
#' SaveVizDE("downregulated_genes.png", plot = p1)
#' @export
SaveVizDE <- function(filename, plot = last_plot(), width = 7, height = 7) {
  y_file_path <- "config/project_config.yml"
  config_file <- yaml.load_file(y_file_path)
  save_dir = config_file$analysis_paths$diff_exp
  ggsave(paste0(save_dir, "/", filename), width = width, height = height)
}





