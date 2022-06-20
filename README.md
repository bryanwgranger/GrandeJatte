# Grande Jatte
## A project manager for scRNAseq with Seurat

This R package is a small organizational tool to be used in conjuction with Seurat when conducting scRNAseq analysis. The package is currently in development and has limited features. It uses a YAML file to set and store configuration options, assisting with creating project directories, saving visualizations, and loading/saving Seurat objects.

## Installation

GrandeJatte currently requires the following packages:

```
ggplot2
yaml
```

To install GrandeJatte, please do so via devtools:


```
devtools::install_github("bryanwgranger/GrandeJatte")
```

## Usage

### First use

When using GrandeJatte for the first time, run the function ```CreateTemplate(home_directory)``` with ```home_directory``` as a directory that will hold all the project directories. This will intialize a default YAML file with options for project directory creation. The file will be saved as ```grandejatte_template.yml``` and will be used by the package in creating and maintaining project directories

### Creating Project Directories

To create a project directory, use the function ```CreateProjectDirectory()``` as such:
```
CreateProjectDirectory(project_title = "sample_project",  # title of the project
						home_directory = "/sample/path",  # home directory that holds all 
														  # projects
						template_file_path = "/sample/path/grandejatte_template.yml")
							# location of the template YAML file 
```

This will create a project directory using information from the template YAML file, and it will create a separate project specific YAML file called ```project_config.yml```. 

### Saving visualizations

Currently, you can save types of vizualizations to specific directories outlined in the ```project_config.yml``` file. Examples:

```
SaveVizQC(filename = "qc_plot.pdf")
SaveVizClustering(filename = "umap_dimplot.pdf")
SaveVizDE(filename = "diff_exp_plot.pdf", width = 12, height = 8)

# A plot parameter can be specified, though the default plot 
# is the last plot via ggplots last_plot():

p1 <- DimPlot(seurat_object)
SaveVizClusterting(filename = "umap_dimplot.pdf", plot = p1, width = 10, height = 10)
```
