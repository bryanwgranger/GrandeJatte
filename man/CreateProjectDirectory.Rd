% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/project_manage_code.R
\name{CreateProjectDirectory}
\alias{CreateProjectDirectory}
\title{Create a project directory}
\usage{
CreateProjectDirectory(project_title, home_directory, template_file_path)
}
\arguments{
\item{project_title}{A string}

\item{home_directory}{A file path}

\item{template_file_path}{A file path to yml file}
}
\description{
Create a project directory
}
\examples{
CreateProjectDirectory("smith", "/Users/asmith01/projects")
CreateProjectDirectory("miller", "~/projects")
}
