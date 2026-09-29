\name{lisa_init_project}
\alias{lisa_init_project}
\title{Create the installed lisaR Quick Start project}
\usage{
lisa_init_project(path)
}
\arguments{
\item{path}{New directory to create.}
}
\value{
Invisibly, the normalized project directory.
}
\description{
Copies a small, synthetic and fully offline project from the installed
package into a new directory. The destination must not already exist.
}
\examples{
project <- lisa_init_project(tempfile("lisa-quick-start-"))
file.exists(file.path(project, "study.yml"))
}
