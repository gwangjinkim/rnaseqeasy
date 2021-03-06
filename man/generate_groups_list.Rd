% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.r
\name{generate_groups_list}
\alias{generate_groups_list}
\title{Group list of paths by change of common substrings automatically.}
\usage{
generate_groups_list(paths, split_ending = c("_", "-"))
}
\arguments{
\item{paths}{the file paths of the symbol counts of each sample}

\item{split_ending}{by which substring the path string should be split and grouped.}
}
\value{
list of grouped file paths named by file name part
}
\description{
Group list of paths.
}
