% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.r
\name{create_meta}
\alias{create_meta}
\title{Generate entire meta file out of information.}
\usage{
create_meta(
  keys,
  db,
  meta_dir,
  project_name,
  source,
  sep = "\\t",
  names_ = NULL
)
}
\arguments{
\item{keys}{list of keys for db. Order determines testing!: First element -> "denom", Second element: "num", rest is ""}

\item{db}{contains key and list of paths (vector of paths) to key files}

\item{meta_dir}{the directory to which meta file should be written to. Name is generated automatically using project_name and source}

\item{project_name}{inserted as signifier into meta filename}

\item{source}{inserted as signifier into meta filename}

\item{sep}{separator in meta file - tab}

\item{names_}{names to be used instead of keys}
}
\value{
list("fpath" = fpath, "core_name" = core_name)
}
\description{
Automatic creation of meta file.
}
