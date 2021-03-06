% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/3-plot_alpha_specific.R
\name{plot_alpha_specific}
\alias{plot_alpha_specific}
\title{Creates a list of alpha diversity plots based on the metrics specified.}
\usage{
plot_alpha_specific(
  phyloseq_object,
  alpha_types = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson",
    "Fisher")
)
}
\arguments{
\item{phyloseq_object}{The phyloseq object containing the relevant dataset.}

\item{alpha_types}{A vector containing the types of alpha diversity metrics to show.}
}
\value{
Returns a list of alpha diversity plots, or raises an error if alpha_types isn't a vector
}
\description{
Creates a list of alpha diversity plots based on the metrics specified.
}
\seealso{
\url{https://www.rdocumentation.org/packages/phyloseq/versions/1.16.2/topics/plot_richness} for the primarily used function.
}
