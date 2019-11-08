\name{canidae}
\alias{canidae}
\docType{data}
\title{canidae
%%   ~~ data name/kind ... ~~
}
\description{ Canidae dataset that can be used as example.
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("canidae")}
\format{
  a set of objects to be used in examples with package mvEvol.
  \describe{
    \item{\code{W}}{Pooled within variance covariance matrix. Can be used as the G matrix in examples.}
    \item{\code{means}}{Matrix containing mean trait values for species. Total of 39 species (lines) and 40 traits (columns).}
    \item{\code{n.s}}{Table containing the sample size for each species. 39 values.}
    \item{\code{tee}}{Dated phylogenetic tree representing the relationship among taxa and node age.}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(canidae)
## maybe str(canidae) ; plot(canidae) ...
}
\keyword{datasets}
