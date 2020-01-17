
\docType{data}
\name{camr}
\alias{camr}
\title{(Data) child mortality in China, cummulative by age and strata of MCHSS.}
\format{A data frame with 720 rows and 6 variables:
\describe{
  \item{agegp}{Age group of children in MCHSS}
  \item{reg}{Region of China (east, middle, or west)}
  \item{res}{rural or urban}
  \item{year}{Year of surveillance}
  \item{qc_nCH9}{Cummulative mortaliy (from 0 to age upper limit) due to pneumonia (per live birth)}
  \item{qc_nCH17}{Cummulative mortaliy (from 0 to age upper limit) due to injury (per live birth)}
  \item{qc_ntotal}{Cummulative mortaliy (from 0 to age upper limit) due to all causes (per live birth)}
    ...
}}
\source{
\url{https://www.sciencedirect.com/science/article/pii/S2214109X16303345}
}
\usage{
camr
}
