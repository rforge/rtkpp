% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/IClusterModel.R, R/global.R,
%   R/ClusterCategorical.R, R/ClusterDiagGaussian.R, R/ClusterGamma.R,
%   R/ClusterMixedData.R, R/ClusterPoisson.R, R/IClusterPredict.R,
%   R/ClusterPredict.R, R/kmm.R, R/kmmMixedData.R
\docType{methods}
\name{summary,IClusterComponent-method}
\alias{summary,IClusterComponent-method}
\alias{summary}
\alias{summary,IClusterModel-method}
\alias{summary,ClusterCategoricalComponent-method}
\alias{summary-ClusterCategoricalComponent,ClusterCategoricalComponent,ClusterCategoricalComponent-method}
\alias{summary,ClusterCategorical-method}
\alias{summary,ClusterDiagGaussian-method}
\alias{summary,ClusterGamma-method}
\alias{summary,ClusterMixedDataModel-method}
\alias{summary,ClusterPoisson-method}
\alias{summary,IClusterPredict-method}
\alias{summary,ClusterPredict-method}
\alias{summary,ClusterPredictMixedData-method}
\alias{summary,KmmModel-method}
\alias{summary,KmmMixedDataModel-method}
\title{Produce summary of a MixAll S4 class.}
\usage{
\S4method{summary}{IClusterComponent}(object, ...)

\S4method{summary}{IClusterModel}(object, ...)

\S4method{summary}{ClusterCategoricalComponent}(object)

\S4method{summary}{ClusterCategorical}(object, ...)

\S4method{summary}{ClusterDiagGaussian}(object, ...)

\S4method{summary}{ClusterGamma}(object, ...)

\S4method{summary}{ClusterMixedDataModel}(object, ...)

\S4method{summary}{ClusterPoisson}(object, ...)

\S4method{summary}{IClusterPredict}(object, ...)

\S4method{summary}{ClusterPredict}(object, ...)

\S4method{summary}{ClusterPredictMixedData}(object, ...)

\S4method{summary}{KmmModel}(object, ...)

\S4method{summary}{KmmMixedDataModel}(object, ...)
}
\arguments{
\item{object}{any cluster model deriving from a \code{\linkS4class{IClusterModel}}
or \code{\linkS4class{IClusterComponent}} object.}

\item{...}{further arguments passed to or from other methods}
}
\value{
NULL. Summaries to standard out.
}
\description{
Produce summary of a MixAll S4 class.
}
