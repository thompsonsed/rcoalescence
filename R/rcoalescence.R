#' rcoalescence: Efficient spatial and non-spatial neutral ecology simulator written in C++
#' @importFrom Rcpp evalCpp
#' @importFrom Rcpp setRcppClass
#' @importFrom methods new
#' @import methods
#' @import RSQLite
#' @import Rcpp
#' @import dplyr
#' @author Sam Thompson
#' @description A convenience wrapper for necsim, a Neutral Ecology Coalescence SIMulator for
#' spatially-explicit neutral models based in C++.
#'
#' rcoalescence allows for large-scale simulations to be performed on fully spatial files, with
#' features to suit a wide variety of use cases.
#'
#' Four classes exist to represent non-spatial/spatial and non-protracted/protracted simulations.
#' See \link[=NeutralTreeSimulation]{here} for a full description.
#'
#' @section Installation:
#' 
#' 
#' \subsection{Installing on Windows}{
#' 
#' \itemize{
#' \item Make sure that you have installed R > 3.6.0 with  
#' \href{https://cran.r-project.org/bin/windows/Rtools/index.html}{Rtools version 3.5 or later} 
#' (available at \url{https://cran.r-project.org/bin/windows/Rtools/index.html}).
#' \item Install the necessary R package requirements using \code{install.packages(c("devtools", "Rcpp", "RSQLite", "dplyr"))}.
#' \item Run \code{devtools::install\_github("thompsonsed/rcoalescence")}. This will also download gdal (if it can't be found already).
#' }
#' }
#' 
#' 
#' \subsection{Installing on macOS}{
#' 
#' \itemize{
#' 	\item Ensure you have R > 3.4.0 installed.
#' 	\item Install gdal >2.4.0. If you have homebrew installed, this can be done using ``brew install gdal``.
#' 	 Check out the [brew site](https://brew.sh/) for the command to install brew.
#' 	\item Make sure xcode is installed ( run \code{xcode-select --install} from the terminal, which is necessary for compiling the package.
#' 	\item Install the necessary R package requirements using \code{install.packages(c("devtools", "Rcpp", "RSQLite", "dplyr"))}.
#' 	\item Run \code{devtools::install\_github("thompsonsed/rcoalescence")}.
#' }
#' }
#' 
#' \subsection{Installing on Linux}{
#' 
#' \itemize{
#' 	\item Ensure you have R > 3.4.0 installed, plus a compiler that supports C++14 or later (e.g. gcc or clang).
#' 	\item Install gdal >2.4.0 from your distribution's repositories (e.g. \code{sudo apt-get install gdal-dev} for Ubuntu), or build from source.
#'  \item Install the necessary R package requirements using \code{install.packages(c("devtools", "Rcpp", "RSQLite", "dplyr"))}.
#'  \item Run \code{devtools::install\_github("thompsonsed/rcoalescence")}.
#' }
#' }
#'                                                                                                         
#' We also suggest a number of other packages contained within the \emph{tidyverse} collection. 
#' All required packages can be installed using \code{install.packages(c("devtools", "Rcpp", "RSQLite", "tidyverse", "rmarkdown", "knitr"))}.
#'
#' \subsection{Development dependencies}{
#'
#' \itemize{
#'     \item R package development prerequisites, which may require additional installation steps,
#'     depending on your platform
#'     \href{https://support.rstudio.com/hc/en-us/articles/200486498}{(see here)}.
#'     \item A C++ compiler that is detectable by R.
#'     \item The devtools package.
#'     \item Gdal library for C++, available \href{http://www.gdal.org/}{here}. Included in rcoalescence for windows.
#'     \item Sqlite3 library for C++, available \href{https://www.sqlite.org/index.html}{here} (comes with most systems by default).
#'
#'     \item Several R packages are dependencies, including Rcpp. These should be installed
#'     automatically.
#' }
#' }
#'
#'
#' @section Installation issues:
#'
#' \subsection{Cannot find C++ compiler}{
#'
#' For Linux systems, the R development package is required. Use
#' \code{sudo apt-get install r-base-dev} for Ubuntu-based systems. For other distributions, check
#' the documentation for installing core software development utilities required for R package
#' development.
#'
#' On MacOS, this requires that XCode is installed properly, and the compiler is detectable by R.
#' See \href{https://support.rstudio.com/hc/en-us/articles/200486498}{here} for additional info.
#'
#' Windows is currently not supported.
#' }
#'
#' \subsection{Missing gdal requirements}{
#'
#' If you get an error about gdal not being detected, check that
#' gdal is properly installed. If the error indicates missing gdal_priv.h/cpl_conv.h, consider
#' compiling gdal directly from source, as many package management systems do not properly install
#' these files. If you know the directory containing gdal, supply it manually to using PKG_LIBS=/dir/
#' during installation.
#' }
#'
#' \subsection{Missing boost requirements}{
#'
#' Boost should not be essential for package functionality. Instead, you can have a C++17 compiler,
#' or depend on the included file system methodsEnsure that boost has been properly compiled.
#' If you still require boost support, try installing the package BH from CRAN - if this
#' works then you should be fine. Alternatively, if you know the directory containing boost, supply
#' it manually to using PKG_LIBS=/dir/ during installation.
#' }
#'
#' \subsection{Still facing issues}{
#'
#' Try downloading the package from github manually and use \code{install_local()}.
#' If you have autotools installed, run \code{autoreconf} and \code{autoconf} from the command line before
#' attempting package installation.
#' }
#'
#' @section Starting a simulation:
#'
#' Methods are contained within the relevant classes \link[=NeutralTreeSimulation]{(see here)} for
#' setting up, running and analysing simulations. The general process is
#' \itemize{
#'     \item Create a new \link{TreeSimulation}, \link{SpatialTreeSimulation},
#'           \link{ProtractedTreeSimulation} or \link{ProtractedSpatialTreeSimulation} object,
#'           depending on use case.
#'     \item Set the simulation parameters and add additional historical maps.
#'     \item Run the simulation (generally the most time-consuming step).
#'     \item Apply speciation rates and other requirements post simulation.
#'     \item Get the species richness from the simulation.
#'     \item Optionally write results to an SQLite database for further manual analysis.
#'     \item If written to an output database, acquire the desired biodiversity metrics, such as
#'     species locations or species abundances.
#'
#' }
#'
#' See examples in the \link[=NeutralTreeSimulation]{relevant classes}.
#'
#' @docType package
#'
"_PACKAGE"