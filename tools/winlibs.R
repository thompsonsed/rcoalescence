# Taken from the rgdal package v 1.4.4
if(getRversion() < "3.3.0") {
  stop("Your version of R is too old. This package requires R-3.3.0 or newer on Windows.")
}

# Download gdal-2.2.0 from rwinlib
VERSION <- commandArgs(TRUE)
if(!file.exists(sprintf("../windows/gdal2-%s/include/gdal/gdal.h", VERSION))){
  if(getRversion() < "3.3.0") setInternet2()
  download.file(sprintf("https://github.com/rwinlib/gdal2/archive/v%s.zip", VERSION), "lib.zip", quiet = TRUE)
  dir.create("../windows", showWarnings = FALSE)
  unzip("lib.zip", exdir = "../windows")
  unlink("lib.zip")
}
