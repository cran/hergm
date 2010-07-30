##########################################################################
#
# Copyright (c) 2009, Michael Schweinberger, Pennsylvania State University
# 
# GPL-3 license
#
##########################################################################

.First.lib <- function(lib, pkg)
  {
  library.dynam("hergm", pkg, lib)
  DESCpath <- file.path(system.file(package="hergm"), "DESCRIPTION")
  info <- read.dcf(DESCpath)
  cat('\nhergm:', info[,"Title"], 
      '\nVersion', info[,"Version"], 'created on', info[,"Date"], '\n')   
  cat(paste("Copyright (c) 2009, Michael Schweinberger, Pennsylvania State University\n", sep = ""))
  cat('To start hergm: enter help(package=\"hergm\")\n')
  cat('To cite hergm: enter citation(\"hergm\")\n')
  }

.Last.lib <- function(libpath)
  {
  library.dynam.unload("hergm",libpath)
  }
