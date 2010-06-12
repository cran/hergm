###########################################################################
# Copyright 2009 Michael Schweinberger                                    #
#                                                                         #
# This file is part of hergm.                                             #
#                                                                         # 
#    hergm is free software: you can redistribute it and/or modify        #
#    it under the terms of the GNU General Public License as published by #
#    the Free Software Foundation, either version 3 of the License, or    #
#    (at your option) any later version.                                  #
#                                                                         # 
#    hergm is distributed in the hope that it will be useful,             #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of       #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        #
#    GNU General Public License for more details.                         #
#                                                                         #
#    You should have received a copy of the GNU General Public License    #
#    along with hergm.  If not, see <http://www.gnu.org/licenses/>.       #
#                                                                         # 
###########################################################################

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
