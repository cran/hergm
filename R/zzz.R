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

.onAttach <- function(lib, pkg)
{
  info <- packageDescription("hergm")
  packageStartupMessage(
  paste('hergm: version ', info$Version, ', created on ', info$Date, '\n',
          "Copyright (c) 2021, Michael Schweinberger, University of Missouri, Columbia\n", 
          "                    with contributions from:\n",
          "                    Mark S. Handcock, University of California - Los Angeles\n",
          "                    Sergii Babkin, Rice University\n",
          "                    Jonathan Stewart, Florida State University\n",
          "                    Duy Vu, University of Melbourne\n",
          "                    Pamela Luna, Baylor College of Medicine\n",
          "                    and selected source files from R package ergm\n",
          "To start hergm: enter help(package=\"hergm\")\n",
          'For license and citation information type citation("hergm").\n',
          '\n',
          'PLEASE NOTE 1: hergm versions >=4 can handle networks with thousands of nodes by using the method = "ml".\n',  
          'PLEASE NOTE 2: hergm versions >=4 use the method = "ml" by default, unless the model includes model terms that cannot be handled by "ml".\n',
          'To reproduce results of hergm versions <4, please use method = "bayes".\n',
          'PLEASE NOTE 3: It is possible that there are issues arising from two recent changes:\n',
          '1. We split hergm into two packages, hergm (unobserved blocks) and mlergm (observed blocks).\n',
          '2. We changed the default method from "bayes" to "ml".\n',
          'If you encounter issues, please contact us, so that we can address them.\n',
           sep="")
 )
}

