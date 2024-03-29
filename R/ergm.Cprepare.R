#  File ergm/R/ergm.Cprepare.R
#  Part of the statnet package, http://statnetproject.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnetproject.org/attribution
#
#  Copyright 2010 the statnet development team
######################################################################
ergm.Cprepare <- function(nw, m) 
{
  # Build an object called Clist that contains all the necessary
  # ingredients to be passed to the C function.
  n <- network.size(nw)
  dir <- is.directed(nw)
  Clist<-list(n=n, dir=dir)
  bip <- nw$gal$bipartite
  if (is.null(bip)) bip <- 0
  Clist$bipartite <- bip
  Clist$ndyads <- n * (n-1) / (2-dir)
  e<-as.matrix.network(nw,matrix.type="edgelist")
  Clist$maxpossibleedges <- min(max(1e+6, 2*nrow(e)), Clist$ndyads)
  if(length(e)==0){
    Clist$nedges<-0
    Clist$heads<-NULL
    Clist$tails<-NULL
  }else{
    if(!is.matrix(e)){e <- matrix(e, ncol=2)}
    Clist$nedges<-dim(e)[1]
    # Ensure that for undirected networks, head<tail.
    if(dir){
      Clist$heads<-e[,1]
      Clist$tails<-e[,2]
    }else{
      Clist$heads<-pmin(e[,1],e[,2])
      Clist$tails<-pmax(e[,1],e[,2])
    }
  }
  mo<-m$terms 
  
  Clist$nterms<-length(mo)
  Clist$nstats<-0
  Clist$fnamestring<-""
  Clist$snamestring<-""
  Clist$inputs<-numeric(0)
  if (Clist$nterms>0) {
    for(i in 1:Clist$nterms) {
      Clist$fnamestring <- paste(Clist$fnamestring, mo[[i]]$name)
      Clist$snamestring <- paste(Clist$snamestring, 
                                 ifelse(is.null(mo[[i]]$soname), "ergm",
                                        mo[[i]]$soname))
      Clist$inputs <- c(Clist$inputs, mo[[i]]$inputs)
      Clist$nstats <- Clist$nstats + mo[[i]]$inputs[2]
    }
  }
  while (substring(Clist$fnamestring, 1, 1)==" ")
    Clist$fnamestring <- substring(Clist$fnamestring, 2)
  while (substring(Clist$snamestring, 1, 1)==" ")
    Clist$snamestring <- substring(Clist$snamestring, 2)

  if("order" %in% names(m)) Clist$order.code <- switch(m$order,
                                                       DissThenForm=1,
                                                       DissAndForm=2,
                                                       FormAndDiss=2,
                                                       FormThenDiss=3,
                                                       FormOnly=4,
                                                       DissOnly=5,
                                                       0)
  
  Clist
}


