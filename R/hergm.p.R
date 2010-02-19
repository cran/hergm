hergm.p <- function(hergm_list) 
{
  p <- .C("P",
  as.integer(hergm_list$Clist$nterms),
  as.integer(hergm_list$d),
  as.double(hergm_list$Clist$inputs), 
  as.double(hergm_list$theta),
  as.integer(hergm_list$Clist$n),
  as.integer(0),
  as.integer(0),
  as.character(hergm_list$Clist$fnamestring),
  as.character(hergm_list$Clist$snamestring),
  p = as.double(choose(hergm_list$Clist$n,2)),
  PACKAGE="hergm")
  p$p
}


