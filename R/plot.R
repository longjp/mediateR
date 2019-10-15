#########
######### PLOTTING AND DISPLAY FUNCTIONS

## returns subgraph (as path matrix) which has causal influence on y
FindSubgraph <- function(params){
  temp <- rowSums(params$path_model[,params$mm_direct!=0,drop=FALSE]!=0)!=0
  xx_to_keep <- temp | params$xx_direct!=0
  mm_to_keep <- params$mm_direct!=0
  return(params$path_model[xx_to_keep,mm_to_keep,drop=FALSE])
}

## creates GraphNEL object which can then
## be plotted
MakeGraphNELObject <- function(path,xx_direct,mm_direct,mm_path=NULL){
  if (!requireNamespace("graph", quietly = TRUE)) {
    stop("Package graph needed for this function to work. Please install it.",
         call. = FALSE)
  }
  if(is.null(mm_path)){
    mm_path <- matrix(0,nrow=length(mm_direct),ncol=length(mm_direct))
  }
  l <- c("h",rownames(path),colnames(path),"y")
  edL <- vector("list",length=length(l))
  names(edL) <- l
  ## draw arrows from h to xxs
  edL[[1]] <- list(edges=2:(nrow(path)+1),weights=rep(1,nrow(path)))
  ## xx arrows
  for(ii in 2:(nrow(path)+1)){
    ix <- which(path[ii-1,]!=0) + nrow(path)+1
    names(ix) <- NULL
    if(xx_direct[ii-1]!=0){
      ix <- c(ix,sum(dim(path))+2)
    }
    edL[[ii]] <- list(edges=ix,weights=rep(1,length(ix)))
  }
  ## mm arrows
  for(ii in (nrow(path)+2):(nrow(path)+ncol(path)+1)){
    jj <- ii - nrow(path) - 1
    if(mm_direct[jj]!=0){
      edL[[ii]] <- list(edges=nrow(path)+ncol(path)+2,weights=1)
    } else {
      edL[[ii]] <- list(edges=c(),weights=c())
    }
    if(length(which(mm_path[,jj]!=0))!=0){
      edL[[ii]]$edges <- c(edL[[ii]]$edges,which(mm_path[,jj]!=0) + nrow(path) + 1)
      edL[[ii]]$weights <- rep(1,length(edL[[ii]]$edges))
    }
  }
  edL[[length(edL)]] <- list(edges=c(),weights=c())
  gR <- graph::graphNEL(nodes=l, edgeL=edL, edgemode="directed")
  return(gR)
}


## makes edge attributes for graphNEL plotting using path_model (mm, xx coeffs) and direct effects
MakeedgeAttrs <- function(path_model,xx_direct,mm_direct,path=NULL){
  ew <- as.vector(path_model)
  names(ew) <- unlist(lapply(colnames(path_model),function(x){paste0(rownames(path_model),"~",x)}))
  ew <- round(ew,2)
  ef <- c(xx_direct,mm_direct)
  names(ef) <- paste0(c(names(xx_direct),names(mm_direct)),"~y")
  ef <- round(ef,2)
  ew <- c(ew,ef)
  if(!is.null(path)){
    enet <- as.vector(t(path))
    names(enet) <- unlist(lapply(colnames(path_model),function(x){paste0(colnames(path_model),"~",x)}))
    enet <- round(enet,2)
    enet <- enet[enet!=0]
    ew <- c(ew,enet)
  }
  eAttrs <- list()
  eAttrs$label <- ew
  return(eAttrs)
}

## find effects which are non-zero in both true and estimates
FindInterestingEffs <- function(eff_comb){
  direct_both0 <- (eff_comb[,1]==0) & (eff_comb[,1]==eff_comb[,2])
  indirect_both0 <- (eff_comb[,3] == 0) & (eff_comb[,3]==eff_comb[,4])
  indirect_both0[is.na(indirect_both0)] <- TRUE
  to_remove <- direct_both0 & indirect_both0
  return(eff_comb[!to_remove,])
}
