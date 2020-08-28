
mod_results <- function(results, nodes){
  if(is.vector(results)){
    a <- as.data.frame(split(results, ceiling(seq_along(results)/nodes)))
    a.m <- t(a)
    a.df <- as.data.frame(t(a))
    rep <- duplicated(
      lapply(1:nrow(a.m), function(y){
        A <- a.m[y, ]
        A[order(A)]
      }))
    a.df.f <- a.df[!rep,]
    return(nrow(a.df.f))
  } else { return(0) }
}



extractCount <- function(results, nodes, uniq){
  if(is.vector(results)){
    a <- as.data.frame(split(results, ceiling(seq_along(results)/nodes)))
    a.m <- t(a)
    a.df <- as.data.frame(t(a))
    rep <- duplicated(
      lapply(1:nrow(a.m), function(y){
        A <- a.m[y, ]
        A[order(A)]
      }))
    a.df.f <- a.df[!rep,]

    count.t <- c()
    for(x in uniq){
      suma <- sum(sum(a.df.f$V1 == x), sum(a.df.f$V2 == x), sum(a.df.f$V3 == x))
      count.t <- c(count.t, suma)
    }
    return(count.t)

  } else {
    all_cero <- rep(0, length(uniq))
    return(all_cero)
  }
}
