# Motif 1
motif1 <- function(pos){

  res <- c()
  uniq.pos <- unique(c(as.character(unique(pos$Source)), as.character(unique(pos$Target))))

  #necessary internal function
  `%!in%` = Negate(`%in%`)
  for(i in uniq.pos){
    l1 <- pos[pos$Source == i | pos$Target == i,]
    for(j in 1:nrow(l1)){
      if(as.character(l1[j,]$Target) != i){
        k <- as.character(l1[j,]$Target)
      } else {k <- as.character(l1[j,]$Source)}
      l2 <- pos[pos$Source == k | pos$Target == k,]
      pott <- unique(c(as.character(l2$Source),as.character(l2$Target)))
      pott <- pott[pott != k & pott != i]
      if( sum(pott %in% c(as.character(l1$Source), as.character(l1$Target))) != 0 ){
        for(m in 1:sum(pott %in% c(as.character(l1$Source), as.character(l1$Target)))){
          res <- c(res, c(i, k, pott[pott %in% c(as.character(l1$Source), as.character(l1$Target))][m]))
        }
      }
    }
  }
  return(res)
}






# Motif 2
motif2 <- function(neg){

  res <- c()
  uniq.neg <- unique(c(as.character(unique(neg$Source)), as.character(unique(neg$Target))))

  #necessary internal function
  `%!in%` = Negate(`%in%`)
  for(i in uniq.neg){
    l1 <- neg[neg$Source == i | neg$Target == i,]
    for(j in 1:nrow(l1)){
      if(as.character(l1[j,]$Target) != i){
        k <- as.character(l1[j,]$Target)
      } else {k <- as.character(l1[j,]$Source)}
      l2 <- neg[neg$Source == k | neg$Target == k,]
      pott <- unique(c(as.character(l2$Source),as.character(l2$Target)))
      pott <- pott[pott != k & pott != i]
      if( sum(pott %in% c(as.character(l1$Source), as.character(l1$Target))) != 0 ){
        for(m in 1:sum(pott %in% c(as.character(l1$Source), as.character(l1$Target)))){
          res <- c(res, c(i, k, pott[pott %in% c(as.character(l1$Source), as.character(l1$Target))][m]))
        }
      } else{NULL}
    }
  }
  return(res)
}




# Motif 3
motif3 <- function(pos, neg){

  res <- c()
  all <- rbind(neg, pos)
  uniq <- unique(c(as.character(unique(all$Source)), as.character(unique(all$Target))))

  #necessary internal function
  `%!in%` = Negate(`%in%`)
  for(i in uniq){
    l1 <- all[all$Source == i | all$Target == i,]
    for(j in 1:nrow(l1)){
      if(as.character(l1[j,]$Target) != i){
        k <- as.character(l1[j,]$Target)
      } else {k <- as.character(l1[j,]$Source)}

      if( l1[j,]$Type == "NEG" ){
        l2 <- pos[pos$Source == k | pos$Target == k,]
        if( nrow(l2) != 0 ){
          pott <- unique(c(as.character(l2$Source),as.character(l2$Target)))
          pott <- pott[pott != k & pott != i]
          l3 <- pos[pos$Source == i | pos$Target == i,]
          if( sum(pott %in% c(as.character(l3$Source), as.character(l3$Target))) != 0 ){
            for(m in 1:sum(pott %in% c(as.character(l3$Source), as.character(l3$Target)))){
              res <- c(res, c(i, k, pott[pott %in% c(as.character(l3$Source), as.character(l3$Target))][m]))
            }
          } else{NULL}
        } else {NULL}
      } else {NULL}
      if( l1[j,]$Type == "POS" ){
        l2 <- neg[neg$Source == k | neg$Target == k,]
        if( nrow(l2) != 0 ){
          pott <- unique(c(as.character(l2$Source),as.character(l2$Target)))
          pott <- pott[pott != k & pott != i]
          l3 <- pos[pos$Source == i | pos$Target == i,]
          if( sum(pott %in% c(as.character(l3$Source), as.character(l3$Target))) != 0 ){
            for(m in 1:sum(pott %in% c(as.character(l3$Source), as.character(l3$Target)))){
              res <- c(res, c(i, k, pott[pott %in% c(as.character(l3$Source), as.character(l3$Target))][m]))
            }
          } else{NULL}
        } else {NULL}
        l2 <- pos[pos$Source == k | pos$Target == k,]
        if( nrow(l2) != 0 ){
          pott <- unique(c(as.character(l2$Source),as.character(l2$Target)))
          pott <- pott[pott != k & pott != i]
          l3 <- neg[neg$Source == i | neg$Target == i,]
          if( sum(pott %in% c(as.character(l3$Source), as.character(l3$Target))) != 0 ){
            for(m in 1:sum(pott %in% c(as.character(l3$Source), as.character(l3$Target)))){
              res <- c(res, c(i, k, pott[pott %in% c(as.character(l3$Source), as.character(l3$Target))][m]))
            }
          } else{NULL}
        } else {NULL}
      } else {NULL}
    }
  }
  return(res)
}




# Motif 4
motif4 <- function(pos, neg){

  res <- c()
  all <- rbind(neg, pos)
  uniq <- unique(c(as.character(unique(all$Source)), as.character(unique(all$Target))))

  #necessary internal function
  `%!in%` = Negate(`%in%`)
  for(i in uniq){
    l1 <- all[all$Source == i | all$Target == i,]
    for(j in 1:nrow(l1)){
      if(as.character(l1[j,]$Target) != i){
        k <- as.character(l1[j,]$Target)
      } else {k <- as.character(l1[j,]$Source)}

      if( l1[j,]$Type == "NEG" ){
        l2 <- pos[pos$Source == k | pos$Target == k,]
        if( nrow(l2) != 0 ){
          pott <- unique(c(as.character(l2$Source),as.character(l2$Target)))
          pott <- pott[pott != k & pott != i]
          l3 <- neg[neg$Source == i | neg$Target == i,]
          if( sum(pott %in% c(as.character(l3$Source), as.character(l3$Target))) != 0 ){
            for(m in 1:sum(pott %in% c(as.character(l3$Source), as.character(l3$Target)))){
              res <- c(res, c(i, k, pott[pott %in% c(as.character(l3$Source), as.character(l3$Target))][m]))
            }
          } else{NULL}
        } else {NULL}

        l2 <- neg[neg$Source == k | neg$Target == k,]
        if( nrow(l2) != 0 ){
          pott <- unique(c(as.character(l2$Source),as.character(l2$Target)))
          pott <- pott[pott != k & pott != i]
          l3 <- pos[pos$Source == i | pos$Target == i,]
          if( sum(pott %in% c(as.character(l3$Source), as.character(l3$Target))) != 0 ){
            for(m in 1:sum(pott %in% c(as.character(l3$Source), as.character(l3$Target)))){
              res <- c(res, c(i, k, pott[pott %in% c(as.character(l3$Source), as.character(l3$Target))][m]))
            }
          } else{NULL}
        } else {NULL}

      } else {NULL}

      if( l1[j,]$Type == "POS" ){
        l2 <- neg[neg$Source == k | neg$Target == k,]
        if( nrow(l2) != 0 ){
          pott <- unique(c(as.character(l2$Source),as.character(l2$Target)))
          pott <- pott[pott != k & pott != i]
          l3 <- neg[neg$Source == i | neg$Target == i,]
          if( sum(pott %in% c(as.character(l3$Source), as.character(l3$Target))) != 0 ){
            for(m in 1:sum(pott %in% c(as.character(l3$Source), as.character(l3$Target)))){
              res <- c(res, c(i, k, pott[pott %in% c(as.character(l3$Source), as.character(l3$Target))][m]))
            }
          } else{NULL}
        } else {NULL}

      } else {NULL}
    }
  }
  return(res)
}






# Motif 5
motif5 <- function(pos, neg){

  res <- c()
  all <- rbind(neg, pos)
  uniq.pos <- unique(c(as.character(unique(pos$Source)), as.character(unique(pos$Target))))

  #necessary internal function
  `%!in%` = Negate(`%in%`)
  for(i in uniq.pos){
    l1 <- pos[pos$Source == i | pos$Target == i,]
    for(j in 1:nrow(l1)){
      if(as.character(l1[j,]$Target) != i){
        k <- as.character(l1[j,]$Target)
      } else {k <- as.character(l1[j,]$Source)}
      l2 <- pos[pos$Source == k | pos$Target == k,]
      pott <- unique(c(as.character(l2$Source),as.character(l2$Target)))
      pott <- pott[pott != k & pott != i]
      l3 <- all[all$Source == i | all$Target == i,]
      if( sum(pott %!in% c(as.character(l3$Source), as.character(l3$Target))) != 0 ){
        for(m in 1:sum(pott %!in% c(as.character(l3$Source), as.character(l3$Target)))){
          res <- c(res, c(i, k, pott[pott %!in% c(as.character(l3$Source), as.character(l3$Target))][m]))
        }
      } else{NULL}
    }
  }
  return(res)
}


# Motif 5 including the ones present on other triangles
motif5.m <- function(pos){

  res <- c()
  uniq.pos <- unique(c(as.character(unique(pos$Source)), as.character(unique(pos$Target))))

  #necessary internal function
  `%!in%` = Negate(`%in%`)
  for(i in uniq.pos){
    l1 <- pos[pos$Source == i | pos$Target == i,]
    for(j in 1:nrow(l1)){
      if(as.character(l1[j,]$Target) != i){
        k <- as.character(l1[j,]$Target)
      } else {k <- as.character(l1[j,]$Source)}
      l2 <- pos[pos$Source == k | pos$Target == k,]
      pott <- unique(c(as.character(l2$Source),as.character(l2$Target)))
      pott <- pott[pott != k & pott != i]
      if( length(pott) != 0 ){
        for(m in 1:length(pott)){
          res <- c(res, c(i, k, pott[m]))
        }
      } else{NULL}
    }
  }
  return(res)
}





# Motif 6
motif6 <- function(pos, neg){

  res <- c()
  all <- rbind(neg, pos)
  uniq.neg <- unique(c(as.character(unique(neg$Source)), as.character(unique(neg$Target))))

  #necessary internal function
  `%!in%` = Negate(`%in%`)
  for(i in uniq.neg){
    l1 <- neg[neg$Source == i | neg$Target == i,]
    for(j in 1:nrow(l1)){
      if(as.character(l1[j,]$Target) != i){
        k <- as.character(l1[j,]$Target)
      } else {k <- as.character(l1[j,]$Source)}
      l2 <- neg[neg$Source == k | neg$Target == k,]
      pott <- unique(c(as.character(l2$Source),as.character(l2$Target)))
      pott <- pott[pott != k & pott != i]
      l3 <- all[all$Source == i | all$Target == i,]
      if( sum(pott %!in% c(as.character(l3$Source), as.character(l3$Target))) != 0 ){
        for(m in 1:sum(pott %!in% c(as.character(l3$Source), as.character(l3$Target)))){
          res <- c(res, c(i, k, pott[pott %!in% c(as.character(l3$Source), as.character(l3$Target))][m]))
        }
      } else{NULL}
    }
  }
  return(res)
}

# Motif 6 including the ones present on other triangles
motif6.m <- function(neg){

  res <- c()
  uniq.neg <- unique(c(as.character(unique(neg$Source)), as.character(unique(neg$Target))))

  #necessary internal function
  `%!in%` = Negate(`%in%`)
  for(i in uniq.neg){
    l1 <- neg[neg$Source == i | neg$Target == i,]
    for(j in 1:nrow(l1)){
      if(as.character(l1[j,]$Target) != i){
        k <- as.character(l1[j,]$Target)
      } else {k <- as.character(l1[j,]$Source)}
      l2 <- neg[neg$Source == k | neg$Target == k,]
      pott <- unique(c(as.character(l2$Source),as.character(l2$Target)))
      pott <- pott[pott != k & pott != i]
      if( length(pott) != 0 ){
        for(m in 1:length(pott)){
          res <- c(res, c(i, k, pott[m]))
        }
      } else{NULL}
    }
  }
  return(res)
}


# Motif 7
motif7 <- function(pos, neg){

  res <- c()
  all <- rbind(neg, pos)
  uniq <- unique(c(as.character(unique(all$Source)), as.character(unique(all$Target))))

  #necessary internal function
  `%!in%` = Negate(`%in%`)
  for(i in uniq){
    l1 <- all[all$Source == i | all$Target == i,]
    for(j in 1:nrow(l1)){
      if(as.character(l1[j,]$Target) != i){
        k <- as.character(l1[j,]$Target)
      } else {k <- as.character(l1[j,]$Source)}

      if( l1[j,]$Type == "NEG" ){
        l2 <- pos[pos$Source == k | pos$Target == k,]
        if( nrow(l2) != 0 ){
          pott <- unique(c(as.character(l2$Source),as.character(l2$Target)))
          pott <- pott[pott != k & pott != i]
          l3 <- all[all$Source == i | all$Target == i,]
          if( sum(pott %!in% c(as.character(l3$Source), as.character(l3$Target))) != 0 ){
            for(m in 1:sum(pott %!in% c(as.character(l3$Source), as.character(l3$Target)))){
              res <- c(res, c(i, k, pott[pott %!in% c(as.character(l3$Source), as.character(l3$Target))][m]))
            }
          } else{NULL}
        } else {NULL}

      } else {NULL}

      if( l1[j,]$Type == "POS" ){
        l2 <- neg[neg$Source == k | neg$Target == k,]
        if( nrow(l2) != 0 ){
          pott <- unique(c(as.character(l2$Source),as.character(l2$Target)))
          pott <- pott[pott != k & pott != i]
          l3 <- all[all$Source == i | all$Target == i,]
          if( sum(pott %!in% c(as.character(l3$Source), as.character(l3$Target))) != 0 ){
            for(m in 1:sum(pott %!in% c(as.character(l3$Source), as.character(l3$Target)))){
              res <- c(res, c(i, k, pott[pott %!in% c(as.character(l3$Source), as.character(l3$Target))][m]))
            }
          } else{NULL}
        } else {NULL}

      } else {NULL}
    }
  }
  return(res)
}



# Motif 7 including the ones present on other triangles
motif7.m <- function(pos, neg){

  res <- c()
  all <- rbind(neg, pos)
  uniq <- unique(c(as.character(unique(all$Source)), as.character(unique(all$Target))))

  #necessary internal function
  `%!in%` = Negate(`%in%`)
  for(i in uniq){
    l1 <- all[all$Source == i | all$Target == i,]
    for(j in 1:nrow(l1)){
      if(as.character(l1[j,]$Target) != i){
        k <- as.character(l1[j,]$Target)
      } else {k <- as.character(l1[j,]$Source)}

      if( l1[j,]$Type == "NEG" ){
        l2 <- pos[pos$Source == k | pos$Target == k,]
        if( nrow(l2) != 0 ){
          pott <- unique(c(as.character(l2$Source),as.character(l2$Target)))
          pott <- pott[pott != k & pott != i]
          if( length(pott) != 0 ){
            for(m in 1:length(pott)){
              res <- c(res, c(i, k, pott[m]))
            }
          } else{NULL}
        } else {NULL}

      } else {NULL}

      if( l1[j,]$Type == "POS" ){
        l2 <- neg[neg$Source == k | neg$Target == k,]
        if( nrow(l2) != 0 ){
          pott <- unique(c(as.character(l2$Source),as.character(l2$Target)))
          pott <- pott[pott != k & pott != i]
          if( length(pott) != 0 ){
            for(m in 1:length(pott)){
              res <- c(res, c(i, k, pott[m]))
            }
          } else{NULL}
        } else {NULL}

      } else {NULL}
    }
  }
  return(res)
}


### 4 nodes motif ###

# Motif 8 - all positive
motif8 <- function(pos){

  res <- c()
  uniq.pos <- unique(c(as.character(unique(pos$Source)), as.character(unique(pos$Target))))

  #necessary internal function
  `%!in%` = Negate(`%in%`)
  for(i in uniq.pos){
    l1 <- pos[pos$Source == i | pos$Target == i,]
    for(j in 1:nrow(l1)){
      if(as.character(l1[j,]$Target) != i){
        k <- as.character(l1[j,]$Target)
      } else {k <- as.character(l1[j,]$Source)}
      l2 <- pos[pos$Source == k | pos$Target == k,]
      if(nrow(l2) != 0){
        for(m in 1:nrow(l2)){
          if(as.character(l2[m,]$Target) != k){
            y <- as.character(l2[m,]$Target)
          } else {y <- as.character(l2[m,]$Source)}
          l3 <- pos[pos$Source == y | pos$Target == y,]
          pott <- unique(c(as.character(l3$Source),as.character(l3$Target)))
          pott <- pott[pott != k & pott != i & pott != y]
        }
        if( sum(pott %in% c(as.character(l1$Source), as.character(l1$Target))) != 0 ){
          for(m in 1:sum(pott %in% c(as.character(l1$Source), as.character(l1$Target)))){
            res <- c(res, c(i, k, y, pott[pott %in% c(as.character(l1$Source), as.character(l1$Target))][m]))
          }
        }
      }

    }
  }
  return(res)
}


# Motif 9 - one negative
motif9 <- function(pos, neg){

  res <- c()
  uniq.neg <- unique(c(as.character(unique(neg$Source)), as.character(unique(neg$Target))))

  #necessary internal function
  `%!in%` = Negate(`%in%`)
  for(i in uniq.neg){
    l1 <- neg[neg$Source == i | neg$Target == i,]
    for(j in 1:nrow(l1)){
      if(as.character(l1[j,]$Target) != i){
        k <- as.character(l1[j,]$Target)
      } else {k <- as.character(l1[j,]$Source)}
      l2 <- pos[pos$Source == k | pos$Target == k,]
      if(nrow(l2) != 0){
        for(m in 1:nrow(l2)){
          if(as.character(l2[m,]$Target) != k){
            y <- as.character(l2[m,]$Target)
          } else {y <- as.character(l2[m,]$Source)}
          l3 <- pos[pos$Source == y | pos$Target == y,]
          pott <- unique(c(as.character(l3$Source),as.character(l3$Target)))
          pott <- pott[pott != k & pott != i & pott != y]
          l4 <- pos[pos$Source == i | pos$Target == i,]
        }
        if( sum(pott %in% c(as.character(l4$Source), as.character(l4$Target))) != 0 ){
          for(m in 1:sum(pott %in% c(as.character(l4$Source), as.character(l4$Target)))){
            res <- c(res, c(i, k, y, pott[pott %in% c(as.character(l4$Source), as.character(l4$Target))][m]))
          }
        }
      }
    }
  }
  return(res)
}


# Motif 10 - half / half
motif10 <- function(pos, neg){

  res <- c()
  uniq.neg <- unique(c(as.character(unique(neg$Source)), as.character(unique(neg$Target))))

  #necessary internal function
  `%!in%` = Negate(`%in%`)
  for(i in uniq.neg){
    l1 <- neg[neg$Source == i | neg$Target == i,]
    for(j in 1:nrow(l1)){
      if(as.character(l1[j,]$Target) != i){
        k <- as.character(l1[j,]$Target)
      } else {k <- as.character(l1[j,]$Source)}
      l2 <- neg[neg$Source == k | neg$Target == k,]
      if(nrow(l2) != 0){
        for(m in 1:nrow(l2)){
          if(as.character(l2[m,]$Target) != k){
            y <- as.character(l2[m,]$Target)
          } else {y <- as.character(l2[m,]$Source)}
          l3 <- pos[pos$Source == y | pos$Target == y,]
          pott <- unique(c(as.character(l3$Source),as.character(l3$Target)))
          pott <- pott[pott != k & pott != i & pott != y]
          l4 <- pos[pos$Source == i | pos$Target == i,]
        }
        if( sum(pott %in% c(as.character(l4$Source), as.character(l4$Target))) != 0 ){
          for(m in 1:sum(pott %in% c(as.character(l4$Source), as.character(l4$Target)))){
            res <- c(res, c(i, k, y, pott[pott %in% c(as.character(l4$Source), as.character(l4$Target))][m]))
          }
        }
      }
    }
  }
  return(res)
}


# Motif 11 - one positive
motif11 <- function(pos, neg){

  res <- c()
  uniq.pos <- unique(c(as.character(unique(pos$Source)), as.character(unique(pos$Target))))

  #necessary internal function
  `%!in%` = Negate(`%in%`)
  for(i in uniq.pos){
    l1 <- pos[pos$Source == i | pos$Target == i,]
    for(j in 1:nrow(l1)){
      if(as.character(l1[j,]$Target) != i){
        k <- as.character(l1[j,]$Target)
      } else {k <- as.character(l1[j,]$Source)}
      l2 <- neg[neg$Source == k | neg$Target == k,]
      if(nrow(l2) != 0){
        for(m in 1:nrow(l2)){
          if(as.character(l2[m,]$Target) != k){
            y <- as.character(l2[m,]$Target)
          } else {y <- as.character(l2[m,]$Source)}
          l3 <- neg[neg$Source == y | neg$Target == y,]
          pott <- unique(c(as.character(l3$Source),as.character(l3$Target)))
          pott <- pott[pott != k & pott != i & pott != y]
          l4 <- neg[neg$Source == i | neg$Target == i,]
        }
        if( sum(pott %in% c(as.character(l4$Source), as.character(l4$Target))) != 0 ){
          for(m in 1:sum(pott %in% c(as.character(l4$Source), as.character(l4$Target)))){
            res <- c(res, c(i, k, y, pott[pott %in% c(as.character(l4$Source), as.character(l4$Target))][m]))
          }
        }
      }
    }
  }
  return(res)
}



# Motif 12 - all negative
motif12 <- function(neg, uniq.neg){

  res <- c()
  uniq.neg <- unique(c(as.character(unique(neg$Source)), as.character(unique(neg$Target))))

  #necessary internal function
  `%!in%` = Negate(`%in%`)
  for(i in uniq.neg){
    l1 <- neg[neg$Source == i | neg$Target == i,]
    for(j in 1:nrow(l1)){
      if(as.character(l1[j,]$Target) != i){
        k <- as.character(l1[j,]$Target)
      } else {k <- as.character(l1[j,]$Source)}
      l2 <- neg[neg$Source == k | neg$Target == k,]
      if(nrow(l2) != 0){
        for(m in 1:nrow(l2)){
          if(as.character(l2[m,]$Target) != k){
            y <- as.character(l2[m,]$Target)
          } else {y <- as.character(l2[m,]$Source)}
          l3 <- neg[neg$Source == y | neg$Target == y,]
          pott <- unique(c(as.character(l3$Source),as.character(l3$Target)))
          pott <- pott[pott != k & pott != i & pott != y]
        }
        if( sum(pott %in% c(as.character(l1$Source), as.character(l1$Target))) != 0 ){
          for(m in 1:sum(pott %in% c(as.character(l1$Source), as.character(l1$Target)))){
            res <- c(res, c(i, k, y, pott[pott %in% c(as.character(l1$Source), as.character(l1$Target))][m]))
          }
        }
      }
    }
  }
  return(res)
}


