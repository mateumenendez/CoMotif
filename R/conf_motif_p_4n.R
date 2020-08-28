#' Random networks confidence interval with paralelization
#'
#' For a given trait, this function ...
#'
#' @param pos Positive interactions. Input consists on a three colums data frame.
#'
#' @param neg Negative interactions. Input consists on a three colums data frame.
#'
#' @param network_name Network name which will be included on the output table.
#'
#' @param num_random_networks Number of random networks generated for calculating the confidence
#' intervals for each motiv.
#'
#' @param out Out format. Choose between "count" (number of motifs) or "normalized"
#' (number of motifs divided by the square of the node number).
#'
#' @param square.motifs TRUE or FALSE to include or not the four node motifs.
#'
#' @param cores Number of cores for paralelization. If missing, the number of available
#' cores will be autimatically detected and the parameter will be stablished as n°cores - 1.
#'
#' @return List of...
#'
#' @examples
#'
#' data(pos)
#' data(neg)
#'
#' conf.motifs.p(pos = pos, neg = neg, network_name = "test_network",
#' num_random_networks = 10, out = "count", cores = 2, square.motifs = TRUE)
#'
#' @export
#'
conf.motifs.p <- function(pos, neg, network_name, num_random_networks, out, cores, square.motifs = FALSE){

  #necesary internal function
  `%!in%` = Negate(`%in%`)
  library(doParallel)
  library(foreach)

  # input modification
  uniq.neg <- unique(c(as.character(unique(neg$Source)), as.character(unique(neg$Target))))
  uniq.pos <- unique(c(as.character(unique(pos$Source)), as.character(unique(pos$Target))))
  all <- rbind(neg, pos)
  uniq <- unique(c(as.character(unique(all$Source)), as.character(unique(all$Target))))


  # get network parameters
  # ## positive network
  # meanD <- c()
  # for(i in uniq.pos){
  #   meanD <- c(meanD, nrow(pos[pos$Source == i | pos$Target == i,]))
  # }
  # pos_mean_degree <- mean(meanD)
  #
  # m <- (pos_mean_degree * length(uniq.pos))/3
  # TT <- length(uniq.pos)^2
  # p.pos <- (m * 2)/TT
  #
  # ## negative network
  # meanD <- c()
  # i <- uniq.neg[1]
  # for(i in uniq.neg){
  #   meanD <- c(meanD, nrow(neg[neg$Source == i | neg$Target == i,]))
  # }
  # neg_mean_degree <- mean(meanD)
  #
  # m <- (neg_mean_degree * length(uniq.neg))/3
  # TT <- length(uniq.neg)^2
  # p.neg <- (m * 2)/TT


  # PARALELIZACION
  if( missing(cores) ){
    cores <-  parallel::detectCores()
  }
  cl <- parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)

  finalMatrix <- foreach::foreach(w = 1:num_random_networks, .combine=rbind,
                         .export = c("motif1", "motif2", "motif3", "motif4", "motif5",
                                     "motif6", "motif7", "motif8", "motif9", "motif10",
                                     "motif11", "motif12")) %dopar% {

    #necessary internal function
    `%!in%` = Negate(`%in%`)

    # output table creation and modification

    if( square.motifs == TRUE ){
      # output table creation and modification
      count <- data.frame(matrix(ncol = 15, nrow = 1))
      colnames(count) <- c("network", "motif_1", "motif_2", "motif_3", "motif_4", "motif_5", "motif_6", "motif_7",
                           "motif_8", "motif_9", "motif_10", "motif_11", "motif_12", "nodes", "links")

    } else {

      # output table creation and modification
      count <- data.frame(matrix(ncol = 10, nrow = 1))
      colnames(count) <- c("network", "motif_1", "motif_2", "motif_3", "motif_4",
                           "motif_5", "motif_6", "motif_7", "nodes", "links")

    }

    count[, colnames(count) == "network"] <- network_name
    class(count$nodes) <- as.numeric(count$nodes)
    count[, colnames(count) == "nodes"] <- length(uniq)
    count[, colnames(count) == "links"] <- nrow(all)

    # ## Comment later ##
    # # simulate networks
    # a <- as.data.frame(simcausal::rnet.gnp(length(uniq.neg), p.neg))
    # a$Source <- as.character(seq(1:nrow(a)))
    # a <- reshape2::melt(a)
    # a <- a[,c(1,3,2)]
    # a$variable <- "NEG"
    # a <- a[stats::complete.cases(a$value),]
    # a$value <- as.character(a$value)
    # colnames(a) <-  c("Source", "Target", "Type")
    # a$Source <- gsub("", "a", a$Source)
    # a$Target <- gsub("", "a", a$Target)
    # neg.s <- a
    #
    # a <- as.data.frame(simcausal::rnet.gnp(length(uniq.pos), p.pos))
    # a$Source <- as.character(seq(1:nrow(a)))
    # a <- reshape2::melt(a)
    # a <- a[,c(1,3,2)]
    # a$variable <- "POS"
    # a <- a[stats::complete.cases(a$value),]
    # a$value <- as.character(a$value)
    # colnames(a) <-  c("Source", "Target", "Type")
    # a$Source <- gsub("", "a", a$Source)
    # a$Target <- gsub("", "a", a$Target)
    # pos.s <- a
    # ##

    pos.length <- nrow(pos)
    neg.length <- nrow(neg)

    pos.s <- data.frame(Source = rep(NA, pos.length), Target = rep(NA, pos.length), Type = rep("POS", pos.length))
    i <- 1
    while(i <= pos.length){
      pair <- sample(x = uniq, size = 2, replace = F)
      tmp1 <- pos.s[pos.s$Source == pair[1] | pos.s$Target == pair[1],]
      if(pair[2] %!in% c(as.character(tmp1$Source), as.character(tmp1$Target))){
        pos.s[i,1] <- pair[1]
        pos.s[i,2] <- pair[2]
        i <- i + 1
        rm(pair, tmp1)
      } else {NULL}
    }
    neg.s <- data.frame(Source = rep(NA, neg.length), Target = rep(NA, neg.length), Type = rep("NEG", neg.length))
    i <- 1
    while(i <= neg.length){
      pair <- sample(x = uniq, size = 2, replace = F)
      temp1 <- pos.s[pos.s$Source == pair[1] | pos.s$Target == pair[1],]
      temp2 <- neg.s[neg.s$Source == pair[1] | neg.s$Target == pair[1],]
      if(pair[2] %!in% c(as.character(temp1$Source), as.character(temp1$Target))){
        if(pair[2] %!in% c(as.character(temp2$Source), as.character(temp2$Target))){
          neg.s[i,1] <- pair[1]
          neg.s[i,2] <- pair[2]
          i <- i + 1
        } else {NULL}
      } else {NULL}
    }
    rm(temp1, temp2, i)


    uniq.neg.s <- unique(c(as.character(unique(neg.s$Source)), as.character(unique(neg.s$Target))))
    uniq.pos.s <- unique(c(as.character(unique(pos.s$Source)), as.character(unique(pos.s$Target))))
    all.s <- rbind(neg.s, pos.s)
    uniq.s <- unique(c(as.character(unique(all.s$Source)), as.character(unique(all.s$Target))))



    # Motif 1
    res <- c()
    results <- motif1(pos = pos.s, uniq.pos = uniq.pos.s, res =  res)

    if(is.vector(results)){
      a <- as.data.frame(split(results, ceiling(seq_along(results)/3)))
      a.m <- t(a)
      a.df <- as.data.frame(t(a))
      rep <- duplicated(
        lapply(1:nrow(a.m), function(y){
          A <- a.m[y, ]
          A[order(A)]
        }))
      a.df.f <- a.df[!rep,]

      if( out == "count" ) {
        count[1 , colnames(count) == "motif_1"] <- nrow(a.df.f)
      } else {NULL}
      if( out == "normalized" ) {
        count[1 , colnames(count) == "motif_1"] <- nrow(a.df.f)/(length(uniq)^2)
      } else {NULL}
    } else { count[1 , colnames(count) == "motif_1"] <- 0 }

    # Motif 2
    res <- c()
    results <- motif2(neg = neg.s, uniq.neg = uniq.neg.s, res = res)

    if(is.vector(results)){
      a <- as.data.frame(split(results, ceiling(seq_along(results)/3)))
      a.m <- t(a)
      a.df <- as.data.frame(t(a))
      rep <- duplicated(
        lapply(1:nrow(a.m), function(y){
          A <- a.m[y, ]
          A[order(A)]
        }))
      a.df.f <- a.df[!rep,]

      if( out == "count" ) {
        count[1 , colnames(count) == "motif_2"] <- nrow(a.df.f)
      } else {NULL}
      if( out == "normalized" ) {
        count[1 , colnames(count) == "motif_2"] <- nrow(a.df.f)/(length(uniq)^2)
      } else {NULL}
    } else { count[1 , colnames(count) == "motif_2"] <- 0 }

    # Motif 3
    res <- c()
    results <- motif3(pos = pos.s, neg = neg.s, all = all.s, uniq = uniq.s, res = res)

    if(is.vector(results)){
      a <- as.data.frame(split(results, ceiling(seq_along(results)/3)))
      a.m <- t(a)
      a.df <- as.data.frame(t(a))
      rep <- duplicated(
        lapply(1:nrow(a.m), function(y){
          A <- a.m[y, ]
          A[order(A)]
        }))
      a.df.f <- a.df[!rep,]

      if( out == "count" ) {
        count[1 , colnames(count) == "motif_3"] <- nrow(a.df.f)
      } else {NULL}
      if( out == "normalized" ) {
        count[1 , colnames(count) == "motif_3"] <- nrow(a.df.f)/(length(uniq)^2)
      } else {NULL}

    } else { count[1 , colnames(count) == "motif_3"] <- 0 }

    # Motif 4
    res <- c()
    results <- motif4(pos = pos.s, neg = neg.s, all = all.s, uniq = uniq.s, res = res)

    if(is.vector(results)){
      a <- as.data.frame(split(results, ceiling(seq_along(results)/3)))
      a.m <- t(a)
      a.df <- as.data.frame(t(a))
      rep <- duplicated(
        lapply(1:nrow(a.m), function(y){
          A <- a.m[y, ]
          A[order(A)]
        }))
      a.df.f <- a.df[!rep,]

      if( out == "count" ) {
        count[1 , colnames(count) == "motif_4"] <- nrow(a.df.f)
      } else {NULL}
      if( out == "normalized" ) {
        count[1 , colnames(count) == "motif_4"] <- nrow(a.df.f)/(length(uniq)^2)
      } else {NULL}
    } else { count[1 , colnames(count) == "motif_4"] <- 0 }

    # Motif 5
    res <- c()
    results <- motif5(pos = pos.s, uniq.pos = uniq.pos.s, all = all.s, res = res)

    if(is.vector(results)){
      a <- as.data.frame(split(results, ceiling(seq_along(results)/3)))
      a.m <- t(a)
      a.df <- as.data.frame(t(a))
      rep <- duplicated(
        lapply(1:nrow(a.m), function(y){
          A <- a.m[y, ]
          A[order(A)]
        }))
      a.df.f <- a.df[!rep,]

      if( out == "count" ) {
        count[1 , colnames(count) == "motif_5"] <- nrow(a.df.f)
      } else {NULL}
      if( out == "normalized" ) {
        count[1 , colnames(count) == "motif_5"] <- nrow(a.df.f)/(length(uniq)^2)
      } else {NULL}

    } else { count[1 , colnames(count) == "motif_5"] <- 0 }

    # Motif 6
    res <- c()
    results <- motif6(neg = neg.s, uniq.neg = uniq.neg.s, all = all.s, res = res)

    if(is.vector(results)){
      a <- as.data.frame(split(results, ceiling(seq_along(results)/3)))
      a.m <- t(a)
      a.df <- as.data.frame(t(a))
      rep <- duplicated(
        lapply(1:nrow(a.m), function(y){
          A <- a.m[y, ]
          A[order(A)]
        }))
      a.df.f <- a.df[!rep,]

      if( out == "count" ) {
        count[1 , colnames(count) == "motif_6"] <- nrow(a.df.f)
      } else {NULL}
      if( out == "normalized" ) {
        count[1 , colnames(count) == "motif_6"] <- nrow(a.df.f)/(length(uniq)^2)
      } else {NULL}
    } else { count[1 , colnames(count) == "motif_6"] <- 0 }

    # Motif 7
    res <- c()
    results <- motif7(pos = pos.s, neg = neg.s, all = all.s, uniq = uniq.s, res = res)

    if(is.vector(results)){
      a <- as.data.frame(split(results, ceiling(seq_along(results)/3)))
      a.m <- t(a)
      a.df <- as.data.frame(t(a))
      rep <- duplicated(
        lapply(1:nrow(a.m), function(y){
          A <- a.m[y, ]
          A[order(A)]
        }))
      a.df.f <- a.df[!rep,]

      if( out == "count" ) {
        count[1 , colnames(count) == "motif_7"] <- nrow(a.df.f)
      } else {NULL}
      if( out == "normalized" ) {
        count[1 , colnames(count) == "motif_7"] <- nrow(a.df.f)/(length(uniq)^2)
      } else {NULL}
    } else { count[1 , colnames(count) == "motif_7"] <- 0 }

    if( square.motifs == TRUE ){

      # Motif 8
      res <- c()
      results <- motif8(pos = pos.s, uniq.pos = uniq.pos.s, res = res)

      if(is.vector(results)){
        a <- as.data.frame(split(results, ceiling(seq_along(results)/4)))
        a.m <- t(a)
        a.df <- as.data.frame(t(a))
        rep <- duplicated(
          lapply(1:nrow(a.m), function(y){
            A <- a.m[y, ]
            A[order(A)]
          }))
        a.df.f <- a.df[!rep,]

        if( out == "count" ) {
          count[1 , colnames(count) == "motif_8"] <- nrow(a.df.f)
        } else {NULL}
        if( out == "normalized" ) {
          count[1 , colnames(count) == "motif_8"] <- nrow(a.df.f)/(length(uniq)^2)
        } else {NULL}
      } else { count[1 , colnames(count) == "motif_8"] <- 0 }

      # Motif 9
      res <- c()
      results <- motif9(pos = pos.s, neg = neg.s, uniq.neg = uniq.neg.s, res = res)

      if(is.vector(results)){
        a <- as.data.frame(split(results, ceiling(seq_along(results)/4)))
        a.m <- t(a)
        a.df <- as.data.frame(t(a))
        rep <- duplicated(
          lapply(1:nrow(a.m), function(y){
            A <- a.m[y, ]
            A[order(A)]
          }))
        a.df.f <- a.df[!rep,]

        if( out == "count" ) {
          count[1 , colnames(count) == "motif_9"] <- nrow(a.df.f)
        } else {NULL}
        if( out == "normalized" ) {
          count[1 , colnames(count) == "motif_9"] <- nrow(a.df.f)/(length(uniq)^2)
        } else {NULL}
      } else { count[1 , colnames(count) == "motif_9"] <- 0 }

      # Motif 10
      res <- c()
      results <- motif10(pos = pos.s, neg = neg.s, uniq.neg = uniq.neg.s, res = res)

      if(is.vector(results)){
        a <- as.data.frame(split(results, ceiling(seq_along(results)/4)))
        a.m <- t(a)
        a.df <- as.data.frame(t(a))
        rep <- duplicated(
          lapply(1:nrow(a.m), function(y){
            A <- a.m[y, ]
            A[order(A)]
          }))
        a.df.f <- a.df[!rep,]

        if( out == "count" ) {
          count[1 , colnames(count) == "motif_10"] <- nrow(a.df.f)
        } else {NULL}
        if( out == "normalized" ) {
          count[1 , colnames(count) == "motif_10"] <- nrow(a.df.f)/(length(uniq)^2)
        } else {NULL}
      } else { count[1 , colnames(count) == "motif_10"] <- 0 }

      # Motif 11
      res <- c()
      results <- motif11(pos = pos.s, neg = neg.s, uniq.pos = uniq.pos.s, res = res)

      if(is.vector(results)){
        a <- as.data.frame(split(results, ceiling(seq_along(results)/4)))
        a.m <- t(a)
        a.df <- as.data.frame(t(a))
        rep <- duplicated(
          lapply(1:nrow(a.m), function(y){
            A <- a.m[y, ]
            A[order(A)]
          }))
        a.df.f <- a.df[!rep,]

        if( out == "count" ) {
          count[1 , colnames(count) == "motif_11"] <- nrow(a.df.f)
        } else {NULL}
        if( out == "normalized" ) {
          count[1 , colnames(count) == "motif_11"] <- nrow(a.df.f)/(length(uniq)^2)
        } else {NULL}
      } else { count[1 , colnames(count) == "motif_11"] <- 0 }

      # Motif 12
      res <- c()
      results <- motif12(neg = neg.s, uniq.neg = uniq.neg.s, res = res)

      if(is.vector(results)){
        a <- as.data.frame(split(results, ceiling(seq_along(results)/4)))
        a.m <- t(a)
        a.df <- as.data.frame(t(a))
        rep <- duplicated(
          lapply(1:nrow(a.m), function(y){
            A <- a.m[y, ]
            A[order(A)]
          }))
        a.df.f <- a.df[!rep,]

        if( out == "count" ) {
          count[1 , colnames(count) == "motif_12"] <- nrow(a.df.f)
        } else {NULL}
        if( out == "normalized" ) {
          count[1 , colnames(count) == "motif_12"] <- nrow(a.df.f)/(length(uniq)^2)
        } else {NULL}
      } else { count[1 , colnames(count) == "motif_12"] <- 0 }

    }

    ## UNA VEZ TODO HA TERMINADO ....

    count #Equivalent to finalMatrix = cbind(finalMatrix, tempMatrix)
  }
  #stop cluster
  parallel::stopCluster(cl)

  return(finalMatrix)

}
