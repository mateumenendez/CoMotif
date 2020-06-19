#' Extract motifs from co-occurrence/co-exclusion networks.
#'
#' Quantify 3 and square 4 nodes mofifs from co-occurrence and co-exclusions networks.
#'
#' @param pos Positive interactions. Input consists on a three colums data frame.
#'
#' @param neg Negative interactions. Input consists on a three colums data frame.
#'
#' @param network_name Network name which will be included on the output table.
#'
#' @param out Out format. Choose between "count" (number of motifs) or "normalized"
#' (number of motifs divided by the square of the node number).
#'
#' @param square.motifs TRUE or FALSE to include or not the four node motifs.
#'
#' @return A list of two data frames...
#'
#' @examples
#'
#' data(pos)
#' data(neg)
#'
#' motifs(pos = pos, neg = neg, network_name = "test_network", out = "count")
#' motifs(pos = pos, neg = neg, network_name = "test_network", out = "normalized", square.motifs = TRUE)
#'
#' @export

motifs <- function(pos, neg, network_name, out, square.motifs = FALSE){
  #necessary internal function
  `%!in%` = Negate(`%in%`)

  # input modification
  uniq.neg <- unique(c(as.character(unique(neg$Source)), as.character(unique(neg$Target))))
  uniq.pos <- unique(c(as.character(unique(pos$Source)), as.character(unique(pos$Target))))
  all <- rbind(neg, pos)
  uniq <- unique(c(as.character(unique(all$Source)), as.character(unique(all$Target))))

  if( square.motifs == TRUE ){
    # output table creation and modification
    count <- data.frame(matrix(ncol = 15, nrow = 1))
    colnames(count) <- c("network", "motif_1", "motif_2", "motif_3", "motif_4", "motif_5", "motif_6", "motif_7",
                         "motif_8", "motif_9", "motif_10", "motif_11", "motif_12", "nodes", "links")

    groups <- as.data.frame(matrix(nrow = length(uniq), ncol = 13, NA))
    colnames(groups) <- c("group", "motif_1", "motif_2", "motif_3",
                          "motif_4", "motif_5", "motif_6", "motif_7",
                          "motif_8", "motif_9", "motif_10", "motif_11", "motif_12")
  } else {

    # output table creation and modification
    count <- data.frame(matrix(ncol = 10, nrow = 1))
    colnames(count) <- c("network", "motif_1", "motif_2", "motif_3", "motif_4",
                               "motif_5", "motif_6", "motif_7", "nodes", "links")

    groups <- as.data.frame(matrix(nrow = length(uniq), ncol = 8, NA))
    colnames(groups) <- c("group", "motif_1", "motif_2", "motif_3",
                          "motif_4", "motif_5", "motif_6", "motif_7")

    }

  count[colnames(count) == "network", 1] <- network_name
  class(count$nodes) <- as.numeric(count$nodes)
  count[1, colnames(count) == "nodes"] <- length(uniq)
  count[1, colnames(count) == "links"] <- nrow(all)

  # tax count

  groups$group <- uniq


  ## MOTIFS ANALYSIS

  # Motif 1
  res <- c()
  results <- motif1(pos = pos, uniq.pos = uniq.pos, res =  res)

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

    for(x in uniq){
      suma <- sum(sum(a.df.f$V1 == x), sum(a.df.f$V2 == x), sum(a.df.f$V3 == x))
      groups[groups$group == x, colnames(groups) == "motif_1"] <- suma
    }
  } else {
    count[1 , colnames(count) == "motif_1"] <- 0
    groups[, colnames(groups) == "motif_1"] <- 0 }


  print("Motif 1 finished")

  # Motif 2
  res <- c()
  results <- motif2(neg = neg, uniq.neg = uniq.neg, res = res)

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

    for(x in uniq){
      suma <- sum(sum(a.df.f$V1 == x), sum(a.df.f$V2 == x), sum(a.df.f$V3 == x))
      groups[groups$group == x, colnames(groups) == "motif_2"] <- suma
    }
  } else {
    count[1 , colnames(count) == "motif_2"] <- 0
    groups[, colnames(groups) == "motif_2"] <- 0 }


  print("Motif 2 finished")

  # Motif 3
  res <- c()
  results <- motif3(pos = pos, neg = neg, all = all, uniq = uniq, res = res)

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

    for(x in uniq){
      suma <- sum(sum(a.df.f$V1 == x), sum(a.df.f$V2 == x), sum(a.df.f$V3 == x))
      groups[groups$group == x, colnames(groups) == "motif_3"] <- suma
    }
  } else {
    count[1 , colnames(count) == "motif_3"] <- 0
    groups[, colnames(groups) == "motif_3"] <- 0 }


  print("Motif 3 finished")

  # Motif 4
  res <- c()
  results <- motif4(pos = pos, neg = neg, all = all, uniq = uniq, res = res)

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

    for(x in uniq){
      suma <- sum(sum(a.df.f$V1 == x), sum(a.df.f$V2 == x), sum(a.df.f$V3 == x))
      groups[groups$group == x, colnames(groups) == "motif_4"] <- suma
    }
  } else {
    count[1 , colnames(count) == "motif_4"] <- 0
    groups[, colnames(groups) == "motif_4"] <- 0 }


  print("Motif 4 finished")

  # Motif 5
  res <- c()
  results <- motif5(pos = pos, uniq.pos = uniq.pos, all = all, res = res)

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

    for(x in uniq){
      suma <- sum(sum(a.df.f$V1 == x), sum(a.df.f$V2 == x), sum(a.df.f$V3 == x))
      groups[groups$group == x, colnames(groups) == "motif_5"] <- suma
    }
  } else {
    count[1 , colnames(count) == "motif_5"] <- 0
    groups[, colnames(groups) == "motif_5"] <- 0 }


  print("Motif 5 finished")

  # Motif 6
  res <- c()
  results <- motif6(neg = neg, uniq.neg = uniq.neg, all = all, res = res)

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

    for(x in uniq){
      suma <- sum(sum(a.df.f$V1 == x), sum(a.df.f$V2 == x), sum(a.df.f$V3 == x))
      groups[groups$group == x, colnames(groups) == "motif_6"] <- suma
    }
  } else {
    count[1 , colnames(count) == "motif_6"] <- 0
    groups[, colnames(groups) == "motif_6"] <- 0 }

  print("Motif 6 finished")


  # Motif 7
  res <- c()
  results <- motif7(pos = pos, neg = neg, all = all, uniq = uniq, res = res)

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

    for(x in uniq){
      suma <- sum(sum(a.df.f$V1 == x), sum(a.df.f$V2 == x), sum(a.df.f$V3 == x))
      groups[groups$group == x, colnames(groups) == "motif_7"] <- suma
    }
  } else {
    count[1 , colnames(count) == "motif_7"] <- 0
    groups[, colnames(groups) == "motif_7"] <- 0 }

  print("Motif 7 finished")




  if( square.motifs == TRUE ){

    # Motif 8
    res <- c()
    results <- motif8(pos = pos, uniq.pos = uniq.pos, res = res)

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

      for(x in uniq){
        suma <- sum(sum(a.df.f$V1 == x), sum(a.df.f$V2 == x), sum(a.df.f$V3 == x))
        groups[groups$group == x, colnames(groups) == "motif_8"] <- suma
      }
    } else {
      count[1 , colnames(count) == "motif_8"] <- 0
      groups[, colnames(groups) == "motif_8"] <- 0 }

    print("Motif 8 finished")


    # Motif 9
    res <- c()
    results <- motif9(pos = pos, neg = neg, uniq.neg = uniq.neg, res = res)

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

      for(x in uniq){
        suma <- sum(sum(a.df.f$V1 == x), sum(a.df.f$V2 == x), sum(a.df.f$V3 == x))
        groups[groups$group == x, colnames(groups) == "motif_9"] <- suma
      }
    } else {
      count[1 , colnames(count) == "motif_9"] <- 0
      groups[, colnames(groups) == "motif_9"] <- 0 }

    print("Motif 9 finished")



    # Motif 10
    res <- c()
    results <- motif10(pos = pos, neg = neg, uniq.neg = uniq.neg, res = res)

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

      for(x in uniq){
        suma <- sum(sum(a.df.f$V1 == x), sum(a.df.f$V2 == x), sum(a.df.f$V3 == x))
        groups[groups$group == x, colnames(groups) == "motif_10"] <- suma
      }
    } else {
      count[1 , colnames(count) == "motif_10"] <- 0
      groups[, colnames(groups) == "motif_10"] <- 0 }

    print("Motif 10 finished")



    # Motif 11
    res <- c()
    results <- motif11(pos = pos, neg = neg, uniq.pos = uniq.pos, res = res)

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

      for(x in uniq){
        suma <- sum(sum(a.df.f$V1 == x), sum(a.df.f$V2 == x), sum(a.df.f$V3 == x))
        groups[groups$group == x, colnames(groups) == "motif_11"] <- suma
      }
    } else {
      count[1 , colnames(count) == "motif_11"] <- 0
      groups[, colnames(groups) == "motif_11"] <- 0 }

    print("Motif 11 finished")


    # Motif 12
    res <- c()
    results <- motif12(neg = neg, uniq.neg = uniq.neg, res = res)

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

      for(x in uniq){
        suma <- sum(sum(a.df.f$V1 == x), sum(a.df.f$V2 == x), sum(a.df.f$V3 == x))
        groups[groups$group == x, colnames(groups) == "motif_12"] <- suma
      }
    } else {
      count[1 , colnames(count) == "motif_12"] <- 0
      groups[, colnames(groups) == "motif_12"] <- 0 }

    print("Motif 12 finished")

  } else {NULL}


  # EXPORT RESULTS

  final.results <- list("count" = count, "groups" = groups)

  return(final.results)
}
