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

motifs <- function(pos, neg, network_name = "Network", square.motifs = FALSE){
  #necessary internal function
  `%!in%` = Negate(`%in%`)

  # CHECKS
  if (missing(pos)) stop("Positive network not provided")
  if (ncol(pos) != 3) stop("Positive network must be on a 3 columns data frame format")
  if (missing(neg)) stop("Negative network not provided")
  if (ncol(neg) != 3) stop("Negative network must be on a 3 columns data frame format")





  # input modification
  colnames(pos) <- c("Source", "Target", "Type")
  pos$Type <- rep("POS", nrow(pos))
  colnames(neg) <- c("Source", "Target", "Type")
  neg$Type <- rep("NEG", nrow(neg))

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

    print("Only 3 nodes motifs will be extracted...")

    }

  count[colnames(count) == "network", 1] <- network_name
  class(count$nodes) <- as.numeric(count$nodes)
  count[1, colnames(count) == "nodes"] <- length(uniq)
  count[1, colnames(count) == "links"] <- nrow(all)

  # tax count

  groups$group <- uniq


  ## MOTIFS ANALYSIS

  # Motif 1
  results <- motif1(pos = pos)
  count[1 , colnames(count) == "motif_1"] <- mod_results(results = results, nodes = 3)
  groups[, colnames(groups) == "motif_1"] <- extractCount(results = results, nodes = 3, uniq = uniq)

  print("Motif 1 finished")

  # Motif 2
  results <- motif2(neg = neg)
  count[1 , colnames(count) == "motif_2"] <- mod_results(results = results, nodes = 3)
  groups[, colnames(groups) == "motif_2"] <- extractCount(results = results, nodes = 3, uniq = uniq)

  print("Motif 2 finished")

  # Motif 3
  results <- motif3(pos = pos, neg = neg)
  count[1 , colnames(count) == "motif_3"] <- mod_results(results = results, nodes = 3)
  groups[, colnames(groups) == "motif_3"] <- extractCount(results = results, nodes = 3, uniq = uniq)

  print("Motif 3 finished")

  # Motif 4
  results <- motif4(pos = pos, neg = neg)
  count[1 , colnames(count) == "motif_4"] <- mod_results(results = results, nodes = 3)
  groups[, colnames(groups) == "motif_4"] <- extractCount(results = results, nodes = 3, uniq = uniq)

  print("Motif 4 finished")

  # Motif 5
  results <- motif5(pos = pos, neg = neg)
  count[1 , colnames(count) == "motif_5"] <- mod_results(results = results, nodes = 3)
  groups[, colnames(groups) == "motif_5"] <- extractCount(results = results, nodes = 3, uniq = uniq)

  print("Motif 5 finished")

  # Motif 6
  results <- motif6(neg = neg, pos = pos)
  count[1 , colnames(count) == "motif_6"] <- mod_results(results = results, nodes = 3)
  groups[, colnames(groups) == "motif_6"] <- extractCount(results = results, nodes = 3, uniq = uniq)

  print("Motif 6 finished")

  # Motif 7
  results <- motif7(pos = pos, neg = neg)
  count[1 , colnames(count) == "motif_7"] <- mod_results(results = results, nodes = 3)
  groups[, colnames(groups) == "motif_7"] <- extractCount(results = results, nodes = 3, uniq = uniq)

  print("Motif 7 finished")


  if( square.motifs == TRUE ){

    # Motif 8
    results <- motif8(pos = pos)
    count[1 , colnames(count) == "motif_8"] <- mod_results(results = results, nodes = 4)
    groups[, colnames(groups) == "motif_8"] <- extractCount(results = results, nodes = 4, uniq = uniq)

    print("Motif 8 finished")

    # Motif 9
    results <- motif9(pos = pos, neg = neg)
    count[1 , colnames(count) == "motif_9"] <- mod_results(results = results, nodes = 4)
    groups[, colnames(groups) == "motif_9"] <- extractCount(results = results, nodes = 4, uniq = uniq)

    print("Motif 9 finished")

    # Motif 10
    results <- motif10(pos = pos, neg = neg)
    count[1 , colnames(count) == "motif_10"] <- mod_results(results = results, nodes = 4)
    groups[, colnames(groups) == "motif_10"] <- extractCount(results = results, nodes = 4, uniq = uniq)

    print("Motif 10 finished")

    # Motif 11
    results <- motif11(pos = pos, neg = neg)
    count[1 , colnames(count) == "motif_11"] <- mod_results(results = results, nodes = 4)
    groups[, colnames(groups) == "motif_11"] <- extractCount(results = results, nodes = 4, uniq = uniq)

    print("Motif 11 finished")

    # Motif 12
    results <- motif12(neg = neg)
    count[1 , colnames(count) == "motif_12"] <- mod_results(results = results, nodes = 4)
    groups[, colnames(groups) == "motif_12"] <- extractCount(results = results, nodes = 4, uniq = uniq)

    print("Motif 12 finished")

  } else {NULL}

  # EXPORT RESULTS

  final.results <- list("count" = count, "groups" = groups)

  return(final.results)
}
