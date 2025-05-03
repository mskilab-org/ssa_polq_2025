#' Find Homologous Regions in a Sequence
#'
#' This function identifies the extent of homology around a given deletion region within a sequence.
#' It compares the bases from two regions: one extending leftwards (5' retained vs 3' deleted) and
#' the other extending rightwards (5' deleted vs 3' retained) from the deletion breakpoint.
#'
#' @param seq A sequence to be analyzed, provided either as a character vector or a DNAString object from Biostrings.
#' @param position An integer indicating the reference position (1-indexed) in the sequence where analysis begins.
#' @param deletion An integer specifying the offset of the deletion relative to the reference position.
#' @param debug Logical flag to enable verbose debugging output. Defaults to FALSE.
#' @param left Logical flag to determine whether to compute homology in the left (5') flanking region. Defaults to TRUE.
#' @param right Logical flag to determine whether to compute homology in the right (3') flanking region. Defaults to TRUE.
#'
#' @details
#' The function operates by iteratively comparing nucleotides in two regions:
#'   - When 'left' is TRUE, bases starting from positions (position - count.x - 1) and (position + deletion - count.x - 1)
#'     are compared until a mismatch is found or the beginning of the sequence is reached.
#'   - When 'right' is TRUE, bases from positions (position + count.y) and (position + deletion + count.y) are compared until
#'     a mismatch is detected or the end of the sequence is reached.
#'
#' The comparison is performed using Biostrings::subseq if seq is a DNAString, and base::substr otherwise.
#'
#' @return
#' A numeric vector containing the count(s) of matching bases. If both left and right homology checks are enabled,
#' the counts for the left region (count.x) are concatenated with the counts for the right region (count.y).
#' If only one direction is enabled, only the corresponding count vector is returned.
#'
#' @examples
#' \dontrun{
#' # Example using a DNAString sequence
#' library(Biostrings)
#' dna_seq <- DNAString("ATGCGATACGCTTACG")
#' position <- 5
#' deletion <- 3
#' homology <- find_homology.list(dna_seq, position, deletion, debug = TRUE)
#' }
#'
#' @export
find_homology.list <- function(seq, position, deletion, debug = F, left = T, right = T) {
  tbl <- data.table(count.x = rep(0, length(seq)),
                    count.y = rep(0, length(seq)))
  if(grepl("DNAString", class(seq))) {
    sub.fn <- Biostrings::subseq
  } else {
    sub.fn <- base::substr
  }
  if(left){
    if(debug) print("5' retained sequence vs 3' deleted sequence")
    if(debug){
      print(paste0("starting 5' retained sequence: ", sub.fn(seq, position - count.x - 1, position - count.x - 1)))
      print(paste0("starting 3' deleted sequence: ", sub.fn(seq, position + deletion - count.x - 1, position + deletion - count.x - 1)))
    }
    #browser()
    while(any(position - tbl$count.x >= 1)){
      index1 <- which((position - tbl$count.x) >= 1)
      if(any((position - tbl$count.x - 1)[index1] == 0))
        break
      eval1 <- (position - tbl$count.x - 1)[index1]
      eval2 <- (position + deletion - tbl$count.x - 1)[index1]
      index2 <- which(sub.fn(seq[index1], eval1 , eval1) == sub.fn(seq[index1], eval2, eval2))
      print(length(index2))
      if(length(index2) == 0)
        break
      tbl[index1[index2], count.x := count.x + 1]
    }    
  }
  if(right){
    if(debug) print("5' deleted sequence vs 3' retained sequence")
    if(debug){
      print(paste0("starting 5' deleted sequence: ", sub.fn(seq, position + deletion + count.y, position + deletion + count.y)))
      print(paste0("starting 3' retained sequence: ", sub.fn(seq, position + count.y, position +  count.y)))
    }
    while(any(position + deletion + tbl$count.y <= nchar(seq))){
      index1 <- which(position + deletion + tbl$count.y <= nchar(seq))
      eval1 <- (position + deletion + tbl$count.y)[index1]
      eval2 <- (position + tbl$count.y)[index1]
      index2 <- which(sub.fn(seq[index1], eval1 , eval1) == sub.fn(seq[index1], eval2, eval2 ))
      print(length(index2))
      if(length(index2) == 0)
        break
      tbl[index1[index2], count.y := count.y + 1]
    }
  }
  if(left && !right) return(tbl$count.x)
  else if(!left && right) return(tbl$count.y)
  return(c(tbl$count.x, tbl$count.y))
}


#' Find Homeology Indices in Sequences
#'
#' This function computes homeology indices by iteratively comparing sub-sequence segments from a given sequence. Depending on the type of the input sequence (e.g., a DNAString object or a character string), the function selects the appropriate sub-sequence extraction, splitting, and reversing routines.
#'
#' The comparison is performed on two sides:
#' \describe{
#'   \item{left (5' retained vs. 3' deleted):}{The function compares positions moving left from a reference point until a predefined deletion length is reached, updating a match-score based on a threshold ratio.}
#'   \item{right (5' deleted vs. 3' retained):}{Similarly, the function compares positions moving right from the reference point, maintaining a match score and a count until the deletion length requirement is fulfilled.}
#' }

#' @param seq A vector of sequences or a single sequence. The function supports both \code{DNAString} objects (from the Biostrings package) and character strings, adapting internal processing accordingly.
#' @param position A numeric vector indicating the reference position(s) in the sequence(s) around which the sub-sequence comparisons are performed.
#' @param deletion A numeric vector specifying the deletion window length for each sequence. This value controls how many positions will be compared.
#' @param debug Logical. If \code{TRUE}, the function prints debug information during execution. Default is \code{FALSE}.
#' @param thresh Numeric. A threshold ratio (default is 0.8) to determine when a sufficient match has been achieved between sub-sequences.
#' @param left Logical. If \code{TRUE}, performs the left side (5' retained vs. 3' deleted) comparison.
#' @param right Logical. If \code{TRUE}, performs the right side (5' deleted vs. 3' retained) comparison.
#'
#' @return Depending on the \code{left} and \code{right} parameters:
#' \itemize{
#'   \item If only \code{left = TRUE} or \code{right = TRUE}, returns a vector of maximum match counts for the specified side.
#'   \item If both are \code{TRUE}, returns a concatenated vector of maximum match counts for the left and right sides.
#' }
#'
#' @details
#' The function initializes a \code{data.table} to track the count of processed positions and the accumulated match score for both sides of the sequence. Iterative loops extend the comparison window until the deletion length is reached, and when the matching ratio (matches/processed count) exceeds the threshold and the compared bases are identical, the current count is recorded as the maximum index for that side.
#'
#' Depending on the sequence type, helper functions for sub-sequence extraction (\code{subseq} or \code{substr}), string splitting, and reversing are selected to ensure correct functionality.
#'
#' @examples
#' \dontrun{
#' library(Biostrings)
#' # For DNAString input:
#' seq <- DNAString("ACGTACGT")
#' pos <- 4
#' del <- 3
#' # Compute homeology indices with debugging enabled.
#' result <- find_homeology.list(seq, pos, del, debug = TRUE)
#' }
#'
#' @export
find_homeology.list <- function(seq, position, deletion, debug = F, thresh = 0.8, left = T, right = T) {
  if (grepl("DNAString", class(seq))) {
    sub.fn <- Biostrings::subseq
    str.split.fn <- base::strsplit
    str.reverse <- Biostrings::reverse
    x.s <- x.t <- y.s <- y.t <- DNAString()
    paste.fn <- xscat
  } else {
    sub.fn <- base::substr
    str.split.fn <- base::strsplit
    str.reverse <- stringi::stri_reverse
    x.s <- x.t <- y.s <- y.t <- ""
    paste.fn <- paste0
  }
  tbl <- data.table(
    count.x = rep(1, length(seq)),
    max.x = rep(0, length(seq)),
    match.x = rep(0, length(seq)),
    max.y = rep(0, length(seq)),
    match.y = rep(0, length(seq)),
    count.y = rep(1, length(seq))
  )

  if (left) {
    if (debug) print("5' retained sequence vs 3' deleted sequence")
    while (any(tbl$count.x < deletion)) {
      i1 <- which(tbl$count.x < deletion)
      print(length(i1))
      e1 <- (position[i1] - tbl$count.x[i1])
      e2 <- (position[i1] + deletion[i1] - tbl$count.x[i1])
      s1 = sub.fn(seq[i1], e1, e1)
      s2 = sub.fn(seq[i1], e2, e2)
      if (debug) print(paste.fn(s1, s2))
      tbl[i1, match.x := match.x + mapply(
        function(x, y) sum(x == y),
        str.split.fn(as.character(s1), ""),
        str.split.fn(as.character(s2), "")
      )]
      if (any(tbl[i1]$match.x / tbl[i1]$count.x >= thresh & s1 == s2)) {
        i2 <- which(tbl[i1]$match.x / tbl[i1]$count.x >= thresh & s1 == s2)
        tbl[i2, max.x := count.x]
      }
      tbl[count.x < deletion, count.x := count.x + 1]
      x.s <- paste.fn(s1, x.s)
      x.t <- paste.fn(s2, x.t)
    }
  }
  if (right) {
    if (debug) print("5' deleted sequence vs 3' retained sequence")
    while (any(tbl$count.y < deletion)) {
      i1 <- which(tbl$count.y < deletion)
      print(length(i1))
      e1 <- position[i1] + tbl$count.y[i1] - 1
      e2 <- position[i1] + tbl$count.y[i1] + deletion[i1] - 1
      s1 = sub.fn(seq[i1], e1, e1)
      s2 = sub.fn(seq[i1], e2, e2)
      if (debug) print(paste0(s1, s2))
      tbl[i1, match.y := match.y + mapply(
        function(x, y) sum(x == y),
        str.split.fn(as.character(s1), ""),
        str.split.fn(as.character(s2), "")
      )]
      if (any(tbl[i1]$match.y / tbl[i1]$count.y >= thresh & s1 == s2)) {
        i2 <- which(tbl[i1]$match.y / tbl[i1]$count.y >= thresh & s1 == s2)
        tbl[i2, max.y := count.y]
      }
      tbl[count.y < deletion, count.y := count.y + 1]
      y.s <- paste.fn(y.s, s1)
      y.t <- paste.fn(y.t, s2)
    }
  }
  if (debug) print(paste0(tbl$max.x, " vs. ", tbl$max.y))
  if (left && !right) {
    return(tbl$max.x)
  } else if (!left && right) {
    return(tbl$max.y)
  }
  return(c(tbl$max.x, tbl$max.y))
}

find_homology.gr <- function(gr, genome, debug = FALSE, left = TRUE, right = TRUE,
                 mc.cores = 1, junction = FALSE) {
  chrom_lengths <- width(genome)
  names(chrom_lengths) <- names(genome)
  
  results <- pbmcapply::pbmclapply(seq_along(gr), function(i) {
  left_match_i <- 0
  right_match_i <- 0


  if (debug)
    message("Processing ", chrom, ":", st, "-", en, " (deletion size = ", n_del, ")")

  if (left) {
    if (st <= n_del) {
    stop("Insufficient left context on ", chrom, " at position ", st, " for deletion of size ", n_del)
    }
    j <- 1
    repeat {
    pos_left <- st - j
    pos_del  <- en - j + 1
    if (pos_left < 1) break

    base_left <- as.character(unlist(genome[GRanges(seqnames = chrom,
                            ranges = IRanges(start = pos_left, end = pos_left))]))
    base_del  <- as.character(unlist(genome[GRanges(seqnames = chrom,
                             ranges = IRanges(start = pos_del, end = pos_del))]))
    if (debug)
      message("  Left j=", j, ": pos_left=", pos_left, " (", base_left, 
          ") vs. pos_del=", pos_del, " (", base_del, ")")

    if (base_left == base_del) {
      left_match_i <- left_match_i + 1
      j <- j + 1
      if (j > n_del) break
    } else {
      break
    }
    }
  }

  if (right) {
    if ((en + n_del) > as.numeric(chrom_lengths[chrom])) {
    stop("Insufficient right context on ", chrom, " at position ", en, 
       " for deletion of size ", n_del)
    }
    j <- 0
    repeat {
    pos_right <- st + j
    pos_adj   <- en + 1 + j
    if (pos_adj > as.numeric(chrom_lengths[chrom])) break

    base_right <- as.character(unlist(genome[GRanges(seqnames = chrom, 
                             ranges = IRanges(start = pos_right, end = pos_right))]))
    base_adj   <- as.character(unlist(genome[GRanges(seqnames = chrom, 
                             ranges = IRanges(start = pos_adj, end = pos_adj))]))
    if (debug)
      message("  Right j=", j, ": pos_right=", pos_right, " (", base_right, 
          ") vs. pos_adj=", pos_adj, " (", base_adj, ")")

    if (base_right == base_adj) {
      right_match_i <- right_match_i + 1
      j <- j + 1
      if (j >= n_del) break
    } else {
      break
    }
    }
  }

  list(left = left_match_i, right = right_match_i)
  }, mc.cores = mc.cores)

  left_match <- sapply(results, function(x) x$left)
  right_match <- sapply(results, function(x) x$right)

  if (left && right)
  return(list(left = left_match, right = right_match))
  else if (left)
  return(left_match)
  else if (right)
  return(right_match)
}

find_homeology.gr <- function(gr, genome, debug = FALSE, thresh = 0.8, left = TRUE, right = TRUE) {
  # gr: GRanges object representing one or more deletion regions.
  # genome: a BSgenome object (or similar) used for sequence look-up.
  # thresh: threshold ratio of matches/count to record homology extension (default 0.8).
  
  # Retrieve chromosome lengths from genome
  chrom_lengths <- seqlengths(genome)
  
  results <- lapply(seq_along(gr), function(i) {
    del <- gr[i]
    chrom <- as.character(seqnames(del))
    st <- start(del)
    en <- end(del)
    n_del <- width(del)
    
    if(debug)
      message("Processing ", chrom, ":", st, "-", en, " (deletion size = ", n_del, ")")
    
    left_match <- 0
    right_match <- 0
    
    # === LEFT HOMEOLOGY SEARCH ===
    if(left) {
      # Initialize counters: we start with count = 1 (first base comparison)
      j <- 1
      match_x <- 0
      count_x <- 0
      max_x <- 0
      while(j <= n_del && (st - j) >= 1) {
        pos_left <- st - j            # base in retained left flank
        pos_del  <- en - j + 1         # corresponding base in deletion region
        # Fetch bases using getSeq():
        base_left <- as.character(getSeq(genome, GRanges(seqnames = chrom,
                                                         ranges = IRanges(start = pos_left, end = pos_left))))
        base_del  <- as.character(getSeq(genome, GRanges(seqnames = chrom,
                                                         ranges = IRanges(start = pos_del, end = pos_del))))
        count_x <- count_x + 1
        if(base_left == base_del) {
          match_x <- match_x + 1
        }
        if(match_x/count_x >= thresh && base_left == base_del) {
          max_x <- j
        } else {
          break
        }
        if(debug)
          message("  Left j=", j, ": pos_left=", pos_left, " (", base_left,
                  ") vs. pos_del=", pos_del, " (", base_del, ") => ratio: ", round(match_x/count_x,3))
        j <- j + 1
      }
      left_match <- max_x
    }
    
    # === RIGHT HOMEOLOGY SEARCH ===
    if(right) {
      j <- 0  # here we start at 0 so that at j = 0 we compare the very first bases
      match_y <- 0
      count_y <- 0
      max_y <- 0
      while(j < n_del && (en + 1 + j) <= as.numeric(chrom_lengths[chrom])) {
        pos_right <- st + j            # retained base in right flank
        pos_adj   <- en + 1 + j          # corresponding base after deletion
        base_right <- as.character(getSeq(genome, GRanges(seqnames = chrom,
                                                          ranges = IRanges(start = pos_right, end = pos_right))))
        base_adj   <- as.character(getSeq(genome, GRanges(seqnames = chrom,
                                                          ranges = IRanges(start = pos_adj, end = pos_adj))))
        count_y <- count_y + 1
        if(base_right == base_adj) {
          match_y <- match_y + 1
        }
        if(match_y/count_y >= thresh && base_right == base_adj) {
          max_y <- j + 1  # j starts at 0 so add 1 for actual count
        } else {
          break
        }
        if(debug)
          message("  Right j=", j, ": pos_right=", pos_right, " (", base_right,
                  ") vs. pos_adj=", pos_adj, " (", base_adj, ") => ratio: ", round(match_y/count_y,3))
        j <- j + 1
      }
      right_match <- max_y
    }
    
    list(left = left_match, right = right_match)
  })
  
  left_vec <- sapply(results, function(x) x$left)
  right_vec <- sapply(results, function(x) x$right)
  
  if(left && right)
    return(list(left = left_vec, right = right_vec))
  else if(left)
    return(left_vec)
  else
    return(right_vec)
}
