
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
#'
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
#'   library(Biostrings)
#'   # For DNAString input:
#'   seq <- DNAString("ACGTACGT")
#'   pos <- 4
#'   del <- 3
#'   # Compute homeology indices with debugging enabled.
#'   result <- find_homeology.list(seq, pos, del, debug = TRUE)
#' }
#'
#' @export
find_homeology.list <- function(seq, position, deletion, debug = F, thresh = 0.8, left = T, right = T) {
  if(grepl("DNAString", class(seq))) {
    sub.fn <- Biostrings::subseq
    str.split.fn <- base::strsplit
    str.reverse <- Biostrings::reverse
    x.s<- x.t<- y.s<- y.t <- DNAString()
    paste.fn <- xscat
  } else {
    sub.fn <- base::substr
    str.split.fn <- base::strsplit
    str.reverse <- stringi::stri_reverse
    x.s<- x.t<- y.s<- y.t <- ""
    paste.fn <- paste0
  }  
  tbl <- data.table(count.x = rep(1, length(seq)),
                    max.x = rep(0, length(seq)),
                    match.x = rep(0, length(seq)),
                    max.y = rep(0, length(seq)),
                    match.y = rep(0, length(seq)),
                    count.y = rep(1, length(seq)))
  
  if(left){
    if(debug) print("5' retained sequence vs 3' deleted sequence")
    while(any(tbl$count.x < deletion)){
      i1 <- which(tbl$count.x < deletion)
      print(length(i1))
      e1 <- (position[i1] - tbl$count.x[i1])
      e2 <- (position[i1] + deletion[i1] - tbl$count.x[i1])
      s1 = sub.fn(seq[i1], e1, e1)
      s2 = sub.fn(seq[i1], e2, e2)
      if(debug) print(paste.fn(s1, s2))
      tbl[i1, match.x := match.x + mapply(function(x,y) sum(x==y),
                                          str.split.fn(as.character(s1),""),
                                          str.split.fn(as.character(s2),"")) ]
      if(any(tbl[i1]$match.x / tbl[i1]$count.x >= thresh & s1 == s2)){
        i2 <- which(tbl[i1]$match.x / tbl[i1]$count.x  >= thresh & s1 == s2)
        tbl[i2, max.x := count.x]
      }
      tbl[count.x < deletion, count.x := count.x + 1]
      x.s <- paste.fn(s1, x.s); x.t <- paste.fn(s2, x.t)
    }
  }
  if(right) {
    if(debug) print("5' deleted sequence vs 3' retained sequence")
    while(any(tbl$count.y < deletion)){
      i1 <- which(tbl$count.y < deletion)
      print(length(i1))
      e1 <- position[i1] + tbl$count.y[i1] - 1
      e2 <- position[i1] + tbl$count.y[i1] + deletion[i1] - 1
      s1 = sub.fn(seq[i1], e1, e1)
      s2 = sub.fn(seq[i1], e2, e2)
      if(debug) print(paste0(s1, s2))
      tbl[i1, match.y := match.y + mapply(function(x,y) sum(x==y),
                                          str.split.fn(as.character(s1),""),
                                          str.split.fn(as.character(s2),""))]
      if(any(tbl[i1]$match.y / tbl[i1]$count.y >= thresh & s1 == s2)){
        i2 <- which(tbl[i1]$match.y / tbl[i1]$count.y >= thresh & s1 == s2)
        tbl[i2, max.y := count.y]
      }
      tbl[count.y < deletion, count.y := count.y + 1]
      y.s <- paste.fn(y.s, s1); y.t <- paste.fn(y.t, s2)
    }
  }
  if(debug) print(paste0(tbl$max.x, " vs. ", tbl$max.y))
  if(left && !right) return(tbl$max.x)
  else if(!left && right) return(tbl$max.y)
  return(c(tbl$max.x, tbl$max.y))
}
