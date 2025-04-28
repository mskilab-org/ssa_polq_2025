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
