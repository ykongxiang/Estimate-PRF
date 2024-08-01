##get all detected PRF period 
extract_sequences <- function(mrna_file) {
  mrna_data <- read_excel(mrna_file)
  mrna_data <- na.omit(mrna_data)  
  mrna_data$FS_start <- as.numeric(mrna_data$FS_start)
  mrna_data$FS_end <- as.numeric(mrna_data$FS_end)
  
  extract_and_reverse_substring <- function(text, start, end, strand) {
    if (strand == '-') {
      text <- as.character(reverseComplement(DNAString(text)))
    }
    # 提取子串
    substring <- substr(text, start - 2, end + 1)
    if (start > end){
      substring <- substr(text,end - 2, start + 1)
    }
    return(substring)
  }
  
  seqs <- sapply(1:nrow(mrna_data), function(i) {
    extract_and_reverse_substring(
      text = mrna_data$Sequence[i],
      start = mrna_data$FS_start[i],
      end = mrna_data$FS_end[i],
      strand = mrna_data$Strand[i]
    )
  })
  return(seqs)
}
#------------------------------------------------------------------------------------------
#calculate all kinds of seqs'frequence and order it
set_data_frame <- function(x){
  string_counts <- table(x)
  total_count <- sum(string_counts)
  string_frequencies <- string_counts / total_count
  data <- data.frame(
    String = names(string_frequencies),
    Frequency = as.numeric(string_frequencies)
  )
  data <- data[order(-data$Frequency), ]
  return(data)
}
#------------------------------------------------------------------------------------------------
#base on data to get first six strings and order it (data should at least contain cols named String and Frequence )
first_six <-function(data){
  first_six$String <- substr(data$String, start = 1, stop = 6)

  first_six <- first_six %>%
    group_by(String) %>%
    summarise(Sum = sum(Frequency))
  return(first_six)
}
#---------------------------------------------------------------------------------------------------------
#base on DNA_seqid to get matched mRNA sequences(when using the blastx data should conain DNA_seqid and the 
#sequence database named all_mrna)
mrna_seq <- function(blastx){
  blastx <- blastx %>%
    mutate(Sequence = all_mrna$Sequence[match(DNA_seqid,all_mrna$DNA_seqid)])
  return(blastx)
}
#-------------------------------------------------------------------------------------------------
#Smith-Waterman 算法打分
SW <- function(x1,x2){
  W1 <- 2
  H <- matrix(0, nrow = length(x2) + 1, ncol = length(x1) + 1)
  for (i in 2:(length(x2) + 1)) {
    for (j in 2:(length(x1) + 1)) {
      s <- ifelse(x1[j-1] == x2[i-1], 3, -3)
      a <- H[i-1, j-1] + s
      H[i, j] <- max(c(a, H[i-1, j] - W1, H[i, j-1] - W1, 0))
    }
  }
  max_H <- max(H)
  n1 <- which(H == max_H, arr.ind = TRUE)[1]
  m1 <- which(H == max_H, arr.ind = TRUE)[2]
  n <- numeric()
  m <- numeric()
  n <- c(n, n1)
  m <- c(m, m1)
  i <- 1
  while (!( n1 == 1 || m1 == 1 )) {
    if (x1[m1-1] == x2[n1-1]) {
      n <- c(n, n1 - 1)
      m <- c(m, m1 - 1)
    } else {
      a <- H[n1-1, m1-1]
      if (a < H[n1, m1-1] && H[n1, m1-1] >= H[n1-1, m1]) {
        n <- c(n, n1)
        m <- c(m, m1 - 1)
      } else if (H[n1-1, m1] > a && H[n1-1, m1] > H[n1, m1-1]) {
        n <- c(n, n1 - 1)
        m <- c(m, m1)
      } else {
        n <- c(n, n1 - 1)
        m <- c(m, m1 - 1)
      }
    }
    i <- i + 1
    n1 <- n[i]
    m1 <- m[i]
  }
  
  n <- rev(n[-1])
  m <- rev(m[-1])
  
  A <- character()
  B <- character()
  for (i in 1:length(n)) {
    if (H[n[i], m[i]] == 0) {
      A <- c(A, "")
      B <- c(B, "")
    } else {
      A <- c(A, x1[m[i]])
      B <- c(B, x2[n[i]])
      if (i < length(n) && n[i+1] == n[i]) B[i] <- "_"
      if (i < length(n) && m[i+1] == m[i]) A[i] <- "_"
    }
  }
  return()
}