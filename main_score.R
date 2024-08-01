library(readxl)
library(tidyverse)
library(dplyr)
library(openxlsx)
oct_blastx <- read_excel("/Volumes/Samsung_T5/移码突变/移码数据/oct_blastx.xlsx")
colnames(oct_blastx) <- c("DNA_seqid", "Pep_seqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qframe", "sframe")
oct_mrna <- read_excel("/Volumes/Samsung_T5/移码突变/移码数据/oct_mrna.xlsx")
oct_prf <- read_excel("/Volumes/Samsung_T5/移码突变/移码数据/oct_prf.xlsx")
oct_mrna$DNA_seqid <- gsub("Contig.*", "", oct_mrna$DNA_seqid)
va_blastx <- read_excel("/Volumes/Samsung_T5/移码突变/移码数据/va_blastx.xlsx")
colnames(va_blastx) <- c("DNA_seqid", "Pep_seqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qframe", "sframe")
va_mrna <- read_excel("/Volumes/Samsung_T5/移码突变/移码数据/va_mrna.xlsx")
va_prf <- read_excel("/Volumes/Samsung_T5/移码突变/移码数据/va_prf.xlsx")

va_prf_mrna <- merge(va_prf, va_mrna, by = "DNA_seqid", all = TRUE)
oct_prf_mrna <- merge(oct_prf, oct_mrna, by ='DNA_seqid', all = TRUE)
#删除contig标注 使DNA信息纯净

na_indices <- which(is.na(oct_prf_mrna$Sequence))
for (i in seq_along(na_indices)) {
  na_index <- na_indices[i]
  for (j in seq(na_index, nrow(oct_prf_mrna))) {
    if (!is.na(oct_prf_mrna$Sequence[j])) {
      oct_prf_mrna$Sequence[na_index] <- oct_prf_mrna$Sequence[j]
      break 
    }
  }
}
na_rows <- which(is.na(oct_prf_mrna$Strand))
oct_prf_mrna <- oct_prf_mrna %>%
  group_by(Sequence) %>%
  mutate(HasNonNA = any(!is.na(Strand))) %>%
  ungroup()
oct_prf_mrna <- oct_prf_mrna %>%
  filter(!(is.na(Strand) & oct_prf_mrna$HasNonNA))

write.xlsx(va_prf_mrna,'/Volumes/Samsung_T5/移码突变/移码数据/va_prf_mrna.xlsx')
write.xlsx(oct_mrna,'/Volumes/Samsung_T5/移码突变/移码数据/oct_mrna.xlsx')

#——————————————————————————————————————————————————————————————————————————————————————————

if (!require("Biostrings")) install.packages("Biostrings", dependencies = TRUE)
library(Biostrings)


#####调用示例：
a <- set_data_frame(extract_sequences('/Volumes/Samsung_T5/移码突变/移码数据/all_prf_data.xlsx'))

# 调用主函数
a <- extract_sequences('/Volumes/Samsung_T5/移码突变/移码数据/all_prf_data.xlsx')
mrna_file <- 'H:/移码突变/移码数据/all_prf_data.xlsx'

mrna_data1 <- read_excel(mrna_file)
mrna_data1 <- na.omit(mrna_data1)



#————————————————————————————————————————————————————————————————————————————————————————————————————
library(Biostrings)

mRNA <- na.omit(read_xlsx('/Volumes/Samsung_T5/移码突变/移码数据/all_prf_data.xlsx'))
# Extract sequences from the BED file
seqs <- list()
orientation <- list()
for (i in 1:nrow(mRNA)) {
  scf <- mRNA$Sequence[i]
  start <- mRNA$FS_start[i]
  end <- mRNA$FS_end[i]
  name <- paste0(mRNA$DNA_seqid, "_", start, "_", end)
  seq <- substr(as.character(scf), start, end + 1)
  strand <- mRNA$Strand[i]
  if (strand == '-'){
    seq <- reverseComplement(DNAString(seq))
  }
  seqs[i]=seq
}



score_frameshift <- function(blast_df) {
  # 计算两段序列的同源性得分
  blast_df$homology_score <- with(blast_df, (pident / 100) * bitscore)
  
  # 计算间隔距离得分
  blast_df$distance_score <- with(blast_df, (7 - (qend - qstart + 1)) / 3)
  
  # 计算总得分
  blast_df$total_score <- with(blast_df, homology_score * distance_score)
  
  # 根据总得分排序并返回结果
  return(blast_df[order(-total_score), ])
}
mrna_seq <- function(blastx){
  blastx <- blastx %>%
    mutate(Sequence = all_mrna$Sequence[match(DNA_seqid,all_mrna$DNA_seqid)])
  return(blastx)
}

char <- 'AAATAA'

test_data <- read.table(system.file("extdata", "test.tab", package = "FScanR"), header=TRUE, sep='\t')
a <- FScanR(test_data)
a <- mrna_seq(a)
a <- extract_sequences(a)
a_score <- smith_waterman(a)
