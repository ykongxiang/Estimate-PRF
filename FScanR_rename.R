sort_blastx <- function(blastx_output,
                        mismatch_cutoff  = 5,
                        evalue_cutoff    = 1e-5,
                        frameDist_cutoff = 10){
    blastx <- blastx_output
    colnames(blastx) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qframe", "sframe")
    blastx <- blastx[complete.cases(blastx) & blastx$evalue <= evalue_cutoff & blastx$mismatch <= mismatch_cutoff,,drop=FALSE]
    #筛选统计学有效数据
    blastx_freq <- table(blastx$qseqid)#计算qseqid列的频率
    blastx_freq_cutoff <- blastx_freq[blastx_freq > 1]
    blastx_cutoff <- blastx[blastx$qseqid %in% names(blastx_freq_cutoff),,drop=FALSE]#保留发生移码的命中行（同一mRNA对应多条蛋白质）
    blastx_cutoff_sort <- blastx_cutoff[order(blastx_cutoff$qseqid, blastx_cutoff$sseqid, blastx_cutoff$qstart), , drop = TRUE]
return(blastx_cutoff_sort)
}

blast_output <- read_excel("H:/移码突变/移码数据/all_blastx.xlsx", 
                           col_types = c("text", "text", "numeric", 
                                         "numeric", "numeric", "numeric", 
                                         "numeric", "numeric", "numeric", 
                                         "numeric", "numeric", "numeric", 
                                         "numeric", "numeric"))

grouped_data <- split(blast_output, blast_output[, "qseqid"])
results <- data.frame()

for (group_name in names(grouped_data)) {
  result <- detect(grouped_data[[group_name]])
  results <- rbind(results, result)
}
for (group_name in names(grouped_data)) {
  group_data <- as.data.frame(grouped_data[[group_name]])
  result <- detect(group_data)
  results <- rbind(results, result)
}

detect <- function(blastx_cutoff_sort){
  #确定命中序列不为无
  for (i in 2:nrow(blastx_cutoff_sort)) {
    mismatch <- blastx_cutoff_sort[i,5]
    qseqid <- blastx_cutoff_sort[i,1]
    sseqid <- blastx_cutoff_sort[i,2]
    qstart <- blastx_cutoff_sort[i,7]
    qend <- blastx_cutoff_sort[i,8]
    sstart <- blastx_cutoff_sort[i,9]
    send <- blastx_cutoff_sort[i,10]
    qframe <- blastx_cutoff_sort[i,13]
    qseqid_last <- blastx_cutoff_sort[i-1,1]
    sseqid_last <- blastx_cutoff_sort[i-1,2]
    qstart_last <- blastx_cutoff_sort[i-1,7]
    qend_last <- blastx_cutoff_sort[i-1,8]
    sstart_last <- blastx_cutoff_sort[i-1,9]
    send_last <- blastx_cutoff_sort[i-1,10]
    qframe_last <- blastx_cutoff_sort[i-1,13]
    strand <- ""
    #数据分类 获取
    if (qseqid == qseqid_last & sseqid == sseqid_last & qframe != qframe_last & qframe * qframe_last > 0) {
      if (qframe > 0 & qframe_last > 0) {
        frameStart <- qend_last
        frameEnd <- qstart
        pepStart <- send_last
        pepEnd <- sstart
        plus_strand <- TRUE
        strand <- '+'
      } else if (qframe < 0 & qframe_last < 0) {
        frameStart <- qstart_last
        frameEnd <- qend
        pepStart <- send
        pepEnd <- sstart_last
        plus_strand <- FALSE
        strand <- '-'
      }
      qDist <- frameEnd - frameStart - 1#计算查询序列长度（核苷酸个数）
      sDist <- pepEnd - pepStart#计算命中序列长度（氨基酸个数）
      FS_type <- qDist + (1 - sDist) * 3#通过比较蛋白质序列和核酸序列分别相映的核苷酸个数确定移码类型
      if (abs(qDist) <= frameDist_cutoff & abs(sDist) <= floor(frameDist_cutoff/3)) { 
        prf_sub <- data.frame(as.character(qseqid), frameStart, frameEnd, as.character(sseqid), send_last + 1, sstart, FS_type, strand)
        prf <- rbind(prf, prf_sub)
      }
      prf <- prf[prf$FS_type < 3 & prf$FS_type > -3,,drop=FALSE]#只保留有意义的移码信息
    }
  }
  if (nrow(prf) > 0) {
    colnames(prf) <- c("DNA_seqid", "FS_start", "FS_end", "Pep_seqid", "Pep_FS_start", "Pep_FS_end", "FS_type", "Strand")#DNA序列号 移码起始位点 结束位点 命中蛋白序列号 蛋白
    prf$loci1 = paste(prf$DNA_seqid, prf$FS_start, sep="_")
    prf$loci2 = paste(prf$DNA_seqid, prf$FS_end, sep="_")
    prf$loci3 = paste(prf$Pep_seqid, prf$Pep_FS_start, sep="_")
    prf$loci4 = paste(prf$Pep_seqid, prf$Pep_FS_end, sep="_")
    #整理有效数据 建立新列
    prf = prf[!duplicated(prf$loci1),,drop=FALSE]
    prf = prf[!duplicated(prf$loci2),,drop=FALSE]
    prf = prf[!duplicated(prf$loci3),,drop=FALSE]
    prf = prf[!duplicated(prf$loci4),,drop=FALSE]
    #删除特征数据 重复序列
    prf = prf[,!colnames(prf) %in% c("loci1", "loci2","loci3", "loci4"),drop=FALSE]
    #返回整洁数据
  }
  return(prf)
}