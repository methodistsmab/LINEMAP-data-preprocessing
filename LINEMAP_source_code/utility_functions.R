
######################3 filter the reads
filter_fastq_stream <- function(
    fq1_in,
    fq1_out,
    fq2_in  = NULL,
    fq2_out = NULL,
    scaffold_seq,
    primer_seq = NULL,
    k = 19,                 # scaffold k-mer 长度：17/19 常用，19更严格更快
    meanQ_min = 34,         # 平均Phred阈值
    drop_N = TRUE,
    chunk = 5e5,            # 每次读多少条（50万比较稳；内存大可1e6）
    verbose_every = 1L
) {
  
  message("===== Step 1: Read filtering =====")
  message("Filtering raw sequencing reads to remove low-quality reads, incomplete alignments, and technical artifacts.")
  
  stopifnot(file.exists(fq1_in))
  if (!is.null(fq2_in)) stopifnot(file.exists(fq2_in))
  if (!is.null(fq2_in)) stopifnot(!is.null(fq2_out))
  
  # ---- 预构建 scaffold k-mer 字典（核心加速）----
  scaffold <- DNAString(scaffold_seq)
  if (length(scaffold) < k) stop("scaffold_seq 长度必须 >= k")
  kmers <- DNAStringSet(vapply(1:(length(scaffold) - k + 1),
                               function(i) as.character(subseq(scaffold, i, width = k)),
                               character(1)))
  kmers <- unique(kmers)
  pd <- PDict(kmers)
  
  primer <- if (!is.null(primer_seq)) DNAString(primer_seq) else NULL
  
  # ---- 清空输出文件（避免 append 叠脏）----
  if (file.exists(fq1_out)) file.remove(fq1_out)
  if (!is.null(fq2_out) && file.exists(fq2_out)) file.remove(fq2_out)
  
  # ---- Streamer ----
  fs1 <- FastqStreamer(fq1_in, n = chunk)
  on.exit(close(fs1), add = TRUE)
  fs2 <- NULL
  if (!is.null(fq2_in)) {
    fs2 <- FastqStreamer(fq2_in, n = chunk)
    on.exit(close(fs2), add = TRUE)
  }
  
  total_in <- 0L
  total_keep <- 0L
  iter <- 0L
  
  repeat {
    iter <- iter + 1L
    
    fq1 <- yield(fs1)
    if (length(fq1) == 0) break
    
    fq2 <- NULL
    if (!is.null(fs2)) {
      fq2 <- yield(fs2)
      if (length(fq2) == 0) stop("R2 读到空，但 R1 还有数据：R1/R2 文件不同步或chunk不一致")
      if (length(fq2) != length(fq1)) stop("R1/R2 chunk 大小不一致：文件可能损坏或不同步")
    }
    
    n <- length(fq1)
    total_in <- total_in + n
    
    # ---- 1) 去 N（可选）----
    keep <- rep(TRUE, n)
    
    if (drop_N) {
      s1 <- sread(fq1)
      keep <- keep & (vcountPattern("N", s1, fixed = TRUE) == 0)
      if (!is.null(fq2)) {
        s2 <- sread(fq2)
        keep <- keep & (vcountPattern("N", s2, fixed = TRUE) == 0)
      }
    } else {
      s1 <- sread(fq1)  # 后面还要用
      if (!is.null(fq2)) s2 <- sread(fq2)
    }
    
    # ---- 2) 平均质量过滤 ----
    q1 <- rowMeans(as(quality(fq1), "matrix"), na.rm = TRUE)
    keep <- keep & (q1 >= meanQ_min)
    
    if (!is.null(fq2)) {
      q2 <- rowMeans(as(quality(fq2), "matrix"), na.rm = TRUE)
      keep <- keep & (q2 >= meanQ_min)
    }
    
    if (!any(keep)) {
      if (verbose_every > 0 && (iter %% verbose_every == 0)) {
        message(sprintf("[iter %d] processed=%d kept=0 (cum kept=%d / %d = %.4f%%)",
                        iter, n, total_keep, total_in, 100*total_keep/total_in))
      }
      next
    }
    
    # ---- 3) scaffold k-mer 命中（任意一个 k-mer 出现即可）----
    # 只在 keep 子集上做，省很多时间

    idx <- which(keep)
    s1k <- sread(fq1)[idx]
    hit_scaf <- colSums(vcountPDict(pd, s1k)) > 0   # ✅这里必须 colSums
    keep[idx] <- hit_scaf
    
    
    if (!any(keep)) {
      if (verbose_every > 0 && (iter %% verbose_every == 0)) {
        message(sprintf("[iter %d] processed=%d kept=0 (cum kept=%d / %d = %.4f%%)",
                        iter, n, total_keep, total_in, 100*total_keep/total_in))
      }
      next
    }
    
    # # ---- 4) primer 命中（可选）----
    # if (!is.null(primer)) {
    #   idx <- which(keep)
    #   s1k <- sread(fq1)[idx]
    #   hit_primer <- vcountPattern(primer, s1k, fixed = TRUE) > 0
    #   keep[idx] <- hit_primer
    # }
    
    # ---- 5) 输出通过的 reads（append 写入 gz）----
    nk <- sum(keep)
    if (nk > 0) {
      writeFastq(fq1[keep], fq1_out, mode = "a", compress = TRUE)
      if (!is.null(fq2)) writeFastq(fq2[keep], fq2_out, mode = "a", compress = TRUE)
      total_keep <- total_keep + nk
    }
    
    if (verbose_every > 0 && (iter %% verbose_every == 0)) {
      message(sprintf("[iter %d] processed=%d kept=%d (cum kept=%d / %d = %.4f%%)",
                      iter, n, nk, total_keep, total_in, 100*total_keep/total_in))
    }
  }
  
  invisible(list(
    total_in = total_in,
    total_keep = total_keep,
    kept_pct = 100 * total_keep / max(1, total_in),
    fq1_out = fq1_out,
    fq2_out = fq2_out
  ))
}


########################################################################################################

append_tsv <- function(df, file) {
  write.table(df, file = file, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = !file.exists(file),
              append = file.exists(file))
}

read.construct = function(id, sequence){
  quality.sequence = 'FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF'
  
  fastqread = matrix('',4*length(id),1)
  for (i in 1 : length(id)){
    fastqread[((i-1)*4+1),1] = id[i]
    fastqread[((i-1)*4+2),1] = sequence[i]
    fastqread[((i-1)*4+3),1] = '+'
    fastqread[((i-1)*4+4),1] = quality.sequence
  }
  return(fastqread)
}



make.config = function(ID,Forward_Reads, Group,GuideRNA, Forward_Primer, Amplicon){
  config = data.frame(ID= '', Barcode = 'Barcode_1', Forward_Reads = '',
                      Reverse_Reads = '', Group = '', Control = 0, guideRNA = '',
                      Forward_Primer = '', Reverse_Primer = '', Direction = 0, Amplicon = '',Donor = '')
  config['ID'] = ID
  config['Forward_Reads'] = Forward_Reads
  config['Group'] = Group
  config['guideRNA'] = GuideRNA
  config['Forward_Primer'] = Forward_Primer
  config['Amplicon'] = Amplicon
  return(config)
  
}


append_alignment_blocks <- function(aes, cell_id, out_file, group_id = NULL,
                                    max_show = Inf, include_score = FALSE) {
  
  pa <- aes@fwdReads[[1]]  # PairwiseAlignments
  pat <- as.character(pattern(pa))
  sub <- as.character(subject(pa))
  sc  <- score(pa)
  
  m <- length(pat)
  
  rc <- tryCatch(aes@readCounts[[1]], error = function(e) NULL)
  if (is.null(rc) || length(rc) != m) rc <- rep(1L, m)
  
  # 按 Count 降序排序（你要的 Count:17 / 6 ... 通常最重要）
  ord <- order(rc, decreasing = TRUE)
  
  n_show <- min(m, max_show)
  ord <- ord[seq_len(n_show)]
  
  con <- file(out_file, open = "a")
  on.exit(close(con), add = TRUE)
  
  # if (!is.null(group_id)) {
  #   cat(sprintf("### Cell: %s  Group: %s ###\n", cell_id, group_id), file = con)
  # }
  
  for (k in seq_along(ord)) {
    j <- ord[k]
    if (include_score) {
      cat(sprintf("ID: %s read_id: %d Count: %d Score: %d\n",
                  cell_id, k, rc[j], sc[j]), file = con)
    } else {
      cat(sprintf("ID: %s read_id: %d Count: %d\n",
                  cell_id, k, rc[j]), file = con)
    }
    cat(pat[j], "\n", file = con)
    cat(sub[j], "\n", file = con)
  }
  
  invisible(m)
}


###############################################

make.align.matrix <- function(cell.read.align) {
  x <- as.character(cell.read.align[, 1])
  
  ## Expect 3-line blocks: header / read / align
  n <- length(x)
  stopifnot(n %% 3 == 0)
  
  header_lines <- x[seq(1, n, by = 3)]
  read_lines   <- x[seq(2, n, by = 3)]
  align_lines  <- x[seq(3, n, by = 3)]
  
  ## Split header lines by whitespace and take first 8 fields
  info_list <- strsplit(header_lines, "\\s+")
  info_mat  <- do.call(rbind, lapply(info_list, function(v) {
    v <- v[v != ""]
    length(v) <- 8          # pad with NA if shorter
    v[1:8]
  }))
  
  ## Build output (same columns as your original function)
  out <- data.frame(
    ID    = info_mat[, 2],
    order = as.numeric(info_mat[, 4]),
    count = as.numeric(info_mat[, 6]),
    score = as.numeric(info_mat[, 8]),
    read  = read_lines,
    align = align_lines,
    stringsAsFactors = FALSE
  )
  
  out
}


fix_leading_align_gaps <- function(read_short,
                                  align_short) {
  if (is.na(read_short) || is.na(align_short)) return(NULL)
  if (nchar(read_short) != nchar(align_short)) return(NULL)
  
  k0 <- nchar(read_short)
  
  # count leading positions where read is '-' and align is not '-'
  k <- 0L
  for (i in seq_len(k0)) {
    a <- substr(align_short,  i, i)
    r <- substr(read_short, i, i)
    if (a == "-" && r != "-") {
      k <- i
    } else {
      break
    }
  }
  
  if (k == 0L) return(NULL)
  
  read_fix  <- substr(read_short,  k + 1L, k0)
  align_fix <- substr(align_short, k + 1L, k0)
  c(read_fix, align_fix, k)
}

##################################

fix_amplican_leading_gap <- function(read_short,
                                     align_short,
                                     ref_seq,
                                     max_prefix = 10) {
  if (is.na(read_short) || is.na(align_short)) return(NULL)
  if (nchar(read_short) != nchar(align_short)) return(NULL)
  
  L <- nchar(align_short)
  if (L < 2) return(NULL)
  
  # split to chars
  a <- strsplit(align_short, "", fixed = TRUE)[[1]]
  r <- strsplit(read_short,  "", fixed = TRUE)[[1]]
  
  # 1) length of leading non-gap prefix (p)
  p <- 0L
  for (i in seq_len(L)) {
    if (a[i] != "-") p <- i else break
  }
  if (p == 0L) return(NULL)              # no leading prefix
  if (p > max_prefix) return(NULL)       # avoid over-aggressive fixes
  if (p >= L) return(NULL)
  
  # 2) length of following gap run (k)
  k <- 0L
  for (i in seq.int(from = p + 1L, to = L)) {
    if (a[i] == "-") k <- k + 1L else break
  }
  
  # 3) insertion sequence in read over the gap-run columns
  ins <- paste0(r[(p + 1L):(p + k)], collapse = "")
  pref <- paste0(a[1L:p], collapse = "")
  
  # need ins at least as long as prefix (it is, but safe)
  if (nchar(ins) < p) return(NULL)
  
  # 4) check: prefix equals tail of insertion
  tail_ins <- substr(ins, nchar(ins) - p + 1L, nchar(ins))
  if (!identical(pref, tail_ins)) return(NULL)
  
  # 5) apply correction: delete those k columns
  keep_idx <- c(seq_len(p), seq.int(from = p + k + 1L, to = L))
  
  read_fix  <- paste0(r[keep_idx], collapse = "")
  align_fix <- paste0(a[keep_idx], collapse = "")
  
  c(read_fix, align_fix, p, k)
  
}


read.separate           <- function(read,
                                    align,
                                    GuideRNA,
                                    primer = "AAGCAGTGGTATCAACGCAGAGTACATGGG",
                                    primer_len = 30) {
  
  GuideRNA = paste0(GuideRNA,'GGG')
  spacer_len <- nchar(GuideRNA)
  
  ## Basic sanity checks
  if (nchar(read) != nchar(align)) return(NULL)
  # if (substr(align, 1, primer_len) != primer) return(NULL)
  
  len <- nchar(align)
  
  ## Find the shortest prefix (starting after primer) whose ungapped length == spacer_len
  ## target.align spans positions [primer_len+1, j]
  target.align <- ""
  
  # find j0: end position where ungapped length == primer_len
  j0 <- NA_integer_
  
  for(j in seq(from = primer_len, to = len)){
    tmp <- substr(align, 1, j)
    if ((nchar(tmp) - stringr::str_count(tmp, "-")) == primer_len) {
      j0 = j
      break
    }
  }
  
  if (is.na(j0)) return(NULL)
  
  for (j in seq(from = j0 + spacer_len, to = len)) {
    tmp <- substr(align, j0+1, j)
    if ((nchar(tmp) - stringr::str_count(tmp, "-")) == spacer_len | j ==len) {
      target.align <- tmp
      break
    }
  }
  
  if (nchar(target.align) == 0) {
    # print('all hgRNA missed') 
    return(NULL)
  }
  
  ## Use the same span to slice `read`
  span_len <- nchar(target.align)
  target.read <- substr(read, j0 + 1, j0 + span_len)
  
  ## Remaining suffix is treated as scaffold (read + align)
  scaffold.align <- substr(align, j0 + span_len + 1, len)
  scaffold.read  <- substr(read,  j0 + span_len + 1, len)
  
  fx <- fix_leading_align_gaps(target.read, target.align)
  if (!is.null(fx)) {
    target.read  <- fx[1]
    target.align <- fx[2]
  }
  
  fx <- fix_amplican_leading_gap(target.read, target.align,
                                 ref_seq = GuideRNA)
  if (!is.null(fx)) {
    target.read  <- fx[1]
    target.align <- fx[2]
  }
  
  
  c(target.read, target.align, scaffold.read, scaffold.align)
}





#######################################

make.mutation.report <- function(target.read,
                                 target.align,
                                 GuideRNA) {
  
  GuideRNA = paste0(GuideRNA,'GGG')
  spacer_len = nchar(GuideRNA) 
  
  ## Return NULL if input is missing
  if (is.na(target.read) || is.na(target.align)) return(NULL)
  
  ## Ensure same length strings (alignment convention)
  if (nchar(target.read) != nchar(target.align)) return(NULL)
  
  a <- strsplit(target.align, "", fixed = TRUE)[[1]]  # reference/alignment string (may contain '-')
  r <- strsplit(target.read,  "", fixed = TRUE)[[1]]  # read string (may contain '-')
  L <- length(a)
  
  ## Base index along the ungapped reference (align) coordinates
  base_idx <- integer(L)
  ref_pos <- 0L
  for (i in seq_len(L)) {
    if (a[i] != "-") {
      ref_pos <- ref_pos + 1L
      base_idx[i] <- ref_pos
    } else {
      base_idx[i] <- ref_pos  # insertion is anchored to previous ref base
    }
  }
  
  ## Helper: get runs of '-' in a character vector
  gap_runs <- function(x) {
    rr <- rle(x == "-")
    ends <- cumsum(rr$lengths)
    starts <- ends - rr$lengths + 1
    which(rr$values) |> (\(idx) {
      if (length(idx) == 0) return(NULL)
      cbind(start = starts[idx], end = ends[idx])
    })()
  }
  
  align_gaps <- gap_runs(a)  # insertions w.r.t reference
  read_gaps  <- gap_runs(r)  # deletions  w.r.t reference
  
  ## Build mutation report
  report <- character(0)
  
  i <- 1L
  while (i <= L) {
    
    ## Case 1: mismatch (both not gaps)
    if (a[i] != "-" && r[i] != "-") {
      if (a[i] != r[i]) {
        report <- c(report, paste0(base_idx[i], " mismatch ", a[i], " ", r[i]))
      }
      i <- i + 1L
      next
    }
    
    ## Case 2: insertion (gap in align/reference, bases in read)
    if (!is.null(align_gaps) && any(align_gaps[, "start"] == i)) {
      row <- which(align_gaps[, "start"] == i)[1]
      s <- align_gaps[row, "start"]; e <- align_gaps[row, "end"]
      ins_ref  <- paste0(a[s:e], collapse = "")
      ins_read <- paste0(r[s:e], collapse = "")
      len_ins  <- e - s + 1L
      anchor   <- base_idx[i]  # anchored to previous ref base
      report <- c(report, paste0(anchor, "<", len_ins, " insert ", ins_ref, " ", ins_read))
      i <- e + 1L
      next
    }
    
    ## Case 3: deletion (gap in read, bases in align/reference)
    if (!is.null(read_gaps) && any(read_gaps[, "start"] == i)) {
      row <- which(read_gaps[, "start"] == i)[1]
      s <- read_gaps[row, "start"]; e <- read_gaps[row, "end"]
      del_ref  <- paste0(a[s:e], collapse = "")
      del_read <- paste0(r[s:e], collapse = "")
      report <- c(report, paste0(base_idx[s], "-", base_idx[e], " del ", del_ref, " ", del_read))
      i <- e + 1L
      next
    }
    
    ## Fallback (should not happen)
    i <- i + 1L
  }
  
  mutation.report <- paste0(paste0(report, ";"), collapse = "")
  
  ## Define a transparent "difference count"
  ## 1) mismatches at non-gap positions
  mismatches <- sum(a != "-" & r != "-" & a != r)
  
  ## 2) insertion length (gaps in align with bases in read)
  ins_len <- if (is.null(align_gaps)) 0L else sum(align_gaps[, "end"] - align_gaps[, "start"] + 1L)
  
  ## 3) deletion length (gaps in read with bases in align)
  del_len <- if (is.null(read_gaps)) 0L else sum(read_gaps[, "end"] - read_gaps[, "start"] + 1L)
  
  number.of.different <- mismatches + ins_len + del_len
  
  ## Handle truncated spacer: if ungapped reference observed < spacer_len, append "notinread"
  ungapped_ref_len <- sum(a != "-")
  if (ungapped_ref_len < spacer_len) {
    missing_from <- ungapped_ref_len + 1L
    missing_to   <- spacer_len
    not_in_read  <- substr(GuideRNA, missing_from, missing_to)
    report2 <- paste0(missing_from, "-", missing_to, " notinread ", not_in_read, " ", strrep("*", missing_to - missing_from + 1L))
    mutation.report <- paste0(mutation.report, report2, ";")
    number.of.different <- number.of.different + (spacer_len - ungapped_ref_len)
  }
  
  list(number.of.different, mutation.report)
}





annotate_mutations <- function(cell.read.matrix,
                               target_id,
                               GuideRNA,
                               verbose_every = 1000) {
  ## ------------------------------------------------------------
  ## Annotate each aligned read with:
  ## 1) spacer vs scaffold split
  ## 2) mutation complexity score and a mutation report string
  ## ------------------------------------------------------------
  
  ## Pre-allocate output columns
  cell.read.matrix$read.short          <- NA_character_
  cell.read.matrix$align.short         <- NA_character_
  cell.read.matrix$scaffold.read       <- NA_character_
  cell.read.matrix$scaffold.align      <- NA_character_
  cell.read.matrix$number.of.different <- NA_real_
  cell.read.matrix$mutation.report     <- NA_character_
  
  for (i in seq_len(nrow(cell.read.matrix))) {
    
    ## Split spacer/scaffold using target-specific spacer length
    res_sep <- read.separate(
      read   = cell.read.matrix[i, "read"],
      align  = cell.read.matrix[i, "align"],
      GuideRNA = GuideRNA
    )
    
    if (!is.null(res_sep) && length(res_sep) == 4) {
      cell.read.matrix$read.short[i]     <- res_sep[1]
      cell.read.matrix$align.short[i]    <- res_sep[2]
      cell.read.matrix$scaffold.read[i]  <- res_sep[3]
      cell.read.matrix$scaffold.align[i] <- res_sep[4]
    } else {
      next
    }
  
    
    ## Generate mutation summary on the spacer region
    res_mut <- make.mutation.report(
      target.read  = cell.read.matrix$read.short[i],
      target.align = cell.read.matrix$align.short[i],
      GuideRNA      = GuideRNA
    )
    
    if (!is.null(res_mut) && length(res_mut) == 2) {
      cell.read.matrix$number.of.different[i] <- res_mut[[1]]
      cell.read.matrix$mutation.report[i]     <- res_mut[[2]]
    }
    
    if (!is.null(verbose_every) && verbose_every > 0 && (i %% verbose_every == 0)) {
      message("Annotated reads: ", i, " / ", nrow(cell.read.matrix),
              " (target ", target_id, ")")
    }
  }
  
  cell.read.matrix
}


##################################################

check.string.different <- function(string1, string2){
  tmp1 = strsplit(string1,'')[[1]]
  tmp2 = strsplit(string2,'')[[1]]
  return(length(which(tmp1 != tmp2)))
}



calculate_mutation_metrics <- function(cell.read.matrix) {
  ## Scaffold length
  cell.read.matrix$len.scaffold <- nchar(cell.read.matrix$scaffold.read)
  
  ## Number of differences in scaffold region
  cell.read.matrix$number.of.different.scaffold <- NA_integer_
  for (i in seq_len(nrow(cell.read.matrix))) {
    cell.read.matrix$number.of.different.scaffold[i] <- 
      check.string.different(
        cell.read.matrix$scaffold.read[i],
        cell.read.matrix$scaffold.align[i]
      )
  }
  
  ## Count number of mutation events
  cell.read.matrix$number.of.mutation <- str_count(
    cell.read.matrix$mutation.report, ";"
  )
  
  ## Compute complexity score
  cell.read.matrix$complexity.score <- 
    cell.read.matrix$len.scaffold -
    cell.read.matrix$number.of.different.scaffold -
    cell.read.matrix$number.of.different -
    cell.read.matrix$number.of.mutation
  
  ## Count specific mutation types from mutation.report (robust)
  cell.read.matrix$number.of.mismatch <-
    stringr::str_count(cell.read.matrix$mutation.report, fixed("mismatch"))
  
  cell.read.matrix$number.of.insertion <-
    stringr::str_count(cell.read.matrix$mutation.report, fixed("insert"))
  
  cell.read.matrix$number.of.deletion <-
    stringr::str_count(cell.read.matrix$mutation.report, fixed("del"))
  
  cell.read.matrix$number.of.notinread <-
    stringr::str_count(cell.read.matrix$mutation.report, fixed("*"))
  
  ## Length-based metrics from aligned sequences
  cell.read.matrix$length.of.insertion <-
    stringr::str_count(cell.read.matrix$align.short, fixed("-"))
  
  cell.read.matrix$length.of.deletion <-
    stringr::str_count(cell.read.matrix$read.short, fixed("-"))
  
  return(cell.read.matrix)
}



#############################################
find.longest.substring <- function(string1, string2.list) {
  vapply(string2.list, function(string2) {
    # 选短的做 a
    if (nchar(string1) <= nchar(string2)) { a <- string1; b <- string2 } else { a <- string2; b <- string1 }
    len1 <- nchar(a)
    
    if (len1 == 0L) return("")
    
    # 从长到短找
    for (L in seq(len1, 1, by = -1)) {
      starts <- 1:(len1 - L + 1)
      for (s in starts) {
        sub <- substr(a, s, s + L - 1)
        if (grepl(sub, b, fixed = TRUE)) return(sub)
      }
    }
    ""  # 理论上不会到这，除非有空字符串
  }, FUN.VALUE = character(1), USE.NAMES = FALSE)
}



add_spacer_similarity_features <- function(cell.read.matrix,
                                           spacer_seq,
                                           prefix = NULL,
                                           verbose_every = 1000) {
  if (is.null(prefix)) prefix <- "spacer"
  
  n <- nrow(cell.read.matrix)
  
  ## pre-allocate
  cell.read.matrix[[paste0("substring.", prefix)]] <- NA_character_
  cell.read.matrix[[paste0("len.substring.", prefix)]] <- NA_integer_
  
  for (i in seq_len(n)) {
    if (!is.null(verbose_every) && i %% verbose_every == 0) {
      message("Computing spacer similarity: ", i, " / ", n)
    }
    
    if (is.na(cell.read.matrix$read.short[i])) next
    
    short.read <- gsub("-", "", cell.read.matrix$read.short[i])
    
    s <- find.longest.substring(spacer_seq, short.read)
    
    cell.read.matrix[[paste0("substring.", prefix)]][i] <- s
    cell.read.matrix[[paste0("len.substring.", prefix)]][i] <- nchar(s)
  }
  
  cell.read.matrix
}


#################################################
compare_cell_read_matrices <- function(mat1,
                                       mat2,
                                       len_gap = 4,
                                       verbose = TRUE) {
  stopifnot(nrow(mat1) == nrow(mat2))
  
  # --- extract vectors ---
  nd1 <- mat1$number.of.different
  nd2 <- mat2$number.of.different
  cs1 <- mat1$complexity.score
  cs2 <- mat2$complexity.score
  ls1 <- mat1$len.substring.spacer
  ls2 <- mat2$len.substring.spacer
  
  ok <- !(is.na(nd1) | is.na(nd2) |
            is.na(cs1) | is.na(cs2) |
            is.na(ls1) | is.na(ls2))
  
  ind1 <- numeric(nrow(mat1))
  ind2 <- numeric(nrow(mat2))
  
  # --- strong rules ---
  strong_A <- ok & (nd1 < nd2) & (cs1 > cs2) & (ls1 > ls2)
  strong_B <- ok & (nd1 > nd2) & (cs1 < cs2) & (ls1 < ls2)
  
  ind1[strong_A] <- 1
  ind2[strong_B] <- 1
  
  # --- weak rules ---
  undecided <- ok & !strong_A & !strong_B
  len_diff  <- ls1 - ls2
  
  ind1[undecided & len_diff >=  len_gap] <- 0.5
  ind2[undecided & len_diff <= -len_gap] <- 0.5
  
  # --- final label ---
  ind <- integer(nrow(mat1))
  ind[ind1 > ind2] <- 1
  ind[ind1 < ind2] <- 2
  
  if (verbose) {
    message("Comparison summary:")
    print(table(ind1, ind2))
    message("Final label counts:")
    print(table(ind))
  }
  
  return(list(
    label = ind,
    ind1 = ind1,
    ind2 = ind2,
    strong_A = strong_A,
    strong_B = strong_B
  ))
}

################################

pick_best_read_per_id <- function(df,
                                  id_col = "ID",
                                  tie_cols = c("complexity.score",
                                               "number.of.mismatch",
                                               "number.of.different",
                                               "count"),
                                  dedup_cols = c("ID", "read.short", "align.short")) {
  
  # 1) 每个 ID：依次按规则过滤（保留最优）
  out <- df %>%
    dplyr::group_by(.data[[id_col]]) %>%
    dplyr::filter(.data[[tie_cols[1]]] == max(.data[[tie_cols[1]]], na.rm = TRUE)) %>%
    dplyr::filter(.data[[tie_cols[2]]] == min(.data[[tie_cols[2]]], na.rm = TRUE)) %>%
    dplyr::filter(.data[[tie_cols[3]]] == min(.data[[tie_cols[3]]], na.rm = TRUE)) %>%
    dplyr::filter(.data[[tie_cols[4]]] == max(.data[[tie_cols[4]]], na.rm = TRUE)) %>%
    dplyr::ungroup()
  
  out <- as.data.frame(out)
  
  # 2) 去重：同一 ID 下完全相同 read/alignment 组合只留一条
  out <- out[!duplicated(out[, dedup_cols]), , drop = FALSE]
  
  # 3) 确保每个 ID 只留一条（如果仍有平局，就保留第一次出现的）
  out <- out[!duplicated(out[, id_col]), , drop = FALSE]
  
  out
}


################################ 

run_step4_cell_level_barcodes <- function(
    base_dir = file.path("result", "mutation_annotation"),
    out_dir  = file.path("result", "cell_barcode"),
    targets  = c("A21", "B25"),
    len_gap  = 4,
    save_csv = TRUE
) {
  stopifnot(length(targets) == 2)
  
  message("===== Step 4: Cell-level barcode identity =====")
  message("Assigning mutation barcodes to individual cells based on read-level mutation annotations.")
  
  # ---- 1) read per-target read-level matrices ----
  cell.read.matrix <- lapply(targets, function(t) {
    readRDS(file.path(base_dir, t, "cell_read_matrix.rds"))
  })
  names(cell.read.matrix) <- targets
  
  m1 <- cell.read.matrix[[targets[1]]]
  m2 <- cell.read.matrix[[targets[2]]]

  print(nrow(m1))
  print(nrow(m2))
  
  # ---- 2) per-read target assignment between the two candidates ----
  res <- compare_cell_read_matrices(m1, m2, len_gap = len_gap)
  ind <- res$label  # 1 = target1 wins, 2 = target2 wins, 0 = undecided
  
  # ---- 3) split reads by assigned target ----
  reads_1 <- m1[ind == 1, , drop = FALSE]
  reads_2 <- m2[ind == 2, , drop = FALSE]
  reads_0_1 <- m1[ind == 0, , drop = FALSE]  # optional: ambiguous reads (A21 view)
  reads_0_2 <- m2[ind == 0, , drop = FALSE]  # optional: ambiguous reads (B25 view)
  
  # ---- 4) per-cell best read selection (cell-level barcode call) ----
  cell_1 <- pick_best_read_per_id(reads_1)
  cell_2 <- pick_best_read_per_id(reads_2)
  
  # ---- 5) build a compact "cell-level barcode table" for downstream Step 5 ----
  # 你可以按需增减字段；这里给出 network 构建最常用的最小集合
  
  cell_barcode <- list(
    A = cell_1,
    B = cell_2
  )
  names(cell_barcode) <- targets
  
  # ---- 6) save outputs into new folder ----
  for (t in targets) {
    dir.create(file.path(out_dir, t), recursive = TRUE, showWarnings = FALSE)
    saveRDS(cell_barcode[[t]], file.path(out_dir, t, "cell_level_barcode.rds"))
  }
  
  saveRDS(list(
    targets = targets,
    assigned_reads = list(
      A = reads_1,
      B = reads_2
    ),
    ambiguous_reads = list(
      A = reads_0_1,
      B = reads_0_2
    ),
    cell_level = cell_barcode,
    assignment = res
  ), file = file.path(out_dir, "step4_result.rds"))
  
  
  # optional: export as csv for quick inspection
  if (save_csv) {
    
    for (t in targets) {
      write.csv(cell_barcode[[t]],
                file.path(out_dir, t, "cell_level_barcode.csv"),
                row.names = FALSE)
    }
    
    # summary counts
    summary_df <- data.frame(
      metric = c("reads_assigned", "cells_called", "reads_ambiguous"),
      A21 = c(nrow(reads_1), nrow(cell_1), nrow(reads_0_1)),
      B25 = c(nrow(reads_2), nrow(cell_2), nrow(reads_0_2))
    )
    colnames(summary_df)[2:3] <- targets
    write.csv(summary_df, file.path(out_dir, "step4_summary.csv"), row.names = FALSE)
  }
  
  invisible(list(
    assigned_reads = setNames(list(reads_1, reads_2), targets),
    cell_level = cell_barcode,
    ambiguous_reads = setNames(list(reads_0_1, reads_0_2), targets),
    assignment = res
  ))
}


################################################# step 5

make.mutation.report.matrix <- function(mutation.report, GuideRNA) {
  
  GuideRNA = paste0(GuideRNA,'GGG')
  spacer_len = nchar(GuideRNA) 
  
  n <- length(mutation.report)
  out <- vector("list", n)
  
  ## row layout:
  ## 1..spacer_len           : reference positions 1..spacer_len
  ## (spacer_len+1)..(2*spacer_len+1) : insertion anchor positions 0..spacer_len
  nrow_mat <- 2 * spacer_len + 1
  
  for (i in seq_len(n)) {
    
    M <- matrix(NA_character_, nrow = nrow_mat, ncol = 5)
    colnames(M) <- c("mutation.type", "location", "number.of.base", "original", "actual")
    
    ## position labels
    M[, "location"] <- c(seq_len(spacer_len), 0:spacer_len)
    M[, "number.of.base"] <- 0
    
    rep_i <- mutation.report[i]
    
    if (is.na(rep_i) || rep_i == "") {
      out[[i]] <- M
      next
    }
    
    events <- strsplit(rep_i, ";", fixed = TRUE)[[1]]
    events <- stringr::str_trim(events)
    events <- events[events != ""]
    
    for (ev in events) {
      tok <- strsplit(ev, " +")[[1]]
      if (length(tok) < 2) next
      
      ## insert: "<anchor><len insert <orig> <actual>"
      if (tok[2] == "insert") {
        a <- strsplit(tok[1], "<", fixed = TRUE)[[1]]
        anchor <- suppressWarnings(as.integer(a[1]))
        inslen <- suppressWarnings(as.integer(a[2]))
        if (is.na(anchor)) anchor <- 0L
        if (is.na(inslen)) next
        
        row <- spacer_len + 1 + anchor  # anchor 0..spacer_len
        if (row >= 1 && row <= nrow_mat) {
          M[row, "mutation.type"] <- "insert"
          M[row, "number.of.base"] <- inslen
          M[row, "original"] <- tok[3]
          M[row, "actual"]   <- tok[4]
        }
        
        ## miss / notinread: "a-b miss ..." or "a-b notinread ..."
      } else if (tok[2] %in% c("del", "notinread")) {
        ab <- suppressWarnings(as.integer(strsplit(tok[1], "-", fixed = TRUE)[[1]]))
        if (length(ab) != 2 || anyNA(ab)) next
        a <- ab[1]; b <- ab[2]
        a <- max(1L, a); b <- min(spacer_len, b)
        if (a > b) next
        
        type <- if (tok[2] == "del") "deletion" else "miss*"
        orig <- tok[3]
        for (k in a:b) {
          M[k, "mutation.type"] <- type
          M[k, "number.of.base"] <- 1
          M[k, "original"] <- substr(orig, k - a + 1, k - a + 1)
          M[k, "actual"]   <- if (type == "deletion") "-" else "*"
        }
        
        ## mismatch: "pos mismatch A T"
      } else if (tok[2] == "mismatch") {
        pos <- suppressWarnings(as.integer(tok[1]))
        if (!is.na(pos) && pos >= 1 && pos <= spacer_len) {
          M[pos, "mutation.type"] <- "mismatch"
          M[pos, "number.of.base"] <- 1
          M[pos, "original"] <- tok[3]
          M[pos, "actual"]   <- tok[4]
        }
      }
    }
    
    out[[i]] <- M
  }
  
  out
}


make.mutation.type.matrix <- function(mutation.report.matrix.list,
                                      ncol_out = NULL,
                                      mismatch_score = 0.5) {
  stopifnot(is.list(mutation.report.matrix.list))
  n <- length(mutation.report.matrix.list)
  if (n == 0L) return(matrix(numeric(0), nrow = 0, ncol = 0))
  
  # 自动确定列数：若没给，就取所有元素中最大行数（更稳）
  if (is.null(ncol_out)) {
    ncol_out <- max(vapply(mutation.report.matrix.list, nrow, integer(1)), na.rm = TRUE)
  }
  
  # 先创建字符矩阵，默认 NA
  type_chr <- matrix(NA_character_, nrow = n, ncol = ncol_out)
  
  # 填充每行：取每个 report matrix 的第 1 列
  for (i in seq_len(n)) {
    x <- mutation.report.matrix.list[[i]]
    if (is.null(x) || nrow(x) == 0L) next
    L <- min(nrow(x), ncol_out)
    type_chr[i, seq_len(L)] <- as.character(x[seq_len(L), 1])
  }
  
  # 映射到数值矩阵
  out <- matrix(0, nrow = n, ncol = ncol_out)  # 默认 0
  
  out[type_chr == "insert"]   <- 1
  out[type_chr %in% c('miss*', 'miss','del', "deletion")]     <- 1
  out[type_chr == "mismatch"] <- mismatch_score
  
  # 其余（包括 NA / 其他字符串）保持 0
  out
}

make.mutation.type.matrix.plus <- function(mutation.report.matrix.list,
                                           ncol_out = NULL,
                                           fill_na = 0) {
  stopifnot(is.list(mutation.report.matrix.list))
  n <- length(mutation.report.matrix.list)
  if (n == 0L) return(matrix(numeric(0), nrow = 0, ncol = 0))
  
  # 自动列数：取所有元素的最大行数（更稳）
  if (is.null(ncol_out)) {
    ncol_out <- max(vapply(mutation.report.matrix.list, nrow, integer(1)), na.rm = TRUE)
  }
  
  # 先建 numeric 矩阵，默认填 fill_na（通常用 0）
  out <- matrix(fill_na, nrow = n, ncol = ncol_out)
  
  for (i in seq_len(n)) {
    x <- mutation.report.matrix.list[[i]]
    if (is.null(x) || nrow(x) == 0L) next
    L <- min(nrow(x), ncol_out)
    
    v <- suppressWarnings(as.numeric(x[seq_len(L), 3]))
    if (!is.null(fill_na)) v[is.na(v)] <- fill_na
    
    out[i, seq_len(L)] <- v
  }
  
  out
}


make.network.direction <- function(data, GuideRNA) {
  
  # ---- containers (avoid rbind growth) ----
  edges <- vector("list", nrow(data))
  e <- 0L
  
  # ---- precompute matrices once ----
  mutation.report.matrix <- make.mutation.report.matrix(data[, "mutation.report"], GuideRNA)
  typeM     <- make.mutation.type.matrix(mutation.report.matrix)
  typePlusM <- make.mutation.type.matrix.plus(mutation.report.matrix)
  
  n <- nrow(data)
  p1 <- ncol(typeM)
  p2 <- ncol(typePlusM)
  
  # helper: j dominates i? (j is ancestor/subset of i)
  # i.e., type[j,] <= type[i,] and typePlus[j,] <= typePlus[i,]
  is_ancestor_of <- function(j_mat, i_vec, dims) {
    # dims = number of columns in matrix
    rowSums(sweep(j_mat, 2, i_vec, `<=`)) == dims
  }
  
  # helper: rows equal?
  is_equal_to <- function(j_mat, i_vec, dims) {
    rowSums(sweep(j_mat, 2, i_vec, `==`)) == dims
  }
  
  # helper: for same-type case, check insert “containment” you implemented
  # (exactly mirror your grepl check)
  insert_contained <- function(i, j) {
    idx <- which(mutation.report.matrix[[j]][, 1] == "insert")
    if (length(idx) == 0L) return(TRUE)  # no insert to compare
    
    j_ins <- mutation.report.matrix[[j]][idx, 5]
    i_ins <- mutation.report.matrix[[i]][idx, 5]
    
    # if i doesn't have corresponding entries, fail safe
    if (length(i_ins) != length(j_ins)) return(FALSE)
    
    for (kk in seq_along(idx)) {
      if (!grepl(j_ins[kk], i_ins[kk], fixed = TRUE)) return(FALSE)
    }
    TRUE
  }
  
  # helper: keep maximal elements among candidates:
  # keep j such that there is no other l in candidates with l >= j (and not equal),
  # i.e., j is not strictly dominated by another candidate.
  keep_maximal <- function(cand_idx) {
    if (length(cand_idx) <= 1L) return(cand_idx)
    
    keep <- rep(TRUE, length(cand_idx))
    for (a in seq_along(cand_idx)) {
      if (!keep[a]) next
      ja <- cand_idx[a]
      
      # compare to all others
      for (b in seq_along(cand_idx)) {
        if (a == b || !keep[a]) next
        jb <- cand_idx[b]
        
        # jb dominates ja and is not equal => drop ja
        if ( all(typeM[jb, ]     >= typeM[ja, ]) &&
             all(typePlusM[jb, ] >= typePlusM[ja, ]) &&
             (any(typeM[jb, ]     > typeM[ja, ]) || any(typePlusM[jb, ] > typePlusM[ja, ])) ) {
          keep[a] <- FALSE
        }
      }
    }
    cand_idx[keep]
  }
  
  # ---- main loop ----
  for (i in seq_len(n)) {
    
    # all candidates except i
    others <- setdiff(seq_len(n), i)
    
    # dominance filter: ancestor candidates
    cand1 <- is_ancestor_of(typeM[others, , drop = FALSE], typeM[i, ], p1)
    cand2 <- is_ancestor_of(typePlusM[others, , drop = FALSE], typePlusM[i, ], p2)
    cand  <- others[cand1 & cand2]
    
    if (length(cand) == 0L) next
    
    # split into A: type different, B: type same but typePlus different
    eq_type     <- is_equal_to(typeM[cand, , drop = FALSE], typeM[i, ], p1)
    eq_typePlus <- is_equal_to(typePlusM[cand, , drop = FALSE], typePlusM[i, ], p2)
    
    A <- cand[!eq_type]                        # your tmp.ind
    B <- cand[ eq_type & !eq_typePlus ]        # your tmp.ind1 candidate set
    
    # for B, apply your insert containment filter
    if (length(B) > 0L) {
      B <- B[ vapply(B, function(j) insert_contained(i, j), logical(1)) ]
    }
    
    if (length(B) > 0L) {
      # choose B with max rowSums(typePlus) (your score1 logic)
      scoreB <- rowSums(typePlusM[B, , drop = FALSE])
      bestB  <- B[scoreB == max(scoreB)]
      
      for (j in bestB) {
        e <- e + 1L
        edges[[e]] <- data.frame(
          from  = data[j, "mutation.report"],
          to    = data[i, "mutation.report"],
          score = data[i, "n"],
          stringsAsFactors = FALSE
        )
      }
      
    } else if (length(A) > 0L) {
      
      # keep only maximal ancestors among A (your tmp.tmp.ind elimination)
      maxA <- keep_maximal(A)
      
      for (j in maxA) {
        e <- e + 1L
        edges[[e]] <- data.frame(
          from  = data[j, "mutation.report"],
          to    = data[i, "mutation.report"],
          score = data[i, "n"],
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  if (e == 0L) {
    return(data.frame(from = character(0), to = character(0), score = numeric(0)))
  }
  do.call(rbind, edges[seq_len(e)])
}



###########################
 
plot_mutation_network <- function(node,

                                  edges,

                                  node_id_col = NULL,

                                  node_size_col = "cell.number",

                                  node_label_col = NULL,

                                  layout = "fr",

                                  label_size_threshold = 30,

                                  seed = 1) {

  node  <- as.data.frame(node)

  edges <- as.data.frame(edges)

  if (is.null(node_id_col)) {

    node_id_col <- if ("node" %in% names(node)) "node" else names(node)[1]

  }

  if (!(node_id_col %in% names(node))) stop("node_id_col not found in node")

  node$id <- as.character(node[[node_id_col]])

  infer_edge_cols <- function(df) {

    cand_pairs <- list(

      c("from", "to"),

      c("parent", "child"),

      c("source", "target"),

      c("src", "dst"),

      c("u", "v")

    )

    for (p in cand_pairs) {

      if (all(p %in% names(df))) return(p)

    }

    names(df)[1:2]

  }

  ec <- infer_edge_cols(edges)

  from_col <- ec[1]; to_col <- ec[2]

  edges$from <- as.character(edges[[from_col]])

  edges$to   <- as.character(edges[[to_col]])

  keep <- edges$from %in% node$id & edges$to %in% node$id

  edges2 <- edges[keep, c("from", "to"), drop = FALSE]

  g <- igraph::graph_from_data_frame(edges2, directed = TRUE, vertices = node)

  # node size

  if (node_size_col %in% names(node)) {

    node[[node_size_col]] <- as.numeric(node[[node_size_col]])

    sz <- node[[node_size_col]]

    sz[is.na(sz)] <- 0

    denom <- max(sz, na.rm = TRUE)

    if (!is.finite(denom) || denom <= 0) denom <- 1

    vsize <- 5 + 20 * (sz / denom)

  } else {

    vsize <- 8

  }

  # labels: top N by size

  if (node_size_col %in% names(node)) {

    sz <- as.numeric(node[[node_size_col]])

    show_ids <- node$id[!is.na(sz) & sz > label_size_threshold]

  } else {

    show_ids <- node$id

  }

  vlab <- rep(NA_character_, igraph::vcount(g))

  names(vlab) <- igraph::V(g)$name

  label_col_to_use <- if (!is.null(node_label_col) && node_label_col %in% names(node)) {

    node_label_col

  } else {

    node_id_col

  }

  lab_map <- setNames(as.character(node[[label_col_to_use]]), node$id)

  vlab[names(vlab) %in% show_ids] <- lab_map[names(vlab)[names(vlab) %in% show_ids]]

  set.seed(seed)

  lay <- switch(tolower(layout),

                "fr"  = igraph::layout_with_fr(g),

                "kk"  = igraph::layout_with_kk(g),

                "dh"  = igraph::layout_with_dh(g),

                "lgl" = igraph::layout_with_lgl(g),

                igraph::layout_with_fr(g))

  plot(g,

       layout = lay,

       vertex.size = vsize,

       vertex.label = vlab,

       vertex.label.cex = 0.6,

       edge.arrow.size = 0.3,

       edge.curved = 0.1)

  invisible(g)

}
 
 
 
 