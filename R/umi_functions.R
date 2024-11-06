

# Trim UMI primers --------------------------------------------------------
trim_umi_primers <- function(fwd, rev, fwd_out, rev_out, for_primer_seq, rev_primer_seq, 
                             for_umi_length=6, rev_umi_length=6,    
                             n = 1e6, qualityType = "Auto", check_paired = TRUE, id.field = NULL, 
                             max_mismatch=0, id.sep = "\\s", compress =TRUE, quiet=FALSE){
  
  ## iterating through forward and reverse files using fastq streaming
  fF <- ShortRead::FastqStreamer(file.path(fwd), n = n)
  on.exit(close(fF))
  fR <- ShortRead::FastqStreamer(file.path(rev), n = n)
  on.exit(close(fR), add=TRUE)
  
  #Check if there were more than 1 primer per sample
  if(any(length(fwd_out), length(rev_out), length(for_primer_seq), length(rev_primer_seq) >0)){
    multi_primer <- TRUE
  }
  
  # Check if number of F and R primers match the number of outfiles
  if (!all.equal(length(fwd_out), length(rev_out), length(for_primer_seq), length(rev_primer_seq))) {
    stop("fwd_out and rev_out must be the same length as for_primer_seq and rev_primer_seq")
  }
  
  #if(!C_isACGT(primer)) stop("Non-ACGT characters detected in primers")
  
  # Delete output files if they already exist
  c(fwd_out, rev_out) %>%
    purrr::walk(remove_if_exists, quiet=quiet)
  
  first=TRUE
  append <- vector("logical", length= length(fwd_out)) 
  remainderF <- ShortRead::ShortReadQ(); remainderR <- ShortRead::ShortReadQ()
  casava <- "Undetermined" #ID field format
  # Setup read tracking
  inseqs <- 0
  outseqs <- vector("numeric", length= length(fwd_out))
  
  while( TRUE ) {
    suppressWarnings(fqF <- ShortRead::yield(fF, qualityType = qualityType))
    suppressWarnings(fqR <- ShortRead::yield(fR, qualityType = qualityType))
    if(length(fqF) == 0 && length(fqR) == 0) { break } # Stop loop if theres no reads left to process
    
    inseqs <- inseqs + length(fqF)
    
    # Make sure that forward and reverse reads are correctly paired
    # Determine the sequence identifier field. Looks for a single 6-colon field (CASAVA 1.8+ id format)
    if(check_paired) {
      if(first) { 
        if(is.null(id.field)) {
          # or a single 4-colon field (earlier format). Fails if it doesn't find such a field.
          id1 <- as.character(ShortRead::id(fqF)[[1]])
          id.fields <- strsplit(id1, id.sep)[[1]]
          ncolon <- sapply(gregexpr(":", id.fields), length)
          ncoltab <- table(ncolon)
          if(max(ncolon) == 6 && ncoltab["6"] == 1) { # CASAVA 1.8+ format
            casava <- "Current"
            id.field <- which(ncolon == 6)
          } else if (max(ncolon) == 4 && ncoltab["4"] == 1) { # CASAVA <=1.7 format
            casava <- "Old"
            id.field <- which(ncolon == 4)
          } else { # Couldn't unambiguously find the seq id field
            stop("Couldn't automatically detect the sequence identifier field in the fastq id string.")
          }
        }
      } else { # !first
        # Prepend the unmatched sequences from the end of previous chunks
        # Need ShortRead::append or the method is not dispatched properly
        fqF <- ShortRead::append(remainderF, fqF)
        fqR <- ShortRead::append(remainderR, fqR)
      }
    } else { # !check_paired
      if(length(fqF) != length(fqR)) stop("Mismatched forward and reverse sequence files: ", length(fqF), ", ", length(fqR), ".")
    }
    
    # Enforce id matching (ASSUMES SAME ORDERING IN F/R, BUT ALLOWS DIFFERENTIAL MEMBERSHIP)
    # Keep the tail of unmatched sequences (could match next chunk)
    if(check_paired) {
      idsF <- sapply(strsplit(as.character(ShortRead::id(fqF)), id.sep), `[`, id.field)
      idsR <- sapply(strsplit(as.character(ShortRead::id(fqR)), id.sep), `[`, id.field)
      if(casava == "Old") { # Drop the index number/pair identifier (i.e. 1=F, 2=R)
        idsF <- sapply(strsplit(idsF, "#"), `[`, 1)
      }
      lastF <- max(c(0,which(idsF %in% idsR)))
      lastR <- max(c(0,which(idsR %in% idsF)))
      if(lastF < length(fqF)) {
        remainderF <- fqF[(lastF+1):length(fqF)]
      } else {
        remainderF <- ShortRead::ShortReadQ() 
      }
      if(lastR < length(fqR)) {
        remainderR <- fqR[(lastR+1):length(fqR)]
      } else {
        remainderR <- ShortRead::ShortReadQ() 
      }
      fqF <- fqF[idsF %in% idsR]
      fqR <- fqR[idsR %in% idsF]
    }
    
    # Demultiplex and trim each primer
    # Only keep reads where primer is detected
    
    for (p in 1:length(for_primer_seq)){
      # Length of primers
      barlenF <- nchar(for_primer_seq[p])
      barlenR <- nchar(rev_primer_seq[p])
      
      # Check for presence of forward primer after umi
      keepF <- Biostrings::neditStartingAt(
        pattern=Biostrings::DNAString(for_primer_seq[p]),
        subject= IRanges::narrow(sread(fqF), 1+for_umi_length, barlenF+for_umi_length),
        starting.at=1,
        with.indels=FALSE,
        fixed=FALSE ) <= max_mismatch
      
      # Check for presence of reverse primer after umi
      keepR <- Biostrings::neditStartingAt(
        pattern= Biostrings::DNAString(rev_primer_seq[p]),
        subject= IRanges::narrow(sread(fqR), 1+rev_umi_length, barlenR+rev_umi_length),
        starting.at=1,
        with.indels=FALSE,
        fixed=FALSE ) <= max_mismatch
      
      # Only keep reads where forward primer is present in F, and reverse in R
      keep <- keepF & keepR
      
      fqF_primer <- suppressWarnings(ShortRead::ShortReadQ(sread=ShortRead::sread(fqF[keep]), 
                                                           quality=Biostrings::quality(Biostrings::quality(fqF[keep])),
                                                           id=ShortRead::id(fqF[keep])))
      
      fqR_primer <- suppressWarnings(ShortRead::ShortReadQ(sread=ShortRead::sread(fqR[keep]), 
                                                           quality=Biostrings::quality(Biostrings::quality(fqR[keep])),
                                                           id=ShortRead::id(fqR[keep])))
      
      # Extract UMIs from left side
      umiF <- as.character(sread(narrow(fqF_primer, start = 1, end = for_umi_length)))
      umiR <- as.character(sread(narrow(fqR_primer, start = 1, end = rev_umi_length)))
      umi_string <- paste0(umiF,"+",umiR)
      
      # Trim primers from left side
      startF <- max(1, barlenF + for_umi_length + 1, na.rm=TRUE)
      startR <- max(1, barlenR + rev_umi_length + 1, na.rm=TRUE)
      
      fqF_primer_trimmed <- narrow(fqF_primer, start = startF, end = NA)
      fqR_primer_trimmed <- narrow(fqR_primer, start = startR, end = NA)
      
      outseqs[p] <- outseqs[p] + length(fqF_primer_trimmed)    
      
      # Update ids for the outputs, including UMI
      fqF_out <- ShortReadQ(
        sread=ShortRead::sread(fqF_primer_trimmed),
        quality=Biostrings::quality(fqF_primer_trimmed), 
        id=Biostrings::BStringSet(paste0(ShortRead::id(fqF_primer_trimmed),":", umi_string )))
      
      fqR_out <- ShortReadQ(
        sread=ShortRead::sread(fqR_primer_trimmed),
        quality=Biostrings::quality(fqR_primer_trimmed), 
        id=Biostrings::BStringSet(paste0(ShortRead::id(fqR_primer_trimmed),":", umi_string )))
      
      #if(!append[p]) {
      #  ShortRead::writeFastq(fqF_out, fwd_out[p], "w", compress = compress)
      #  ShortRead::writeFastq(fqR_out, rev_out[p], "w", compress = compress)
      #  append[p]=TRUE
      #  first=FALSE
      #} else {
      ShortRead::writeFastq(fqF_out, fwd_out[p], "a", compress = compress)
      ShortRead::writeFastq(fqR_out, rev_out[p], "a", compress = compress)
      #}
    }
  }
  
  if(!quiet) {
    outperc <- purrr::map(outseqs, ~{
      round(.x * 100 / inseqs, 1)
    }) %>%
      unlist()
    outperc <- paste(" (", outperc, "%),", sep="")
    message("Read in ", inseqs, " paired-sequences, output ", paste(" ", outseqs, ",", sep=""), " ", outperc, " primer-trimmed paired-sequences.", sep="")
  }
  
  if(sum(outseqs)==0) {
    message(paste("No reads remaining for:", fwd, "and", rev))
    file.remove(fwd_out)
    file.remove(rev_out)
  }
  out <- data.frame(for_primer_seq = for_primer_seq,
                    rev_primer_seq = rev_primer_seq,
                    trimmed_input = inseqs,
                    trimmed_output = outseqs)
  return(out)
}

# Merge paired end reads --------------------------------------------------
merge_paired_reads <- function(fwd, rev, out, 
                               prefer=1, maxMismatch=0, minOverlap=20, trimOverhang=TRUE,
                               n = 1e6, qualityType = "Auto", quiet=FALSE, compress=TRUE){
  
  ## iterating through forward and reverse files using fastq streaming
  fF <- ShortRead::FastqStreamer(file.path(fwd), n = n)
  on.exit(close(fF))
  fR <- ShortRead::FastqStreamer(file.path(rev), n = n)
  on.exit(close(fR), add=TRUE)
  
  # Delete output files if they already exist
  out %>%
    purrr::walk(remove_if_exists, quiet=quiet)
  
  first=TRUE
  append <- vector("logical", length= length(fwd_out)) 
  remainderF <- ShortRead::ShortReadQ(); remainderR <- ShortRead::ShortReadQ()
  casava <- "Undetermined" #ID field format
  # Setup read tracking
  inseqs <- 0
  outseqs <- vector("numeric", length= length(fwd_out))
  
  while( TRUE ) {
    suppressWarnings(fqF <- ShortRead::yield(fF, qualityType = qualityType))
    suppressWarnings(fqR <- ShortRead::yield(fR, qualityType = qualityType))
    if(length(fqF) == 0 && length(fqR) == 0) { break } # Stop loop if theres no reads left to process
    
    inseqs <- inseqs + length(fqF)
    
    fqF_in <- suppressWarnings(ShortRead::ShortReadQ(sread=ShortRead::sread(fqF), 
                                                     quality=Biostrings::quality(Biostrings::quality(fqF)),
                                                     id=ShortRead::id(fqF)))
    
    fqR_in <- suppressWarnings(ShortRead::ShortReadQ(sread=ShortRead::sread(fqR), 
                                                     quality=Biostrings::quality(Biostrings::quality(fqR)),
                                                     id=ShortRead::id(fqR)))
    # Sequences for alignment
    pairdf <- data.frame(forward=unname(as.character(ShortRead::sread(fqF_in))), 
                         reverse=rc(unname(as.character(ShortRead::sread(fqR_in)))), 
                         id=as.character(ShortRead::id(fqR_in)))
    
    # Unique pairs
    ups <- unique(pairdf[,1:2])
    ups$sequence <- ""
    keep <- !is.na(ups$forward) & !is.na(ups$reverse)
    ups <- ups[keep, ]
    
    # Adjusting align params to prioritize zero-mismatch merges
    tmp <- getDadaOpt(c("MATCH", "MISMATCH", "GAP_PENALTY"))
    if(maxMismatch==0) {
      setDadaOpt(MATCH=1L, MISMATCH=-64L, GAP_PENALTY=-64L)
    } else {
      setDadaOpt(MATCH=1L, MISMATCH=-8L, GAP_PENALTY=-8L)
    }
    
    #Needleman–Wunsch alignment of each UNIQUE pair
    alvecs <- mapply(function(x,y) nwalign(x,y,band=-1), ups$forward, ups$reverse, SIMPLIFY=FALSE)
    
    # Return dada2 parameters to previous after alignment
    setDadaOpt(tmp)
    
    # Evaluate the alignments
    outs <- t(sapply(alvecs, function(x) dada2:::C_eval_pair(x[1], x[2])))
    ups$nmatch <- outs[,1]
    ups$nmismatch <- outs[,2]
    ups$nindel <- outs[,3]
    ups$prefer <- 1 # Change htis to prefer the one with higher quality! would be better if it was done for each base
    maxMismatch=0
    minOverlap=20
    trimOverhang=TRUE
    ups$accept <- (ups$nmatch >= minOverlap) & ((ups$nmismatch + ups$nindel) <= maxMismatch)
    
    # Make the merged sequence
    ups$sequence <- mapply(dada2:::C_pair_consensus, sapply(alvecs,`[`,1), sapply(alvecs,`[`,2), ups$prefer, trimOverhang);
    
    # merge back into paired table
    joint_table <- pairdf %>%
      dplyr::left_join(ups, by = c("forward", "reverse"))
    
    write_csv(joint_table, "read_merging_qc.csv")
    
    accepted_table <- joint_table %>%
      dplyr::filter(accept)
    # Create output file - QUALITY SCORES ARE LOST HERE - NEED BETTER CONSENSUS FUNCTION!
    fq_out <- ShortRead( 
      sread=Biostrings::DNAStringSet(accepted_table$sequence),
      #quality=Biostrings::quality(fqF_in), 
      id=Biostrings::BStringSet(accepted_table$id))
    
    ShortRead::writeFasta(fq_out, out, "a", compress = compress)
    
    outseqs <- outseqs + length(fq_out)   
  }
  
  if(!quiet) {
    outperc <- purrr::map(outseqs, ~{
      round(.x * 100 / inseqs, 1)
    }) %>%
      unlist()
    outperc <- paste(" (", outperc, "%),", sep="")
    message("Read in ", inseqs, " paired-sequences, output ", paste(" ", outseqs, ",", sep=""), " ", outperc, " merged sequences.", sep="")
  }
  
  if(sum(outseqs)==0) {
    message(paste("No reads remaining for:", fwd, "and", rev))
    file.remove(fwd_out)
    file.remove(rev_out)
  }
  out <- data.frame(merged_input = inseqs,
                    merged_output = outseqs)
  return(out)
}

