#!/usr/bin/env Rscript

# Processes a directory of fa files. Proceeds to extend them, blast them, 
# align them, cluster the results, and find the edges of the resulting consensus,
# if available. Then iterates the process, resulting in a folder "final" with the
# consensus obtained.

# Requirements: bedtools, samtools, blast, and a modified version of make_align_from_blast.sh
# from Anna Protasio original TE_Manannot scripts.

# The resulting consensus should afterwards be classified, and their TE_AID plots
# generated to perform the final revision of them.

# Uncomment this to run in slurm if you need to install packages in your own folder.

local({r <- getOption("repos")
r["CRAN"] <- "https://cran.r-project.org"
options(repos=r)
})
.libPaths( c( "~/R/x86_64-pc-linux-gnu-library/3.6" , .libPaths() ) )

# To make the sampling reproducible.
set.seed(1234) 

# Get path of the script, to use as base path also for the other in bash.
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)

# Libraries
if (!require("data.table", quietly = TRUE))
  install.packages('data.table')
if (!require("stringi", quietly = TRUE))
  install.packages('stringi')
if (!require("getopt", quietly = TRUE))
  install.packages('getopt')
if (!require("dbscan", quietly = TRUE))
  install.packages('dbscan')
if (!require("ape", quietly = TRUE))
  install.packages("ape")
if (!require("bioseq", quietly = TRUE)) 
  install.packages('bioseq') 
if (!require("keypress", quietly = TRUE)) 
  install.packages('keypress')           

# Get parameters
spec = matrix(
  c(
    'input',                'i', 1, "character",
    'output',               'o', 1, "character",
    'genome',               'g', 1, "character",
    'max_sequences',        'x', 1, "integer",
    'help',                 'h', 0, "logical",
    'min_plurality',        'p', 1, "integer",
    'min_saturation',       's', 1, "integer",
    'cluster_factor',       'c', 1, "integer",
    'min_cluster',          'm', 1, "integer",
    'group_outliers',       'l', 0, "logical",
    'interactive',          'n', 0, "logical",
    'top_rounds',           'r', 1, "integer",
    'min_seqs_per_cluster', 'q', 1, "integer",
    'extend',               'e', 1, "integer",
    'identity',             'd', 1, "double"
  ),
  byrow = TRUE,
  ncol = 4
)
opt = getopt(spec)

if (!is.null(opt$help)) {
  cat(getopt(spec, usage = TRUE))
  q(status = 1)
}
if (is.null(opt$input)) {
  print("-input missing")
  q(status = 1)
}
if (is.null(opt$genome)) {
  print("-genome missing")
  q(status = 1)
}
if (is.null(opt$output)) {
  print("-output missing")
  q(status = 1)
}

# Max number of sequences to sample
if ( is.null(opt$max_sequences)) {
  opt$max_sequences = 200 
}

# Min percentage of bases forming the plurality.
if ( is.null(opt$min_plurality )) {
  opt$min_plurality  = 40 
}

# Min value of saturation corrected by density to find the edges.
if ( is.null(opt$min_saturation)) {
  opt$min_saturation <- 80 
}

# Min cluster size has to be at least >= max(# seqs / cluster_factor, min_cluster / 2)
if ( is.null(opt$cluster_factor)) {
  opt$cluster_factor <- 10 
}

# Min number of total seqs to try to cluster.
if ( is.null(opt$min_cluster)) {
  opt$min_cluster <- 8 
}

# Threshold to define if the edges found are too close to the end of the alignment
# and hence whether it needs to be extended or not.
if ( is.null(opt$end_threshold)) {
  opt$end_threshold <- 100 
}

# If FALSE tries to remove outliers from the clusters and produces a cluster with them.
# Consensus might be a little better, but we lose sequences for the alignments.
if ( is.null(opt$group_outliers)) {
  opt$group_outliers <- TRUE 
}   

# To run in an interactive way. Generates plots and shows each alignment in Aliview
if ( is.null(opt$interactive)) {
  opt$interactive <- FALSE 
}   

# Max number of rounds to iterate.
if ( is.null(opt$top_rounds)) {
  opt$top_rounds <- 6 
} 

# Bases to extend
if ( is.null(opt$extend)) {
  opt$extend <- 1000 
} 

# Percentage of identity for blast.
if ( is.null(opt$identity)) {
  opt$identity <- 0.9
} 

# Read a maf file, store it on a data.table and include on the first column 
# its name to carry it on the rest of functions.

read_maf <- function(maf_file) {
  name <- gsub(".*/", "", maf_file)
  name <- gsub("\\..*$", "", name)
  print(paste("[+++] Processing:",name))
  # Avoid processing a maf already processed.
  if ((file.size(maf_file) != 0L) & 
      !(file.exists(paste0(opt$output, "/", name, "_alt_1.con.fa")) |
      (file.exists(paste0(opt$output, "/", name, "_alt_0.con.fa"))) |
      (file.exists(paste0(opt$output, "/", name, ".con.fa"))))) { 
    maf_seqs <-
    fread(
      cmd = paste(
        "awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) { print \"\\n\"$0} else {printf $0}}}' ",
        maf_file,
        "| awk 'NR%2==0 {print $1}'"
      ),
      header = F,
      col.names = c("seqs")
    )
    
  # Convert string to chars, one for each row. Rows are bases, cols are sequences.
  data_m <- stri_extract_all(maf_seqs$seqs, regex = ".")
  r <- data.table(rbindlist(list(data_m)))
  
  # Store the name of the maf as name of the first column to pass it between functions.
  colnames(r)[1] <- name
  return(r)
  }
  else {
    return(NULL)
  }   
}

# Divide the maf seqs in clusters using dbscan

seq_clus <- function(maf) {
  if (is.null(maf)) {
    return(NULL)
  }
  else {
    name <- colnames(maf)[1]
    d <- dist.dna(as.DNAbin(t(maf)), pairwise.deletion = TRUE)
    kd <- mean(d) # Average Kimura distance of all the maf 
    print(paste("[+++] Kimura Distance of full alignment:",kd))
    
    if (ncol(maf) > opt$min_cluster) { # Don't try to cluster if less than min_cluster sequences
      print(paste("[+++] Proceeding to cluster ", name))
      if (ncol(maf) > opt$max_sequences) { # If more than max_sequences get a sample of them
        maf <- sample(maf, opt$max_sequences)
        # Recalculate d if we are using a sample.
        d <- dist.dna(as.DNAbin(t(maf)), pairwise.deletion = TRUE)
        print(paste("[+++] Sampling",opt$max_sequences, "sequences from the alignment."))
      }
      num_seqs <- ncol(maf)
      # Calculating mp value for DBSCAN based on the clustering factor and the minimum size of cluster
      mp <- ncol(maf) %/% opt$cluster_factor
      mp <- max(mp,opt$min_cluster%/%2) 

      print(paste("[+++] Mpoints = ", mp))
      
      # Calculating the eps for DBSCAN as the one that produces more "stable" clusters.
      # We define "stable as clusters that remain for a distance of at least 0.04 (3 eps "steps")
      # The idea is to ignore clusters that are happen just at specific values as outliers.

      points <- c()
      # Convert NAs to maximum distance.
      d[is.na(d)] <- 1 
      # min eps granularity of 0.2. Could be made a parameter. 
      for (eps in (seq(0.0, 1, .02))) { 
        candidate <- sum(unique(dbscan(d, eps = eps, minPts = mp)$cluster) != 0)
        points <- append(points,candidate)
      }
      
      print(points)
      # 3 means the cluster must extend at least over 3 eps candidate steps ("stable cluster") 
      # Could be made a parameter.
      nclusters <- max(as.numeric(names(table(points)[table(points)>3]))) 
      if (!opt$group_outliers) {
        eps <- (match(nclusters,points)-1)*.02 
      } else {  # If this version is used the 0 clusters are probably just junk and should be discarded. WE CAN DO IT HERE (new option to keep them)
        eps <- (51-(match(nclusters,rev(points))))*.02 # highest point with the top clusters (We get more reads on the clusters)
      }

      clusters <- dbscan(d, eps = eps, minPts = mp)$cluster
      print(paste("[+++] Eps:",eps))
      if (length(table(points))>1) { # Show clusters if there are more than one (or there are outliers)
          print("[+++] Detected clusters:")
          print(clusters)
      } else {
        print("[+++] No clusters detected.")
      }

      mafs <- list()
      for (i in unique(clusters)) {
        if (i != 0) {
          mafs[[i]] <- data.table(maf[, c(clusters == i), with = FALSE])
          colnames(mafs[[i]])[1] <- paste0(name, "_alt_", i)
        }
      }  
      # The 0 case needs to be handled different as the list is 1 indexed, 
      # to keep the order. But only if more than 1 seq 
      # (it does not follow the MinPoints) And only if more than 1 cluster (LTRs)
      # eps = 0 is the case in which no cluster is stable enough. Can happen on small
      # alignments. In that case we keep the whole alignment.
      if ((0 %in% unique(clusters)) & (sum(clusters==0)>1) & max(clusters>1) | eps == 0) { 
        mafs[[max(unique(clusters))+1]] <- data.table(maf[, c(clusters == 0), with = FALSE])
        colnames(mafs[[max(unique(clusters))+1]])[1] <- paste0(name, "_alt_0")
      }  
    }  
    else {
      print(paste("[+++] No need to cluster ", name))
      mafs <- list() 
      mafs[[1]] <- data.table(maf)
      colnames(mafs[[1]])[1] <- paste0(name, "_alt_1")
    }
    # The average KD is passed as the name of the first cluster.
    names(mafs)[1] <- round(kd,3) 
    return(mafs) 
  }
}

# Processes a list of maf alignments, there can be several for the same maf file 
# if more than one cluster have been found.
# As earlier, columns are sequences, rows are positions (bases)

process_maf <- function(list_r) {
  if (is.null(list_r)) {
    return(NULL)
  }
  else {
  stats <- list()
  kd <- names(list_r)[1]
  for (i in seq(1, length(list_r), 1)) {
    r <- list_r[[i]]
    seqs <- ncol(r)
    cols <- colnames(r)
    name <- colnames(r)[1]
    # Calculate number of empty places and create an index of rows (bases)
    r[, res := rowSums(r == "-")]
    r[, ID := .I]
    # remove fully empty positions (can happen after clustering)
    r <- r[res != seqs] 
    l <- nrow(r)
    r[, resA := rowSums(r == "a")]
    r[, resC := rowSums(r == "c")]
    r[, resG := rowSums(r == "g")]
    r[, resT := rowSums(r == "t")]
    r$mb <- apply(r[, c("resA", "resC", "resG", "resT")], 1, max) # Majority base
    # "Saturation" is the ratio of bases forming the consensus compared to the total
    # number of sequences.
    r[,satur := mb/seqs]
    # Remove empty positions and positions with just one base
    r <- r[res+1 < seqs] 
    # Get majority base.
    r$base <- ""  
    r[r$resA==r$mb]$base <- "a"
    r[r$resC==r$mb]$base <- "c"
    r[r$resG==r$mb]$base <- "g"
    r[r$resT==r$mb]$base <- "t"
    
    # Define edges
    # lm and rm are equivalent to saturation, but requiring that the previous (lm) or
    # next (rm) position forms the consensus too. 
    r[,lm:=0]
    for (i in seq(2,nrow(r),1)) {  
      r[i]$lm <- sum(r[i,1:seqs]==r[i]$base & r[i-1,1:seqs]==r[i-1]$base, na.rm = TRUE)/seqs
    }
    r[,rm:=shift(lm,-1)]
    m1 <- mean(r[lm>0.5]$lm)
    sd1 <- sd(r[lm>0.5]$lm)
    cutpoint1 <- m1-sd1*2
    # Select edges (lm and rm could be unified as the max of both)
    min_r <- r[lm >= cutpoint1 | rm >= cutpoint1] 
    # Correct edges for low number of seqs (need more confidence -> more consecutive bases)
    # Take the edges that match too this rule
    confidence_needed <- 20%/%seqs+1
    v <- min_r[1:confidence_needed]$rm >= cutpoint1
    while (sum(v)<confidence_needed) {
      min_r <- min_r[(confidence_needed-(match(FALSE,rev(v)))+2):nrow(min_r)]
      v <- min_r[1:confidence_needed]$rm >= cutpoint1
    }
    v <- min_r[(nrow(min_r)-confidence_needed+1):nrow(min_r)]$lm >= cutpoint1
    while (sum(v)<confidence_needed) {
      min_r <- min_r[1:(nrow(min_r)-confidence_needed+(match(FALSE,v))-1)]
      v <- min_r[(nrow(min_r)-confidence_needed+1):nrow(min_r)]$lm >= cutpoint1
    }
    
    # Get the final set of bases that we will use to build the consensus.
    minmin_r <- r[ID>= min(min_r$ID) & ID<=max(min_r$ID)]
    # kimura distance of the final alignment BUT STILL WITH INTERNAL LOW COUNT BASES THAT WILL NOT BE PART OF THE CONSENSUS.
    kdc <- mean(dist.dna(as.DNAbin(t(minmin_r[,1:seqs]))), pairwise.deletion = TRUE) 
    print(paste0("[+++] Kimura distance for ", name, ": ",kdc))
    
    # Build the consensus with a simple plurality model. It can be made as sophisticated as needed.
    minmin_r <- minmin_r[satur>(opt$min_plurality/100)]
    s <- paste0(minmin_r$base, collapse = "")
    # Alternative 
    #    s <- paste0(seqinr::consensus(t(as.matrix(minmin_r[,1:seqs])), method = "IUPAC"),collapse="")
    
    # Draw plots and open cluster alignments in AliView when in interactive mode.
    if (opt$interactive) { 
      png(paste0("plot_",name,".png"))
      plot(r$ID,r$lm, main = name , xlab=paste("left =", min(minmin_r$ID)," - right =",max(minmin_r$ID)))
      clip(1,min(minmin_r$ID),0,1)
      abline(h= cutpoint1, col="blue")
      clip(max(minmin_r$ID),max(r$ID),0,1)
      abline(h= cutpoint1, col="blue")
      clip(0,max(r$ID),0,1)
      abline(v =min(minmin_r$ID), col="red")
      abline(v =max(minmin_r$ID), col="green")
      clip(min(minmin_r$ID),max(minmin_r$ID),0,1)
      abline(h =opt$min_plurality/100, col="black")
      dev.off()
      aliview(as.DNAbin(t(minmin_r[,1:(ncol(minmin_r)-9)])))
      # aliview(as.DNAbin(t(r[,1:(ncol(minmin_r)-9)])))
      print("Press a key to continue")
  #    xxx <- keypress()
   #   invisible(readline(prompt="Press [enter] to continue"))
    }

    # Generating stats for each processed consensus. 
    min_l <- nrow(minmin_r)
    if (is.na(kd)) {kd <- Inf}
    stats <- rbindlist(list(stats,  
                            data.table(seq = name, 
                                              seqs_for_consensus = seqs,
                                              alig_size = l,
                                              cons_size = nchar(s),  
                                              start_con = minmin_r[1]$ID,
                                              end_con= minmin_r[min_l]$ID,
                                              surplus_left = nrow(r[ID<minmin_r[1]$ID]),
                                              surplus_right = nrow(r[ID>minmin_r[min_l]$ID]),
                                              dust_left = 1-(sum(r[ID<minmin_r[1]$ID]$res)/(seqs*nrow(r[ID<minmin_r[1]$ID]))),  
                                              dust_right = 1-(sum(r[ID>minmin_r[min_l]$ID]$res)/(seqs*nrow(r[ID>minmin_r[min_l]$ID]))),
                                              maf_kd = kd,
                                              cluster_kd = kdc,
                                              end_l = nrow(r[ID<minmin_r[1]$ID]) - (sum(r[ID<minmin_r[1]$ID]$res)/seqs) > opt$end_threshold,
                                              end_r = nrow(r[ID>minmin_r[min_l]$ID]) -(sum(r[ID>minmin_r[min_l]$ID]$res)/seqs) > opt$end_threshold)))
    # Write fasta for final result
    if (stats[nrow(stats)]$end_l & stats[nrow(stats)]$end_r) {
      fileConn <- file(paste0(opt$output, "/final/", name, ".con.fa"))
      writeLines(c(paste0(">", name), s), fileConn)
      close(fileConn)
    } else {
      end <- FALSE
      fileConn <- file(paste0(opt$output,"/to_extend_",round, "/", name, ".con.fa"))
      writeLines(c(paste0(">", name), s), fileConn)
      close(fileConn)
    }
    }
  return(stats)
  }
}

###
### Main body
###

round <- 0 # Counter of extension iterations
end <- FALSE # TRUE if no more sequences to extend.
dir.create(file.path(".", opt$output), showWarnings = FALSE)
dir.create(file.path(".", paste0(opt$output,"/final")), showWarnings = FALSE)

while (round < opt$top_rounds & !(end) ) {
round <- round + 1
dir.create(file.path(".", paste0(opt$output,"/mafs_",round)), showWarnings = FALSE)
dir.create(file.path(".", paste0(opt$output,"/to_extend_",round)), showWarnings = FALSE)

print(paste0("[+++] STARTING ROUND ",round))
print(paste0("[+++] Sequences folder: ",opt$input))
print(paste0("[+++] Mafs folder: ",opt$output,"/mafs_",round))

# Fasta table for the special first case (all should try to be extended)
if (round==1) {
  fasta.table <- fread(cmd = paste0("cat ",opt$input,"/*.fa |  awk '$0 ~ \">\" {if (NR > 1) {print c \"\tFALSE\tFALSE\";} c=0;printf substr($0,2,100) \"\t\"; } $0 !~ \">\" {c+=length($0);} END { printf c \"\tFALSE\tFALSE\";}'"))
}
# Read fa files and generate their maf. Uses a modified version of make_align_from_blast from TE_AID
filenames <-
  list.files(opt$input, pattern = "*.fa", full.names = TRUE)
setwd(paste0(opt$output,"/mafs_",round))
for (i in 1:length(filenames)) {
  system(command = paste0("../../",script.basename,"/make_align_from_blast_alt.sh ",
                          opt$genome," ", "../../",filenames[i]," ", fasta.table[i,2]*opt$identity, " ", opt$extend," ",
                          fasta.table[i,3]," ",fasta.table[i,4]))
}
# Delete potential empty files on the maf folder. 
system(command ="find . -size 0 -print -delete")

setwd("../..")

filenames <-
  list.files(paste0(opt$output,"/mafs_",round), pattern = "*.maf.*", full.names = TRUE)

# Divide in sets of alignments to avoid overflow processing them. Sould be a variable
lf <- split(filenames, ceiling(seq_along(filenames)/20))
stats_all <- list()

# Main loop for round: reads maf, generates a list of sub mafs, and finally 
# writes the consensus for each one and returns a list of stats about them.
for (i in lf) {
  stats <- lapply(i, function(x) {process_maf(seq_clus(read_maf(x)))})
  stats <- rbindlist(stats)
  stats_all <- rbindlist(list(stats_all,stats))
}

# Write the final stats
fwrite(stats_all, paste0(opt$output,"/stats_round_", round), sep = "\t", append = TRUE) 
print(paste("[+++] Round,", round, "completed.", nrow(stats_all),"consensus sequences generated." ))

# Dealing with the special case of overlapping opposite edges.
stats[,short:=gsub("#.*$","",seq)]
left <- stats[!(end_l) & end_r]
right <- stats[!(end_r) & end_l]
candidates <- merge(left,right, by ="short", allow.cartesian=TRUE)
merged_list <- c()
if (nrow(candidates)!=0) {
  for (name in seq(1,nrow(candidates),1)) {
    merged_name <- paste0(candidates[name]$short,"_merged")
    if (!(merged_name %in% merged_list)) {
    if (candidates[name]$end_r.x) {
      sleft <- candidates[name]$seq.y
      sright <- candidates[name]$seq.x
    } else {
      sleft <- candidates[name]$seq.x
      sright <- candidates[name]$seq.y
    }
    seq_l <- fread(cmd = paste0("awk '(NR)%2==0' ",opt$output,"/to_extend_",round,"/",sleft,".con.fa"),header = F, col.names=c("sequence"))
    seq_r <- fread(cmd = paste0("awk '(NR)%2==0' ",opt$output,"/to_extend_",round,"/",sright,".con.fa"),header = F, col.names=c("sequence"))
    merged <- FALSE
    for (i in seq(50,(min(nchar(seq_l),nchar(seq_r))-1),1)) {
      if (mapply(function(x,y) sum(x!=y),
                 strsplit(substr(seq_r,1,i),""),
                 strsplit(substr(seq_l,nchar(seq_l)-i+1,nchar(seq_l)),"")) < 5) { # Replace with KD between both
        print(paste("Merging edges for",candidates[name]$short,"- Overlap =", i))
        merged_seq <- paste0(seq_l, substr(seq_r,i+1,nchar(seq_r)))
        clean_fasta <- cbind(c(t(data.table(paste0(">",merged_name),merged_seq))))
        fwrite(clean_fasta,paste0(opt$output,"/final/", merged_name,"-",name, ".con.fa"), sep=" ", col.names= FALSE)
        unlink(paste0(opt$output,"/to_extend_",round,"/",sleft,".con.fa"))
        unlink(paste0(opt$output,"/to_extend_",round,"/",sright,".con.fa"))
        stats_all[seq==paste0(sleft,".con.fa")] <- NULL
        stats_all[seq==paste0(sright,".con.fa")] <- NULL
        merged <- TRUE
      }
      if (merged) { 
        merged_list <- c(merged_list,merged_name)
        break
      }
    }
  }
  }
}

fasta.table <- stats_all[,c(1,4,13,14)]
fasta.table <- fasta.table[!(end_l==TRUE & end_r == TRUE)] 
opt$input <- paste0(opt$output,"/to_extend_",round)
if (nrow(fasta.table)==0) { end <- TRUE}

}

