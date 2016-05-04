#!/usr/bin/env Rscript

##########################################################################################
# MSKCC DMP
# Run Facet
# https://github.com/mskcc/facets
# author: Raghu Chandra Mohan, Ronak H Shah
# Date: May 2 2016
# Input: Counts File
# Output: Copy Number Files
##########################################################################################

LocationOfThisScript = function()
  # Function LocationOfThisScript returns the location of this .R script (may be needed to source other files in same dir)
{
  this.file = NULL
  # This file may be 'sourced'
  for (i in-(1:sys.nframe())) {
    if (identical(sys.function(i), base::source))
      this.file = (normalizePath(sys.frame(i)$ofile))
  }
  
  if (!is.null(this.file))
    return(dirname(this.file))
  
  # But it may also be called from the command line
  cmd.args = commandArgs(trailingOnly = FALSE)
  cmd.args.trailing = commandArgs(trailingOnly = TRUE)
  cmd.args = cmd.args[seq.int(from = 1,
                              length.out = length(cmd.args) - length(cmd.args.trailing))]
  res = gsub("^(?:--file=(.*)|.*)$", "\\1", cmd.args)
  
  # If multiple --file arguments are given, R uses the last one
  res = tail(res[res != ""], 1)
  if (0 < length(res))
    return(dirname(res))
  
  # Both are not the case. Maybe we are in an R GUI?
  return(NULL)
}

####################################################
### Run analysis
####################################################
####################################################
if (!interactive()){
  pkgs = c('argparse', 'facets')
  junk <-
    lapply(pkgs, function(p) {
      suppressPackageStartupMessages(require(p, character.only = T))
    })
  rm(junk)
  currwd <- getwd()
  path = LocationOfThisScript()
  source(paste(path, "funcs.R", sep = "/"))
  cc <- function(...) {
    paste(..., sep = '_')
  }
  ### Parse input arguments
  parser = ArgumentParser()
  parser$add_argument("-c", "--counts_file", nargs = 1, help = "Paired Counts File")
  parser$add_argument("-p", "--project_name", nargs = 1, default="facetProject", help = "Project Name")
  parser$add_argument("-t", "--tumor_name", nargs = 1, default="Tumor", help = "Tumor Sample Name")
  parser$add_argument("-n", "--normal_name", nargs = 1, default = "Normal", help = "Normal Sample Name")
  args = parser$parse_args()
  ### Get File, Base and Pool Information
  FILE = args$counts_file
  tumorName = args$tumor_name
  normalName = args$normal_name
  projectName = args$project_name
  BASE = basename(FILE)
  POOL = basename(dirname(FILE))
  ### Print Facet Version Information
  buildData = installed.packages()["facets", ]
  cat("#Module Info\n")
  for (fi in c("Package", "LibPath", "Version", "Built")) {
    cat("#", paste(fi, ":", sep = ""), buildData[fi], "\n")
  }
  version = buildData["Version"]
  cat("\n")
  ### Get variables for facet input
  BASE = gsub(".dat.*", "", BASE)
  #sampleNames = gsub(".*recal_", "", strsplit(BASE, "____")[[1]])
  #tumorName = sampleNames[1]
  #normalName = sampleNames[2]
  #projectName = paste("MSK-IMPACT", POOL, sep = "-")
  TAG = paste("facets",
              projectName,
              tumorName,
              normalName,
              "cval",
              "150",
              sep = "__")
  TAG1 = cc(projectName, tumorName, normalName)
  chromLevels = c(1:22, "X")
  
  ### Run Analysis
  dat = preProcSample(
    FILE,
    snp.nbhd = 250,
    cval = 100,
    chromlevels = chromLevels,
    ndepth = 35,
    hetscale = TRUE,
    unmatched = FALSE,
    het.thresh = 0.25
  )
  out = procSample(dat,
                   cval = 100,
                   min.nhet = 15,
                   dipLogR = NULL)
  fit = emcncf(out, min.nhet = 15)
  
  
  pdf(
    file = cc(TAG, "BiSeg.pdf"),
    height = 10,
    width = 8
  )
  plotSample(
    out,
    chromlevels = chromLevels,
    emfit = fit,
    plot.type = "em",
    sname = paste(POOL, tumorName, normalName, sep = "_")
  )
  #text(-.08,-.08,paste(projectName,"[",tumorName,normalName,"]","cval =150"),xpd=T,pos=4)
  dev.off()
  out$IGV = formatSegmentOutput(out, TAG1)
  save(out,
       fit,
       file = cc(TAG, ".Rdata"),
       compress = T)
  
  ff = cc(TAG, ".out")
  cat("# TAG =", TAG1, "\n", file = ff)
  cat("# Version =",
      version,
      "\n",
      file = ff,
      append = T)
  cat("# Input =",
      basename(FILE),
      "\n",
      file = ff,
      append = T)
  cat("# snp.nbhd =",
      "250",
      "\n",
      file = ff,
      append = T)
  cat("# cval =", "100", "\n", file = ff, append = T)
  cat("# min.nhet =", "15", "\n", file = ff, append = T)
  cat("# Project =",
      projectName,
      "\n",
      file = ff,
      append = T)
  cat("# Tumor =",
      tumorName,
      "\n",
      file = ff,
      append = T)
  cat("# Normal =",
      normalName,
      "\n",
      file = ff,
      append = T)
  cat("# Purity =",
      fit$purity,
      "\n",
      file = ff,
      append = T)
  cat("# Ploidy =",
      fit$ploidy,
      "\n",
      file = ff,
      append = T)
  cat("# dipLogR =",
      fit$dipLogR,
      "\n",
      file = ff,
      append = T)
  cat("# dipt =", fit$dipt, "\n", file = ff, append = T)
  cat("# loglik =",
      fit$loglik,
      "\n",
      file = ff,
      append = T)
  cat("# emflag =",
      fit$emflags,
      "\n",
      file = ff,
      append = T)
  cat("# flags =",
      out$flags,
      "\n",
      file = ff,
      append = T)
  
  write.table(cbind(out$IGV[, 1:4], fit$cncf[, 2:ncol(fit$cncf)]),
              cc(TAG, "cncf.txt"), row.names = F)
  write.table(
    out$IGV,
    file = cc(TAG, '.seg'),
    row.names = F,
    quote = F,
    sep = "\t"
  )
  
}


#  pdf(file=cc(TAG,".pdf"),height=10,width=8)

#plotSampleCNCF.custom(out$jointseg,out$out,fit,
#                      main=paste(projectName,"[",tumorName,normalName,"]","cval =",CVAL))
# plotSampleCNCF(out,fit)

# dev.off()



##################################################################################
