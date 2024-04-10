#Dependencies
#core algorithms
#1.  data.table
#2.  foreach
#3.  Rfast
#4.  Rfast2
#5.  igraph
#6.  gstat
#7.  collapse
#8.  fastmatch

#parallel processing and reporting
#9.  future
#10. doFuture
#11. progressr

#Plotting
#12: ggplot2
#13. ggrepel
#14. GGally
#15. cowplot
#16. patchwork
#17. ggforce


#===========================================================================
# Internal functions
#===========================================================================

# suppress print, message, warnings etc
sup.all <- function(x) {
  capture.output(out <- suppressMessages(suppressWarnings(x)))
  return(out)
}

# verbose function, prints message if verbose option is set to TRUE
vrb <- function(m) {
  if (parent.frame(n = sys.nframe() - 1)$verbose) { # access verbose setting from parent function
    cat(m)
  }
}

# check for appropriate plR input type
class.check <- function(plR.input, req.class, env) {
  if (is.null(plR.input)) {
    stop("Required plR.input object is not present", call. = FALSE)
  } else {
    plR.class <- plR.input$plR.class
    mss <- paste0("plR.input object is not valid\nplR.class = ",
                  plR.class, " is not among the following permissible classes: ", 
                  paste(req.class, collapse = " or "))
    if (is.null(plR.class)) {
      stop(mss, call. = FALSE)
    } else {
      if (plR.class %in% req.class) { # return output to parent environment
        list2env(plR.input, envir = env)
      } else {
        stop(mss, call. = FALSE)
      }
    }
  }
}

# check for appropriate plR input type for plotting functions
class.check0 <- function(plR.input, env) {
  if (is.null(plR.input)) {
    stop("Required plR.input object is not present", call. = FALSE)
  } else {
    plR.class <- plR.input$plR.class
    input <- NULL
    if (is.null(plR.class)) {
      stop("Required plR.input object is not present", call. = FALSE)
    } else {
      plr.class <- strsplit(plR.class, "\\|")[[1]] # retrieve user parameters
      if (plr.class[3] == 2) {
        input <- c(input, "pruned")
      } else {
        if (!is.null(plR.input$vario.summary) | plr.class[2] == 2) {
          if (!is.null(plR.input$vario.summary)) {
            input <- c(input, "variogram")
          }
          if (plr.class[2] == 2) {
            input <- c(input, "rescaled")
          }
        } else {
          if (plr.class[1] %in% 2:3) {
            if (!is.null(plR.input$perm.mat)) {
              input <- c(input, "matched")
            }
          }
        }
      }
    }
    if (is.null(input)) {
      stop(paste0("plR.input object is not valid\n", 
                  "Must be matched, rescaled, pruned, or include ",
                  "fitted variogram summary"), call. = FALSE)
    } else {
      list2env(plR.input, envir = env) # unpack plR object
      list2env(list(input = input, plr.class = plr.class), envir = env)
    }
  }
}

# plR file reader
file.load <- function(input.path, group, bp.to.cM, env) {
  # check valid file paths
  if (!dir.exists(input.path)) {
    stop("Invalid path to parent directory\n", call. = FALSE)
  }
  
  ll <- list.files(input.path)
  if (!is.null(group)) {
    ll <- ll[grep(group, ll)]
  }
  
  if (length(ll) == 0) {
    stop("No input files detected, check path to parent directory",
         call. = FALSE)
  }
  
  req.files <- c("ObjInfo", "SetInfo", "SetObj")
  if (bp.to.cM) {
    opt.files <- c("VarInfo", "RecRate")
  } else {
    opt.files <- "VarInfo"
  }
  EMSS <- MSS <- NULL
  for (i in list(req.files, opt.files)) {
    check.f <- sapply(i, function(x) grepl(x, ll))
    nT <- Rfast::colsums(check.f)
    if (any(nT != 0)) {
      if (any(nT > 1)) {
        fg1 <- paste(i[which(nT > 1)], collapse = ", ")
        mss0 <- paste0("Multiple copies of optional file",
                       ifelse(length(fg1) > 1, "s ", " "),
                       fg1, " detected\n")
        EMSS <- c(EMSS, mss0)
      }
      if (any(nT == 0)) {
        f0 <- paste(i[which(nT == 0)], collapse = ", ")
        mss0 <- paste0("[Note: optional file", 
                       ifelse(length(f0) == 2, "s ", " "),
                       f0, " not detected]\n")
        if (length(i) == 3) {
          EMSS <- c(EMSS, mss0)
        } else {
          MSS <- c(MSS, mss0)
        }
      }
    }
  }
  
  if (!is.null(EMSS)) {
    stop(paste0(EMSS), 
         "Check path and ensure that file names are correct", call. = FALSE)
  }
  
  # load data
  vrb(paste0("Loading required files ", 
             paste(req.files, collapse = ", "), "\n"))
  sh.pr <- parent.frame()$verbose
  for (j in 1:3) {
    i <- req.files[j]
    I <- c("obj.info", "set.info", "set.obj")[j]
    file.in <- file.path(input.path, ll[grep(i, ll)])
    assign(x = I, envir = env, 
           value = data.table::fread(file.in, showProgress = sh.pr))
  }
  
  cA <- Rfast::colAny(check.f)
  if (sum(cA) > 0){
    vrb(paste0("Loading optional file",
               ifelse(sum(cA) == 1, " ", "s "),
               paste(opt.files[cA], collapse = ", "), "\n\n"))
    if (!is.null(MSS)) {
      vrb(paste0(MSS, "\n"))
    } 
    for (j in 1:length(opt.files)) {
      i <- opt.files[j]
      I <- c("var.info", "rec.rate")[j]
      if (cA[j]) {
        file.in <- file.path(input.path, ll[grep(i, ll)])
        assign(x = I, envir = env, 
               value = data.table::fread(file.in, showProgress = sh.pr))
      } else {
        assign(x = I, envir = env, value = NULL)
      }
    }
  } else {
    vrb(paste0(MSS, "\n"))
  }
}

# check file completeness and create contiguous integer values for core parameters
file.check <- function(obj.info, set.info = NULL, set.obj = NULL,
                       var.info = NULL, rec.rate = NULL, update.indices = TRUE, 
                       env) {
  
  ##--------------------------------------------------------------------------##
  # check that appropriate columns are present in each input
  ##--------------------------------------------------------------------------##
  
  if (!data.table::is.data.table(obj.info)) { # convert inputs to data.table objects
    data.table::setDT(obj.info)
  }
  oi.copy <- data.table::copy(obj.info)
  c.oi <- c("objID", "objStat", "startpos", "endpos")
  chk.oi <- all(c.oi %in% names(oi.copy))
  if (chk.oi) { # check numeric columns (chr class not constrained)
    chk.oi <- all(apply(oi.copy[, .SD, .SDcols = c.oi], 2, is.numeric))
  }
  
  if(!data.table::is.data.table(set.obj)) { # convert inputs to data.table objects
    data.table::setDT(set.obj) 
  }
  so.copy <- data.table::copy(set.obj)
  c.so <- c("setID", "objID")
  chk.so <- all(c.so %in% names(so.copy))
  if (chk.so) {
    chk.so <- all(apply(so.copy[, .SD, .SDcols = c.so], 2, is.numeric))
  }
  
  if(!data.table::is.data.table(set.info)) { # convert inputs to data.table objects
    data.table::setDT(set.info)
  }
  si.copy <- data.table::copy(set.info)
  chk.si <- "setID" %in% names(si.copy)
  if (chk.si) {
    chk.si <- is.numeric(si.copy$setID)
  }

  # var.info file
  if (!is.null(var.info)) {
    if(!data.table::is.data.table(var.info)) { # convert inputs to data.table objects
      data.table::setDT(var.info)
    }
    c.vi <- c("pos", "value")
    chk.vi <- all(c.vi %in% names(var.info))
    if (chk.vi) {
      chk.vi <- all(apply(var.info[, .SD, .SDcols = c.vi], 2, is.numeric))
    }
  } else {
    chk.vi <- TRUE
  }
  
  # rec.rate file
  if (!is.null(rec.rate)) {
    if(!data.table::is.data.table(rec.rate)) { # convert inputs to data.table objects
      data.table::setDT(rec.rate)
    }
    c.rr <- c("startpos", "endpos", "rec.rate")
    chk.rr <- all(c.rr %in% names(rec.rate))
    if (chk.rr) {
      chk.rr <- all(apply(rec.rate[, .SD, .SDcols = c.rr], 2, is.numeric))
      if (chk.rr) {
        chk.rr <- class(rec.rate$chr) %in% c("numeric", "integer", "character")
      }
    }
  } else {
    chk.rr <- TRUE
  }
  
  mW <- which(!c(chk.oi, chk.si, chk.so, chk.vi, chk.rr))
  if (length(mW) > 0) {
    miss.obj <- c("obj.info", "set.info", "set.obj", "var.info", "rec.rate")[mW]
    stop("At least one column label missing or incorrect class for ", 
         paste0(miss.obj, collapse = " and "), 
         "\nSee ?polylinkR::plR.read for required columns", call. = FALSE)
  }
  
  ##--------------------------------------------------------------------------##
  # check for shared gene, gene set, chromosome labels
  ##--------------------------------------------------------------------------##
  
  sxv <- si.copy$setID %in% so.copy$setID
  xsv <- so.copy$setID %in% si.copy$setID
  if(!all(sxv) | !all(xsv)) {
    stop("Gene set labels [setID] differ between set.info and set.obj",
         "\nPlease ensure common gene set labels across input files", 
         call. = FALSE)
  }
  
  # check for shared gene and gene set labels
  xov <- so.copy$objID %in% oi.copy$objID
  if(!all(xov)) {
    stop("Some set.obj gene labels [objID] are not found in obj.info",
         "\nPlease ensure set.obj gene labels are present in obj.info", 
         call. = FALSE)
  }

  # check for shared chromosome labels
  oc <- unique(oi.copy$chr)
  ocn <- length(oc)
  
  if (exists("var.info")) {
    if (!is.null(var.info)) { # variant file
      vT <- TRUE
      vc <- unique(var.info$chr)
      wovc <- which(oc %in% vc)
      if (length(wovc) < ocn) {
        e.mss1 <- c("obj.info chromosomes ", paste0(oc[wovc], collapse = ", "),  
                   " are absent from var.info\n")
      } else {
        e.mss1 <- ""
      }
    } else {
      vT <- FALSE
      wovc <- vc <- oc
      e.mss1 <- ""
    } 
  } else {
    vT <- FALSE
    wovc <- vc <- oc
    e.mss1 <- ""
  }
  
  if (exists("rec.rate")) {
    if (!is.null(rec.rate)) { # recombination rate file
      rc <- unique(rec.rate$chr)
      worc <- which(oc %in% rc)
      if (length(worc) < ocn) {
        e.mss2 <- c("obj.info chromosomes ", paste0(oc[worc], collapse = ", "),  
                   " are absent from rec.rate\n")
      } else {
        e.mss2 <- ""
      }
    } else {
      worc <- rc <- oc
      e.mss2 <- ""
    }
  } else {
    worc <- rc <- oc
    e.mss2 <- ""
  }
  
  if (length(wovc) < ocn | length(worc) < ocn) {
    stop(paste(e.mss1, e.mss2, 
               "Please ensure that all obj.info chromosomes are present in",
               "optional files"), call. = F)
  } else {
    if (length(vc) > ocn) { # removing var.info chromosomes not in obj.info
      var.info <- var.info[chr %in% oc]
    }
    if (length(rc) > ocn) { # removing rec.rate chromosomes not in obj.info
      rec.rate <- rec.rate[chr %in% oc]
    }
  }

  # check consistent distance units among obj.info and var.info files
  # for all chromosomes, ratio of lengths should be > 0.5 and < 2; i.e. abs(log2(ratio)) < 1
  if (vT) {
    v.chr.d <- var.info[, diff(range(pos)), by = chr]$V1
    o.chr.d <- oi.copy[, max(endpos) - min(startpos), by = chr]$V1
    if (all(abs(log2(v.chr.d / o.chr.d)) > 1)) {
      stop("Inconsistent distance units among obj.info and var.info\n",
           "Check distance units in input files", call. = FALSE)
    }
  }
  
  ##--------------------------------------------------------------------------##
  # create integer indices for genes, gene sets, and chromosomes
  ##--------------------------------------------------------------------------##
  
  if (update.indices) {
    # ensure ojbIDs are contiguous and obj.info is ordered by objIDs
    oi.copy[, objID.orig := objID]
    oi.copy <- oi.copy[order(objID.orig)]
    oi.copy[, objID := 1:.N]
    
    # create unique integer for each chromosome
    oi.copy[, chr.orig := chr]
    oi.copy[, chr := as.integer(factor(chr.orig), 
                                levels = sort(unique(chr.orig)))]
    
    if (!("midpos" %in% names(oi.copy))) { # add midpos for genes
      oi.copy[, midpos := (startpos + endpos) / 2]
    }
    
    n.genes <- oi.copy[, .N] # gene count
    n.chr <- oi.copy[, data.table::uniqueN(chr)] # chromosome count
    
    if (vT) {
      if (!is.integer(var.info$chr)) {
        var.info[, chr := as.integer(factor(chr), 
                                     levels = sort(unique(chr)))]
      }
    }
    
    # ensure setIDs are contiguous and set.info is ordered by setIDs
    si.copy[, setID.orig := setID]
    si.copy <- si.copy[order(setID.orig)]
    si.copy[, setID := 1:.N]
    
    # regenerate set x object table
    SO.tmp <- data.table::merge.data.table(so.copy[, .(setID.orig = setID, 
                                                       objID.orig = objID)], 
                                           oi.copy[, .(objID, objID.orig)], 
                                           by="objID.orig")
    
    so.copy <- data.table::merge.data.table(SO.tmp, 
                                            si.copy[, .(setID, setID.orig)], 
                                            by="setID.orig")
    
    so.copy <- so.copy[order(setID, objID)][, .(setID, objID, 
                                                setID.orig, objID.orig)]
    
    set.genes <- so.copy[, sort(unique(objID))] # count of genes in gene sets
    n.set.genes <-  length(set.genes) # total genes occurring in at least 1 gene set
    n.sets <- si.copy[, .N] # gene set count
    
    list2env(list(n.genes = n.genes, set.genes = set.genes, 
                  n.set.genes = n.set.genes, n.sets = n.sets, n.chr = n.chr),
             envir = env)
  }
  
  # return output to parent environment
  list2env(list(oi.copy = oi.copy, si.copy = si.copy, so.copy = so.copy),
           envir = env)
}

# reconstitute original parameter labels
file.reset <- function(oi.copy = NULL, si.copy = NULL, so.copy = NULL) {
  oi.copy[, objID := objID.orig]
  oi.copy[, objID.orig := NULL]
  oi.copy[, chr := chr.orig]
  oi.copy[, chr.orig := NULL]
  
  si.copy[, setID := setID.orig]
  si.copy[, setID.orig := NULL]
  
  so.copy[, setID := setID.orig]
  so.copy[, setID.orig := NULL]
  so.copy[, objID := objID.orig]
  so.copy[, objID.orig := NULL]
}

# Set up parallel processing
pp.params <- function(n.cores, fut.plan) {
  pMSS <- NULL # catch messages
  pWMSS <- NULL # catch warnings
  
  # detect available cores
  N.cores <- future::availableCores()
  
  # check n.cores option
  nc.mss <- paste0("n.cores incorrectly specified\n",
                   "must be either positive integer, 'default' or 'max'")
  if (n.cores %in% c("default", "Default", "DEFAULT")) {
    nc <- N.cores - 1
  } else {
    if (n.cores %in% c("max", "Max", "MAX")) {
      nc <- N.cores
    } else {
      if (is.numeric(n.cores)) {
        nc <- as.integer(n.cores) # force integer
        if (n.cores <= 0) {
          stop(nc.mss, call. = FALSE)
        } else {
          if (n.cores > N.cores) {
            nc <- N.cores
            ncNc.mss <- paste0("n.cores exceeds available no. of cores [", 
                               N.cores, "], setting n.cores to max")
          }
        }
      } else {
        stop(nc.mss, call. = FALSE)
      }
    }
  }
  
  # check fut.plan option
  valid.fp <- c("default", "Default", "DEFAULT",
                "sequential", "multisession", "multicore", "cluster")
  if (!fut.plan %in% valid.fp) {
    stop("fut.plan incorrectly specified\nMust be one of: ",
         paste(valid.fp[c(1, 4:7)], collapse = ", "), 
         call. = FALSE)
  } else {
    if (fut.plan == "sequential") { # user has chosen sequential processing
      if (n.cores %in% c("default", "Default", "DEFAULT")) {
        nc <- 1
      }
      if (nc > 1) {
        nc <- 1
        pWMSS <- c(pWMSS, 
                   paste0("Mismatched parallel processing parameters: ",
                          "fut.plan = sequential but n.cores > 1, ",
                          "setting n.cores = 1\n"))
      } else {
        pMSS <- c(pMSS, "Sequential processing chosen [fut.plan = sequential]\n")
      }
      pMSS <- c(pMSS, "All processes will utilise a single core\n") 
    } else {
      if (fut.plan %in% c("multisession", "multicore", "cluster")) {
        if (n.cores == 1) {
          pWMSS <- c(pWMSS, 
                     paste0("Mismatched parameters: fut.plan = ", fut.plan,
                            " but n.cores = 1, setting fut.plan = sequential\n"))
          fut.plan <- "sequential"
          pMSS <- c(pMSS, "All processes will utilise a single core\n") 
        } else {
          mc.mss <- paste0("Parallel processing chosen [fut.plan = ", 
                           fut.plan, "]")
          if (N.cores == 1) {
            pWMSS <- c(pWMSS, 
                       paste0(mc.mss, " but only single core available\n",
                             "All processes will utilise a single core\n"))
            nc <- 1
          } else {
            if (fut.plan == "multicore") { # check for multicore support
              mcT <- future::supportsMulticore() 
              if(!mcT){
                fut.plan <- "multisession"
                pWMSS <- c(pWMSS, 
                           paste0(mc.mss, 
                                  "but multicore processing not supported\n",
                                  "Using multisession processing instead\n"))
              }
            }
            if (exists("ncNc.mss")) {
              pWMSS <- c(pWMSS, ncNc.mss)
            }
            pMSS <- c(pMSS, 
                      paste0(mc.mss, "\nSpecific processes will utilise ", 
                            nc, " cores\n"))
          }
        }
      }
    }
  }

  return(list(n.cores = nc, fut.plan = fut.plan, pWMSS = pWMSS, pMSS = pMSS))
}

# screen argument settings of core functions
arg.check <- function(plr.func, env) {
  # retrieve function arguments
  list2env(mget(ls(envir = env), envir = env), envir = environment()) 
  
  # check for violations in verbose setting
  if (!is.logical(verbose) | length(verbose) > 1) {
    stop("verbose is not logical\nRespecify or use default", call. = FALSE)
  }
  
  MSS <- NULL # collect messages
  WMSS <- NULL # collect warnings
  
  # perform function specific tests
  if (plr.func == "read") { # check polylinkR::plR.read arguments
    if (!is.character(input.path)) {
      stop("Invalid input.path [non-character string]", call. = FALSE)
    }
    
    e.mss <- paste0("merge.set.prop incorrectly specified [",
                    "must be numeric value between 0 and 1]\n",
                    "Respecify or use default")
    if(!is.numeric(merge.set.prop)) {
      stop(e.mss, call. = FALSE)
    } else {
      if (merge.set.prop > 1 | merge.set.prop <= 0) {
        stop(e.mss, call. = FALSE)
      }
    }
    
    e.mss <- paste0("bp.to.cM is not a logical value\n",
                    "Respecify or use default")
    if (!is.logical(bp.to.cM)) {
      stop("bp.to.cM is not a logical value\nRespecify or use default")
    } else {
      if (bp.to.cM) {
        mf.names <- c("Haldane", "Kosambi", "Carter-Falconer", "Morgan")
        fmss <- c("map.fun is not valid\nOptions: ",
                  paste(paste(1:4, "="), mf.names, collapse = ", "))
        if (!is.numeric(map.fun) | length(map.fun) > 1) {
          stop(fmss, call. = FALSE)
          if (!(map.fun %in% 1:4)) {
            stop(fmss, call. = FALSE)
          }
        }
      }
    }
    
    if (!is.numeric(min.set.n) | length(min.set.n) > 1) {
      stop("min.set.n incorrectly specified ",
           "[must be numeric value between 2 and Inf]\n",
           "Respecify or use default", call. = FALSE)
    }
    
    if (!is.numeric(max.set.n) | length(max.set.n) > 1) {
      stop("max.set.n incorrectly specified [",
           "must be numeric value between 10 and Inf]\n",
           "Respecify or use default", call. = FALSE)
    }
    
    # check set size filtering criteria
    ss.ch <- c(min.set.n < 2, max.set.n < 10, min.set.n > max.set.n)
    mss.ch <- c("min.set.n < 2", "max.set.n < 10", "min.set.n > max.set.n")
    if (any(ss.ch)) {
      W <- mss.ch[which(ss.ch)]
      stop("Set size filters incorrectly specified: ",
           paste(W, collapse = " and "),
           "\nRespecify or use default", call. = FALSE)
    }
    
    if (!is.logical(rem.ident.genes) | length(rem.ident.genes) > 1){
      stop("rem.ident.genes must be logical\nRespecify or use default", 
           call. = FALSE)
    }
    
    # check gene set and gene filters
    if (!is.null(set.in) & !is.null(set.out)) {
      stop("Use either set.in or set.out option, not both", 
           call. = FALSE)
    }
    
    if (!is.null(obj.in) & !is.null(obj.out)) {
      stop("Use either obj.in or obj.out option, not both", 
           call. = FALSE)
    }
    
    ss.i <- foreach::foreach(I = list(set.in, set.out, obj.in, obj.out),
                             i = iterators::icount(), .combine = c) %do% {
      if (is.null(I)) {
        si <- FALSE
      } else {
        si <- !is.numeric(I) | !is.vector(I)
      }
      si
    }
    
    mss.ch <- c("set.in", "set.out", "obj.in", "obj.out")
    if (any(ss.i)) {
      W <- which(ss.i)
      stop(ifelse(W %in% 1:2, "Set ", ""), 
           ifelse(W %in% 1:2 & W %in% 3:4, "and ", ""),
           ifelse(W %in% 3:4, ifelse(W %in% 1:2, "Gene ", "gene "), ""),
           "filters incorrectly specified\n",
           "Must be integer vectors with specifying genes or gene sets that",
           "should be kept [set.in or obj.in] or removed [set.out or obj.out]",
           "\nRespecify or use default", call. = FALSE)
    }
  } else { # other plR functions
    # check seed
    if (!is.null(user.seed)) {
      if (!is.numeric(user.seed)) {
        stop("Incorrect random seed specification:\n",
             "Must be integer value x such that abs(x) <= ", 
             .Machine$integer.max, "\nRespecify or use default", call. = FALSE)
      } else {
        if (!is.integer(user.seed)) {
          if (user.seed %% 1 > 0) {
            seed <- as.integer(round(user.seed, 0))
            WMSS <- c(WMSS, 
                      paste0("user.seed is non-integer\n", 
                             "Converted user.seed to nearest integer value [", 
                             seed, "]\n"))
          } else {
            seed <- as.integer(user.seed)
          }
        } else {
          seed <- user.seed
        }
      }
    } else {
      # generate random seed (useful to reproduce results if seed not set)
      seed <- sample.int(.Machine$integer.max, 1) * sample(c(-1, 1), 1)
    }
    
    # check tail fit options for null distributions
    e.mss <- paste0("tail.fit incorrectly specified\n",
                    "Must be integer value either 0, 1, or 2\n",
                    "[0 = no tail fitting, 1 = max likelihood, 2 = method of moments]")
    n.mss <- paste0("n.tail incorrectly specified\n",
                    "Must be integer value between 10 and 100")
    if (!is.numeric(tail.fit)) {
      stop(e.mss, call. = FALSE)
    } else {
      if (!tail.fit %in% 0:2) {
        stop(e.mss, call. = FALSE)
      } else {
        if (tail.fit == 0) {
          fit.meth <- NULL
        } else {
          fit.meth <- ifelse(tail.fit == 1, "mle", "pwm")
          if (!is.numeric(tail.fit)) {
            stop(n.mss, call. = FALSE)
          } else {
            if (n.tail < 10 | n.tail > 100) {
              stop(n.mss, call. = FALSE)
            } else {
              n.tail <- as.integer(n.tail)
              if (n.tail %% 10 != 0) {
                n.tail <- n.tail %/% 10 * 10 # make value divisible by 10
                WMSS <- c(WMSS, 
                          paste0("n.tail must be divisible by 10\n",
                                 "n.tail reset to ", n.tail, "\n"))
              }
            }
          }
        }
      }
    }
    
    if (plr.func == "permute") { # check polylinkR::plR.permute arguments
      # check for appropriate matching choices
      if (all(perm.type %in% 1:3) & length(perm.type) == 1) {
        if (perm.type == 1) {
          # matching not chosen
          if (!is.null(match.wt)) {
            stop("perm.type = 1 but user has provided matching weights\n",
                 "[perm.type = 1 and match.wt = NULL will result in random sampling]\n",
                 "[perm.type = 2 or 3 and match.wt = NULL will result in",
                 " matched sampling using inbuilt algorithm]\n",
                 "[perm.type = 2 or 3 and match.wt != NULL will result in",
                 " user-defined matched sampling]", call. = F)
          }
        } else { # matching chosen
          if (!is.null(match.fun)) {
            if (!is.function(match.fun)) {
              stop("match.fun is not an R function", call. = FALSE)
            }
          }
          match.z <- ifelse(perm.type == 2, FALSE, TRUE)
          assign(x = "match.z", envir = env, value = match.z)
          # check for appropriate covariate columns
          cov.cols <- grep("Cov|COV", names(oi.copy))
          assign(x = "cov.cols", envir = env, value = cov.cols)
          e.mss <- c("Appropriately labelled covariate columns were not detected ", 
                     "in oi.copy\nLabels must contain 'Cov' or 'COV' followed by ",
                     "integer from 1 to no. of covariates\n")
          if (length(cov.cols) == 0) { 
            stop(e.mss, call. = FALSE) # reset flag to perform standard permutations
          } else {
            cov.names <- names(oi.copy[, .SD, .SDcols = cov.cols])
            cov.int <- gsub("Cov|COV", "", cov.names)
            if (!all(cov.int == 1:length(cov.cols))) {
              stop(e.mss, call. = FALSE)
            }
            cov.numeric <- oi.copy[, sapply(.SD, is.numeric), .SDcols = cov.cols]
            wC <- which(!cov.numeric)
            if (length(wC) > 0) {
              stop("Coverage column", ifelse(length(wC) == 1, " ", "s "),
                   paste0(cov.names[wC], collapse = ", "), " are non-numeric: ",
                   call. = FALSE)
            }
          }
          if (!is.null(match.wt)) { # check suitability of user supplied matching weights
            if (!is.matrix(match.wt)) {
              stop("User supplied match.wt object is not a matrix\n",
                   "Respecify or use default", call. = FALSE)
            } else {
              if (!is.null(obj.info.merged)) {
                req.n.obj <- max(max(obj.info.merged$objID), 
                                 max(oi.copy$objID.orig))
              } else {
                req.n.obj <- max(oi.copy$objID.orig)
              }
              if (nrow(match.wt) != req.n.obj | ncol(match.wt) != req.n.obj) {
                stop("Dimensions of match.wt matrix not correct\n",
                     "# of rows and columns of match.wt must equal #", 
                     " of genes in original input file [", req.n.obj, 
                     "]\nRespecify or use default", call. = FALSE)
              }
            }
            if (n.set.genes < n.genes) { # resize matrix
              assign(x = "match.wt", envir = env, value = match.wt[set.genes, ])
            }
          } else { # check MD scaling options
            if (is.null(match.fun)) {
              e.mss <- c("match.scale must be non-negative real number\n", 
                         "Respecify or use default")
              if (!is.numeric(match.scale) | length(match.scale) != 1) {
                stop(e.mss, call. = FALSE)
              } else {
                if (match.scale < 0) {
                  stop(e.mss, call. = FALSE)
                }
                if (match.scale == 0) {
                  stop("match.scale = 0; all sampling probabilities will be equal\n",
                       "Set match = F if equal sampling probabilities are required", 
                       call. = FALSE)
                }
              }
            } else {
              if (!is.function(match.fun)) {
                stop("match.fun must be an R function\nRespecify or use default", 
                     call. = FALSE)
              }
            }
            if (!is.numeric(match.cap) | length(match.cap) != 1) {
              stop("match.cap is not numeric\nRespecify or use default", 
                   call. = FALSE)
            } else {
              if (match.cap == 1 / (n.genes - 1)) {
                stop("match.cap = 1 / (n.genes - 1); ",
                     "all sampling probabilities are equal\n",
                     "Set match = F if equal sampling probabilities are required", 
                     call. = FALSE)
              } else {
                if (match.cap < 1 / (n.genes - 1)) {
                  stop("match.cap must be >= 1 / (genes - 1); ", 
                       "probability cap is lower than equal sampling weights\n",
                       "Respecify or use default\n",
                       call. = FALSE)
                }
                if (match.cap > 1) {
                  stop("match.cap must be <= 1\nRespecify or use default ",
                       "[Set match.cap = 1 if no upper bound required]", 
                       call. = FALSE)
                }
              }
            }
          }
        }
      } else { 
        stop("perm.type must be either 1, 2 or 3\n",
             "[1 = random, 2 = unstandardised matched, 3 = standardised matched]",
             call. = FALSE)
      }
      
      # permuation parameters
      if (!is.numeric(n.perm) | n.perm < 0 | length(n.perm) != 1) {
        stop("n.perm is not numeric or < 0\nRespecify or use default", 
             call. = FALSE)
      } else {
        assign(x = "n.perm", envir = env, value = as.integer(n.perm))
        if (n.perm < 1e4){
          WMSS <- c(WMSS, 
                    paste0("Permutation size (", n.perm, ") is less than 10k; ", 
                           "resetting to 10k\n"))
          n.perm <- 10000L
          assign(x = "n.perm", envir = env, value = n.perm)
        }
        if (n.perm %% 100 != 0){
          n.perm <- round(n.perm / 100, 0) * 100
          assign(x = "n.perm", envir = env, value = n.perm)
          WMSS <- c(WMSS,
                    paste0("Number of permuations not divisible by 100\nUsing ", 
                           n.perm, " iterations\n"))
        }
      }
    } else {
      if (plr.func == "rescale") { # check polylinkR::plR.rescale arguments
        use.var <- !is.null(var.info)
        assign(x = "use.var", envir = env, value = use.var)
        # check for violations in logical settings
        lp <- c("rescale", "inc.nugget")
        clp <- sapply(lp, function(x) is.logical(get(x)))
        if (!all(clp)) {
          stop(paste(lp[!clp], collapse = " & "),
               ifelse(sum(!clp) == 1, " is", " are"), " incorrectly specified\n",
               "Use logical value or default", call. = FALSE)
        }
        
        if (!is.null(ac)) { # check user provided autocovariance object
          exp.col <- c("objID.A", "objID.B", "midpos.dist", "covar")
          emm <- paste0("ac must be a dataframe with the following columns: '",
                        paste(exp.col, collapse = "', '"), "'",
                        "\nAND include a separate row for all gene pairs with ",
                        "non-zero autocovariance\n",
                        "[NB: appropriate ac object can be generated using ",
                        "polylinkR::plr.rescale; see manual]")
          if (!is.data.frame(ac)) {
            stop(emm, call. = FALSE)
          } else {
            data.table::setDT(ac) # ensure data.table format
            if (!all(exp.col %in% colnames(ac))) {
              stop(emm, call. = FALSE)
            } else {
              # check that gene variance is included
              ac.diag <- ac[objID.A == objID.B]$covar
              if (length(ac.diag) != n.genes | any(ac.diag <= 0)) {
                stop("Variance for all genes must be included and be postive value",
                     "[NB: appropriate ac object can be generated using ",
                     "polylinkR::plr.rescale; see manual]", 
                     call. = FALSE)
              }
            }
          }
        } else { # check for violations in parameter settings for variogram estimation
          if (use.var) {
            if (var.info[, .N] < n.genes) {
              WMSS <- c(WMSS, 
                        "Fewer variants in var.info than genes in obj.info\n")
            }
          }
          
          # check distances 
          # (obj.info and var.info already checked for similarity, so only need to check obj.info)
          if (!any(dist %in% c("check", "bp", "cM")) | length(dist) > 1) {
            stop("dist must be one of the three following options:\n",
                 "1) 'bp' (base pairs; physical distance)\n",
                 "2) 'cM' (centiMorgans; recombination-based distance\n",
                 "3) 'check' (default: distance information is used to guess",
                 " genomic distance type)", 
                 call. = FALSE)
          } else {
            oi.l <- oi.copy[, .(DIFF = max(endpos) - min(startpos)), by = chr]
            t1 <- oi.l[, median(DIFF)] # minimum chromosome length
            oi.d <- oi.copy[, .(DIFF = diff(sort(startpos))), by = chr]
            t2 <- oi.d[DIFF > 0, median(DIFF)] # minimum inter-gene distance
            if (dist == "check") {
              dist <- ifelse(t1 < 1e3 | t2 < 1, "cM", "bp") # guess units
              assign(x = "dist", envir = env, value = dist)
              MSS <- c(MSS, paste0("Assuming distance units are ", dist, "\n"))
            } else {
              if (dist == "cM") {
                if (t1 >= 1e3 | t2 >= 1) {
                  mss <- paste0("dist = cM but",
                                ifelse(t1 >= 1e3, " median chromosome length >= 1e3", ""),
                                ifelse(t1 >= 1e3 & t2 >= 1, " and ", " "),
                                ifelse(t2 >= 1, "median intergene distance >= 1", ""),
                                       "\nMay be physical distances [i.e. base pairs] rather than ",
                                       "recombination-based distances [i.e. centiMorgans]\n")
                  if (t1 >= 1e3 & t2 >= 1) {
                    stop(mss, "Check distance units in input files", 
                         call. = FALSE) 
                  } else {
                    WMSS <- c(WMSS, mss)
                  }
                }
              } else {
                if (t1 < 1e3 | t2 < 1) {
                  mss <- paste0("dist = bp but",
                                ifelse(t1 < 1e3, " median chromosome length < 1e3", ""),
                                ifelse(t1 < 1e3 & t2 < 1, " and ", " "),
                                ifelse(t2 < 1, "median intergene distance < 1", ""),
                                "\nMay be recombination-based distances [i.e. centiMorgans] ",
                                "rather than physical distances [i.e. base pairs]\n")
                  if (t1 < 1e3 & t2 < 1) {
                    stop(mss, "Check distance units in input files", 
                         call. = FALSE) 
                  } else {
                    WMSS <- c(WMSS, mss)
                  }
                }
              }
            }
          }
          
          # check kappa
          if (!is.numeric(kappa) | length(kappa) != 2) {
            stop("kappa must be a numeric vector with two entries specifying",
                 " upper and lower boundaries\nReset or use default\n",
                 call. = FALSE)
          } else {
            kappa <- sort(kappa) # make sure kappa is in ascending order
            assign(x = "kappa", envir = env, value = kappa)
            if (kappa[1] < 0) {
              stop("kappa values must be postive\nReset or use default\n",
                   call. = FALSE)
            } else {
              if (any(kappa > 1e5)) {
                WMSS <- c(WMSS,
                          paste0("kappa values > 1e5:",
                                 "Large kappa values can prolong optimisation\n"))
              }
            }
          }
          
          # check minimum correlation coefficient (min.rho)
          mr.mss <- c("min.rho incorrectly entered\n",
                      "Must be numerical value between 0 and 0.1")
          if (!is.numeric(min.rho) | length(min.rho) > 1) {
            stop(mr.mss, call. = FALSE)
          } else {
            if (min.rho < 0 | min.rho > 0.1) {
              stop(mr.mss, call. = FALSE)
            }
          }
          
          # check that variogram fitting method is supported
          fit.range <- c(1, 2, 6, 7)
          if (!is.numeric(vg.fit) | length(vg.fit) > 1) {
            stop("vg.fit incorrectly specified\n",
                 "Use default value or number in range: ",
                 paste(fit.range, collapse = ", "),
                 "\n[See ?gstat::fit.variogram for more information]", 
                 call. = FALSE)
          } else {
            if (!(vg.fit %in% fit.range)) {
              stop("User provided vg.fit value (", vg.fit, 
                   ") is not an option for gstat::fit.variogram\n",
                   "Valid options: ", paste(fit.range, collapse = ", "), 
                   "; default = 7\n[See ?gstat::fit.variogram for more information]",
                   call. = FALSE)
            }
          }
          
          #bin size, block size, maximum variogram distance
          EMSS <- function(x, mss) {
            paste0("User specified value for ", x, " is incorrect\n",
                   mss, "Either re-enter value or use default")
          }
          
          # check maximum variogram distance
          if (max.link != "default") {
            min.maxl <- ifelse(dist == "bp", 1e4, 0.01)
            if (!is.numeric(max.link) | block.size < min.maxl) {
              stop(EMSS(x = "max.link", 
                        mss = paste0("max.link must be numeric value >= ",
                                     min.maxl, dist, "\n")), call. = FALSE)
            }
          }
          
          #check block sizes (no. genes or variants per block)
          if (block.size != "default") {
            min.bs <- ifelse(use.var, 1e3, 100)
            ii <- ifelse(use.var, " variants\n", " genes\n")
            emss <- EMSS(x = "block.size", 
                         mss = paste0("block.size must be numeric value >= ",
                                      min.bs, ii))
            if (!is.numeric(block.size) | length(block.size) > 1){
              stop(emss, call. = FALSE)
            } else {
              if (block.size < min.bs) {
                stop(emss, call. = FALSE)
              }
            }
            # check that block size is smaller than no. variants / genes on least populated genome
            max.n.chr <- ifelse(use.var, 
                                max(var.info[, .N, by = chr]$N),
                                max(oi.copy[, .N, by = chr]$N))
            if (block.size > max.n.chr) {
              WMSS <- c(WMSS, 
                        paste0("block.size exceeds no. of genes on chromosome with",
                               " least no. of genes\nWhole chromosomes will be used\n"))
              assign("block.size", envir = env, value = Inf)
            }
          }
          
          # check bin sizes
          if (bin.size != "default") {
            min.bin <- ifelse(dist == "cM", 1e-4, 1e2)
            emss <- EMSS(x = "bin.size", 
                         mss = paste0("bin.size must be numeric value >= ",
                                     min.bin, dist, "\n"))
            if (!is.numeric(bin.size) | length(bin.size) > 1) {
              stop(emss, call. = FALSE)
            } else { # sufficient no. of distinct gene bins
              big.bin <- ifelse(dist == "cM", 0.1, 1e5)
              if (bin.size > big.bin) {
                WMSS <- c(WMSS, 
                          paste0("bin.size is > ", paste(big.bin, dist),
                                 "\nLarge bins coarsen inference ",
                                 "and lead to poor variogram estimation\n"))
              } else {
                if (bin.size < min.bin) {
                  stop(emss, call. = FALSE)
                } 
              }
            }
          }
        }
        
        # create logical capturing user rescaling choice
        assign(x = "rescaled", envir = env, 
               value = strsplit(plR.class, "\\|")[[1]][2] == 2)
      } else { # check polylinkR::plR.prune arguments
        hlt <- "Get min.set.n value used in polylinkR::plr.read"
        if (is.null(min.set.n)) {
          stop("min.set.n needs to be specified\n", hlt, call. = FALSE)
        } else {
          (!is.numeric(min.set.n) | length(min.set.n) > 1)
          if (min.set.n > min(si.copy$N) | min.set.n < 2) {
            stop("min.set.n incorrectly specified\n", hlt, call. = FALSE)
          }
        }
        
        n.perm <- ncol(setScore.exp)
        assign(x = "n.perm", envir = env, value = n.perm)
        e.mss <- paste0("n.fdr incorrectly specified\n",
                        "n.fdr must be integer value >= 100 and <=", 
                        n.perm / 20, "[i.e. n.perm / 20]")
        if (!is.numeric(n.fdr) | length(n.fdr) > 1) {
          stop(e.mss, call. = FALSE)
        } else {
          if (n.fdr < 100 | n.perm / n.fdr < 20) {
            stop(e.mss, call. = FALSE)
          }
        }
        
        if (!pr.method %in% 1:2 | length(pr.method) > 1) {
          stop("pr.method incorrectly specified\n",
               "Use 1 for regression and 2 for standard\n",
               "[Standard pruning not possible for rescaled gene sets]", 
               call. = FALSE)
        } else {
          if (pr.method == 1) { # use regression pruning
            pr.type <- "reg"
          } else {
            if (rescaled) {
              stop("pr.method incorrectly specified\n",
                   "Standard pruning not possible for rescaled gene sets", 
                   call. = FALSE)
            } else {  # use gene removal pruning
              pr.type <- "rem"
            }
          }
        }
        
        lp <- c("est.pi0", "fast.fdr")
        clp <- sapply(lp, function(x) is.logical(get(x)))
        if (!all(clp)) {
          stop(paste(lp[!clp], collapse = " & "), 
               ifelse(sum(!clp) == 1, " is", " are"), " incorrectly specified\n",
               "Use logical value or default", call. = FALSE)
        }
      }
    }
    #set up parallel back end
    pp.check <- pp.params(n.cores = n.cores, fut.plan = fut.plan)
    list2env(pp.check, envir = environment()) # return objects to function environment
    list2env(list(seed = seed, n.tail = n.tail, fit.meth = fit.meth, 
                  n.cores = n.cores, fut.plan = fut.plan,
                  pMSS = pMSS, pWMSS = pWMSS), envir = env)
  }
  # return objects to parent environment
  list2env(list(MSS = MSS, WMSS = WMSS), envir = env)
}

# screen argument settings of plotting functions
arg.check0 <- function(plr.func, env) {
  list2env(mget(ls(envir = env), envir = env), envir = environment()) 
  MSS <- NULL # collect messages
  WMSS <- NULL # collect warnings
  
  # check appropriate variable settings
  if (!dir.exists(plot.path)) {
    stop("plot.path does not exist\n", call. = FALSE)
  }
  
  if (!is.null(plot.name)) {
    if (!is.character(plot.name)) {
      stop("plot.name must be character string\n", call. = FALSE)
    }
  }
  
  # check common logical parameters 
  lp <- c("keep.summary", "verbose", "cov.mat")
  clp <- sapply(lp, function(x) is.logical(get(x)))
  if (!all(clp)) {
    stop(paste(lp[!clp], collapse = " & "), 
         ifelse(sum(!clp) == 1, " is", " are"), " incorrectly specified\n", 
         "Logical input required [i.e. TRUE or FALSE]\n", call. = FALSE)
  }
  
  # check facet specification
  if (is.numeric(max.facets)) {
    if (max.facets < 1) {
      stop("max.facets must be >= 1", call. = FALSE)
    }
    max.facets <- round(max.facets, 0)
  } else {
    stop("max.facets must be genetic", call. = FALSE)
  }
  
  #check that highlighted sets are present
  if (!is.null(highlight.sets)) {
    emss <- paste0("highlight.sets incorrectly specified\n",
                   "Must be numeric vector where each entry is a set ID\n")
    if (!is.numeric(highlight.sets)) {
      stop(emss, call. = FALSE)
    } else {
      if (length(highlight.sets) == 0) {
        stop(emss, call. = FALSE)
      } else {
        hs.check <- highlight.sets %in% si.copy$setID.orig
        if (!all(hs.check)) {
          stop("Some gene sets [", 
               paste(highlight.sets[!hs.check], collapse = ", "),
               "] in highlight.sets are not present in set.info\n", 
               call. = FALSE)
        }
      }
    }
  }
  
  if ("matched" %in% input) { # check polylinkR::check.match arguments
    # check that covariate columns are present
    w.cov <- names(oi.copy)[grep("Cov|COV", names(oi.copy))]
    if (length(w.cov) > 0){
      assign(x = "w.cov", envir = env, value = w.cov)
      if (!is.null(cov.labels)) {
        if (length(cov.labels) != length(w.cov)) {
          stop(paste0("length of cov.labels [", length(cov.labels),
                      "] != number of covariate columns in obj.info [" , 
                      length(w.cov), "]\n"), call. = F)
        }
      }
    } else {
      stop("No covariate columns detected in obj.info\n", call. = F)
    }

    n.perm <- perm.param$n.perm
    if (!is.logical(fast.plot)) {
      stop("fast.plot incorrectly specified\n", 
           "Logical input required\n", call. = FALSE)
    } else { 
      if (fast.plot & n.perm > 1e4) { # implement fast plotting
        n.perm <- 1e4
        perm.mat <- perm.mat[1:n.perm, ]
        assign(x = "perm.mat", envir = env, value = perm.mat)
      }
    }
    assign(x = "n.perm", envir = env, value = n.perm)
  } else {
    if ("variogram" %in% input) { # check polylinkR::check.ac.rs arguments
      lp <- c("log.gamma", "log.dist")
      clp <- sapply(lp, function(x) is.logical(get(x)))
      if (!all(clp)) {
        stop(paste(lp[!clp], collapse = " & "),
             ifelse(sum(!clp) == 1, " is", " are"), " incorrectly specified\n", 
             "Logical input required [i.e. TRUE or FALSE]\n", call. = FALSE)
      }
    }
  }
  
  # optional covariate biplots
  if (cov.mat) {
    if (!("matched" %in% input)) {
      w.cov <- names(oi.copy)[grep("Cov|COV", names(oi.copy))]
    }
    # check covariate labels
    if (length(w.cov) > 0){
      assign(x = "w.cov", envir = env, value = w.cov)
      if (!is.null(cov.labels)) {
        if (length(cov.labels) != length(w.cov)) {
          stop(paste0("length of cov.labels [", length(cov.labels),
                      "] != number of covariate columns in obj.info [" , 
                      length(w.cov), "]\n"), call. = F)
        }
      }
    } else {
      WMSS <- c(WMSS, 
                paste0("cov.mat = T but no covariate columns detected in ",
                       "obj.info\nCovariance biplots will not be generated\n"))
      assign(x = "cov.mat", envir = env, value = FALSE)
    }
  }
  
  # return objects to parent environment
  list2env(list(MSS = MSS, WMSS = WMSS), envir = env)
}

# convert physical distances to genetic distances
bp2cM <- function(rec.rate, obj.info, var.info, map.fun, env) {
  # check that rec.rate and obj.info file are in physical units
  if (any(obj.info$startpos %% 1 != 0) | any(obj.info$endpos %% 1 != 0)) {
    t1 <- TRUE
  } else {
    t1 <- FALSE
  }
  
  if (any(rec.rate$startpos %% 1 != 0) | any(rec.rate$endpos %% 1 != 0)) {
    t2 <- TRUE
  } else {
    t2 <- FALSE
  }
  
  if (t1 | t2) {
    mss <- paste0("bp.to.cM = T but non-integer distance values detected in ",
                  paste0(c("obj.info", "rec.rate")[c(t1, t2)], 
                         collapse = " and "),
                  "\nConverstion will be performed but set bp.to.cM = F if ",
                  "distances are already genetic units (i.e. centiMorgans)")
    WMSS <- c(WMSS, mss)
  }

  # get appropriate map function
  if (map.fun == 1) mf <- function(rr) qtl:::mf.h(rr * 1e6) 
  if (map.fun == 2) mf <- function(rr) qtl:::mf.k(rr * 1e6)
  if (map.fun == 3) mf <- function(rr) qtl:::mf.cf(rr * 1e6)
  if (map.fun == 4) mf <- function(rr) qtl:::mf.m(rr * 1e6)
  
  # compute genetic distances
  mf.names <- c("Haldane","Kosambi","Carter-Falconer","Morgan")
  vrb(paste("Performing physical to genetic distances coversion using the", 
            mf.names[map.fun], "map function\n"))
  
  rec.rate[, `:=` (START = startpos, END = endpos)]
  rec.rate[, `:=` (startpos = NULL, endpos = NULL)]
  rec.rate[, pos.l := END - START]
  rec.rate[, cM.l := mf(rr = rec.rate) * (END - START) / 100, by = chr]
  rec.rate[, endcM := cumsum(cM.l), by = chr]
  rec.rate[, startcM := c(0, endcM[-.N]), by = chr]
  data.table::setkey(rec.rate, chr, START, END)
  
  # check that sufficient coverage on each genome
  rr.range <- rec.rate[, .(MIN = min(START), MAX = max(END)), by = chr]
  oi.range <- obj.info[, .(MIN = min(startpos), MAX = max(endpos)), by = chr]
  rr <- merge(oi.range, rr.range, by = "chr")
  wR <- rr[, which(MIN.x < MIN.y | MAX.x > MAX.y)]
  if (length(wR) > 0) {
    warning("rec.rate boundaries do not span all obj.info gene boundaries ",
            "for chromosome", ifelse(length(wR) > 1, "s ", " "),
            paste(wR, collapse = ", "), "\nGenes lying outside of ",
            "rec.rate chromosome boundaries will be lost")
  }
  
  # impute obj.info file
  #start pos
  vrb("Converting gene boundaries from bp to cM for obj.info\n")
  obj.info[, `:=` (startpos.orig = startpos, endpos.orig = endpos)]
  oi.names <- names(obj.info)
  obj.info <- data.table::foverlaps(obj.info, rec.rate, type = "within",
                                    mult = "first",
                                    by.x = c("chr", "startpos", "startpos.orig"))
  obj.info[, startpos := startcM + cM.l * ((startpos - START) / pos.l)]
  obj.info <- obj.info[, .SD, .SDcols = oi.names]
  
  #end pos
  obj.info <- data.table::foverlaps(obj.info, rec.rate,  type = "within",
                                    mult = "first",
                                    by.x = c("chr", "endpos", "endpos.orig"))
  obj.info[, endpos := startcM + cM.l * ((endpos - START) / pos.l)]
  obj.info <- obj.info[, .SD, .SDcols = oi.names]
  assign(x = "oi.copy", envir = env, value = obj.info)
  
  if (!is.null(var.info)) {
    vi.range <- var.info[, .(MIN = min(pos), MAX = max(pos)), by = chr]
    rr <- merge(vi.range, rr.range, by = "chr")
    wR <- rr[, which(MIN.x < MIN.y | MAX.x > MAX.y)]
    if (length(wR) > 0) {
      warning("rec.rate boundaries do not span all var.info gene boundaries ",
              "for chromosome", ifelse(length(wR) > 1, "s\n", "\n"),
              paste(wR, collapse = ", "), "\nVariants lying outside of ",
              "rec.rate chromosome boundaries will be lost")
    }
    vrb("Converting gene boundaries from bp to cM for var.info\n")
    # impute var.info file
    var.info[, pos.orig := pos]
    var.info <- data.table::foverlaps(var.info, rec.rate, type = "within",
                                      mult = "first",
                                      by.x = c("chr", "pos", "pos.orig"))
    var.info[, pos := startcM + cM.l * ((pos - START) / pos.l)]
    var.info <- var.info[, .(chr, pos, value, pos.orig)]
    assign(x = "var.info", envir = env, value = var.info)
  }
  vrb("\n")
}

# function to determine final set merging results
final.merge <- function(I, J) { 
  # find sets that have been merged multiple times
  mgm <- merge(I, J, by.x = "setID.merge", by.y = "setID", all.y = TRUE)
  K <- mgm[, .(setID = ifelse(is.na(setID), setID.merge, setID),
               TEMP = setID.merge.y)]
  out <- merge(I, K, by = "setID", all.x = TRUE)
  out[, .(setID, setID.merge = ifelse(is.na(TEMP), setID.merge, TEMP))]
}

# function that determines subnetwork of genes sets sharing min.sim proportion genes
get.subnetworks <- function(SO, min.sim) {
  PxG <- Matrix::sparseMatrix(i = SO$setID, j = SO$objID, x = 1L)
  
  # compute number of shared genes
  sim.mat <- as.matrix(Matrix::tcrossprod(PxG))
  
  # create adjacency matrix
  adj.mat <- (sim.mat / diag(sim.mat)) >= min.sim
  
  # all gene sets with > min.sim genes in both cases
  igg <- igraph::graph_from_adjacency_matrix(adj.mat + t(adj.mat) == 2)
  
  # find connected components
  sgg <- igraph::components(igg)$membership
  
  # keep sets with genes
  as.numeric(as.factor(sgg[diag(sim.mat) != 0]))
}

# merge gene sets based proportion of shared genes
merge.sim.sets <- function(SI, SO, min.sim, env) {
  # initialise
  R <- 0
  SI[, setG := get.subnetworks(SO = SO, min.sim = min.sim)]
  SO.0 <- data.table::copy(SO)
  merged.sets <- list()
  
  #repeat until all gene sets fulfill gene sharing requirements 
  while(SI[, max(setG) < .N]){
    R <- R + 1
    SO.0 <- merge(SO.0, SI[, .(setID, setG)], by = "setID")
    
    #rename merged sets by set with most genes, or sort alphanumerically if ties
    so <- SO.0[, .(setG = unique(setG), .N), by = setID]
    so.max <- so[, .(N.merged = .N, setID.merge = setID[which.max(N)]), 
                 by = setG]
    set.merge <- merge(so[, .(setG, setID)], so.max, by = "setG")
    merged.sets[[R]] <- set.merge[N.merged > 1, .(setID, setID.merge)] 
    
    #create new files
    so.temp <- merge(set.merge[!duplicated(setG), .(setG, setID.merge)], 
                     SO.0, by = "setG")
    so.temp <- so.temp[, .(objID = unique(objID)), by = setID.merge] # remove gene dups
    so.temp[, setID := setID.merge]
    so.temp[, setID.merge := NULL]
    SO.0 <- so.temp[order(setID, objID), .(setID, objID)]
    
    SI <- SI[setID %in% SO.0$setID]
    SI[, N := SO.0[, .N, by = setID][order(setID)]$N]
    
    #check that new sets fulfill merging requirements
    SI[, setG := get.subnetworks(SO = SO.0, min.sim = min.sim)]
  }
  
  #prepare output
  if (R == 0) {
    SO.merged <- NA
  } else {
    SO <- SO.0
    SO.merged <- Reduce(final.merge, merged.sets)
    SI[setID %in% SO.merged$setID, setName := paste0(setName, "!")]
  }
  SI[, setG := NULL]
  
  list2env(list(so.copy = SO, si.copy = SI, set.info.merged = SO.merged, R = R),
           envir = env)
}

# record cumulative run time at specific points
get.time <- function(start.time) {
  tt <- ceiling(as.numeric(difftime(Sys.time(), start.time, units = "secs")))
  H <- tt %/% 3600
  rr <- tt %% 3600
  if(rr > 0) {
    M <- rr %/% 60
    rr2 <- rr %% 60
    if(rr2 > 0) {
      S <- round(rr2)
    } else {
      S <- 0
    }
  } else {
    M <- 0
    S <- 0
  }
  return(paste0(H, "h ", M, "m ", S, "s"))
}

# calculate Mahalanobis distances for all genes in gene sets vs all genes
# uses Cholesky factorisation for efficiency
# core algorithm from: https://stats.stackexchange.com/questions/65705/pairwise-mahalanobis-distances
chol.md <- function(X) {
  X <- as.matrix(X)
  if(ncol(X) == 1) {
    mahal.dist <- as.matrix(dist(X))
  } else {
    cX <- Rfast::cova(X)
    # check for singularity
    qrX <- qr(cX)
    if (qrX$rank < ncol(X)) {
      col.keep <- qrX$pivot[1:qrX$rank]
      cX <- cX[col.keep, col.keep]
      col.rem <- setdiff(1:ncol(X), col.keep)
      warning("Covariate matrix is singular\nRedundant covariates [", 
              paste0("Cov", col.rem, collapse = ", "),
              "] will not be used in Mahalanobis Distance calculations", 
           call. = FALSE, immediate. = T)
    } else {
      decX <- chol(cX)
      tmp <- forwardsolve(t(decX), t(X))
      mahal.dist <- as.matrix(dist(t(tmp)))
    }
  }
  return(mahal.dist)
}

# generate matched gene weights [VERBOSE?]
get.obj.mw <- function(md, set.genes, match.cap, match.scale = NULL, 
                       match.fun = NULL, env) {
  # check user supplied objects
  if (!is.null(match.fun)) {
    emss <- "Inappropriate user-defined scaling function\n"
    md.trans <- match.fun(md)
    
    if (!is.matrix(md.trans)) {
      stop(emss, "Returned object is not a matrix.", call. = FALSE)
    } else {
      if (any(dim(md.trans) == dim(md))) {
        stop(emss, "Returned object does not preserve", nrow(md), 
             "row by", ncol(md), "column matrix structure.", 
             call. = FALSE)
      }
    }
  } else { # use default (inverse scaling)
    md.trans <- exp(-match.scale * md)
    if (length(set.genes) < ncol(md)) {
      md.trans[cbind(1:length(set.genes), set.genes)] <- 0
    } else {
      diag(md.trans) <- 0
    }
  }
  
  # scale weights to sum to 1
  md.trans <- md.trans / Rfast::rowsums(md.trans) 
  
  # Optional: apply probability cap
  if (match.cap < 1) {
    wGT0 <- which(Rfast::rowAny(md.trans > match.cap))
    if (length(wGT0) > 0){
      mdT <- md.trans[wGT0, ]
      wGT <- which(mdT >= match.cap, arr.ind = TRUE)
      while (nrow(wGT) > 0) {
        mdT[wGT] <- NA
        ns <- 1 - Rfast::rowTrue(is.na(mdT)) * match.cap
        ds <- Rfast::rowsums(mdT, na.rm = TRUE)
        rs <- ns / ds
        mdT <- mdT * rs
        wGT <- which(mdT > match.cap, arr.ind = TRUE)
      }
      
      # create output matrix
      mdT[is.na(mdT)] <- match.cap
      md.trans[wGT0, ] <- mdT
      mw.mss <- paste0(" capped at ", match.cap, " (", length(wGT0), 
                       " set genes impacted)")
    } else {
      mw.mss <- vrb(paste0("; match.cap not applied (no weight > ", 
                           match.cap, ")"))
    }
  } else {
    mw.mss <- ""
  }
  assign(x = "mw.mss", envir = env, value = mw.mss)
  
  return(md.trans)
}

# generate weighted sample
get.wt.samp <- function(mw, ng, np) {
  wt.s <- collapse::dapply(mw, MARGIN = 1, FUN = function(P) {
    sample.int(ng, np, replace = TRUE, prob = P)
  })
  c(t(wt.s))
}

# create genomic blocks for autocovariance estimation
ac.blocks <- function(score, block.size) {
  #split chromosomes
  out <- split(score, f = score$chr)
  if (!is.infinite(block.size)) {
    half.bs <- block.size / 2
    out <- foreach::foreach(ss = out, .combine = c) %do% {
      if (nrow(ss) >= block.size + half.bs){
        data.table::setkey(ss, I)
        OID <- ss$I
        b1.seq <- seq(min(OID), max(OID) + block.size, block.size) - 1
        b1 <- split(OID, f = cut(OID, b1.seq))
        b2.seq <- seq(min(OID) + half.bs, max(OID) + block.size, block.size) - 1
        b2 <- split(OID, f = cut(OID, b2.seq))
        
        #remove partial blocks
        n1 <- length(b1)
        n2 <- length(b2)
        if (all(b1[[n1]] %in% b2[[n2]]) | all(b2[[n2]] %in% b1[[n1]])) {
          if (length(b1[[n1]]) > length(b2[[n2]])) {
            b2[[n2]] <- NULL
          } else {
            b1[[n1]] <- NULL
          }
        }
        
        blocks <- list()
        blocks[seq(1, 2 * length(b1), by = 2)] <- b1
        blocks[seq(2, 2 * length(b2), by = 2)] <- b2

        #combine small blocks
        threeq.bs <- 0.75 * block.size
        n0 <- length(blocks)
        if (length(blocks[[n0]]) < threeq.bs) {
          blocks[[n0 - 1]] <- sort(unique(c(blocks[[n0 - 1]], blocks[[n0]])))
          blocks[[n0]] <- NULL
        }
        
        OUT <- foreach::foreach(bl = blocks) %do% ss[.(bl)]
      } else {
        OUT <- list(ss)
      }
      OUT
    }
  }
  return(out)
}

# fit variogram function
fit.vg <- function(x, v.now, vg.fit) {
  bf <- gstat::fit.variogram(object = v.now, fit.method = vg.fit,
                             model = gstat::vgm(, "Mat", kappa = x, 
                                                nugget = NA))
  assign(x = "bf", envir = parent.frame(n = 3), value = bf)
  attr(bf, "SSErr")
}

# estimate autocovariance across pairs of genes
vg.est <- function(ss.in, max.link, bin.size, vg.fit, kappa) {
  # estimate empirical variogram
  v.now <- gstat::variogram(OS ~ 1, locations = ~ midpos + midpos.dummy, 
                            data = ss.in, cressie = TRUE, cross = FALSE,
                            width = bin.size, pseudo = 0, cutoff = max.link)
  
  # hill climbing algorithm for best fitting 
  sup.all(optimize(f = fit.vg, interval = kappa, v.now = v.now, 
                   vg.fit = vg.fit))

  return(list(bf, v.now))
}

# create parts of autocovariance object
make.ac <- function(p2i, vgm.est, nug, MD) {
  pair.dist <- c(as.matrix(dist(p2i$midpos)))
  GxG <- expand.grid(1:nrow(p2i), 1:nrow(p2i))
  
  if (!nug) { # set nugget to zero (only use spatial variance)
    vgm.est[1, 2] <- 0
  }
  
  #generate covariance estimates for each gene pair separated by < min.dist
  v.est <- gstat::variogramLine(vgm.est, dist_vector = pair.dist, 
                                maxdist = MD, covariance = TRUE)
  
  out <- data.table::data.table(objID.A = p2i$I[GxG[, 2]], 
                                objID.B = p2i$I[GxG[, 1]],
                                midpos.dist = pair.dist,
                                covar = v.est$gamma)
  
  return(out)
}

# function to estimate design effect and rescale set scores
ss.rescale <- function(pm, AC, OID, SID) {
  # generate data.table of all unique ac values for each set
  pmu <- unique(pm)
  ACx <- AC[.(pmu)][.(pmu), on = "B"]
  ACy <- data.table::data.table(C = SID, X = pm[OID])[, .N, by = c("C", "X")]
  ACa <- ACx[ACy, on = .(A = X), allow.cartesian = TRUE]
  ACb <- ACa[ACy, on = .(B = X, C == C)]
  
  ACc <- rbind(ACb, ACb[A != B, .(A = B, B = A, covar, C, N = i.N, i.N = N)])
  
  # sum diagonal and full components, allowing for duplication
  X <- ACc[A == B, sum(N * covar), by = C]$V1
  Y <- ACc[, sum(N * i.N * covar), by = C]$V1
  
  # return ratio
  X / Y 
}

# calculate the maximum empirical p value with lower boundary <= 1/n
est.tail.p <- function(x0, X, gdp.th, n.tail, fit.meth) {
  check <- TRUE
  N <- length(X)
  while (check & gdp.th >= n.tail) { # decrease exceedences if no convergence
    x.th.prop <- gdp.th / N
    x.th <- quantile(X, 1 - x.th.prop) # threshold value
    gdp.fit <- try(fExtremes::gpdFit(x = X[X >= x.th], #exceedences
                                     u = x.th, fit = fit.meth)@fit$par.ests,
                   silent = T)
    check <- "try-error" %in% class(gdp.fit)
    if (check) gdp.th <- gdp.th - 10
  }
  gpd.p <- fExtremes::pgpd(q = x0 - x.th, xi = gdp.fit[1], 
                           beta = gdp.fit[2], lower.tail = FALSE)[[1]]
  gpd.p * x.th.prop
}

# correct set scores for shared genes using linear regression
prune.sets.reg <- function(SSE, SSO, p0, SOC, mss, CM, n.tail, fit.meth) {
  # keep track top gene sets and their genes
  n.sets <- length(SSO)
  si.now <- matrix(0, nrow = n.sets, ncol = 6)
  set.rem <- 1:n.sets
  n.set.rem <- length(set.rem)
  n.obj.rem <- data.table::uniqueN(SOC$X)
  np.update <- cbind(SOC[, .N, by = A], p = p0)
  data.table::setkey(np.update, A)
  N.tail <- n.tail + 1
  Np <- ncol(SSE) + 1
  K <- 0
  
  while (n.set.rem > 1) { # continue while more than one valid set remains
    K <- K + 1 # increment

    ts.now <- np.update[p == min(p)]
    if (nrow(ts.now) > 1) { # use largest standardised set score in case of ties
      s0.sd <- sqrt(CM[.(ts.now$A, ts.now$A)]$C)
      ts.now <- ts.now[which.max(SSO[ts.now$A] / s0.sd)]
    }
    
    topSet <- ts.now$A # top gene set ID
    sso.now <- SSO[topSet] # top gene set score
    
    #determine genes and sets to update
    covXY <- CM[.(topSet)] # covariance with top scoring gene set
    
    # update step
    if (nrow(covXY) > 1) { # only update if correlated gene sets exist
      # estimate regression coefficient
      cov.xy <- covXY$C
      var.x <- covXY[.(topSet), on = "B"]$C
      b <- cov.xy / var.x 
      
      # update gene set scores
      ss.u <- covXY$B # rows to modify in set score matrices
      sse.now <- SSE[topSet, ]
      sse.new <- SSE[ss.u, ] - tcrossprod(b, sse.now)
      sso.new <- SSO[ss.u] - b * sso.now
      SSE[ss.u, ] <- sse.new
      SSO[ss.u] <- sso.new
      
      # update p values
      p.update <- Rfast::rowTrue(sse.new >= sso.new) + 1
      if (!is.null(fit.meth) & min(p.update) <= N.tail) { # estimate tail probability using GPD
        for (p.i in which(p.update <= N.tail)) {
          s.i <- ss.u[p.i]
          p.update[p.i] <- est.tail.p(x0 = SSO[s.i], X = SSE[s.i, ], 
                                      gdp.th = 200, fit.meth = fit.meth, 
                                      n.tail = n.tail) * Np # put on same scale as empirical p values
        }
      }
      np.update[.(ss.u), p := p.update]
      
      # update covariance matrix
      cm.ij <- CM[data.table::CJ(c(ss.u), c(ss.u)), nomatch = NULL]
      cov.ix <- cov.xy[fastmatch::fmatch(cm.ij$A, ss.u)]
      cov.jx <- cov.xy[fastmatch::fmatch(cm.ij$B, ss.u)]
      cm.ij[, C := C - cov.ix * cov.jx / var.x] # decrement covariance
      CM[.(cm.ij), C := i.C]
    }
    
    # update looping variables and remove unused genes and sets
    obj.out <- SOC[.(topSet)]$X
    n.new <- SOC[.(obj.out), on = "X"][, .N, by = A] # no. of shared genes
    np.update[.(n.new), N := N - i.N] # update gene set size
    np.update <- np.update[N >= mss]
    set.in <- np.update$A # sets with sufficient genes
    set.out <- setdiff(set.rem, set.in) # remove gene sets with insufficient genes
    set.rem <- set.in # update remaining gene sets
    CM <- CM[!.(set.out)][!.(set.out), on = "B"]
    SOC <- SOC[!.(set.out)][!.(obj.out), on = "X"]
    n.set.rem <- length(set.rem)
    n.obj.rem <- n.obj.rem - length(obj.out)
    
    # keep track of successive top sets
    si.now[topSet, ] <- c(K, ts.now$N, n.obj.rem, n.set.rem, sso.now, ts.now$p)
  }
  
  si.now
}

# correct set scores for shared genes by gene remvoal
# only possible for unscaled gene sets
prune.sets.rem <- function(SSE, SSO, p0, SOC, mss, PxG, PI, PSM, 
                           n.tail, fit.meth) {
  # keep track top gene sets and their genes
  n.sets <- n.set.rem <- length(SSO)
  n.obj.rem <- ncol(PSM)
  si.now <- matrix(0, nrow = n.sets, ncol = 6)
  set.rem <- 1:n.sets
  np.update <- cbind(SOC[, .N, by = A], p = p0)
  data.table::setkey(np.update, A)
  N.tail <- n.tail + 1
  Np <- ncol(SSE) + 1
  K <- 0
  
  while (n.set.rem > 1) { # continue while more than one valid set remains
    K <- K + 1 # increment
    
    # determine most significant gene set
    ts.now <- np.update[p == min(p)]
    if (nrow(ts.now) > 1) { # use largest standardised set score in case of ties
      s0 <- SSE[ts.now$A, ]
      s0.hat <- Rfast::rowmeans(s0)
      s0.sd <- Rfast::rowVars(s0, std = TRUE)
      ts.now <- ts.now[which.max((SSO[ts.now$A] - s0.hat) / s0.sd)]
    }
    
    topSet <- ts.now$A # top gene set ID
    ts.n <- ts.now$N
    sso.now <- SSO[topSet] / ts.n # top gene set score

    #determine genes and sets to update
    nXY <- PI[.(topSet)] # shared genes with top scoring gene set
    obj.out <- SOC[.(topSet)]$X # genes to remove
    
    if (nrow(nXY) > 1) {  # update gene set scores
      ss.u <- nXY$B # rows to modify in set score matrices
      PxG.now <- PxG[ss.u, obj.out]
      ss.dec <- as.matrix(Matrix::tcrossprod(PxG.now, PSM[, obj.out]))
      sse.new <- SSE[ss.u, ] - ss.dec[, -1] # revised expected scores
      sso.new <- SSO[ss.u] - ss.dec[, 1] # revised observed scores
      SSE[ss.u, ] <- sse.new
      SSO[ss.u] <- sso.new
      
      # update p values
      p.update <- Rfast::rowTrue(sse.new >= sso.new) + 1
      if (!is.null(fit.meth) & min(p.update) <= N.tail) { # estimate tail probability using GPD
        for (p.i in which(p.update <= N.tail)) {
          s.i <- ss.u[p.i]
          p.update[p.i] <- est.tail.p(x0 = SSO[s.i], X = SSE[s.i, ], 
                                      gdp.th = 200, fit.meth = fit.meth, 
                                      n.tail = n.tail) * Np
        }
      }
      np.update[.(ss.u), p := p.update]
      
      # update count matrix
      PxN.now <- as.matrix(Matrix::tcrossprod(PxG.now))
      PxN.gt0 <- which(PxN.now > 0, arr.ind = TRUE) # non-zero elements
      PI.update <- data.table::data.table(A = ss.u[PxN.gt0[, 2]], 
                                          B = ss.u[PxN.gt0[, 1]],
                                          N = PxN.now[PxN.gt0])
      PI[.(PI.update), N := N - i.N]
    }
    
    # remove unused genes and sets
    np.update[.(nXY[, .(A = B, N)]), N := N - i.N]
    np.update <- np.update[N >= mss]
    set.in <- np.update$A # sets with sufficient genes
    set.out <- setdiff(set.rem, set.in) # remove gene sets with insufficient genes
    set.rem <- set.in # update remaining gene sets
    PI <- PI[!.(set.out)][!.(set.out), on = "B"][N > 0] # remove unused gene sets
    SOC <- SOC[!.(set.out)][!.(obj.out), on = "X"] # update genes in sets
    n.set.rem <- length(set.in)
    n.obj.rem <- n.obj.rem - ts.n
    
    # keep track of successive top sets
    si.now[topSet, ] <- c(K, ts.n, n.obj.rem, n.set.rem, sso.now, ts.now$p)
  }
  
  si.now
}

# compute number of p values falling within specified intervals
fdr.bin.counts <- function(O, E, p.bins, n.fdr) {
  O[, FDR.bin := cut(setScore.pr.p, p.bins)]
  rp <- O[, .(OBS = .N), by = FDR.bin]
  
  n.rs <- nrow(O) / n.fdr # expected count rescaling factor
  E[, FDR.bin := cut(setScore.pr.p, p.bins)]
  e.bin <- E[, .N, by = c("FDR.rep", "FDR.bin")]
  e.bin[, N.tot := sum(N), by = FDR.rep]
  vp <- e.bin[!is.na(FDR.bin), .(EXP = sum(N / N.tot) * n.rs), by = FDR.bin]
  
  #merge and fill in empty bins
  out <- merge(rp, vp, by = "FDR.bin", all = TRUE)
  out[is.na(OBS), OBS := 0] 
  out[is.na(EXP), EXP := 0]
  return(out)
}

# Matern covariance function
# from https://scikit-gstat.readthedocs.io/en/latest/reference/models.html, using full range
matern <- function(kappa, nugget, range, psill, x) {
  kx <- x / range
  A <- 1 / (2 ^ (kappa - 1) * gamma(kappa))
  B <- kx ^ kappa
  C <- besselK(kx, nu = kappa)
  return(nugget + psill * (1 - A * B * C))
}

# determine optimal facet arrangement
facet.opt <- function(OPT, fac.per.row = 3:6) {
  if (OPT <= min(fac.per.row)) {
    dim1 <- 1
    dim2 <- OPT
  } else {
    rm1 <- OPT %% fac.per.row
    if (all(rm1 != 0)) { # get no. redundant facets
      fac.tot <- ceiling(OPT / fac.per.row) * fac.per.row 
      rm1 <- fac.tot - OPT
    }
    dim1 <- fac.per.row[rm1 == min(rm1)]
    if (length(dim1) > 1) {
      dim1 <- dim1[which.min(abs(dim1 - sqrt(OPT)))]
    }
    dim2 <- ceiling(OPT / dim1)
  }
  out <- sort(c(dim1, dim2))
  return(list(dim1 = out[1], dim2 = out[2]))
}


#===========================================================================
# Function: plR.read(input.path, min.set.n, max.set.n, si.path, oi.path, so.path, 
#                    merge.set.prop, rem.ident.genes, obj.in, obj.out, set.in, set.out, 
#                    group, verbose)
#
# Read in all required files (obj.info, set.info, set.obj)
#
# - in.path        : character; path to directory with input files
#                    default = getwd(), current working director
#                    Will search for the following files:
#                    - objInfo: tab separated file with fields (those in square brackets are optional):
#                      objID (integer), [objName (character)], objStat (numeric), chr (numeric), startpos (numeric), endpos (numeric), [Cov1 (numeric)], [Cov2 (numeric)], [Cov (numeric) ...]
#                    - setInfo: tab separated file with fields (those in square brackets are optional):
#                      setID (integer), [setName (character)], [N (integer], [setSource (character)]
#                    - setObj : tab separated file with fields:
#                      setID (integer), objID (integer)
#                    - VarInfo : tab separated file with fields:
#                      chr (numeric), pos (numeric), value (numeric)
#                    - RecRate : tab separated file with fields:
#                      chr (numeric), startpos (numeric), endpos (numeric), rec.rate (numeric)
# - min.set.n      : integer; exclude gene sets with size below min.set.n
#                    default = 2; range [2, max.set.n)
# - max.set.n      : integer; exclude gene sets with size above max.set.n
#                    default = Inf; range [min.set.n, Inf)
# - merge.set.prop : numeric; minimum proportion of shared genes used to concatenate gene sets
#                    default = 0.95; range (0, 1]
# - rem.ident.genes: logical; should genes with identical start and end positions be removed?
#                    default = FALSE
# - obj.in         : integer vector; objID of genes to retain
#                    default = NULL
# - obj.out        : integer vector; objID of genes to remove
#                    default = NULL
# - set.in         : integer vector; setID of gene sets to retain
#                    default = NULL
# - set.out        : integer vector; setID of gene sets to remove
#                    default = NULL
# - group          : character; if in.path option is used, label specifying which set of files to load 
#                    (only needed if multiple file sets are present)
#                    default = NULL (assumes that only one set of files is present)
# - verbose.       : logical, should progress reporting be enabled?
#                  : default = TRUE
#===========================================================================

plR.read <- function(input.path = getwd(), min.set.n = 2, max.set.n = Inf,
                     group = NULL, obj.in = NULL, obj.out = NULL, set.in = NULL, 
                     set.out = NULL, merge.set.prop = 0.95, rem.ident.genes = FALSE,
                     bp.to.cM = FALSE, map.fun = 1, verbose = TRUE) {
  
  startT <- Sys.time()
  suppressMessages(library(foreach))
  suppressMessages(library(data.table))
  
  # check arguments
  arg.check(plr.func = "read", env = environment())
  
  #----------------------------------------------------------------------------#
  ## Read and check plR input files
  #----------------------------------------------------------------------------#
  
  file.load(input.path = input.path, group = group, bp.to.cM = bp.to.cM, 
            env = environment())

  # run internal file checks
  file.check(set.info = set.info, obj.info = obj.info, set.obj = set.obj,
             var.info = var.info, rec.rate = rec.rate, update.indices = FALSE,
             env = environment())
  
  #----------------------------------------------------------------------------#
  ## Filtering data
  #----------------------------------------------------------------------------#
  
  #NOTE: removal of sets will not remove orphan genes from gene list
  # since these can still be used for permutations
  # however, gene removal can result in removal of sets 
  # if all genes are removed, or when gene sets fall below a specific size
  km.ch <- list(obj.in, obj.out, set.in, set.out) 
  w.ch <- sapply(km.ch, is.numeric)
  if (any(w.ch)) {
    W.ch <- which(w.ch)
    nor0 <- nrow(oi.copy)
    sor0 <- nrow(si.copy)
    nsr.d <- nor.d <- 0
    if (1 %in% W.ch) {
      oi.copy <- oi.copy[objID %in% obj.in]
      nor <- nor0 - nrow(oi.copy)
      vrb(paste0("obj.in option enacted: ", nor, " gene", 
                 ifelse(nor == 1, "", "s"), " removed\n"))
      nor.d <- length(obj.in) - nrow(oi.copy)
    }
    if (2 %in% W.ch) {
      oi.copy <- oi.copy[!(objID %in% obj.out)]
      nor <- nor0 - nrow(oi.copy)
      vrb(paste0("obj.out option enacted: ", nor, " gene", 
                 ifelse(nor == 1, "", "s"), " removed\n"))
      nor.d <- length(obj.out) - nor
    }
    if (nor.d > 0) {
      vrb(paste0("[NB: ", nor.d, " gene ID", 
                 ifelse(nor.d == 1, " was", "s were"),
                 " not detected in obj.info]\n"))
    }
    # update gene sets
    si.copy <- si.copy[setID %in% so.copy$setID]
    
    if (3 %in% W.ch) {
      si.copy <- si.copy[setID %in% set.in]
      nsr <- sor0 - nrow(si.copy)
      vrb(paste0("set.in option enacted: ", nsr, " gene set ID", 
                 ifelse(nsr == 1, "", "s"), " removed\n"))
      nsr.d <- length(set.in) - nrow(si.copy)
    }
    if (4 %in% W.ch) {
      si.copy <- si.copy[!(setID %in% set.out)]
      nsr <- sor0 - nrow(si.copy)
      vrb(paste0("set.out option enacted: ", nsr, " gene set", 
                 ifelse(nsr == 1, "", "s"), " removed\n"))
      nsr.d <- length(set.out) - nsr
    }
    if (nsr.d > 0) {
      vrb(paste0("[NB: ", nsr.d, " gene set ID", 
                 ifelse(nsr.d == 1, " was", "s were"),
                 " not detected in set.info]\n"))
    }
    
    # check that at least one set remains
    if (nrow(si.copy) == 0) {
      stop("All gene sets removed\nCheck filtering criteria", call. = FALSE)
    }
    
    # update set.obj
    so.copy <- so.copy[setID %in% si.copy$setID][objID %in% oi.copy$objID]
    
    vrb("\n")
  }
  
  # check set size filtering criteria
  maxN <- max(so.copy[, .N, by = setID]$N)
  if (min.set.n > maxN) {
    stop("min.set.n exceeds no. genes in largest set [N = ", maxN, 
         "] after filtering\nCheck filtering criteria", call. = FALSE)
  }
  
  
  #----------------------------------------------------------------------------#
  ## Cleaning data
  #----------------------------------------------------------------------------#
  
  # clean var.info (can be slow for large files)
  if (!is.null(var.info)) {
    data.table::setkey(var.info, chr, pos)
    var.dup <- var.info[, .N, by = c("chr", "pos")][N > 1, .(chr, pos)]
    n.dup <- nrow(var.dup)
    vrb(paste0(nrow(var.dup), " duplicated variants detected\n"))
    if (nrow(var.dup) > 1) {
      vv <- var.info[.(var.dup)]
      var.info <- rbind(vv[duplicated(chr, pos)], # take single instance of each duplicated variant
                        var.info[!(var.dup)])[order(chr, pos)]
      vrb(paste0(nrow(var.info), 
                 " variants remaining after duplicate removal\n\n"))
    }
  }
  
  # remove genes with no information
  no.scores <- oi.copy[is.na(objStat), objID]
  nsl <- length(no.scores)
  if (nsl > 0) {
    oi.copy <- oi.copy[!(objID %in% no.scores)]
    so.copy <- so.copy[!(objID %in% no.scores)]
    vrb(paste0("Removed ", nsl, " gene", ifelse(nsl == 1, " ", "s "), 
                   "with no information\n\n"))
  }
  
  # check for duplicated genes
  obj.dups <- oi.copy[, .(.N, objID), by = c("chr", "startpos", "endpos")][N > 1]
  
  if (nrow(obj.dups) > 0) {
    vrb(paste0("Found ", nrow(obj.dups), " genes that have an identical ",
               "genomic position with other genes\n"))
    if (rem.ident.genes) { # Optional: remove genes with identical start and end positions on same chromosome
      # privilege retention of genes in gene sets (or appear in most sets)
      obj.dups[, ID := as.numeric(as.factor(paste(chr, startpos, endpos)))]
      obj.dups <- merge(obj.dups[, .(objID, ID)], 
                        so.copy[objID %in% obj.dups$objID, .N, by = objID], 
                        by = "objID", all.x = TRUE) # enumerate number of times gene occurs in sets
      obj.dups[is.na(N), N := 0]
      
      obj.dups[, objID.keep := objID[which.max(N)], by = ID]
      obj.keep <- unique(obj.dups$objID.keep)
      obj.rem <- setdiff(obj.dups$objID, obj.keep)
      
      oi.copy <- oi.copy[!(objID %in% obj.rem)]
      so.copy <- so.copy[!(objID %in% obj.rem)] 
      
      obj.info.merged <- obj.dups[, .(objID, objID.keep, N.sets = N)]
      vrb(paste0(length(obj.keep), " unique genes remain after merging\n",
                 "[merge info available in oi.copy.merge]\n\n"))
    } else {
      obj.info.merged <- NULL
      vrb("rem.ident.genes = F, identical genes will be retained\n\n")
    }  
  } else {
    obj.info.merged <- NULL
  }
  
  # check that all genes are unique
  so.check <- so.copy[, .N, by = c("setID", "objID")][setID %in% setID[N > 1]]
  so.n <- data.table::uniqueN(so.check$setID)
  if (so.n > 0) {
    vrb(paste0("Indentified ", so.n, " gene set", ifelse(so.n == 1, "", "s"),
               " with redundant genes, forcing all gene sets to have unique objIDs\n",
               "[info available in repeat.objID]\n\n"))
    repeat.objID <- so.check
    so.copy <- so.copy[, .(objID = unique(objID)), by = setID] # recreate
  } else {
    repeat.objID <- NULL
  }
  
  # update set counts
  si.copy[, N := so.copy[, .N, by = setID][order(setID)]$N]
  
  # merge similar sets
  n.sets.orig <- nrow(si.copy)
  merge.sim.sets(SI = data.table::copy(si.copy), 
                 SO = data.table::copy(so.copy), 
                 min.sim = merge.set.prop, env = environment())
  
  if (R == 0) {
    vrb(paste("No gene sets >", merge.set.prop, "similarity\n\n"))
  } else {
    vrb(paste0("Performed ", R, ifelse(R == 1, " round", " rounds"), 
               " of gene set merging\n"))
    vrb(paste0(nrow(set.info.merged), " gene sets with >= ", 
               merge.set.prop * 100, "% similarity were merged into ",  
               data.table::uniqueN(set.info.merged$setID.merge), " gene set",
               ifelse(n.sets.orig - nrow(si.copy) == 1, "\n", "s\n"),
               "[merge info available in set.info.merge]\n\n"))
  }
  
  # Optional: remove gene sets that have too many/too few genes (and orphan genes)
  if (min.set.n > 1 | !is.infinite(max.set.n)){
    fmss <- paste0("with >= ", min.set.n, " and <= ", max.set.n, " genes")
    s.out <- si.copy[N < min.set.n | N > max.set.n]$setID
    if (length(s.out) == nrow(si.copy)) {
      stop("All gene sets removed\nCheck filtering criteria", call. = FALSE)
    } else {
      # update objects
      removed.sets <- si.copy[setID %in% s.out] # keep track of removed sets
      si.copy <- si.copy[!(setID %in% s.out)]
      so.copy <- so.copy[!(setID %in% s.out)]
    }
  } else {
    fmss <- ""
    removed.sets <- NULL
  }
  
  #----------------------------------------------------------------------------#
  # Optional: convert genome coordinates from physical (bp) to genetic (cM)
  #----------------------------------------------------------------------------#
  
  if (bp.to.cM) {
    if (is.null(rec.rate)) {
      warning("bp.to.cM = T but no RecRate file provided\n",
              "Physical to genetic distance coversion not performed")
    } else {
      bp2cM(obj.info = oi.copy, var.info = var.info, 
            rec.rate = rec.rate, map.fun = map.fun, 
            env = environment())
      rm(rec.rate); gc()
    }
  }
  
  #----------------------------------------------------------------------------#
  ## Prepare input for polylinkR
  #----------------------------------------------------------------------------#
  
  read.args <- list(min.set.n =min.set.n, max.set.n = max.set.n,
                    obj.in = obj.in, obj.out = obj.out, set.in = set.in, 
                    set.out = set.out, merge.set.prop = merge.set.prop, 
                    rem.ident.genes = rem.ident.genes, bp.to.cM = bp.to.cM, 
                    map.fun = map.fun)
    
  vrb(paste0("Loaded ", nrow(si.copy), " gene sets ", fmss, " and ",
             nrow(oi.copy), " unique genes in total\n"))
  
  obj.nx <- nrow(oi.copy) - data.table::uniqueN(so.copy$objID)
  if(obj.nx > 0) {
    vrb(paste0("NB: ", obj.nx, " genes [", 
               ceiling(obj.nx / nrow(oi.copy) * 100),
               "% of total genes] are not in any gene set ",
               "but are retained for computing null and autocovariance\n\n"))
  } else {
    vrb("\n")
  }
  
  return(list(set.info = si.copy, obj.info = oi.copy, set.obj = so.copy, 
              var.info = var.info, repeat.objID = repeat.objID, 
              obj.info.merged = obj.info.merged, set.info.merged = set.info.merged, 
              removed.sets = removed.sets, read.args = read.args, 
              plR.class = "0|0|0"))
}


#===============================================================================#
# Function: plR.permute(plR.input, set.info, obj.info, set.obj, match, n.perm, 
#                       match.wt, match.z, match.cap, match.fun, match.scale, user.seed,
#                       verbose, fut.plan, n.cores)
#
# Generate permuted gene scores using matching or standard permutations
#
# - plR.input      : R object; output of polylinkR::plR.input
#.                  default = NULL; must be provided if obj.info / set.info / set.obj not provided
# - obj.info       : data.frame or data.table with columns [optional]:
#                      objID (integer), [objName (character)], objStat (numeric), chr (numeric), start (numeric), end (numeric), [Cov1 (numeric)], [Cov2 (numeric)], [Cov (numeric) ...]
#.                   default = NULL; must be provided if plR.input not provided
# - set.info.      : data.frame or data.table with columns [optional]:
#                      setID (integer), [setName (character)], [N (integer)], [setSource (character)]
#.                   default = NULL; must be provided if plR.input not provided
# - set.obj.       : data.frame or data.table with columns [optional]:
#                      setID (integer), objID (integer)
#.                   default = NULL; must be provided if plR.input not provided
# - perm.type      : integer; determine the type of permutation to be performed
#                    default = 1; options: 1 = random, 2 = unstandardised matched, 3 standardised matched
# - n.perm         : integer; number of permutations used to generate null distributions
#                    default = 10000L; range [10000L, Inf) -- increases in increments of 100
# - match.wt       : numeric matrix; user provided matching weights
#                    default = NULL
# - match.cap      : numeric; apply maximum resampling probability for any single gene 
#                    default = 0.05; range = (1 / (n.genes - 1), 1]
# - match.scale    : numeric; value of exponent in soft max function used to convert Mahalanobis Distances to weights
#                    default = 2; range = (0, Inf]
# - match.fun      : function; user specified function to convert Mahalanobis Distances to weights
#                    default = NULL
# - tail.fit       : integer; GPD fitting procedure for empirical p values below threshold value (i.e. n.tail)
#                    default = 1; options: 0 = no fitting, 1 = maximum likelihood, 2 = probability weighted moments
# - n.tail         : integer; threshold no. of exceedences (null set scores >= observed set score) to enact GPD tail fitting (if tail.fit != 0)
#                    default = 10; range = [10, 100] in increments of 10
# - user.seed      : integer: value used to set random seed
#                    default = NULL; range: [-.Machine$integer.max, .Machine$integer.max]
# - verbose.       : logical, should progress reporting be enabled?
#                  : default = TRUE
# - n.cores        : integer; number of cores to utilise in parallelisable operations
#                    default = max no. cores - 1; range = [1, max no. cores]
# - fut.plan       : character; parallel processing option from the future package
#                    default = "multisession"; other options = "multicore" or "sequential"
#===============================================================================#

plR.permute <- function(plR.input = NULL, perm.type = 1, n.perm = 10000L, 
                        match.wt = NULL, match.cap = 0.05, match.scale = 2, 
                        match.fun = NULL, tail.fit = 1, n.tail = 10, 
                        user.seed = NULL, verbose = TRUE, n.cores = "default", 
                        fut.plan = "multisession") {
  
  ##==========================================================================##
  ##PART 1: clean data and run checks
  ##==========================================================================##
  
  startT <- Sys.time()
  suppressMessages(library(foreach))
  suppressMessages(library(data.table))
  suppressMessages(library(doFuture))
  
  # run preliminary file inspection
  class.check(plR.input = plR.input, req.class = "0|0|0", env = environment())
  rm(plR.input); gc()
  file.check(set.info = set.info, obj.info = obj.info, set.obj = set.obj, 
             var.info = var.info, env = environment())
  
  # argument check
  arg.check(plr.func = "permute", env = environment())
  
  # checks finished, announce function
  vrb("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")
  vrb("      Running plR.permute: Generating random set scores\n")                
  vrb("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n")
  
  # reporting
  if (!is.null(WMSS)) warning(WMSS, "\n", immediate. = T, call. = F)
  if (!is.null(MSS)) vrb(paste(paste0(MSS, collapse = ""), "\n"))
  if (!is.null(pWMSS)) warning(pWMSS, "\n", immediate. = T, call. = F)
  if (!is.null(pMSS)) vrb(paste(paste0(pMSS, collapse = ""), "\n"))
  
  # set random seed
  set.seed(seed)
  
  # set up reporting
  if (verbose) {
    progressr::handlers(list(
      progressr::handler_progress(
        format = "[:bar] :percent in :elapsed ETA: :eta",
      )
    ))
  } else {
    progressr::handlers(progressr::handler_void())
  }
  
  
  ##==========================================================================##
  ##PART 2: permute genes
  ##==========================================================================##
  
  oS <- oi.copy$objStat
  set.oi <- oi.copy[objID %in% so.copy$objID] # genes occurring in sets
  set.oi[, objID.i := 1:.N] # add new gene ID index for fast set score calculation
  
  #--------------------------------------------------------------------------#
  ## randomly sample scores based on user-specified sampling scheme
  #--------------------------------------------------------------------------#
  
  if (perm.type == 1) { # equal weights assigned to all genes
    vrb("Random sampling chosen: ")
    match.wt <- matrix(1 / (n.genes - 1), nrow = n.set.genes, ncol = n.genes)
    match.wt[cbind(1:n.set.genes, set.genes)] <- 0
  } else { # calculate Mahalanobis distances to weight sampling
    vrb("Matched sampling chosen: ")
    if (is.null(match.wt)) { # no user supplied matching object
      vrb(paste0(length(cov.cols), " covariate columns identified\n",
                 "Calculating ", ifelse(match.z, "standardised ", ""), 
                 "Mahalanobis distances ", 
                 ifelse(n.set.genes == n.genes, 
                        paste0("for all gene pairs [n = ", n.genes, " ^ 2]\n"), 
                        paste0("from each set gene to all genes [n = ",
                               n.set.genes, " x ", n.genes, "]\n"))))
      
      # calculate Mahalanobis distance for each gene based on covariate values
      if (match.z) { # standardise covariate values
        m.dist <- chol.md(X = oi.copy[, lapply(.SD, scale), .SDcols = cov.cols])
      } else { # use raw covariate values
        m.dist <- chol.md(X = oi.copy[, .SD, .SDcols = cov.cols])
      }
      
      # exclude non set genes
      if (n.set.genes < n.genes) {
        m.dist <- m.dist[set.genes, ]
      }
      
      if (is.null(match.fun)) { # calculate weights using softmax function
        match.wt <- get.obj.mw(md = m.dist, set.genes = set.genes, 
                               match.scale = match.scale, 
                               match.cap = match.cap, env = environment())
        f.choice <- "softmax rescaled MD-based"
      } else { # calculate weights using user defined function
        match.wt <- get.obj.mw(md = m.dist, set.genes = set.genes, 
                               match.fun = match.fun, 
                               match.cap = match.cap, env = environment())
        f.choice <- "user-defined rescaling of MD-based"
      }
      
      rm(m.dist); gc() # clean up unused large objects
    } else {
      match.wt <- match.wt[oi.copy[objID %in% set.genes]$objID.orig, ]
      mw.rs <- Rfast::rowsums(match.wt)
      if (any(zapsmall(mw.rs) != 1)) { # convert to probability weights
        match.wt <- match.wt / mw.rs
      }
      f.choice <- "user-defined"
    }
  }
  
  # generate permuted scores based on weighted sampling
  vrb(paste0("Generating ", n.perm, 
             ifelse(perm.type == 1, " random", " matched"),
             " samples for each set gene\n"))
  if (perm.type > 1) {
    vrb(paste0("[applying ", f.choice, " sampling weights", mw.mss, "]\n"))
  }
  
  rrl <- seq(5, 25, 5) * 10 # between 50 to 250 genes passed to each core
  mc <- rrl[which.max(((n.set.genes / n.cores) / 5) %/% rrl == 0)] # try ~ 5 iterations per core
  ch.n <- ceiling(n.set.genes / mc)
  
  if (fut.plan == "sequential") {
    future::plan("sequential", gc = TRUE)
  } else {
    future::plan(fut.plan, workers = n.cores, gc = TRUE)
    options(future.globals.maxSize = Inf)
  }
  
  progressr::with_progress({
    prog <- progressr::progressor(steps = ch.n)
    p.l <- foreach::foreach(m1 = itertools::isplitRows(match.wt, 
                                                       chunkSize = mc),
                            .options.future = list(packages = "collapse", 
                                                   seed = TRUE)) %dofuture% {
      pm.now <- get.wt.samp(mw = m1, ng = n.genes, np = n.perm)
      prog()
      pm.now
    }
  })

  # extract position and score matrices
  p.l <- unlist(p.l)
  perm.mat <- matrix(p.l, nrow = n.perm)
  perm.score.mat <- matrix(oS[p.l], nrow = n.perm)
  
  rm(p.l); gc() # clean up unused large objects
  
  
  ##==========================================================================##
  ##PART 3: calculate set scores
  ##==========================================================================##
  
  #generate set by gene indicator matrix for fast set score computation
  vrb(paste0("\nCalculating observed scores and [", n.perm,
             "] expected scores for each gene set\n"))
  
  PxG.info <- merge(so.copy[, .(objID, setID)], 
                    set.oi[, .(objID, objID.i)], by="objID")

  PxG <- Matrix::sparseMatrix(i = PxG.info$setID, j = PxG.info$objID.i, 
                              x = 1L, dims = c(n.sets, n.set.genes),
                              check = FALSE)
  
  # calculate set scores
  ss.N <- si.copy$N
  ss.exp <- as.matrix(Matrix::tcrossprod(PxG, perm.score.mat)) / ss.N
  ss.obs <- as.vector(PxG %*% oS[set.genes]) / ss.N
  
  vrb("Calculating p values for each gene set\n")
  p.raw <- (Rfast::rowTrue(ss.exp >= c(ss.obs)) + 1) / (n.perm + 1)
  
  min.p <- (n.tail + 1) / (n.perm + 1)
  if (!is.null(fit.meth) & min(p.raw) <= min.p) { # use GPD to estimate small p values
    vrb(paste0("[using GDP to estimate p values <= ", 
               formatC(min.p, format = "e", digits = log10(n.perm)), 
               "; i.e. (n.tail + 1) / (n.perm + 1)]\n"))
    wP <- which(p.raw <= min.p)
    for (i in wP) {
      p.raw[i] <- est.tail.p(x0 = ss.obs[i], X = ss.exp[i, ], gdp.th = 200, 
                             n.tail = n.tail, fit.meth = fit.meth)
    }
  }
  p.raw[p.raw == 0] <- 2e-16
  
  
  ##==========================================================================##
  ##PART 4: create output
  ##==========================================================================##
  
  # add estimated values to set.info
  si.copy[, setScore.raw := ss.obs]
  si.copy[, setScore.raw.p := p.raw]
  
  # recreate original file formats
  file.reset(oi.copy = oi.copy, si.copy = si.copy, so.copy = so.copy)
  
  # set object class
  plr.class <-  strsplit(plR.class, "\\|")[[1]]
  plr.class[1] <- perm.type
  plr.class[2:3] <- 1
  plr.class <- paste(plr.class, collapse = "|")
  
  # record permutation arguments
  perm.param <- list(n.perm = n.perm, perm.type = perm.type, 
                     match.cap = match.cap, match.scale = match.scale, 
                     match.fun = match.fun)

  vrb("\nCompleted permutations\n")
  vrb(paste0(" *** Run time: ", get.time(start.time = startT), " *** \n\n"))
  
  return(list(set.info = si.copy, obj.info = oi.copy, set.obj = so.copy, 
              var.info = var.info, match.wt = match.wt, perm.mat = perm.mat, 
              setScore.exp = ss.exp, perm.param = perm.param, seed = seed, 
              min.set.n = read.args$min.set.n, plR.class = plr.class))
}


#===============================================================================#
# Function: plR.rescale(plR.input, set.info, obj.info, set.obj, match, n.perm, 
#                       match.wt, match.z, match.cap, match.fun, match.scale,
#                       user.seed, verbose, fut.plan, n.cores)
#
# Computes genomic autocorrelation and (optionally) rescales gene scores
#
# - plR.input.     : R object; output of polylinkR::plR.permute
#                    default = NULL; must be provided if obj.info not provided, required for set score rescaling
# - obj.info       : data.frame or data.table with columns [optional]:
#                      objID (integer), [objName (character)], objStat (numeric), chr (numeric), start (numeric), end (numeric), [Cov1 (numeric)], [Cov2 (numeric)], [Cov (numeric) ...]
#                    default = NULL; must be provided if plR.input not provided, only autocovariance estimation possible
# - ac             : numeric matrix; user provided matrix of intergene autocovariances
#                    default = NULL
# - dist           : character; type of genomic distance measure either "bp" for physical (basepairs) or "cm" for genetic (centimorgans)
#                    default = NULL
# - max.link       : numeric; maximum genomic interval (base pairs) used to evaluate intergene autocovariance
#                    default = 5e6 (bp) or 5 (cm); range = [5e5, Inf) (bp) or [0.5, Inf) (cm)
# - bin.size       : numeric; size of genomic intervals (in base pairs or centimorgans) used in variogram estimation
#                    default = 1e4 (bp) or 0.01 (cm), range [1e2, max.link / 50] (bp) or [1e-4, max.link / 50] (cm)
# - block.size     : integer; number of contiguous genes on chromosome used for independent variogram estimation
#                    default = Inf (i.e. whole chromosomes), range [100L, Inf)
# - kappa          : numeric vector; range of kappa values for Matern model used in variogram fitting
#                    default = (1e-3, 100)
# - min.rho        : numeric; type of genomic distance measure either "bp" for physical (basepairs) or "cm" for genetic (centimorgans)
#                    default = NULL
# - inc.nugget     : logical; should the nugget (non-spatial covariance) be used in covariance estimation?
#                    default = FALSE
# - vg.fit         : integer; determines method employed to fit chosen vg.model(s) to estimated autocovariances
#                    default = 7; see ?gstat::fit.variogram() for full list of options and further details 
# - tail.fit       : integer; GPD fitting procedure for empirical p values below threshold value (i.e. n.tail)
#                    default = 1; options: 0 = no fitting, 1 = maximum likelihood, 2 = probability weighted moments
# - n.tail         : integer; threshold no. of exceedences (null set scores >= observed set score) to enact GPD tail fitting (if tail.fit != 0)
#                    default = 10; range = [10, 100] in increments of 10
# - user.seed      : integer: value used to set random seed
#                    default = NULL; range: [-.Machine$integer.max, .Machine$integer.max]
# - verbose.       : logical, should progress reporting be enabled?
#                  : default = TRUE
# - n.cores        : integer; number of cores to utilise in parallelisable operations
#                    default = max no. cores - 1; range = [1, max no. cores]
# - fut.plan       : character; parallel processing option from the future package
#                    default = "multisession"; other options = "multicore" or "sequential"
#===============================================================================#

plR.rescale <- function(plR.input = NULL, rescale = TRUE, dist = "check", 
                        ac = NULL, max.link = "default", bin.size = "default", 
                        block.size = "default", vg.fit = 7, kappa = c(1e-3, 1e2), 
                        min.rho = 1e-3, inc.nugget = FALSE, tail.fit = 1, 
                        n.tail = 10, user.seed = NULL, verbose = TRUE, 
                        n.cores = "default", fut.plan = "multisession") {
  
  ##==========================================================================##
  ##PART 1: clean data and run checks
  ##==========================================================================##
  
  startT <- Sys.time()
  suppressMessages(library(foreach))
  suppressMessages(library(data.table))
  suppressMessages(library(doFuture))
  
  # run preliminary file inspection
  rc <- paste0(1:3,  "|1|1")
  class.check(plR.input = plR.input, env = environment(),  req.class = rc)
  rm(plR.input); gc()
  file.check(set.info = set.info, obj.info = obj.info, set.obj = set.obj, 
             var.info = var.info, env = environment())
  
  # argument check
  arg.check(plr.func = "rescale", env = environment())
  
  #checks passed, announce function
  if (verbose) {
    lc <- "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"
    if (rescale) {
      if (is.null(ac)) {
        vrb(paste0(lc, 
                   "      Running plR.rescale: estimate intergene autocovariance &\n",
                   "                      rescale gene set scores\n",
                   lc, "\n"))
      } else {
        vrb(paste0(lc,
                   "         Running plR.rescale: rescale gene set scores\n",
                   lc, "\n"))
      }
    } else {
      vrb(paste0(lc,
                 "      Running plR.rescale: estimate intergene autocovariance\n",
                 lc, "\n"))              
    }
  }
  
  # reporting
  if (!is.null(WMSS)) warning(WMSS, "\n", immediate. = T, call. = F)
  if (!is.null(MSS)) vrb(paste(paste0(MSS, collapse = ""), "\n"))
  if (!is.null(pWMSS)) warning(pWMSS, "\n", immediate. = T, call. = F)
  if (!is.null(pMSS)) vrb(paste(paste0(pMSS, collapse = ""), "\n"))
  
  # set random seed
  set.seed(seed)
  
  # set up reporting
  if (verbose) {
    progressr::handlers(list(
      progressr::handler_progress(
        format = "[:bar] :percent in :elapsed ETA: :eta",
      )
    ))
  } else {
    progressr::handlers(progressr::handler_void())
  }
  
  
  ##==========================================================================##
  ##PART 2: Calculate genomic autocovariance
  ##==========================================================================##
  
  if (is.null(ac)) {
    #create inputs for autocovariance estimation functions
    s.u <- ifelse(dist == "cM", 1e6, 1)
    s.d <- 1 / s.u
    if (use.var) {
      score.max <- max(var.info$value)
      rs <- ifelse(score.max > 1e2, floor(log10(score.max)), 0) # scaling factor
      score <- var.info[, .(I = 1:.N, chr, 
                            midpos = pos * s.u, # rescale distances to improve inference
                            midpos.dummy = 0, OS = value * 10 ^ -rs)] # rescale gene scores to avoid integer overflow
    } else {
      score.max <- max(oi.copy$objStat)
      rs <- ifelse(score.max > 1e2, floor(log10(score.max)), 0) # scaling factor
      score <- oi.copy[, .(I = objID, chr, 
                           midpos = midpos * s.u, # rescale distances to improve inference
                           midpos.dummy = 0, OS = objStat * 10 ^ -rs)] # rescale gene scores to avoid integer overflow
      
      #jitter genes with same mean position
      data.table::setkey(score, chr, midpos)
      pos.dups <- score[, .N, by = c("chr", "midpos")][N > 1]
      nn <- sum(pos.dups$N)
      mm <- runif(nn, 0.1, 0.5)
      score[.(pos.dups[, .(chr, midpos)]), midpos := midpos + mm]
    }
    
    # determine maximum variogram distance
    if (max.link == "default") {
      max.link <- ifelse(dist == "cM", 2, 2e6)
    }
    ML <- max.link * s.u # rescale for gstat inference
    
    # determine block sizes
    if (block.size == "default") {
      l <- ifelse(use.var, 1e3, 100)
      bs.l <- ifelse(use.var, max(1e3, var.info[, .N] / 200), 
                     max(100, n.genes / 50))
      block.size <- ceiling(bs.l / l) * l
    }
    
    # create blocks
    vario.block <- ac.blocks(score = score, block.size = block.size)
    
    # determine bin size
    if (bin.size == "default") {
      sn <- ifelse(use.var, 1e2, 10)
      l <- seq(1e3, 1e6, 1e3)
      # determine maximum inter- variant / gene distance in which specific # of variants / genes is reached
      rn <- sapply(vario.block, function(x) {
        x <- sort(x$midpos)
        sx <- sort(diff(x))
        sx[sn]
      })
      base.l <- max(rn)
      BS <- l[which.min(base.l %/% l)]
      bin.size <- BS * s.d # rescale for gstat inference
    } else {
      BS <- bin.size * s.u # rescale for gstat inference
    }

    # report settings
    vrb(paste0(ifelse(use.var, "var.info file detected: f", "F"),
               "itting variagrams for"))
                      
    if (is.infinite(block.size)) {
      vrb(" each chromosome\n")
    } else {
      sc.type <- ifelse(use.var, " variant", " gene") 
      vrb(paste0(" successive blocks of ", block.size, sc.type, 
                 "s sliding at ",
                 floor(block.size / 2), sc.type, " intervals\n"))
    }
    
    vrb(paste0("Estimating covariance for ", 
               ifelse(use.var, "gene" , "variant"), " pairs binned every ", 
               bin.size, dist, " up to maximum of ", max.link, dist, " per ", 
               ifelse(is.infinite(block.size), "chromosome\n", "block\n")))
    
    if (max.link / bin.size > 2000) {
      warning("max.link / bin.size > 2000\n",
              "Variogram model fitting is slow for large numbers of bins", 
              immediate. = TRUE, call. = FALSE)
    }
    
    if (fut.plan == "sequential") {
      future::plan("sequential", gc = TRUE)
    } else {
      future::plan(fut.plan, workers = n.cores, gc = TRUE)
      options(future.globals.maxSize = Inf)
    }
    
    # calculate best fitting variogram for each block
    ff.exp <- list(packages = c("gstat", "data.table"), seed = TRUE)
    progressr::with_progress({
      prog <- progressr::progressor(steps = length(vario.block))
      block.fit <- foreach::foreach(ss = vario.block, 
                                    .options.future = ff.exp) %dofuture% {
        # get best fitting variogram models
        bf.bl <- vg.est(ss.in = ss, max.link = ML, bin.size = BS, 
                        vg.fit = vg.fit, kappa = kappa)
        prog()
        bf.bl
      }
    })
    
    vrb("Variogram fitting complete, checking parameter estimates\n")
    # simplify variogram blocks for export
    vb0 <- foreach::foreach(ss = vario.block, kk = iterators::icount(),  
                   .combine = rbind) %do% {
      ss[, .(Block = kk, chr = chr[1], startpos = min(midpos) * s.d,
             endpos = max(midpos) * s.d)]
    }
    
    # ensure contiguous blocks
    # replace endpos with one bp less than first position of next Block
    for (i in 0:1) {
      w <- vb0[, which(Block %% 2 == i)]
      vb0[w, endpos2 := c((startpos - s.d)[-1], 0), by = chr]
    }
    vb0[, endpos := Rfast::rowMaxs(cbind(endpos, endpos2), value = TRUE)]
    vb0[, endpos2 := NULL]
    min.max1 <- vb0[, .(min(Block), max(Block)), by = chr]
    min.max2 <- oi.copy[, .(min(startpos), max(endpos)), by = chr]
    vb0[Block %in% min.max1$V1, startpos := min.max2$V1]
    vb0[Block %in% min.max1$V2, endpos := min.max2$V2]
    
    # generate gene Blocks if variants used for autocovariance estimation
    if (use.var) {
      vb0[, N.poly := sapply(vario.block, nrow)]
      vb0[, Var.poly := sapply(vario.block, function(x) var(x$OS))]
      data.table::setkey(vb0, chr, startpos, endpos)
      vb1 <- data.table::foverlaps(oi.copy[, .(I = objID, chr, 
                                               startpos = midpos, 
                                               endpos = midpos, 
                                               OS = objStat)], 
                                   vb0, type = "within", nomatch = NULL)
      vario.block <- split(vb1[, .(I, chr, midpos = i.startpos * s.u, OS)], 
                           f = factor(vb1$Block))
    }
    
    vb0[, N.genes := sapply(vario.block, nrow)]
    vb0[, First.gene := sapply(vario.block, function(x) min(x$I))]
    vb0[, Last.gene := sapply(vario.block, function(x) max(x$I))]
    vb0[, Var.gene := sapply(vario.block, function(x) var(x$OS))]

    # add best fitting model info
    vario.fit <- lapply(block.fit, "[[", 1)
    vb0[, Nugget := sapply(vario.fit, "[", 1, 2)]
    vb0[, pSill := sapply(vario.fit, "[", 2, 2)]
    vb0[, Range := sapply(vario.fit, "[", 2, 3) * s.d]
    vb0[, Kappa := sapply(vario.fit, "[", 2, 4)]
    vb0[, StdErr := sapply(vario.fit, attr, "SSErr")]
    vb0[, Singular := sapply(vario.fit, attr, "singular")]
    
    vario.param <- list(use.var = use.var, dist = dist, max.link = max.link, 
                        bin.size = bin.size, block.size = block.size, 
                        vg.fit = vg.fit, kappa = kappa, min.rho = min.rho, 
                        inc.nugget = inc.nugget)
    
    # create empirical variogram output
    vario.emp <- lapply(block.fit, "[[", 2)
    
    if (dist == "cM") { # rescale genetic distance
      for (i in 1:length(vario.emp)) {
        vario.emp[[i]]$dist <- vario.emp[[i]]$dist * s.d
      }
    }
    
    # construct autocovariance object
    if (any(vb0$Range < 0) | any(vb0$Singular)) { # check proper variogram estimation
      if (any(vb0$Singular)) {
        nfb <- sum(vb0$Singular)
        message("Error: Matern covariance function fitting did not converge for ", 
                nfb, ifelse(is.infinite(block.size), " chromosome", " block"),
                ifelse(nfb > 1, "s", ""))
      }
      if (vb0[, any(Range < 0 & !Singular)]) {
        nfb <- vb0[Range < 0 & !Singular, .N]
        message(ifelse(any(vb0$Singular), "Also n", "Error: N"),
                "egative range values detected for ",
                nfb, ifelse(is.infinite(block.size), " chromosome", " block"),
                ifelse(nfb > 1, "s", ""))
      }
      message("Check variogram fit and paramaters [variogram fitting ",
              "can be visualised using polylinkR::check.ac.rs]")
      rescale <- FALSE
      AC <- TRUE
      ac <- NULL
      ac.ok <- FALSE
    } else { # variogram estimation OK
      vrb("Calculating autocovariance for all gene pairs\n")
      if (!inc.nugget) {
        vrb(paste0("[inc.nugget = FALSE: ",
                   "Only spatial covariance component will be calculated]\n"))
      }
      
      # make autocovariance object
      ac.tmp <- foreach::foreach(i = vb0[N.genes > 0]$Block, # only use blocks with genes
                                 .combine = rbind) %do% {
        bl.ac <- make.ac(p2i = vario.block[[i]], vgm.est = vario.fit[[i]], 
                         nug = inc.nugget, MD = max.link)
        cbind(Block = i, bl.ac)
      }
      
      # compute correlations & sample weights
      ac.tmp[, `:=` (rho = covar / max(covar), N = .N), by = Block] 
      data.table::setkey(ac.tmp, objID.A, objID.B)
      
      # compute weighted means of duplicated gene pairs
      ac.N <- ac.tmp[, .N, by = c("objID.A", "objID.B")]
      ac.tmp0 <- ac.tmp[.(ac.N[N > 1])][, SW := N / sum(N), 
                                        by = c("objID.A", "objID.B")]
      ac.tmp0 <- ac.tmp0[, .(covar = sum(covar * SW), rho = sum(rho * SW)), 
                         by = c("objID.A", "objID.B", "midpos.dist")]
      ac <- rbind(ac.tmp[.(ac.N[N == 1]), 
                         .(objID.A, objID.B, midpos.dist, covar, rho)],
                  ac.tmp0)
  
      # remove gene pairs with correlations < min.rho
      ac <- ac[order(objID.A, objID.B)][rho > min.rho]
      if (dist == "cM") {
        ac[, midpos.dist := midpos.dist * s.d] # rescale cM
      }
      AC <- FALSE
      ac.ok <- TRUE
      rm(vario.block, block.fit, ac.tmp, ac.tmp0); gc() # clean up large objects
    }
  } else {
    AC <- TRUE
    ac.ok <- TRUE
    vb0 <- NULL
    vario.emp <- NULL
    vario.param <- NULL
    vario.fit <- NULL
  }
  
  
  ##==========================================================================##
  ##PART 3: rescale set scores
  ##==========================================================================##
  
  if (rescale) {
    if (verbose) {
      vrb(paste0("Rescaling gene set scores using ", 
                 ifelse(AC, "user-provided", "estimated"),
                 " autocovariance matrix:\n"))
    }
    # create input data.table and simplified covariance matrix
    set.oi <- oi.copy[objID %in% so.copy$objID] # genes occurring in sets
    set.oi[, objID.i := 1:.N] # add new gene ID index for fast set score calculation
    soi <- merge(so.copy[, .(setID, objID)], 
                 set.oi[, .(objID, objID.i)], 
                 by="objID")
    soi <- soi[order(setID, objID), .(SID=setID, OID=objID.i)]
    ac.i <- ac[, .(A = objID.A, B = objID.B, covar)]
    data.table::setkey(ac.i, A, B)
    
    # compute rescaling factor for observed set scores
    ss.obs.rs <- ss.rescale(pm = set.genes, AC = ac.i, 
                            OID = soi$OID, SID = soi$SID)
    
    # compute rescaling factor for expected set scores
    n.perm <- nrow(perm.mat)
    vrb(paste0("Computing rescaling coefficients for ", n.perm, " sets of ",
               n.sets, " permuted gene set scores\n"))

    if (fut.plan == "sequential") {
      future::plan("sequential", gc = TRUE)
    } else {
      future::plan(fut.plan, workers = n.cores, gc = TRUE)
      options(future.globals.maxSize = Inf)
    }
    ac.chunk <- 100 # chunk size
    
    progressr::with_progress({
      prog <- progressr::progressor(steps = n.perm / ac.chunk)
      x.r <- foreach::foreach(pc = itertools::isplitRows(perm.mat, 
                                                         chunkSize = ac.chunk),
                              .options.future = list(packages = "data.table")) %dofuture% {
        out <- foreach::foreach(pm = itertools::isplitRows(pc, chunkSize = 1)) %do% {
          ss.rescale(pm = c(pm), AC = ac.i, OID = soi$OID, SID = soi$SID)
        }
        prog()
        do.call(cbind, out)
      }
    })
    
    ss.exp.rs <- do.call(cbind, x.r) # round values, will set to integer
    rm(x.r, perm.mat); gc() # clean up
    
    ss.exp <- setScore.exp * ss.exp.rs # rescale expected scores
    ss.obs <- si.copy$setScore.raw * ss.obs.rs # rescale observed scores
    
    vrb("\nCalculating p values for each gene set\n")
    p.rs <- (Rfast::rowTrue(ss.exp >= ss.obs) + 1) / (n.perm + 1)
    
    min.p <- (n.tail + 1) / (n.perm + 1)
    if (!is.null(fit.meth) & min(p.rs) <= min.p) { # use GPD to estimate small p values
      vrb(paste0("[using GDP to estimate p values <= ", 
                 formatC(min.p, format = "e", digits = log10(n.perm)), 
                 "; i.e. (n.tail + 1) / (n.perm + 1)]\n"))
      wP <- which(p.rs <= min.p)
      for (i in wP) {
        p.rs[i] <- est.tail.p(x0 = ss.obs[i], X = ss.exp[i, ], gdp.th = 200, 
                              n.tail = n.tail, fit.meth = fit.meth)
      }
    }
    p.rs[p.rs == 0] <- 2e-16 # minimum p value
  }

  
  ##==========================================================================##
  ##PART 4: create output
  ##==========================================================================##
  
  ex.mss <- paste0("\nCompleted ", 
                   ifelse(AC, "", "estimating intergene autocovariance"),
                   ifelse(!AC & rescale, " and ", ""),
                   ifelse(rescale, "rescaling set scores", ""), "\n")
  
  if (rescale) { 
    # add estimated values to set.info
    si.copy[, setScore.rs := ss.obs]
    si.copy[, setScore.rs.p := p.rs]
    
    #update plR class
    plr.class <- strsplit(plR.class, "\\|")[[1]]
    plr.class[2] <- 2
    plR.class <- paste0(plr.class, collapse = "|")
  } else { # just return autocovariance objects
    ss.exp <- NULL
    ss.exp.rs <- NULL
    match.wt <- NULL
  }
  
  # reconstitute original parameters
  file.reset(oi.copy = oi.copy, si.copy = si.copy, so.copy = so.copy)
  
  if (ac.ok) {
    vrb(ex.mss)
    vrb(paste0(" *** Run time: ", get.time(start.time = startT), " *** \n\n"))
  }
  
  return(list(obj.info = oi.copy, set.info = si.copy, set.obj = so.copy,
              setScore.exp = ss.exp, setScore.rs = ss.exp.rs, ac = ac, 
              vario.emp = vario.emp, vario.param = vario.param, 
              vario.fit = vario.fit, vario.summary = vb0, match.wt = match.wt, 
              seed = seed, min.set.n = min.set.n, plR.class = plR.class))
}


#===============================================================================#
# Function: plR.prune(plR.input, n.fdr, fast.fdr, est.pi0, pr.method, wCI.alpha,
#.                    min.set.n, user.seed, verbose, fut.plan, n.cores)
#
# Uses regression to rescale gene score due to impact of shared genes and linkage
#
# - plR.input.     : R object; either output of polylinkR::plR.permute or polylinkR::plR.rescale
#                    required object
# - n.fdr          : integer; no. of times to replicate pruning procedure to estimate FDR-corrected p (q) values using histogram method
#                    default = 300; range = [100, Inf]
# - fast.fdr       : logical; should fast pruning be enacted by using only 10000 permutations for null distribution?
#                    default = TRUE
# - est.pi0        : logical; should pi0 be estimated in FDR correction?
#                    default = TRUE; if FALSE pi0 is set to 1
# - pr.method.     : integer; specifies pruning methodology
#                    default = 1; options 1 = regression pruning, 2 = gene removal pruning (not available for rescaled set scores)
# - min.set.n     : integer; minimum gene set size; required only if files not loaded with polylinkR::plR.input
#                    default: NULL
# - tail.fit       : integer; GPD fitting procedure for empirical p values below threshold value (i.e. n.tail)
#                    default = 1; options: 0 = no fitting, 1 = maximum likelihood, 2 = probability weighted moments
# - n.tail         : integer; threshold no. of exceedences (null set scores >= observed set score) to enact GPD tail fitting (if tail.fit != 0)
#                    default = 10; range = [10, 100] in increments of 10
# - user.seed      : integer: value used to set random seed
#                    default = NULL; range: [-.Machine$integer.max, .Machine$integer.max]
# - verbose        : logical, should progress reporting be enabled?
#                  : default = TRUE
# - n.cores        : integer; number of cores to utilise in parallelisable operations
#                    default = max no. cores - 1; range = [1, max no. cores]
# - fut.plan       : character; parallel processing option from the future package
#                    default = "multisession"; other options = "multicore" or "sequential"
#===============================================================================#

plR.prune <- function(plR.input, n.fdr = 300, fast.fdr = TRUE, pr.method = 1, 
                      est.pi0 = TRUE, tail.fit = 1, n.tail = 10, 
                      user.seed = NULL, verbose = TRUE, n.cores = "default", 
                      fut.plan = "multisession") {
  
  
  ##==========================================================================##
  ##PART 1: clean data and run checks
  ##==========================================================================##
  
  startT <- Sys.time()
  suppressMessages(library(foreach))
  suppressMessages(library(data.table))
  suppressMessages(library(doFuture))
  
  # run preliminary file inspection
  rc <- paste0(rep(1:3,  each = 2), "|", rep(1:2, 3), "|", 1)
  class.check(plR.input = plR.input, env = environment(), req.class = rc)
  rm(plR.input); gc()
  file.check(set.info = set.info, obj.info = obj.info, set.obj = set.obj, 
             var.info = var.info, env = environment())
  
  # argument check
  arg.check(plr.func = "prune", env = environment())
  
  # checks completed, announce function
  vrb("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")
  vrb("          plR.prune: correct gene set scores for shared genes\n")
  vrb("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n")
  
  # reporting
  if (!is.null(WMSS)) warning(WMSS, "\n", immediate. = T, call. = F)
  if (!is.null(MSS)) vrb(paste(paste0(MSS, collapse = ""), "\n"))
  if (!is.null(pWMSS)) warning(pWMSS, "\n", immediate. = T, call. = F)
  if (!is.null(pMSS)) vrb(paste(paste0(pMSS, collapse = ""), "\n"))
  
  # set random seed
  set.seed(seed)
  
  # set up reporting
  if (verbose) {
    progressr::handlers(list(
      progressr::handler_progress(
        format = "[:bar] :percent in :elapsed ETA: :eta",
      )
    ))
  } else {
    progressr::handlers(progressr::handler_void())
  }
  
  
  ##==========================================================================##
  ##PART 2: generate FDR sets and empirical prepare objects for pruning step
  ##==========================================================================##
  
  #generate random scores for q-value adjustment
  vrb("Preparing objects for pruning step:\n")
  
  oS <- oi.copy$objStat
  set.oi <- oi.copy[objID %in% so.copy$objID] #genes occurring in sets
  set.oi[, objID.i := 1:.N] # add new gene ID index for fast set score calculation
  soi <- merge(so.copy[, .(setID, objID)], set.oi[, .(objID, objID.i)], 
               by="objID")
  soi <- soi[order(setID, objID), .(SID = setID, OID = objID.i)]
  PxG.info <- merge(so.copy[, .(objID, setID)], 
                    set.oi[, .(objID, objID.i)], by="objID")
  PxG <- Matrix::sparseMatrix(i = PxG.info$setID, j = PxG.info$objID.i, 
                              x = 1L, dims = c(n.sets, n.set.genes),
                              check = FALSE)
  PxN <- as.matrix(Matrix::tcrossprod(PxG)) # shared gene matrix

  # replicate sampling scheme used to generate null
  plr.class <- strsplit(plR.class, "\\|")[[1]]
  vrb(paste0("Generating ", n.fdr, 
             ifelse(plr.class[1] == 3, " standardised ", " "),
             ifelse(plr.class[1] %in% 2:3, "matched", "random"),
             " gene permutations\n"))
  
  rrl <- seq(5, 25, 5) * 10 # between 50 to 250 genes passed to each core
  mc <- rrl[which.max(((n.set.genes / n.cores) / 5) %/% rrl == 0)] # try ~ 5 iterations per core
  ch.n <- ceiling(n.set.genes / mc)
  
  if (fut.plan == "sequential") {
    future::plan("sequential", gc = TRUE)
  } else {
    future::plan(fut.plan, workers = n.cores, gc = TRUE)
    options(future.globals.maxSize = Inf)
  }
  
  progressr::with_progress({
    prog <- progressr::progressor(steps = ch.n)
    p.l <- foreach::foreach(m1 = itertools::isplitRows(match.wt, 
                                                       chunkSize = mc),
                            .options.future = list(packages = "collapse",
                                                   seed = TRUE)) %dofuture% {
      pm.now <- get.wt.samp(mw = m1, ng = n.genes, np = n.fdr)
      prog()
      pm.now
    }
  })
  
  #extract position and score matrices
  p.l <- unlist(p.l)
  fdr.mat  <- matrix(p.l, nrow = n.fdr)
  fdr.score.mat <- matrix(oS[p.l], nrow = n.fdr)
  rm(match.wt, p.l); gc()
  
  vrb("Calculating raw gene set scores and p values\n")
  
  #create FDR matrix
  FDR <- as.matrix(Matrix::tcrossprod(PxG, fdr.score.mat)) / si.copy$N
  
  rescaled <- ifelse(plr.class[2] == 2, TRUE, FALSE)
  if (rescaled) { # regenerate unscaled set scores
    a.p <- Rfast::rowRanks(FDR, descending = TRUE, method = "max")
    b.p <- Rfast::rowRanks(cbind(FDR, setScore.exp / setScore.rs), 
                           descending = TRUE, method = "max")[, 1:n.fdr]
  } else {
    a.p <- Rfast::rowRanks(FDR, descending = TRUE, method = "max")
    b.p <- Rfast::rowRanks(cbind(FDR, setScore.exp), descending = TRUE,
                           method = "max")[, 1:n.fdr]
  }
  fdr.raw.p <- (b.p - a.p + 1) / (n.perm + 1)
  
  min.p <- (n.tail + 1) / (n.perm + 1) # threshold empirical p value for GPD estimation
  if (!is.null(fit.meth) & min(fdr.raw.p) <= min.p) { # use GPD to estimate small p values
    vrb(paste0("Using GDP to estimate raw p values <= ", 
               formatC(min.p, format = "e", digits = floor(log10(n.perm))), 
               " [i.e. (n.tail + 1) / (n.perm + 1)] \n"))
    
    wP <- as.data.frame(which(fdr.raw.p <= min.p, arr.ind = TRUE))
    wP <- data.table::setDT(wP); data.table::setkey(wP, row)
    for (p.row in unique(wP$row)) {
      wp <- as.matrix(wP[.(p.row)])
      if (rescaled) { # regenerate unscaled set scores
        fdr.raw.p[wp] <- est.tail.p(x0 = FDR[wp], X = setScore.exp[p.row, ], 
                                    gdp.th = 200, n.tail = n.tail, 
                                    fit.meth = fit.meth)
      } else {
        fdr.raw.p[wp] <- est.tail.p(x0 = FDR[wp], 
                                    X = setScore.exp[p.row, ] / setScore.rs[p.row, ], 
                                    gdp.th = 200, n.tail = n.tail, 
                                    fit.meth = fit.meth)
      }
    }
    fdr.raw.p[which(fdr.raw.p == 0, arr.ind = TRUE)] <- 2e-16
  }
  
  #summarise results for FDR gene sets
  FDR.pr <- list(FDR.rep = rep(1:n.fdr, each = n.sets), 
                 setID = rep(1:n.sets, n.fdr),
                 setScore.raw = c(FDR), 
                 setScore.raw.p = c(fdr.raw.p))
  data.table::setDT(FDR.pr)
  obj2rm <- c("fdr.raw.p", "a.p", "b.p")
  
  # OPTIONAL: rescale FDR set scores if input set scores have been rescaled
  if (rescaled) {
    vrb("Calculating rescaled FDR gene set scores\n")
    
    # recreate autocovariance matrix
    AC.i <- ac[covar > 0, .(A = objID.A, B = objID.B, covar)]
    data.table::setkey(AC.i, A, B)
    
    #create matrix with set scores for FDR calculations 
    if (fut.plan == "sequential") {
      future::plan("sequential", gc = TRUE)
    } else {
      future::plan(fut.plan, workers = n.cores, gc = TRUE)
      options(future.globals.maxSize = Inf)
    }
    
    progressr::with_progress({
      prog <- progressr::progressor(steps = n.fdr)
      f.r <- foreach::foreach(.options.future = list(packages = "data.table"),
                              fm = itertools::isplitRows(fdr.mat, 
                                                         chunkSize = 1)) %dofuture% {
        out <- ss.rescale(pm = c(fm), AC = AC.i, OID = soi$OID, SID = soi$SID)
        prog()
        out
      }
    })
    
    rs <- do.call(cbind, f.r)
    FDR <- FDR * rs # rescale FDR gene set scores
    
    vrb(paste0("Calculating rescaled p values\n"))
    a.p <- Rfast::rowRanks(FDR, descending = TRUE, method = "max")
    b.p <- Rfast::rowRanks(cbind(FDR, setScore.exp), descending = TRUE,
                           method = "max")[, 1:n.fdr]
    fdr.rs.p <- (b.p - a.p + 1) / (n.perm + 1)
    
    if (!is.null(fit.meth) & min(fdr.rs.p) <= min.p) { # use GPD to estimate small p values
      vrb(paste0("Using GDP to estimate rescaled p values <= ", 
                 formatC(min.p, format = "e", digits = floor(log10(n.perm))), 
                 " [i.e. (n.tail + 1) / (n.perm + 1)] \n"))
      
      wP <- as.data.frame(which(fdr.rs.p <= min.p, arr.ind = TRUE))
      wP <- data.table::setDT(wP); data.table::setkey(wP, row)
      for (p.row in unique(wP$row)) {
        wp <- as.matrix(wP[.(p.row)])
        fdr.rs.p[wp] <- est.tail.p(x0 = FDR[wp], X = setScore.exp[p.row, ], 
                                   gdp.th = 200, n.tail = n.tail, 
                                   fit.meth = fit.meth)
      }
      fdr.rs.p[which(fdr.rs.p == 0, arr.ind = TRUE)] <- 2e-16
    }
    
    
    FDR.pr <- cbind(FDR.pr, setScore.rs = c(FDR), setScore.rs.p = c(fdr.rs.p))
    obj2rm <- c(obj2rm, "rs", "fdr.rs.p")
  }
  
  rm(list = obj2rm); gc()
  
  #--------------------------------------------------------------------------#
  ## Prepare objects needed for pruning step
  #--------------------------------------------------------------------------#
  
  # determine number of permutations to use for pruning step
  if (fast.fdr) {
    if (n.perm == 1e4) { # no downsampling required
      fast.fdr <- F
      n.perm.pr <- n.perm
      SSE <- setScore.exp
    } else { # downsampling required
      n.perm.pr <- 1e4
      sub.sample <- sample.int(n.perm, n.perm.pr)
      perm.mat <- perm.mat[sub.sample, ] # downsample perm.mat
      SSE <- setScore.exp[, sub.sample] # downsample permuted gene set scores
    }
  } else { # no downsampling required
    n.perm.pr <- n.perm
    SSE <- setScore.exp
  }
  
  # determine gene sets with non zero covariance and shared genes
  PxN.gt0 <- which(PxN > 0, arr.ind = TRUE)
  
  #create table to track removed genes and gene sets
  sso <- paste0("setScore.", ifelse(rescaled, "rs", "raw"))
  soc.index <- so.copy[, .(A = setID, X = as.numeric(as.factor(objID)))]
  data.table::setkey(soc.index, A, X)
  
  if (pr.method == 1) { # generate covariance matrix
    vrb("Calculating empirical gene set covariance\n")
    
    # centre set scores
    ss.hat <- Rfast::rowmeans(setScore.exp)
    setScore.exp <- setScore.exp - ss.hat
    SSE <- SSE - ss.hat
    SSO <- si.copy[, get(sso)] - ss.hat
    FDR <- FDR - ss.hat # centre FDR gene set scores
    
    # only calculate non-zero covariance elements
    sse.keep <- data.table::as.data.table(PxN.gt0)[row >= col]
    data.table::setkey(sse.keep, col)
    
    if (fut.plan == "sequential") {
      future::plan("sequential", gc = TRUE)
    } else {
      future::plan(fut.plan, workers = n.cores, gc = TRUE)
      options(future.globals.maxSize = Inf)
    }
    
    progressr::with_progress({
      prog <- progressr::progressor(steps = n.sets)
      c.mat <- foreach::foreach(.options.future = list(packages = "data.table"),
                                i = 1:n.sets) %dofuture% {
        cm.i <- setScore.exp[sse.keep[.(i)]$row, ] %*% setScore.exp[i, ]
        prog()
        cm.i
      }
    })
    
    #create data table comprising gene sets with non-zero covariance
    C.mat <- cbind(sse.keep[, .(A = col, B = row)], 
                   C = unlist(c.mat) / (n.perm - 1))
    C.mat <- rbind(C.mat, C.mat[A != B][, .(A = B, B = A, C)]) # upper triangle
    data.table::setkey(C.mat, A, B)
  
    rm(c.mat, setScore.exp); gc()
  } else {
    vrb("Regenerating permuted gene scores\n")
    
    #update sums rather than means
    SSE <- SSE * si.copy$N
    SSO <- si.copy[, setScore.raw * N]
    FDR <- FDR * si.copy$N
    
    #generate gene statistic matrices
    oS <- oi.copy$objStat
    perm.score.mat <- matrix(oS[perm.mat], nrow = n.perm.pr)
    FDR.score.mat <- matrix(oS[fdr.mat], nrow = n.fdr)
    
    # create data table comprising gene sets with non-zero shared genes
    PI <- list(A = PxN.gt0[, 2], B = PxN.gt0[, 1], 
               N = c(PxN[as.matrix(PxN.gt0)]))
    data.table::setDT(PI)
    data.table::setkey(PI, A, B)
    C.mat <- NULL
    
    rm(perm.mat, setScore.exp); gc()
  }
  
  
  ##==========================================================================##
  ##PART 3: Prune scores and use permuted datasets to estimate FDR
  ##==========================================================================##
  
  vrb("\nPruning step:")
  
  if (pr.method == 1) {
    vrb(" using regression pruning method")
  } else {
    vrb(" using gene removal pruning method")
  }
  
  if (fast.fdr & n.perm > 1e4) { # speed up fdr calculations by using subset of permutations
    vrb(paste0("\n[Fast pruning chosen, only first 10000 permuted gene set scores", 
               " will be used for pruning]"))
    
  }
  
  vrb("\nPruning observed gene set scores\n")
  
  p.init <- si.copy[, get(paste0(sso, ".p"))] * (n.perm.pr + 1)
  
  if (pr.method == 1) {
    SSO.pr <- prune.sets.reg(SSE = SSE, SSO = SSO, p0 = p.init,
                             SOC = data.table::copy(soc.index),
                             mss = min.set.n, CM = data.table::copy(C.mat),
                             n.tail = n.tail, fit.meth = fit.meth)
  } else {
    SSO.pr <- prune.sets.rem(SSE = SSE, SSO = SSO, p0 = p.init,
                             SOC = data.table::copy(soc.index),
                             mss = min.set.n, PxG = PxG, 
                             PI = data.table::copy(PI), 
                             PSM = rbind(oS[set.genes], perm.score.mat),
                             n.tail = n.tail, fit.meth = fit.meth)
  }
  
  c.names <- c("Rank", "N.pr", "Tot.obj.rem",
               "Tot.set.rem", "setScore.pr", "setScore.pr.p")
  colnames(SSO.pr) <- c.names
  SSO.pr <- data.table::as.data.table(SSO.pr)
  SSO.pr <- cbind(setID = 1:n.sets, SSO.pr)
  SSO.pr[Rank == 0, (c.names) := list(NA, NA, NA, NA, NA, NA)]
  if (pr.method == 1) { # recentre set scores
    SSO.pr[, setScore.pr := setScore.pr + ss.hat] 
  }
  SSO.pr[, setScore.pr.p := setScore.pr.p / (n.perm.pr + 1)]
  SSO.pr[setScore.pr.p == 0, setScore.pr.p := 2e-16]
  
  # FDR set scores
  vrb(paste0("Pruning FDR gene set scores [", n.fdr, " replicates]\n"))
  
  if (fut.plan == "sequential") {
    future::plan("sequential", gc = TRUE)
  } else {
    future::plan(fut.plan, workers = n.cores, gc = TRUE)
    options(future.globals.maxSize = Inf)
  }
  
  fdr.p.init <- split(FDR.pr[, get(paste0(sso, ".p"))] * (n.perm.pr + 1), 
                      FDR.pr$FDR.rep)
  fp <- c("data.table", "Rfast", "Matrix")
  if (pr.method == 1) {
    progressr::with_progress({
      prog <- progressr::progressor(steps = n.fdr)
      ps <- foreach::foreach(p.init = fdr.p.init, 
                             fdr.i = itertools::isplitCols(FDR, chunkSize = 1),
                             .options.future = list(packages = fp,
                                                    seed = TRUE)) %dofuture% {
        ps.i <- prune.sets.reg(SSE = SSE, SSO = as.vector(fdr.i), p0 = p.init,
                               SOC = data.table::copy(soc.index),
                               mss = min.set.n, CM = data.table::copy(C.mat),
                               n.tail = n.tail, fit.meth = fit.meth)
        prog()
        ps.i
      }
    })
  } else {
    perm.score.mat <- rbind(FDR.score.mat, perm.score.mat) # combine fdr and permutated gene sets
    j.s <- 1:n.perm.pr + n.fdr
    progressr::with_progress({
      prog <- progressr::progressor(steps = n.fdr)
      ps <- foreach::foreach(p.init = fdr.p.init, j = iterators::icount(),
                             fdr.i = itertools::isplitCols(FDR, chunkSize = 1),
                             .options.future = list(packages = fp,
                                                    seed = TRUE)) %dofuture% {
        ps.i <- prune.sets.rem(SSE = SSE, SSO = as.vector(fdr.i), p0 = p.init, 
                               SOC = data.table::copy(soc.index),
                               mss = min.set.n, PxG = PxG, 
                               PI = data.table::copy(PI), 
                               PSM = perm.score.mat[c(j, j.s), ],
                               n.tail = n.tail, fit.meth = fit.meth)
        prog()
        ps.i
      }
    })
  }
  
  coln <- c(colnames(FDR.pr), c.names)
  FDR.pr <- cbind(FDR.pr, data.table::as.data.table(do.call(rbind, ps)))
  data.table::setnames(FDR.pr, coln)
  FDR.pr[Rank == 0, (c.names) := list(NA, NA, NA, NA, NA, NA)]
  if (pr.method == 1) { # recentre set scores
    FDR.pr[, setScore.pr := (setScore.pr + rep(ss.hat, n.fdr))]
  }
  FDR.pr[, setScore.pr.p := setScore.pr.p / (n.perm.pr + 1)]
  FDR.pr[setScore.pr.p == 0, setScore.pr.p := 2e-16]
  rm(ps); gc()
  
  
  ##==========================================================================##
  ##PART 4: Compute adjusted p-values and generate output
  ##==========================================================================##
  
  vrb("\nComputing FDR corrected p values [q values] and preparing output\n")
  
  #--------------------------------------------------------------------------#
  ## estimate pi0 using histogram method
  #--------------------------------------------------------------------------#
  
  pi0 <- 1 # set initial value
  n.fdr.sets <- nrow(SSO.pr[!is.na(Rank)]) # number of pruned p-values evaluated
  if (est.pi0) { # estimate pi0
    n.cuts <- min(floor(n.fdr.sets), 100) + 1 # number of bins up to max 100
    pi0.p.bins <- seq(0, 1, length.out = n.cuts)
    
    # create percentile bins and calculate observed and expected values
    pi0.dt <- fdr.bin.counts(O = SSO.pr[!is.na(Rank)], 
                             E = FDR.pr[!is.na(Rank)], 
                             p.bins = pi0.p.bins, n.fdr = n.fdr)
    
    # determine obs < exp up to first failure and evaluate probabilities in prior bins
    tolerance <- 10 ^ -4 # difference between successive iterations of piO estimation
    obs.excess <- pi0.dt[, which(cumsum(OBS < EXP) == 1) - 1][1]
    if(obs.excess > 0) { 
      track.pi0 <- pi0
      k <- 0; TT <- TRUE
      pi0.dt[, EXP.orig := EXP] 
      while (isTRUE(TT)) { # run to convergence
        k <- k + 1
        true.pos <- pi0.dt[1:obs.excess, sum(OBS - EXP) / n.fdr.sets]
        pi0 <- 1 - true.pos
        pi0.dt[, EXP.scaled := pi0 * EXP.orig]
        obs.excess <- pi0.dt[, which(cumsum(OBS < EXP) == 1) - 1][1]
        track.pi0 <- c(track.pi0, pi0)
        TT <- abs(track.pi0[k] - pi0) > tolerance # test convergence
      }
    }
  }
  
  #--------------------------------------------------------------------------#
  ## calculate q values
  #--------------------------------------------------------------------------#
  
  # bins to evaluate expected number of true nulls
  emp.p.bins <- SSO.pr[, c(0, sort(unique(setScore.pr.p)))] 
  
  q.dt <- fdr.bin.counts(O = SSO.pr[!is.na(Rank)], 
                         E = FDR.pr[!is.na(Rank)], 
                         p.bins = emp.p.bins, n.fdr = n.fdr)
  q.dt[, pi0 := pi0]
  q.dt[, VP := cumsum(EXP)]
  q.dt[, RP := cumsum(OBS)]
  q.dt[, FDR := Rfast::rowMins(cbind(1, pi0 * VP / RP), value = TRUE)]
  q.dt[FDR == 0, FDR := 1 / (n.fdr.sets + 1)]
  
  # compute q values from FDR values
  FDR.vect <- q.dt$FDR
  M <- min(FDR.vect)
  I <- max(which(FDR.vect == M))
  K <- 1
  q.vect <- c(rep(M, I), rep(NA, length(FDR.vect) - I))
  while(I < length(q.vect)) { # reset FDRs for larger p values with smaller FDR
    K <- I + 1
    M <- min(FDR.vect[K:length(FDR.vect)])
    I <- max(which(FDR.vect == M))
    q.vect[K:I] <- M
  }
  
  q.dt[, q := q.vect]
  SSO.pr[, FDR.bin := cut(setScore.pr.p, emp.p.bins)]
  SSO.pr <- merge(SSO.pr, 
                  q.dt[, .(FDR.bin, setScore.pr.q = q)], 
                  by = "FDR.bin")
  
  #--------------------------------------------------------------------------#
  ## compile output
  #--------------------------------------------------------------------------#
  
  SSO.pr[, FDR.bin := NULL]
  si.copy <- merge(si.copy, SSO.pr, by = "setID", all.x = TRUE)
  
  # return original format
  file.reset(oi.copy = oi.copy, si.copy = si.copy, so.copy = so.copy)
  
  #update plR class
  plr.class[3] <- 2
  plr.class <- paste0(plr.class, collapse = "|")
  
  # reconstitute original parameters
  prune.param <- list(n.fdr = n.fdr, fast.fdr = fast.fdr, 
                      pr.method = pr.method, est.pi0 = est.pi0)
  
  vrb("\nCompleted pruning and FDR correction\n")
  vrb(paste0(" *** Run time: ", get.time(start.time = startT), " *** \n\n"))
  
  return(list(set.info = si.copy, obj.info = oi.copy, set.obj = so.copy,
              c.mat = C.mat, q.dt = q.dt, FDR.pr = FDR.pr, 
              prune.param = prune.param, seed = seed, plR.class = plr.class))
}


#===============================================================================#
# Function: plR.plot(plR.input, plot.path, plot.name,  keep.summary, verbose)
#
# Visual inspection of changes in p values due to pruning and/or rescaling
#
# - plR.input      : R object; output of polylinkR::plR.prune.
#                    required object
# - plot.path      : character; pathway to output folder
#                    default = NULL; current working directory will be used
# - plot.name      : character; name(s) to be appended to default output file name
#                    default = ""; default file name will be used
# - cov.labels     : character vector; vector of covariate labels to use in plots [only used for plR.permute output]
#                    default = NULL; standard labels used
# - highlight.sets : integer vector; vector of setIDs to highlight in plots [only used for plR.permute output]
#                    default: NULL; no sets highlighted
# - cov.mat.       : logical; should all pairwise covariate relationships be plotted? [only used for plR.permute output]
#                    default TRUE
# - fast.plot      : logical; should only the first 10000 permutations be used to generate plots? [only iused for plR.permute output]
#                    default TRUE
# - log.dist       : logical; should logarithmic scaling be used for genomic distances? [only used for plR.rescale output]
#                    default = FALSE
# - log.gamma      : logical; should logarithmic scaling be used for autocovariance measure? [only used for plR.rescale output]
#                    default = FALSE
# - max.facets     : numeric; maximum number of facets to print per page
#                    default = 24
# - keep.summary   : logical; should summarised results be returned?
#                    default = FALSE
# - verbose        : logical, should progress reporting be enabled?
#                  : default = TRUE
#===============================================================================#

plR.plot <- function(plR.input, plot.path = getwd(), plot.name = NULL, 
                     highlight.sets = NULL, keep.summary = FALSE, cov.mat = FALSE, 
                     cov.labels = NULL, max.facets = 24, fast.plot = TRUE, 
                     log.dist = TRUE, log.gamma = TRUE, verbose = TRUE) {
  
  #--------------------------------------------------------------------------#
  ## unpack data and perform basic checks
  #--------------------------------------------------------------------------#
  
  suppressMessages(library(foreach))
  suppressMessages(library(data.table))
  suppressMessages(library(ggplot2))
  suppressMessages(library(ggforce))
  
  # run preliminary file inspection
  class.check0(plR.input = plR.input, env = environment())
  rm(plR.input); gc()
  
  # file check
  file.check(set.info = set.info, obj.info = obj.info, set.obj = set.obj, 
             var.info = var.info, env = environment())
  
  # argument check
  arg.check0(plr.func = input, env = environment())
  
  # reporting
  if (!is.null(WMSS)) warning(WMSS, "\n", immediate. = T, call. = F)
  if (!is.null(MSS)) vrb(paste(paste0(MSS, collapse = ""), "\n"))
  
  if ("matched" %in% input) {
    #--------------------------------------------------------------------------#
    ## Evaluate matching
    #--------------------------------------------------------------------------#

    vrb("Calculating matching statistics\n")
    
    cv.vals <- oi.copy[, .SD, .SDcols = w.cov] # covariate matrix
    set.oi <- oi.copy[objID %in% so.copy$objID] #genes occurring in sets
    set.oi[, objID.i := 1:.N] # add new gene ID index for fast set score calculation
    
    # observed covariate values by gene
    obj.Mo <- cv.vals[set.oi$objID]
    
    mc <- n.perm / 100
    ch.n <- n.perm / mc
    grp <- collapse::GRP(rep(1:n.set.genes, each = mc))
    fp <- c("data.table", "collapse")
    prog <- progressr::progressor(steps = ch.n)
    p.l <- foreach::foreach(P = itertools::isplitRows(perm.mat, chunkSize = mc),
                            .options.future = list(packages = fp,
                                                   seed = TRUE)) %dofuture% {
     collapse::fmean(cv.vals[c(P), ], grp)
    }

    obj.Me <- (foreach::foreach(x = p.l, .combine = `+`) %do% as.matrix(x)) / ch.n
    obj.Me <- as.data.frame(obj.Me)
    data.table::setDT(obj.Me)
    data.table::setnames(obj.Me, colnames(obj.Mo))
    rm(p.l); gc()
    
    # observed covariate values by gene set
    soi <- merge(set.oi[, .SD, .SDcols = c("objID", "objID.i", w.cov)], 
                 so.copy, by = "objID")[order(setID, objID.i)]
    PxG <- Matrix::sparseMatrix(i = soi$setID, j = soi$objID.i, x = 1L,
                                dims = c(n.sets, n.set.genes))
    
    # mean covariate values by gene set
    set.Mo <- as.matrix(PxG %*% as.matrix(obj.Mo) / si.copy$N)
    set.Me <- as.matrix(PxG %*% as.matrix(obj.Me) / si.copy$N)
    
    vrb("\nGenerating plots\n")
    
    # prepare output
    gg.obj <- merge(
      data.table::melt(data.table::data.table(objID = set.oi$objID.orig, obj.Me), 
                       id.vars = "objID", variable.name = "Covariate", 
                       value.name = "Matched"),
      data.table::melt(data.table::data.table(objID = set.oi$objID.orig, obj.Mo), 
                       id.vars = "objID", variable.name = "Covariate", 
                       value.name = "Observed"),
      by = c("objID", "Covariate"))
    
    gg.set <- merge(
      data.table::melt(data.table::data.table(setID = si.copy$setID.orig, set.Me), 
                       id.vars = "setID", variable.name = "Covariate", 
                       value.name = "Matched"),
      data.table::melt(data.table::data.table(setID = si.copy$setID.orig, set.Mo), 
                       id.vars = "setID", variable.name = "Covariate", 
                       value.name = "Observed"),
      by = c("setID", "Covariate"))
    
    #create repeller
    gg.set[, LAB := ifelse(setID %in% highlight.sets, setID, "")]
    
    #create labeller
    if (is.null(cov.labels)) {
      f.nm <- w.cov
    } else {
      f.nm <- cov.labels
    }
    names(f.nm) <- w.cov
    
    # set plot dimensions
    n.cov <- length(w.cov)
    plot.dims <- facet.opt(OPT = min(c(n.cov, max.facets)))
    
    # set plot paths and names
    if (is.null(perm.param)) {
      lab <- "user.weights"
    } else {
      list2env(perm.param, envir = environment())
      lab <- "soft.max.md.weights"
      if (!is.null(match.fun)) {
        lab <- "user.scaled.md.weights"
      } else {
        lab <- c(lab, paste0("scale.", match.scale))
        if (match.cap < 1) {
          lab <- c(lab, paste0("cap.", match.cap))
        }
      }
      if (plr.class[1] == 3) {
        lab <- gsub("md", "md.z", lab)
      }
    }
    lab <- paste0(lab, collapse = "_")
    out.name1 <- paste(plot.name, lab, "obj.pdf",  sep = "_")
    out.name2 <- gsub("obj", "set", out.name1)
    
    # create plots
    vrb("Generating expected vs observed covariate scores")
    vrb(ifelse(length(highlight.sets) > 0, " and marking highlighted sets\n", 
               "\n"))
    n.pages <- ceiling(n.cov / max.facets)
    
    for (i in 1:n.pages) {
      gg1 <- ggplot(gg.obj, aes(x = Matched, y = Observed)) +
        geom_point(size = 0.25, col = "gray25", alpha = 0.5) +
        ggforce::facet_wrap_paginate( ~ Covariate, scales = "free", 
                                     nrow = plot.dims$dim1, 
                                     ncol = plot.dims$dim2, page = i,
                                     labeller = labeller(Covariate = f.nm)) +
        geom_abline(slope = 1, intercept = 0, col = "orange", lty = 2) +
        geom_smooth(linewidth = 0.5, method = "lm", col = "orange") +
        labs(x = "Matched [mean]") +
        theme_bw()  + 
        theme(strip.text.x = element_text(margin = margin(0.1, 0, 0.1, 0, "cm")))
      
      gg2 <- ggplot(gg.set, aes(x = Matched, y = Observed)) +
        geom_point(data = gg.set[LAB == ""],
                   size = 0.25, col = "gray25", alpha = 1/3) + 
        ggforce::facet_wrap_paginate( ~ Covariate, scales = "free", 
                                     nrow = plot.dims$dim1, 
                                     ncol = plot.dims$dim2, page = i,
                                     labeller = labeller(Covariate = f.nm)) +
        geom_abline(slope = 1, intercept = 0, col = "orange", lty = 2) +
        geom_smooth(linewidth = 0.5, method = "lm", col = "orange")  +
        labs(x = "Matched [mean]") +
        theme_bw() + 
        theme(strip.text.x = element_text(margin = margin(0.1, 0, 0.1, 0, "cm")))
      
      if (!is.null(highlight.sets)) {
        library(ggrepel)
        gg2 <- gg2 + 
          geom_text_repel(aes(label = LAB), size = 2, col = "red", 
                          box.padding = 0.5, max.overlaps = Inf,
                          segment.size = 0.2) +
          geom_point(data = gg.set[LAB != ""], size = 0.8, color = "red",
                     alpha = 2/3)
      }
      
      if (n.pages > 1) {
        out.name1 <- gsub(".pdf", paste0("_pg", i, ".pdf"), out.name1)
        out.name2 <- gsub(".pdf", paste0("_pg", i, ".pdf"), out.name2)
      }
      
      pH <- plot.dims$dim1 * 6
      pW <- (plot.dims$dim2 * 6) * 1.025
      
      # save plots
      suppressMessages(ggsave(gg1, filename = file.path(plot.path, out.name1),
                              width = pW, height = pH, units = "cm"))
      
      suppressMessages(ggsave(gg2, filename = file.path(plot.path, out.name2),
                              width = pW, height = pH, units = "cm"))
    }
    
    if (!is.null(highlight.sets)) {
      vrb("Generating covariate distributions for highlighted gene sets\n")
      # visualise covariate distributions for highlighted sets
      hs <- soi[setID.orig %in% highlight.sets]
      hs <- hs[, .(setID = setID.orig, objID = objID.i)]
      hs[, objID.i := as.numeric(as.factor(objID))]
      hs[, setID.i := as.numeric(as.factor(setID))]
      
      # calculate set values for each permutation
      nS <- length(highlight.sets)
      nO <- data.table::uniqueN(hs$objID.i)
      pm.now <- perm.mat[, sort(unique(hs$objID))]
      PxG <- Matrix::sparseMatrix(i = hs$setID.i, j = hs$objID.i, 
                                  x = 1L, dims = c(nS, nO))
      sid <- si.copy[setID.orig %in% highlight.sets][order(setID)]
      sN <- sid$N
      ssc <- foreach::foreach(cc = names(cv.vals)) %do% {
        cv.now <- matrix(cv.vals[, get(cc)][c(pm.now)], nrow = n.perm)
        as.matrix(Matrix::tcrossprod(PxG, cv.now) / sN)
      }
      set.cv <- foreach::foreach(x = ssc) %do% {
        c(t(x))
      }
      set.cv <- c(list(rep(sid$setID.orig, each = n.perm)), set.cv)
      data.table::setDT(set.cv)
      data.table::setnames(set.cv, c("setID", w.cov))
      set.cv <- data.table::melt(set.cv, id.vars = "setID", 
                                 variable.name = "Covariate")
      
      set.cv.mean <- set.cv[, .(value = mean(value)), 
                            by = c("Covariate", "setID")]
      
      all.cv.mean <- list(Covariate = w.cov, value = sapply(cv.vals, mean))
      data.table::setDT(all.cv.mean)
      
      # derive empirical p-values
      s.m <- gg.set[setID %in% highlight.sets]
      emp.quant <- foreach::foreach(i = 1:ncol(cv.vals)) %do% {
        EXP <- ssc[[i]]
        OBS <- s.m[Covariate == paste0("Cov", i)]$Observed
        (Rfast::rowTrue(EXP <= OBS) + 1) / (n.perm + 1)
      }
      
      s.m <- cbind(s.m[order(Covariate, setID)], 
                   Quantile = unlist(emp.quant))
      
      # set plot dimensions
      n.hs <- length(highlight.sets)
      pH <- plot.dims$dim1 * (2 + ((n.hs - 1) * 0.5))
      pW <- (plot.dims$dim2 * 6) * 1.025
      
      out.name3 <- gsub("obj", "covariate.distribution", out.name1) # plot name
      
      # create plot
      for (i in 1:n.pages) {
        gg3 <- ggplot(set.cv) +
          geom_violin(aes(x = value, y = factor(setID)), 
                      linewidth = 0.125, scale = "width") +
          geom_point(data = s.m, size = 2, shape = 23, alpha = 0.75,
                     aes(x = Observed, y = factor(setID), fill = Quantile)) + 
          geom_point(data = set.cv.mean, size = 3, shape = 3, alpha = 0.75,
                     aes(x = value, y = factor(setID))) + 
          geom_vline(data = all.cv.mean, aes(xintercept = value),
                     lty = 2, linewidth = 0.5) +
          scale_shape_manual(values = c(3, 4)) + 
          scale_fill_gradient2(low = "blue", mid = "gray95", high = "red",
                               midpoint = 0.5) +
          ggforce::facet_wrap_paginate( ~ Covariate, scales = "free_x", 
                                       nrow = plot.dims$dim1, 
                                       ncol = plot.dims$dim2, page = i,
                                       labeller = labeller(Covariate = f.nm)) +
          labs(x = "Covariate value", y = "setID") +
          theme_bw() +
          theme(strip.text.x = element_text(margin = margin(0.1, 0, 0.1, 0, "cm")))
        
        if (n.pages > 1) {
          out.name3 <- gsub(".pdf", paste0("_pg", i, ".pdf"), out.name3)
        }
        
        suppressMessages(ggsave(gg3, filename = file.path(plot.path, out.name3),
                                width = pW, height = pH, units = "cm"))
      }
      
      # calculate objstat standardised residuals
      vrb("Generating standardised gene score distributions for highlighted gene sets\n")
      match.wt <- match.wt / Rfast::rowsums(match.wt)
      os <- oi.copy$objStat
      obj.mat <-  list(os = oi.copy[set.genes]$objStat,
                       os.hat = c(match.wt %*% os), # weighted mean object stat
                       os.sq.hat = c(match.wt %*% (os ^ 2)), # weighted mean object stat
                       os.d.eff = 1 / (1 - Rfast::rowsums(match.wt ^ 2))) # sampling design effect
      data.table::setDT(obj.mat)
      obj.mat[, os.var:= (os.sq.hat - os.hat ^ 2) * os.d.eff] # object stat variance]
      obj.mat[, resid := (os - os.hat) / sqrt(os.var)]
      
      hs.obj <- cbind(setID = soi[setID.orig %in% highlight.sets]$setID.orig,
                      obj.mat[soi[setID.orig %in% highlight.sets]$objID.i])
      
      hs.obj.m <- hs.obj[, .(m = mean(resid)), by = setID]
      obj.mat.m <- obj.mat[, .(m = mean(resid))]
      
      n.bins <- 30
      f.seq <- obj.mat[, seq(min(resid), max(resid), length.out = n.bins)]
      f.bins <- (f.seq[-1] + f.seq[-length(f.seq)]) / 2
      n.exp <- table(cut(obj.mat$resid, f.seq, include.lowest = T))
      gg.f <- foreach(hs.i = highlight.sets, .combine = rbind) %do% {
        n.obs <- table(cut(hs.obj[setID == hs.i]$resid, f.seq, 
                           include.lowest = T))
        out <- data.table::data.table(setID = hs.i, Bins = f.bins,
                                      n.obs = c(n.obs), n.exp = c(n.exp))
      }
      
      gg.f[, f.obs := n.obs / sum(n.obs), by = setID]
      gg.f[, f.exp := n.exp / sum(n.exp), by = setID]
      
      plot.dims <- facet.opt(OPT = min(c(n.hs, max.facets)))
      pH <- plot.dims$dim1 * 6
      pW <- (plot.dims$dim2 * 6) * 1.025
      out.name4 <- gsub("obj", "objstat.distribution", out.name1)
      n.pages <- ceiling(n.hs / max.facets)
      s.n <- gg.f[, sum(n.obs), by = setID]
      s.nm <- s.n[, paste0("setID = ", setID, " N = ", V1)]
      names(s.nm) <- s.n$setID
      
      for (i in 1:n.pages) {
        gg4 <- ggplot(gg.f) +
          geom_bar(aes(x = Bins, y = f.exp, group = setID), stat = "identity", 
                   linewidth = 0.01, alpha = 0.25, col = "black") +
          geom_bar(aes(x = Bins, y = f.obs, group = setID), stat = "identity", 
                   linewidth = 0.01, alpha = 0.25, col = "red", fill = "red") +
          geom_vline(data = obj.mat.m, aes(xintercept = m), linetype = 2, 
                     linewidth = 0.75) +
          geom_vline(data = hs.obj.m, aes(xintercept = m), linetype = 2, 
                     linewidth = 0.75, col = "red") +
          #geom_text(aes(label)) +
          labs(x = "Gene score residual", y = "Density") +
          scale_x_continuous(expand = expansion(mult = 0.01)) +
          scale_y_continuous(expand = expansion(mult = 0.02)) +
          ggforce::facet_wrap_paginate(setID ~ ., nrow = plot.dims$dim1, 
                                       ncol = plot.dims$dim2, page = i,
                                       labeller = labeller(setID = s.nm)) +
          theme_bw() +
          theme(strip.text.x = element_text(margin = margin(0.1, 0, 0.1, 0, "cm")))
        
        if (n.pages > 1) {
          out.name4 <- gsub(".pdf", paste0("_pg", i, ".pdf"), out.name4)
        }
        
        # save plots
        suppressMessages(ggsave(gg4, filename = file.path(plot.path, out.name4),
                                width = pW, height = pH, units = "cm"))
      }
    }
    
    if (keep.summary) {
      gg.set[, LAB := NULL]
      summ.list <- list(obj.gg = gg.obj, set.gg = gg.set)
    }
  } else {
    if (any(c("variogram", "rescaled") %in% input)) {
      if ("variogram" %in% input){
        #-----------------------------------------------------------------------#
        ## Evaluate variogram fitting
        #-----------------------------------------------------------------------#
        vrb("Generating variogram fit plots\n")
        
        # get variogram parameters
        list2env(vario.param, envir = environment())
        
        # create variogram fit objects
        dist.unit <- ifelse(dist == "cM", 1e-6, 1)
        if (log.dist) { #create logarithmic distances if required
          dist.seq <- 10 ^ seq(log10(dist.unit), log10(max.link), length.out = 1000)
        } else { # use standard distances
          dist.seq <- seq(dist.unit, max.link, length.out = 1000)
        }
        
        vario.fit <- foreach::foreach(vs = itertools::isplitRows(vario.summary, 
                                                                 chunkSize = 1),
                                      i = iterators::icount(), 
                                      .combine = rbind) %do% {
          if (!vs$Singular) { # omit singular model fits
            vs.m <- vs[, matern(kappa = Kappa, psill = pSill, nugget = Nugget, 
                                range = Range, x = dist.seq)]
          } else {
            vs.m <- NA
          }
          data.table::data.table(Block = i, dist = dist.seq, 
                                 gamma = vs.m, StdErr = log10(vs$StdErr))
        }
        
        if (keep.summary) { # add column with fitted gamma values
          vario.emp <- foreach::foreach(ve = vario.emp, i = iterators::icount(),
                                        vs = itertools::isplitRows(vario.summary, 
                                                                   chunkSize = 1),
                                        .combine = rbind) %do% {
            if (!vs$Singular) { # omit singular model fits
              vs.m <- vs[, matern(kappa = Kappa, psill = pSill, nugget = Nugget, 
                                  range = Range, x = ve$dist)]
            } else {
              vs.m <- NA
            }                   
            data.table::data.table(Block = i, ve[, 1:3], gamma.fit = vs.m)
          }
        } else {
          vario.emp <- foreach::foreach(ve = vario.emp, i = iterators::icount(),
                                        .combine = rbind) %do% {
            data.table::data.table(Block = i, ve[, 1:3])
          }
        }
        
        # scale physical distances
        DIST <- ifelse(dist == "bp", "Mbp", "cM")
        vs <- data.table::copy(vario.summary)
        if (dist == "bp") {
          vs[, startpos := startpos * 1e-6]
          vs[, endpos := endpos * 1e-6]
          vario.fit[, dist := dist * 1e-6]
        }
        
        # create plots
        x.tr <- ifelse(log.dist, "log10", "identity") # log transform axes if requested
        y.tr <- ifelse(log.gamma, "log10", "identity") # log transform axes if requested
        x.max <- max.link
        x.min <- 10 ^ floor(log10(vario.emp[, min(dist)]))
        y.max <- vario.emp[, max(gamma)]
        y.min <- 10 ^ floor(log10(vario.emp[, min(gamma)]))
        
        # create block labels
        LABS <- vs[, paste0(chr, ":", round(startpos, 2), "-", round(endpos, 2), 
                            " SE = ", formatC(StdErr, format = "e", digits = 2))]
        names(LABS) <- vs$Block
        
        # set up plot structure
        n.blocks <- vario.emp[, uniqueN(Block)]
        plot.dims <- facet.opt(OPT = min(n.blocks, max.facets))
        n.pages <- ceiling(n.blocks / max.facets)
        
        for (i in 1:n.pages) {
          gg0 <- ggplot(data = vario.emp) +
            geom_point(aes(x = dist, y = gamma, size = np), alpha = 0.5) +
            geom_line(data = vario.fit, linewidth = 0.75,
                      aes(x = dist, y = gamma, col = StdErr)) +
            geom_point(data = vs[Nugget >= y.min],
                       aes(x = x.min, y = Nugget), size = 2, pch = 4) +
            geom_hline(data = vs, linewidth = 0.75,
                       aes(yintercept = Var.gene), col = "blue", lty = 2) +
            scale_x_continuous(expand = expansion(mult = 0.02), trans = x.tr) +
            scale_y_continuous(expand = expansion(mult = 0.01), trans = y.tr) +
            scale_color_gradient2(low = "yellow", mid = "orange", high = "red",
                                  midpoint = mean(vs[, log10(StdErr)])) +
            coord_cartesian(xlim = c(x.min, x.max),
                            ylim = c(y.min, y.max)) +
            ggforce::facet_wrap_paginate(~ Block, ncol = plot.dims$dim2, 
                                         nrow = plot.dims$dim1, page = i,
                                         labeller = as_labeller(LABS)) +
            labs(x = paste0("Intergene distance [", DIST, "]"),
                 y = "Autocovariance", col = "log10(SE)",
                 size = paste("#", ifelse(use.var, "poly.", "gene"), "pairs")) +
            theme_bw() +
            theme(strip.text.x = element_text(margin = margin(0.1, 0, 0.1, 0, "cm")))
          
          if (log.dist) {
            gg0 <- gg0 +
              annotation_logticks(short = unit(0.075, "cm"), 
                                  mid = unit(0.1, "cm"),
                                  long = unit(0.25, "cm"),
                                  sides = "b", size = 0.25)
          }
          
          if (log.gamma) {
            gg0 <- gg0 +
              annotation_logticks(short = unit(0.075, "cm"), 
                                  mid = unit(0.1, "cm"),
                                  long = unit(0.25, "cm"),
                                  sides = "l", size = 0.25)
          }
          
          out.name1 <- paste0(plot.name, "_variogram_fit.pdf")
          if (n.pages > 1) {
            out.name1 <- gsub(".pdf", paste0("_pg", i, ".pdf"), out.name1)
          }
          
          #plot dimensions
          pH <- plot.dims$dim1 * 6
          pW <- (plot.dims$dim2 * 6) * 1.025 
          
          suppressMessages(ggsave(gg0, filename = file.path(plot.path, out.name1),
                                  width = pW, height = pH, units = "cm"))
        }
      } else {
        vario.emp <- NULL
      }
      
      if ("rescaled" %in% input) {
        #--------------------------------------------------------------------------#
        ## Evaluate set score rescaling
        #--------------------------------------------------------------------------#
        vrb("Generating set score rescaling plots\n")
        
        rs.gg <- cbind(si.copy[, .(setID = setID.orig, N, 
                                   rs.obs = setScore.rs / setScore.raw)], 
                       rs.exp = Rfast::rowmeans(setScore.rs))
        ss1 <- seq(0, 1, length = 31)
        rs.gg2 <- list((ss1[-1] + ss1[-length(ss1)]) / 2,
                       table(cut(rs.gg$rs.obs, ss1, include.lowest = T)),
                       table(cut(rs.gg$rs.exp, ss1, include.lowest = T)))
        data.table::setDT(rs.gg2)
        data.table::setnames(rs.gg2, c("Mids", "Obs", "Exp"))
        
        #create repeller
        rs.gg[, LAB := ifelse(setID %in% highlight.sets, setID, "")]
        
        gg1 <- ggplot(rs.gg2) +
          geom_bar(aes(x = Mids, y = Obs, fill = "Obs"), stat = "identity", 
                   linewidth = 0.01, alpha = 0.25, col = "red") +
          geom_bar(aes(x = Mids, y = Exp, fill = "Exp"), stat = "identity", 
                   linewidth = 0.01, alpha = 0.25, col = "black") +
          scale_x_continuous(expand = expansion(mult = 0.01)) +
          scale_y_continuous(expand = expansion(mult = 0.01)) +
          labs(x = "Rescaling coeff.", y = "Frequency") +
          scale_fill_manual(name = "", breaks = c("Exp", "Obs"),
                              values = c("Obs" = "red", "Exp" = "black")) +
          theme_bw() +
          theme(legend.title = element_blank(),
                legend.position = c(0.01, 0.99),
                legend.justification = c("left", "top"),
                legend.direction = "horizontal",
                legend.background = element_rect(fill = "white", 
                                                 linewidth = 0.25, 
                                                 color = "black"),
                legend.key.height = unit(0.75, "line"),
                legend.text = element_text(size = 5)) +
        guides(color = guide_legend(override.aes = list(size = 0.5)))
        
        gg2 <- ggplot(rs.gg, aes(x = rs.exp, y = rs.obs)) +
          geom_point(data = rs.gg[LAB == ""], alpha = 0.25, size = 0.75) +
          geom_abline(slope = 1, intercept = 0, linetype = 2) +
          geom_smooth(method = lm, linewidth = 0.5) + 
          labs(x = "Exp. rescal. coeff.", y = "Obs. rescal. coeff.") +
          scale_x_continuous(expand = expansion(mult = 0.01)) +
          scale_y_continuous(expand = expansion(mult = 0.01)) +
          theme_bw()
        
        gg3 <- ggplot(rs.gg, aes(x = N, y = rs.exp)) +
          geom_point(data = rs.gg[LAB == ""], alpha = 0.25, size = 0.75) +
          geom_smooth(method = lm, linewidth = 0.5) + 
          labs(x = "Gene set size", y = "Exp. rescal. coeff.") +
          scale_x_continuous(expand = expansion(mult = 0.01)) +
          scale_y_continuous(expand = expansion(mult = 0.01)) +
          theme_bw()
        
        gg4 <- ggplot(rs.gg, aes(x = N, y = rs.obs)) +
          geom_point(data = rs.gg[LAB == ""], alpha = 0.25, size = 0.75) +
          geom_smooth(method = lm, linewidth = 0.5) + 
          labs(x = "Gene set size", y = "Obs. rescal. coeff.") +
          scale_x_continuous(expand = expansion(mult = 0.01)) +
          scale_y_continuous(expand = expansion(mult = 0.01)) +
          theme_bw()
        
        if (!is.null(highlight.sets)) {
          library(ggrepel)
          gg2 <- gg2 + 
            ggrepel::geom_text_repel(aes(label = LAB), size = 2, col = "red", 
                                     box.padding = 0.5, max.overlaps = Inf,
                                     segment.size = 0.2) +
            geom_point(data = rs.gg[LAB != ""], size = 1.25, color = "red",
                       alpha = 2/3)
          
          gg3 <- gg3 + 
            ggrepel::geom_text_repel(aes(label = LAB), size = 2, col = "red", 
                                    box.padding = 0.5, max.overlaps = Inf,
                                    segment.size = 0.2) +
            geom_point(data = rs.gg[LAB != ""], size = 1.25, color = "red",
                       alpha = 2/3)
          
          gg4 <- gg4 + 
            ggrepel::geom_text_repel(aes(label = LAB), size = 2, col = "red", 
                                     box.padding = 0.5, max.overlaps = Inf,
                                     segment.size = 0.2) +
            geom_point(data = rs.gg[LAB != ""], size = 1.25, color = "red",
                       alpha = 2/3)
        }
        
        gg.out <- patchwork::wrap_plots(list(gg1, gg2, gg3, gg4), ncol = 2)
        
        out.name2 <- paste0(plot.name, "_set.score_rescaling.pdf")
        suppressMessages(ggsave(gg.out, filename = file.path(plot.path, out.name2),
                                width = 16, height = 16, units = "cm"))
        
      } else {
        rs.gg <- NULL
      }
      
      if (keep.summary) { # return summary data if requested
        return(list(vario.gg = vario.emp, rs.gg = rs.gg))
      }
    } else {
      #--------------------------------------------------------------------------#
      ## Evaluate pruning
      #--------------------------------------------------------------------------#
      vrb("Generating pruning plots\n")
      
      rescaled <- plr.class[2] == 2
      n.breaks <- min(100,  nrow(si.copy) %/% 5)
      nn <- si.copy[, which(!is.na(Rank))]
      
      #distributions
      if (rescaled) {
        seq.l <- 1:3
        K <- c(4, 7, 8, 2, 3, 6) 
        J <- list(1:2, c(1, 3), 2:3)
      } else {
        seq.l <- c(1, 3)
        K <- c(3, 2)
        J <- list(c(1, 3))
      }
      
      library(ggplot2)
      gg.x <- list()
      for(I in seq.l) {
        i <- c("raw", "rs", "pr")[I]
        labs <- paste0("setScore.", i)
        si.now <- si.copy[, .(name = c("Raw", "Rescaled", "Pruned")[I],
                              x = get(labs))][!is.na(x)]
        
        col.now <- ifelse(i == "pr", "orange", "gray50")
        gg0 <- ggplot(si.now, aes(x = x)) +
          geom_histogram(bins = n.breaks, fill = col.now)
        
        if (i != "pr") {
          gg0 <- gg0 +
            geom_histogram(data = si.now[nn], fill = "orange", bins = n.breaks) +
            geom_vline(xintercept = mean(si.now[nn]$x), linetype = 2, 
                       col = "red")
        }
        
        if (i == "raw") {
          gg0 <- gg0 +
            facet_wrap(~ name, strip.position = "top")
        }
        
        if (i == "pr") {
          gg0 <- gg0 +
            facet_wrap(~ name, strip.position = "right")
        }
        
        col.now <- ifelse(i == "pr", "red", "black")
        gg0 <- gg0 + 
          geom_vline(xintercept = mean(si.now$x), linetype = 2, col = col.now) +
          theme_bw() +
          theme(axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                axis.text.x = element_text(size = 7.5),
                axis.text.y = element_text(size = 7.5))
        
        X <- ifelse(rescaled, I * 4 - 3, I * 1.5 - 0.5)
        gg.x[[X]] <- gg0
      }
      
      si.copy[, `q values` := cut(setScore.pr.q, c(0, 0.01, 0.05, 0.20, 0.5, 1))]
      si.copy[is.na(`q values`), `q values` := "Pruned"]
      sh.x <- c(21:23, 3:4, 20)
      sz.x <-  c(1.5, 1.5, 1.25, 1, 1, 1)
      cl.x <- c("pink", "pink", "pink", "black", "black", "black")
      al.x <- c(8, 8, 8, 5, 5, 5) / 10
      names(sh.x) <- names(cl.x) <- names(al.x) <- si.copy[, levels(`q values`)]
      
      #pairwise set score scatterplots
      k <- 0
      for (i in c("", ".p")) {
        for (j in J) {
          k <- k + 1
          stats <- c("raw", "rs", "pr")[j]
          nm <- c("Raw", "Rescaled", "Pruned")[j]
          labs <- paste0("setScore.", stats, i)
          si.now <- si.copy[, .(setID = setID.orig, name = nm[1], name2 = nm[2],
                                x = get(labs[1]), y = get(labs[2]), `q values`)]
          si.now[, LAB := ifelse(setID %in% highlight.sets, setID, "")]
          si.now <- si.now[!is.na(y)] # remove NAs
          
          gg0 <- ggplot(si.now, aes(x = x, y = y)) +
            geom_point(aes(shape = `q values`, fill = `q values`,
                           alpha = `q values`, size = `q values`), stroke = 0.4) +
            scale_shape_manual(values = sh.x, drop = FALSE) +
            scale_size_manual(values = sz.x, drop = FALSE) +
            scale_fill_manual(values = cl.x, drop = FALSE) +
            scale_alpha_manual(values = al.x, drop = FALSE) +
            geom_abline(intercept = 0, slope = 1, col = "red", 
                        linewidth = 0.5, alpha = 0.5) +
            theme_bw() +
            theme(axis.title.x = element_blank(),
                  axis.title.y = element_blank(),
                  axis.text.x = element_text(size = 7.5),
                  axis.text.y = element_text(size = 7.5))
          
          if (i == "") { # add lines at origin
            if (any(si.now$x < 0)) {
              gg0 <- gg0 +
                geom_vline(xintercept = 0, linetype = 2, linewidth = 0.5)
            }
            if (any(si.now$y < 0)) {
              gg0 <- gg0 +
                geom_hline(yintercept = 0, linetype = 2, linewidth = 0.5)
            }
          } else { # log scale for p values
            gg0 <- gg0 +
              scale_x_log10() +
              scale_y_log10() +
              annotation_logticks(short = unit(0.05, "cm"),
                                  mid = unit(0.01, "cm"),
                                  long = unit(0.15, "cm"), size = 0.25) +
              coord_flip()
          }
          
          if (rescaled) {
            if (K[k] == 2) {
              gg0 <- gg0 +
                facet_wrap(~ name2, strip.position = "top")
            }
            if (K[k] == 6) {
              gg0 <- gg0 +
                facet_wrap(~ name, strip.position = "right")
            }
            if (K[k] == 3) {
              gg0 <- gg0 +
                facet_grid(name ~ name2)
            }
          } else {
            if (K[k] == 2) {
              gg0 <- gg0 +
                facet_grid(name ~ name2)
            }
          }
          
          if (!is.null(highlight.sets)) {
            library(ggrepel)
            gg0 <- gg0 + 
              ggrepel::geom_text_repel(aes(label = LAB), size = 2, col = "red", 
                                       box.padding = 0.5, max.overlaps = Inf,
                                       segment.size = 0.2)
          }
          
          gg.x[[K[k]]] <- gg0
        }
      }  
      
      gg1 <- patchwork::wrap_plots(gg.x, ncol = length(seq.l)) + 
        patchwork::plot_layout(guides = 'collect')
      
      if (rescaled) {
        sfx <- paste(plr.class, collapse = "_")
      } else {
        sfx <- paste(plr.class[-2], collapse = "_")
      }
      out.name1 <- paste0("set.score_distributions_", sfx, ".pdf")
      if (plot.name != "") {
        out.name1 <- paste0(plot.name, "_", out.name1)
      }
      
      ph <- length(seq.l) * 5
      pw <- ph * ifelse(rescaled, 4/3, 7/5)
      suppressMessages(ggsave(gg1, filename = file.path(plot.path, out.name1),
                              width = pw, height = ph, units = "cm"))
      
      #p setscore p value distributions
      q.vect <- c(1, 5, 25, 50, 75, 95, 99) / 100
      fdr.n <- FDR.pr[!is.na(Rank), .N, by = Rank]
      fdr.n[, Rank.i := ifelse(N == max(N), Rank, min(Rank[N < max(N)]))]
      FDR.pr <- merge(FDR.pr, fdr.n, by = "Rank")
      fdr.summ <- FDR.pr[!is.na(Rank), .(Quant = q.vect, 
                                         Value = quantile(setScore.pr.p, q.vect)), 
                         by = Rank.i]
      fdr.summ <- data.table::dcast(fdr.summ, Rank.i ~ Quant,
                                    value.var = "Value")
      
      gg2a <- ggplot(fdr.summ, aes(y = Rank.i)) +
        geom_path(aes(x = `0.5`), size = 0.5) + 
        geom_ribbon(aes(xmin = `0.01`, xmax = `0.99`), 
                    fill = "grey50", alpha = 0.1) + 
        geom_ribbon(aes(xmin = `0.05`, xmax = `0.95`), 
                    fill = "grey50", alpha = 0.2) + 
        geom_ribbon(aes(xmin = `0.25`, xmax = `0.75`), 
                    fill = "grey50", alpha = 0.2) + 
        scale_x_continuous(expand = expansion(mult = c(0.005, 0.005))) +
        scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
        labs(x = "p value", y = "Rank") +
        theme_bw()
      
      n.bins <- min(100, si.copy[!is.na(Rank), .N],
                    FDR.pr[!is.na(Rank), .N, by = FDR.rep][, floor(mean(N))])
      bins <- seq(0, 1, length.out = n.bins + 1)
      fdr.x <- FDR.pr[!is.na(Rank), hist(setScore.pr.p, breaks = bins, 
                                         plot = FALSE)[c(3:4)], by = FDR.rep]
      fdr.x <- fdr.x[, .(m = mean(density),
                         q05 = quantile(density, 0.05),
                         q25 = quantile(density, 0.25),
                         q50 = quantile(density, 0.50), 
                         q75 = quantile(density, 0.75),
                         q95 = quantile(density, 0.95)), by = mids]
      si.y <- si.copy[!is.na(Rank), hist(setScore.pr.p, breaks = bins, 
                                         plot = FALSE)[c(3:4)]]
      
      # p value distribution
      gg2b <- ggplot(fdr.x, aes(x = mids)) +
        geom_point(aes(y = q50 / sum(q50)), shape = 3, size = 0.75) + 
        geom_point(aes(y = m / sum(m)), size = 0.75) + 
        geom_line(aes(y = m / sum(m)), linewidth = 0.25) + 
        geom_line(data = si.y,  aes(y = density / sum(density)), 
                  col = "red", linewidth = 0.25) + 
        geom_point(data = si.y, aes(y = density / sum(density)), 
                   col = "red") +
        geom_ribbon(aes(ymin = q05 / sum(q50), ymax = q95 / sum(q50)), 
                    fill = "grey50", alpha = 0.25) + 
        geom_ribbon(aes(ymin = q25 / sum(q50), ymax = q75 / sum(q50)), 
                    fill = "grey50", alpha = 0.25) + 
        scale_x_continuous(expand = expansion(mult = c(0.005, 0.005))) +
        scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
        labs(x = "p value", y = "Density") +
        theme_bw()
      
      # FDR plot
      gg2c <- ggplot(q.dt) +
        geom_point(aes(x = 1:nrow(q.dt), y = FDR)) +
        geom_line(aes(x = 1:nrow(q.dt), y = q), col = "red") +
        labs(x = "Rank", y = "q / FDR") +
        scale_x_continuous(expand = expansion(mult = c(0.005, 0.005))) +
        scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
        theme_bw()
      
      p.scale <- 0.5 * (2 + (nrow(fdr.summ) - 1) / 100)
      gg2 <- patchwork::wrap_plots(list(gg2a, gg2b, gg2c), ncol = 1) + 
        patchwork::plot_layout(heights = c(p.scale, 1, 1))
      
      out.name2 <- gsub("distributions", "q.summary", out.name1)
      suppressMessages(ggsave(gg2, filename = file.path(plot.path, out.name2),
                              width = 20, height = 15, units = "cm"))
      
      if (keep.summary) { # return summary data if requested
        summ.list <- list(FDR.gg = fdr.x)
      }
    }
  }
  
  if (cov.mat) {
    #--------------------------------------------------------------------------#
    # plot pairwise relationships
    #--------------------------------------------------------------------------#
    # create covariate matrix
    vrb("Generating pairwise relationships for covariates and gene scores\n")
    if (input != "matched") {
      n.cov <- length(w.cov)
      #create labeller
      if (is.null(cov.labels)) {
        f.nm <- w.cov
      } else {
        f.nm <- cov.labels
      }
      names(f.nm) <- w.cov
      cv.vals <- cbind(oi.copy$objStat, oi.copy[, .SD, .SDcols = w.cov]) 
    } else {
      cv.vals <- cbind(oi.copy$objStat, cv.vals)
    }
    data.table::setnames(cv.vals, c("Gene score", f.nm))
    
    my_gg <- function(data, mapping){
      ggplot(data = data, mapping = mapping) +
        geom_hex(bins = 40) +
        geom_smooth(color = "orange", linewidth = 0.75) +
        geom_smooth(method = lm, color = "red", linewidth = 0.75)
    }
    
    gg5 <- GGally::ggpairs(cv.vals, progress = FALSE,
                           lower = list(continuous = my_gg)) +
      theme_bw()
    
    out.name <- paste0(plot.name, "_covariate_biplots.pdf")
    pD <- (n.cov + 1) * 2.5
    suppressMessages(ggsave(gg5, filename = file.path(plot.path, out.name),
                            width = pD, height = pD, units = "cm"))
  }
  
  # return data summaries used to make objects
  if (keep.summary) {
    vrb("keep.summary = TRUE; returning plotting objects")
    return(summ.list)
  }
}


