## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

suppressPackageStartupMessages({
  library(epitopes)
  library(seqinr)})

## ---- eval=FALSE--------------------------------------------------------------
#  library(epitopes)
#  
#  # Download the full IEDB export into a given folder "XYZ"
#  epitopes::get_IEDB(save_folder = "XYZ")
#  
#  # Extract only the linear B-cell epitopes from the IEDB export and save the
#  # result to folder "ABC"
#  epitopes <- epitopes::get_LBCE(data_folder = "XYZ",
#                                 ncpus = parallel::detectCores() - 2,
#                                 save_folder = "ABC")
#  
#  # Retrieve the relevant proteins and save them into the same folder "ABC"
#  proteins <- epitopes::get_proteins(uids = unique(epitopes$protein_id),
#                                     save_folder = "ABC")
#  
#  # Retrieve a taxonomy list for all pathogens and save it into the same folder "ABC"
#  taxonomy <- epitopes::get_taxonomy(uids = unique(epitopes$sourceOrg_id),
#                                     save_folder = "ABC")

## ---- eval = FALSE------------------------------------------------------------
#  # Join proteins into epitope dataset
#  jdf <- epitopes::prepare_join_df(epitopes, proteins,
#                                   min_epit        = 8,
#                                   max_epit        = 25,
#                                   only_exact      = FALSE,
#                                   pos.mismatch.rm = "all",
#                                   set.positive    = "mode")

## ---- eval = FALSE------------------------------------------------------------
#  OrgID   <- 6282 # Pathogen: O. volvulus
#  hostIDs <- 9606 # Host: humans
#  
#  # Filter data
#  jdf <- epitopes::filter_epitopes(jdf,
#                                   orgIDs   = OrgID,
#                                   hostIDs  = hostIDs,
#                                   tax_list = taxonomy)
#  

## ---- eval = FALSE------------------------------------------------------------
#  # Prepare windowed representation
#  wdf <- epitopes::make_window_df(jdf, ncpus = parallel::detectCores() - 2)
#  
#  # Calculate features
#  wdf <- epitopes::calc_features(wdf, max.N = 2, ncpus = parallel::detectCores() - 2)
#  

## ---- eval = FALSE------------------------------------------------------------
#  library(seqinr)
#  
#  prots <- proteins[proteins$UID %in% unique(jdf$protein_id), ]
#  
#  if(!dir.exists("./BLASTp")) dir.create("./BLASTp")
#  fn <- "./BLASTp/prots.fasta"
#  seqinr::write.fasta(as.list(prots$TSeq_sequence), names = prots$UID,
#                      file.out = fn, as.string = TRUE, nbchar = 100000)
#  
#  # Run BLASTp to determine protein similarities
#  system(paste0("makeblastdb -in ", fn, " -dbtype prot"))
#  system(paste0("blastp -query ", fn, " -db ", fn, " -seg no ",
#                "-outfmt '6 qseqid sseqid length qlen ",
#                "slen nident pident qcovhsp mismatch gaps qstart ",
#                "qend sstart send evalue score' > ", fn, "-BLAST"))
#  
#  # Split data into training/test/final validation sets, splitting by protein ID
#  # and using BLASTp to minimise similarities across splits
#  splits <- epitopes::split_epitope_data(wdf, split_level = "prot",
#                                         split_perc = c(75, 25),
#                                         split_names = c("01_training", "02_holdout"),
#                                         save_folder = "./splits",
#                                         blast_file = "./BLASTp/prots.fasta-BLAST",
#                                         coverage_threshold = 60, identity_threshold = 60)

## ---- eval=FALSE--------------------------------------------------------------
#  # Fit model using ranger's standard parameters
#  my.model <- epitopes::fit_model(data.train = splits[[1]]$wdf,
#                                  data.test  = splits[[2]]$wdf,
#                                  rnd.seed   = 1234,
#                                  ncpus      = parallel::detectCores() - 2)
#  
#  # Calculate performance
#  my_perf <- epitopes::calc_performance(truth = splits[[2]]$wdf$Class,
#                                        pred  = my.model$rf_class,
#                                        prob  = my.model$rf_probs,
#                                        ret.as.list = TRUE)

## ---- eval=FALSE--------------------------------------------------------------
#  # Plot ROC curve
#  plot(my_perf$fpr, my_perf$tpr, type = "p", pch = 20, las = 1,
#       xlab = "FPR", ylab = "TPR",
#       main = "ROC curve for O. volvulus predictor",
#       sub = paste("AUC = ", signif(my_perf$auc, digits = 3)))
#  
#  print(unlist(my_perf[c(5:12)]))

## ---- include=FALSE-----------------------------------------------------------
# this comes from the actual simulation
perf_vec <- c(0.7624003, 0.6491852, 0.7012283, 0.7167080, 
              0.7305359, 0.4147487, 0.7079694, 0.7783327)
names(perf_vec) <- c("sens", "spec", "ppv", "npv", 
                     "f1", "mcc", "accuracy", "auc")
print(perf_vec)

## ---- eval = FALSE------------------------------------------------------------
#  # Get the first protein from the hold-out set as an example:
#  myprot <- proteins[proteins$UID == splits[[2]]$wdf$Info_protein_id[1]]
#  
#  # Make windowed representation:
#  myprot_w <- epitopes::make_window_df(myprot,
#                                       window_size = nchar(splits[[2]]$wdf$Info_window_seq[1]))
#  myprot_w <- epitopes::calc_features(myprot_w, max.N = 2,
#                                      ncpus = parallel::detectCores() - 2)
#  
#  # Get predictions from the model
#  myprot_w <- as.data.frame(myprot_w)[, which(!grepl("Info_", names(myprot_w)))]
#  preds    <- stats::predict(my.model$rf_model,
#                             data = myprot_w)

## ---- eval=FALSE--------------------------------------------------------------
#  # Smooth predictions (to remove very short positive sequences)
#  myclass <- epitopes::smooth_predictions(as.numeric(preds$predictions[, 2] > 0.5),
#                                          minsize = 8)
#  
#  plot(preds$predictions[, 2], type = "l", lwd = .5, las = 1,
#       main = paste("Predictions for protein", myprot$UID[1]),
#       xlab = "Protein position",
#       ylab = "Probability / Prediction",
#       ylim = c(0, 1.05))
#  points(myclass, type = "p", pch = 20, col = myclass + 1)
#  points(myclass, type = "l", lwd = .3, col = "#77777777")

