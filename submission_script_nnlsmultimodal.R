##################################################################################################
### PLEASE only edit the program function between YOUR CODE BEGINS/ENDS HERE                   ###
##################################################################################################

#' The function to predict the A matrix
#' In the provided example, we use basic non-negative least squares (package "nnls"), which consists in minimizing the error term $||Mix - Ref \times Prop||^2$ with only positive entries in the prop matrix.
#'
#' @param mix a matrix of bulk samples (columns) and features (rows)
#' @param ref a matrix of pure cell types (columns) and features (rows)
#' @param ... other parameters that will be ignored
#' 
#' @return the predicted A matrix
#' @examples
#' 
program = function(mix_rna=NULL, 
                   ref_bulkRNA=NULL, 
                   mix_met=NULL,
                   ref_met=NULL,
                   ref_scRNA=NULL) {
  ##
  ## YOUR CODE BEGINS HERE
  ##

  # # Pseudo bulk reference from scRNAseq
  # ref_sc_peng <- as.matrix(ref_scRNA$ref_sc_peng$counts)
  # colnames(ref_sc_peng) <- ref_scRNA$ref_sc_peng$metadata$cell_type
  # i1 <- which(colnames(ref_sc_peng)=="immune")
  # i2 <- which(colnames(ref_sc_peng)=="endo")
  # i3 <- which(colnames(ref_sc_peng)=="fibro")
  # i4 <- which(colnames(ref_sc_peng)=="classic")
  # i5 <- which(colnames(ref_sc_peng)=="basal")
  # ref_immune <- ref_sc_peng[,i1]
  # ref_endo <- ref_sc_peng[,i2]
  # ref_fibro <- ref_sc_peng[,i3]
  # ref_classic <- ref_sc_peng[,i4]
  # ref_basal <- ref_sc_peng[,i5]
  # # 
  # ref_sc_constrained <- cbind(rowS(ref_immune), 
  #                             rowMeans(ref_endo), 
  #                             rowMeans(ref_fibro), 
  #                             rowMeans(ref_classic), 
  #                             rowMeans(ref_basal))
  # colnames(ref_sc_constrained) = c("immune", "endo", "fibro", "classic", "basal")

  # Function for filtering
  calculate_fold_change <- function(x) {
    sorted_values <- sort(x, decreasing = TRUE) # Sort expression values for the gene
    if (length(sorted_values) > 1) {
      return(sorted_values[1] / sorted_values[2]) # Fold change: top value / second top value
    } else {
      return(NA) # Handle cases where there are less than 2 cell types
    }
  }

  # Prefiltering on ref_bulkRNA
  # Identify genes with highest FC between first and second most expressed cell types. 
  # Calculate the fold change between the top two cell types for each gene
  # Apply the function to each row of the matrix
  fold_changes <- apply(ref_bulkRNA, 1, calculate_fold_change)

  # Create a dataframe of genes and their fold changes
  fold_change_df <- data.frame(
    Gene = rownames(ref_bulkRNA),
    FoldChange = fold_changes,
      log2FC = log2(fold_changes+1)
  )

  # Filter for non-NA fold changes and sort by fold change in descending order
  top_genes <- fold_change_df[!is.na(fold_change_df$FoldChange), ]
  top_genes <- top_genes[order(-top_genes$FoldChange), ]

  # Select the top 500 genes
  top_genes_FC <- head(top_genes, 300)

  genelist_top300FC = top_genes_FC$Gene
  new_ref_bulkRNA <- ref_bulkRNA[genelist_top300FC,]

  # Prefiltering on met

  # Apply the function to each row of the matrix
  met_fold_changes <- apply(ref_met, 1, calculate_fold_change)

  # Create a dataframe of sites and their fold changes
  met_fold_change_df <- data.frame(
      Site = rownames(ref_met),
      FoldChange = met_fold_changes,
      met_log2FC = log2(met_fold_changes+1),
      rank_FC = rank(-met_fold_changes)
  )

  # Filter for non-NA fold changes and sort by fold change in descending order
  top_sites <- met_fold_change_df[!is.na(met_fold_change_df$FoldChange), ]
  top_sites <- top_sites[order(-top_sites$FoldChange), ]

  # Select the top 500 sites
  top_sites_FC <- head(top_sites, 1000)

  sitelist_top1000FC = top_sites_FC$Site
  new_ref_met <- ref_met[sitelist_top1000FC,]
  
  # # Instlaling EpiDISH for the RCP function
  # BiocManager::install("EpiDISH")
  # library(EpiDISH)
  
  if ( !( is.null(x = mix_rna) ) ) {

    # prop_rna = t(epidish(beta.m = mix_rna, ref.m = ref_bulkRNA, method = "RPC")$estF)
    
    idx_feat = intersect(rownames(mix_rna), rownames(new_ref_bulkRNA))
    mix_rna = mix_rna[idx_feat,]

  # Add a tranformation to mix_rna  

    new_ref_bulkRNA = new_ref_bulkRNA[idx_feat,]
    
    prop_rna = apply(mix_rna, 2, function(b, A) {
      tmp_prop = nnls::nnls(b=b, A=A)$x
      tmp_prop = tmp_prop / sum(tmp_prop) # Sum To One
      return(tmp_prop)
    }, A=new_ref_bulkRNA)  
    rownames(prop_rna) = colnames(new_ref_bulkRNA)

    # # Test a random solution
    # # Set a seed for reproducibility
    # set.seed(4)

    # nrows <- dim(ref_bulkRNA)[2]
    # ncols <- dim(mix_rna)[2]

    # # Create a dataframe with random values between 0 and 1 for each cell
    # prop_rna <- matrix(runif(nrows * ncols), nrow = nrows, ncol = ncols)

    # rownames(prop_rna) <- colnames(ref_bulkRNA)
    # colnames(prop_rna) <- colnames(mix_rna)
    
  }

  # return(prop_rna)

  # print(prop_rna)

  
  # ## we compute the estimation of A for the methylation data set :
  if ( !( is.null(mix_met) ) ) {

    # prop_met = t(epidish(beta.m = mix_met, ref.m = ref_met, method = "RPC")$estF)

    # colnames(prop_met) = colnames(ref_met)
    # print(colnames(ref_met))
    
    idx_feat = intersect(rownames(mix_met), rownames(new_ref_met))
    mix_met = mix_met[idx_feat,]
 
    # Linear methods
    new_ref_met = new_ref_met[idx_feat,]
    
    prop_met = apply(mix_met, 2, function(b, A) {
      tmp_prop = nnls::nnls(b=b, A=A)$x
      tmp_prop = tmp_prop / sum(tmp_prop) # Sum To One
      return(tmp_prop)
    }, A=new_ref_met)  
    rownames(prop_met) = colnames(new_ref_met)

    return(prop_met)

    }
      

  # ## we compute the mean of all the estimated A matrices as the final A matrix :
  
  
  if ( !is.null(x = mix_met) ) {
    if ( !is.null(x = mix_rna) ) {
      
      stopifnot( identical(x = dim(prop_rna), y = dim(prop_met)) ) # stop if not same number of cell type or samples
      ## we have an estimation of prop based on the the methylation and transcriptome data sets
      prop <- (prop_rna + prop_met) / 2
      
    } else {
      prop <- prop_met
    }
  } else {
    prop <- prop_rna
  }
  
  if (any(colSums(prop) != 1)) { # Sum To One 
    prop <- sapply(
      1:ncol(prop), 
      function(col_i) { prop[,col_i] / sum(prop[,col_i]) }
    )
  }

  return(prop)

  ##
  ## YOUR CODE ENDS HERE
  ##
}


##############################################################
### Generate a prediction file /!\ DO NOT CHANGE THIS PART ###
##############################################################

validate_pred <- function(pred, nb_samples = ncol(mix_rna) , nb_cells= ncol(ref_rna),col_names = colnames(ref_met) ){

  error_status = 0   # 0 means no errors, 1 means "Fatal errors" and 2 means "Warning"
  error_informations = ''

  ## Ensure that all sum ofcells proportion approximately equal 1
  if (!all(sapply(colSums(pred), function(x) isTRUE(all.equal(x, 1) )))) {
    msg = "The prediction matrix does not respect the laws of proportions: the sum of each columns should be equal to 1\n"
    error_informations = paste(error_informations,msg)
    error_status = 2
  }

  ##Ensure that the prediction have the correct names ! 
  if(! setequal(rownames(pred),col_names) ){
    msg = paste0(    "The row names in the prediction matrix should match: ", toString(col_names),"\n")
    error_informations = paste(error_informations,msg)
    error_status = 2
  }

  ## Ensure that the prediction return the correct number of samples and  number of cells. 
  if (nrow(pred) != nb_cells  | ncol(pred) != nb_samples)  {
    msg= paste0('The prediction matrix has the dimention: ',toString(dim(pred))," whereas the dimention: ",toString(c(nb_cells,nb_samples))," is expected\n"   )
    error_informations = paste(error_informations,msg)
    error_status = 1
  }

  if(error_status == 1){
    # The error is blocking and should therefor stop the execution. 
    # tryCatch(message("hello\n"), message=function(e){cat("goodbye\n")})  use this here ? 
    stop(error_informations)
  }
  if(error_status == 2){
    print("Warning: ")
    warning(error_informations)
  }  
}

dir_name = paste0("data",.Platform$file.sep)
dataset_list = list.files(dir_name,pattern="mixes*")

reference_data <- readRDS(file =  paste0(dir_name, "reference_pdac.rds"))


predi_list = list()
for (dataset_name in dataset_list){

  print(paste0("generating prediction for dataset:",toString(dataset_name) ))

  mixes_data <- readRDS(file = paste0(dir_name, dataset_name))

  if ("mix_rna" %in% names(mixes_data)) {
    mix_rna = mixes_data$mix_rna
  } else {
    mix_rna = mixes_data
  }
  if ("mix_met" %in% names(mixes_data)) {
    mix_met = mixes_data$mix_met  
  } else {
    mix_met = NULL
  }

  if ("ref_bulkRNA" %in% names(reference_data)) {
    ref_bulkRNA = reference_data$ref_bulkRNA
  } else {
    ref_bulkRNA = reference_data
  }
  if ("ref_met" %in% names(reference_data)) {
    ref_met = reference_data$ref_met  
  } else {
    ref_met = NULL
  }
  if ("ref_scRNA" %in% names(reference_data)) {
    ref_scRNA = reference_data$ref_scRNA  
  } else {
    ref_scRNA = NULL
  }

  # we use the previously defined function 'program' to estimate A :
  pred_prop <- program(mix_rna, ref_bulkRNA, mix_met=mix_met, ref_met=ref_met, ref_scRNA=ref_scRNA)
  validate_pred(pred_prop,nb_samples = ncol(mix_rna),nb_cells = ncol(ref_bulkRNA),col_names = colnames(ref_met))
  predi_list[[dataset_name]] = pred_prop

}


##############################################################
### Check the prediction /!\ DO NOT CHANGE THIS PART ###
##############################################################


###############################
### Code submission mode


print("")
for (package in c("zip") ) {
  if ( !{ package %in% installed.packages( ) } ) {
        print(x = paste("Installation of ", package, sep = "") )
        install.packages(
            pkgs = "zip"
          , repos = "https://cloud.r-project.org"
        )
    } 
}


# we generate a zip file with the 'program' source code

if ( !dir.exists(paths = "submissions") ) {
    dir.create(path = "submissions")
}

# we save the source code as a R file named 'program.R' :
dump(
    list = c("program")
    # list = new_functions
  , file = paste0("submissions", .Platform$file.sep, "program.R")
)

date_suffix = format(x = Sys.time( ), format = "%Y_%m_%d_%H_%M_%S")



zip_program <- paste0("submissions", .Platform$file.sep, "program_", date_suffix, ".zip")
zip::zip(zipfile= zip_program
  , files   = paste0("submissions", .Platform$file.sep, "program.R")
  , mode = "cherry-pick")

if(dir.exists("attachement")) {
  zip::zip_append(
      zipfile = zip_program
      , files= paste0("attachement", .Platform$file.sep)
    , mode = "cherry-pick"
  )
}

zip::zip_list(zip_program)
print(x = zip_program)




# # we create the associated zip file :
# zip_program <- paste0("submissions", .Platform$file.sep, "program_", date_suffix, ".zip")
# zip::zip(zipfile= zip_program
#                 , files= paste0("submissions", .Platform$file.sep, "program.R")
#                 , mode = "cherry-pick"
#                 )

# zip::zip_list(zip_program)
# print(x = zip_program)

###############################
### Result submission mode  

#  Generate a zip file with the prediction
if ( !dir.exists(paths = "submissions") ) {
    dir.create(path = "submissions")
}

prediction_name = "prediction.rds"

## we save the estimated A matrix as a rds file named 'results.rds' :
saveRDS(
object = predi_list
, file   = paste0("submissions", .Platform$file.sep, prediction_name)) 

# write_rds(pred_prop, file = "prediction_hugo.rds")

## we create the associated zip file :
zip_results <- paste0("submissions", .Platform$file.sep, "results_", date_suffix, ".zip")
zip::zipr(
         zipfile = zip_results
       , files   = paste0("submissions", .Platform$file.sep, c(prediction_name) )
     )
print(x = zip_results)

sessionInfo( )

###############################################################
### How to submit the zip file? /!\ DO NOT CHANGE THIS PART ###
###############################################################
#
# The code above generates the files *`r zip_program`*  and *`r zip_results`*  (the 1st one for code submission, the 2nd one for result submission).
#
# Submit the zip submission file on the challenge in the `My Submission` tab, fill the metadata, select the task you want to submit to and upload your submission files
#
# On the codalab challenge web page, The *STATUS* become :
#   - Submitting
#   - Submitted
#   - Running
#   - Finished
#
# When it’s finished :
#   - You refresh the page and click on the green button 'add to leaderboard' to see your score
#   - If enable, details for report could be downloaded by clicking on your submission
#   - Some execution logs are available to check and or download.
#   - Metadata are editable when you click on your submission
#   - Leader board is updated in the `Results` tab.
#

