#' Reorder Components in a Lucid Model
#' 
#' This function reorders various components within a Lucid model based on the specified model type and order. Only work for K = 2.
#' 
#' @param lucidus_fit The fitted Lucid model object.
#' @param lucid_model A character string indicating the type of Lucid model ("early", "parallel", or "serial_sub_with_sub_before_it_reordered").
#' @param order An optional parameter specifying the order for reordering.
#' 
#' @return The modified Lucid model object after reordering components based on the specified model type and order.
#' 
#' @details
#' The function is designed to reorder different components of a Lucid model based on the specified model type.
#' For "early" model type, it reorders components such as res_Beta, res_Mu, res_Sigma, res_Gamma, and inclusion.p based on the specified order.
#' For "parallel" model type, it performs specific reordering operations on components based on the given order.
#' For "serial_sub_with_sub_before_it_reordered" model type, it reorders components in a particular manner and returns the modified object.
#' 
#' @examples
#' # Example usage:
#' # reorder_lucid(lucidus_fit, "early", order = 2)
#' 
#' @export
#' 
reorder_lucid <- function(lucidus_fit,lucid_model,
                          order = NULL) {
  if (lucid_model == "early"){
    # reorder beta
    GtoX <- lucidus_fit$res_Beta
    lucidus_fit$res_Beta <- (-1)^(order - 1) * GtoX
    
    # reorder mu
    XtoZ <- lucidus_fit$res_Mu
    x <- c(1, 2) # order of clusters
    if(order == 2) {
      x <- c(2, 1)
      XtoZ = XtoZ[x, ]
    } else{
      XtoZ = XtoZ[x, ]
    }
    lucidus_fit$res_Mu = XtoZ
    
    # reorder Sigma
    Sigma <- lucidus_fit$res_Sigma
    
    if(order == 2) {
      Sigma_original <- Sigma
      Sigma[[1]] = Sigma_original[[2]]
      Sigma[[2]] = Sigma_original[[1]]
    } 
    
    lucidus_fit$res_Sigma = Sigma
    
    
    # reorder gamma
    XtoY <- lucidus_fit$res_Gamma$beta
    XtoY[1] <- XtoY[1] + sum(XtoY[-1] * (order - 1)) # reference level using the new reference
    XtoY[-1] <- (-1)^(order - 1) * XtoY[-1] # if reference = 2, flip the estimates
    lucidus_fit$res_Gamma$beta <- XtoY
    #lucidus_fit$res_Gamma$fit$coefficients <- XtoY
    
    # reorder inclusion.p
    IP <- lucidus_fit$inclusion.p
    
    if(order == 2) {
      IP_original <- IP
      IP[,1] = IP_original[,2]
      IP[,2] = IP_original[,1]
    } 
    lucidus_fit$inclusion.p = IP
    
    # return the object using the new reference
    return(lucidus_fit)
    
  }else if (lucid_model == "parallel"){
    if(is.null(order)) {
      warning("no reference specified, return the original model")
      return(lucidus_fit)
    }
    reference <- order
    n_omic <- length(order)
    
    # reorder beta
    GtoX <- lucidus_fit$res_Beta$Beta
    lucidus_fit$res_Beta$Beta <- lapply(1:n_omic, function(i) {
      (-1)^(reference[i] - 1) * GtoX[[i]] # if reference = 1, no changes; 
      # if reference = 2, flip the reference and negate the estimates
      # intercept is wrong actually, need to change, but we might not use it
    })
    
    
    # reorder mu
    XtoZ <- lucidus_fit$res_Mu
    lucidus_fit$res_Mu <- lapply(1:n_omic, function(i) {
      x <- c(1, 2) # order of clusters
      if(reference[i] == 2) {
        x <- c(2, 1)
        XtoZ[[i]][, x]
      } else{
        XtoZ[[i]][, x]
      }
    }) 
    
    # reorder Sigma
    Sigma <- lucidus_fit$res_Sigma
    lucidus_fit$res_Sigma <- lapply(1:n_omic, function(i) {
      x <- c(1, 2) # order of clusters
      if(reference[i] == 2) {
        x <- c(2, 1)
        Sigma[[i]][, , x]
      } else{
        Sigma[[i]][, , x]
      }
    }) 
    
    # reorder gamma
    XtoY <- lucidus_fit$res_Gamma$fit$coefficients[1:4]
    XtoY[1] <- XtoY[1] + sum(XtoY[-1] * (reference - 1)) # reference level using the new reference
    XtoY[-1] <- (-1)^(reference - 1) * XtoY[-1] # if reference = 2, flip the estimates
    lucidus_fit$res_Gamma$Gamma$mu <- XtoY
    lucidus_fit$res_Gamma$fit$coefficients <- XtoY
    
    # reorder ip
    ip <- lucidus_fit$inclusion.p
    lucidus_fit$inclusion.p <- lapply(1:n_omic, function(i) {
      x <- c(1, 2) # order of clusters
      if(reference[i] == 2) {
        x <- c(2, 1)
        ip[[i]][, x]
      } else{
        ip[[i]][, x]
      }
    }) 
    
    # return the object using the new reference
    return(lucidus_fit)
    
  }else if (lucid_model == "serial_sub_with_sub_before_it_reordered"){
    #to reorder a regular submodel for LUCID, just use early
    #Hoever, after reordering, the next sub model's beta need to be ordered..
    #the previous layer is reordered, leading to the change of beta...
    
    # reorder beta
    GtoX <- lucidus_fit$res_Beta
    lucidus_fit$res_Beta <- (-1)^(order - 1) * GtoX
    
    
    # return the object using the new reference
    return(lucidus_fit)
  }
  
}


library(venn)
#' Generate Venn Diagrams for Feature Comparison
#' 
#' This function generates a Venn diagram comparing features obtained from different stages 
#' (Early, Intermediate, Late) across three LUCID models: early, parallel, and serial.
#' 
#' @param fit_early Fitted LUCID model for the early stage.
#' @param fit_parallel Fitted LUCID model for the parallel stage.
#' @param fit_serial Fitted LUCID model for the serial stage.
#' @param pattern A pattern used for feature selection.
#' @param omics The type of omics data (e.g., genomics, proteomics) being compared.
#' 
#' @details
#' The function extracts features from each stage of the LUCID models based on the provided pattern 
#' and generates a Venn diagram to visualize the overlap or unique features among the stages.
#' 
#' @examples
#' # Example usage:
#' # venn_feature(fit_early, fit_parallel, fit_serial, "pattern", "omics_type")
#' 
#' @export
venn_feature <- function(fit_early, fit_parallel, fit_serial, pattern, omics) {
  z_early = fit_early$var.names[[2]]
  z_early = z_early[grep(pattern, z_early)]
  
  z_parallel = unlist(fit_parallel$var.names$Znames)
  z_parallel = z_parallel[grep(pattern, z_parallel)]
  
  z_serial = unlist(fit_serial$var.names$Znames)
  z_serial = z_serial[grep(pattern, z_serial)]
  
  list_feature = list(z_early , z_parallel, z_serial)
  names(list_feature) = c("Early","Intermediate","Late")
  venn(list_feature,zcolor = "style",ilabels = TRUE,plotsize = 8,box = FALSE,sncs= 0.4)
  
}

#' Generate Venn Diagrams for High-Risk Group Comparison
#' 
#' This function generates a Venn diagram comparing high-risk groups across different omics layers 
#' within the three LUCID models: early, parallel, and serial.
#' 
#' @param fit_early Fitted LUCID model for the early stage.
#' @param fit_parallel Fitted LUCID model for the parallel stage.
#' @param fit_serial Fitted LUCID model for the serial stage.
#' @param omics_layer_index Index of the omics layer for comparison.
#' 
#' @details
#' The function identifies high-risk clusters within each omics layer of the LUCID models 
#' and generates a Venn diagram to visualize the overlap or unique clusters among the stages.
#' 
#' @examples
#' # Example usage:
#' # venn_cluster(fit_early, fit_parallel, fit_serial, omics_layer_index)
#' 
#' @export
venn_cluster <- function(fit_early, fit_parallel, fit_serial, omics_layer_index) {
  ip = as.data.frame(fit_early$inclusion.p)
  cluster_early = ifelse(ip$V1 > ip$V2, 0, 1)
  
  ip = as.data.frame(fit_parallel$inclusion.p[[omics_layer_index]])
  cluster_parallel = ifelse(ip$V1 > ip$V2, 0, 1)
  
  ip = as.data.frame(fit_serial$inclusion.p[[omics_layer_index]])
  cluster_serial = ifelse(ip$V1 > ip$V2, 0, 1)
  
  name_e <- rep(NA, length(cluster_early))
  name_e[cluster_early == 1] <- which(cluster_early == 1)
  name_e = name_e[complete.cases(name_e)]
  
  name_p_m <- rep(NA, length(cluster_parallel))
  name_p_m[cluster_parallel == 1] <- which(cluster_parallel == 1)
  name_p_m = name_p_m[complete.cases(name_p_m)]
  
  name_s_m <- rep(NA, length(cluster_serial))
  name_s_m[cluster_serial == 1] <- which(cluster_serial == 1)
  name_s_m = name_s_m[complete.cases(name_s_m)]
  
  
  list_r = list(name_e, name_p_m, name_s_m)
  names(list_r) = c("Early","Intermediate","Late")
  venn(list_r,zcolor = "style",ilabels = TRUE,plotsize = 8,box = FALSE,sncs= 0.4)
  
}


#' Generate Venn Diagrams for High-Risk Group Comparison Within a Model
#' 
#' This function generates a Venn diagram comparing high-risk groups within a single LUCID model 
#' across different omics layers.
#' 
#' @param fit Fitted LUCID model.
#' @param omics_list_name List of omics layer names.
#' 
#' @details
#' The function identifies high-risk clusters within each omics layer of the provided LUCID model 
#' and generates a Venn diagram to visualize the overlap or unique clusters among the layers.
#' 
#' @examples
#' # Example usage:
#' # venn_cluster_within_model(fit, omics_list_name)
#' 
#' @export
venn_cluster_within_model <- function(fit, omics_list_name) {
  layer_num = length(omics_list_name)
  cluster_list = rep(list(0), layer_num)
  name_list = rep(list(0), layer_num)
  for (i in 1:layer_num){
    ip = as.data.frame(fit$inclusion.p[[i]])
    cluster_list[[i]] = ifelse(ip$V1 > ip$V2, 0, 1)
    name_list[[i]] <- rep(NA, length(cluster_list[[i]]))
    name_list[[i]][cluster_list[[i]] == 1] <- which(cluster_list[[i]] == 1)
    name_list[[i]] = name_list[[i]][complete.cases(name_list[[i]])]
  }
  
  
  list_r = name_list
  names(list_r) = omics_list_name
  venn(list_r,zcolor = "style",ilabels = TRUE,plotsize = 8,box = FALSE,sncs= 0.4)
}


extract_X_assignment_auxi <- function(lucid_model, type) {
  ip = lucid_model$inclusion.p
  if (type == "early"){
    n = nrow(lucid_model$Z)
    x <- sapply(1:n, function(x) return(nnet::which.is.max(ip[x, ])))
    x <- x - 1 #reference cluster is 0
  }else if (type == "parallel"){
    nOmics = length(lucid_model$K)
    n = lucid_mode$N
    x <- vector(mode = "list", length = nOmics)
    for (i in 1:nOmics) {
      x[[i]] <- sapply(1:n, function(x) return(nnet::which.is.max(ip[[i]][x, ])))
      x[[i]] <- pred_X[[i]] - 1
    }
  }
  return(x)
}

