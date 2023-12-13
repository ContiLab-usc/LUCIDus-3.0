library(purrr)
library(dplyr)

#' Plot Omics Profiles from LUCID Models
#' 
#' This function generates bar plots to visualize omics profiles obtained from different stages 
#' (Early, Intermediate, Late) within LUCID models: early, parallel, or serial, based on the integration type.
#' 
#' @param fit Fitted LUCID model object.
#' @param integration_type Type of integration ("early", "intermediate", "late").
#' @param pattern_list List of patterns used for feature selection in each omics layer.
#' @param omics_list List of omics data types (e.g., genomics, proteomics) being compared.
#' 
#' @return A bar plot representing omics profiles across different stages and omics layers.
#' 
#' @details
#' The function extracts omics profiles from different stages of LUCID models based on the integration type 
#' and generates bar plots to display the profiles of features selected by provided patterns in each omics layer.
#' 
#' @examples
#' # Example usage:
#' # plot_omics_profiles(fit, "early", pattern_list, omics_list)
#' 
#' @export

plot_omics_profiles <- function(fit, integration_type, pattern_list, omics_list) {
  list_num = lapply(1:length(omics_list), function(x) paste(x))
  color_p = c("#2fa4da", "#A77E69", "#e7b6c1","yellow", "green")
  if(integration_type == "early"){
    M_mean = as.data.frame(fit$res_Mu)
    M_mean$cluster = as.factor(1:2)
    # Reshape the data
    M_mean_melt <- M_mean %>% 
      pivot_longer(cols = -cluster, names_to = "variable", values_to = "value")
    
    M_mean_melt <- M_mean_melt %>% 
      mutate(cluster = ifelse(cluster == 2, "High Risk", "Low Risk"))
    
    # add color label for omics layer
    M_mean_melt = M_mean_melt %>%
      mutate(color_label <- map2(pattern_list, omics_list, 
                                 ~ case_when(
                                   str_detect(variable, .x) ~ .y,
                                   TRUE ~ NA_character_  # Default condition if no pattern matches
                                 )
      ) %>%
        reduce(coalesce) %>%
        tibble(color_label = .)
      
      ) 
    
    
    
    
    fig <- ggplot(M_mean_melt, 
                  aes(fill = color_label, y = value, x = variable)) +
      geom_bar(position="dodge", stat="identity") +
      ggtitle("Omics profiles for the two latent clusters") +
      facet_grid(rows = vars(cluster), cols = vars(color_label),  scales = "free_x", space = "free") +
      theme(legend.position="none") +
      geom_hline(yintercept = 0) +
      xlab("") +
      theme(text = element_text(size=10),
            axis.text.x = element_text(angle = 90, vjust = 1,
                                       hjust = 1),
            plot.margin = margin(10, 10, 10, 80),
            panel.background = element_rect(fill="white"), 
            strip.background = element_rect(fill = "white"),
            axis.line.x = element_line(color = "black"),
            axis.line.y = element_line(color = "black"),) +
      scale_fill_manual(values = color_p[1:length(omics_list)])
    
  } else if(integration_type == "intermediate"){
    M_mean = map_dfr(fit$res_Mu, ~ as_tibble(.x, rownames = "variable")) %>%
      mutate(omic_layer = map2(pattern_list, omics_list, 
                               ~ case_when(
                                 str_detect(variable, .x) ~ .y,
                                 TRUE ~ NA_character_  # Default condition if no pattern matches
                               )
      ) %>%
        reduce(coalesce))
    
    # Reorder results because mirna order is reversed
    M_mean1 <- M_mean %>% 
      #left_join(meta_df, by = c("variable" = "ftr_name")) %>%
      # mutate(`Low Risk`  =  if_else(omic_layer == "RNA", V2, V1), 
      #`High Risk` =  if_else(omic_layer == "RNA", V1, V2)) %>%
      mutate(`Low Risk`  =  V1, 
             `High Risk` =  V2) %>%
      dplyr::select(-c("V1", "V2"))
    
    # Pivot longer for figure 
    M_mean_l <- M_mean1 %>% 
      pivot_longer(cols = c(`Low Risk`, `High Risk`),
                   names_to = "cluster",
                   values_to = "value")
    
    # add color label for omics layer
    M_mean2 = M_mean_l %>%
      mutate(color_label = as.character(as.numeric(factor(omic_layer, levels = unique(omic_layer)))), 
             low_high = if_else(str_detect(cluster, "Low"), 0,1),
             omic = if_else(omic_layer == "miRNA", 
                            "mir",
                            str_sub(omic_layer, end = 1) %>% toupper()),
             omic_cluster = str_c(omic, low_high))
    
    # Filter only the top ## differential expressed features 
    M_mean2_top <- M_mean2 %>% 
      group_by(variable) %>% 
      filter(abs(value) == max(abs(value))) %>% 
      ungroup() %>% 
      arrange(max(abs(value))) %>% 
      group_by(omic_layer) %>% 
      slice_head(n=12) %>%
      ungroup()
    
    # Plots top 12 features
    fig <- ggplot(M_mean2  %>% filter(variable %in% M_mean2_top$variable),
                  aes(fill = color_label, y = value, x = variable)) +
      geom_bar(position="dodge", stat="identity") +
      ggtitle("Omics profiles for 2 latent clusters - Lucid in Parallel") +
      facet_grid(rows = vars(cluster),
                 cols = vars(omic_layer), scales = "free_x", space = "free") +
      theme(legend.position="none") +
      geom_hline(yintercept = 0) +
      xlab("") +
      theme(text = element_text(size=10),
            axis.text.x = element_text(angle = 90, vjust = 1,
                                       hjust = 1),
            plot.margin = margin(10, 10, 10, 80),
            panel.background = element_rect(fill="white"),
            strip.background = element_rect(fill = "white"),
            axis.line.x = element_line(color = "black"),
            axis.line.y = element_line(color = "black"),) +
      scale_fill_manual(values = color_p[1:length(omics_list)])
  }
  
  if(integration_type == "late"){
    M_mean = map_dfr(fit$submodel, ~ as_tibble(t(.x$res_Mu), rownames = "variable")) %>%
      mutate(omic_layer = map2(pattern_list, omics_list, 
                               ~ case_when(
                                 str_detect(variable, .x) ~ .y,
                                 TRUE ~ NA_character_  # Default condition if no pattern matches
                               )
      ) %>%
        reduce(coalesce))
    
    # Reorder results because mirna order is reversed
    M_mean1 <- M_mean %>% 
      #left_join(meta_df, by = c("variable" = "ftr_name")) %>%
      mutate(`Low Risk`  =  V1, 
             `High Risk` =  V2) %>%
      dplyr::select(-c("V1", "V2"))
    
    # Pivot longer for figure 
    M_mean_l <- M_mean1 %>% 
      pivot_longer(cols = c(`Low Risk`, `High Risk`),
                   names_to = "cluster",
                   values_to = "value")
    
    # add color label for omics layer
    M_mean2 = M_mean_l %>%
      mutate(color_label = as.character(as.numeric(factor(omic_layer, levels = unique(omic_layer)))),
             low_high = if_else(str_detect(cluster, "Low"), 0,1),
             omic = if_else(omic_layer == "miRNA", 
                            "mir",
                            str_sub(omic_layer, end = 1) %>% toupper()),
             omic_cluster = str_c(omic, low_high))
    
    # Filter only the top ## differential expressed features 
    M_mean2_top <- M_mean2 %>% 
      group_by(variable) %>% 
      filter(abs(value) == max(abs(value))) %>% 
      ungroup() %>% 
      arrange(max(abs(value))) %>% 
      group_by(omic_layer) %>% 
      slice_head(n=12) %>%
      ungroup()
    
    # Plots top 12 features
    fig <- ggplot(M_mean2  %>% filter(variable %in% M_mean2_top$variable),
                  aes(fill = color_label, y = value, x = variable)) +
      geom_bar(position="dodge", stat="identity") +
      ggtitle("Omics profiles for 2 latent clusters - Lucid in Parallel") +
      facet_grid(rows = vars(cluster),
                 cols = vars(omic_layer), scales = "free_x", space = "free") +
      theme(legend.position="none") +
      geom_hline(yintercept = 0) +
      xlab("") +
      theme(text = element_text(size=10),
            axis.text.x = element_text(angle = 90, vjust = 1,
                                       hjust = 1),
            plot.margin = margin(10, 10, 10, 80),
            panel.background = element_rect(fill="white"),
            strip.background = element_rect(fill = "white"),
            axis.line.x = element_line(color = "black"),
            axis.line.y = element_line(color = "black"),) +
      scale_fill_manual(values = color_p[1:length(omics_list)])
  }
  
  return(fig)
}




