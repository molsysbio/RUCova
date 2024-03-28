
filter <- dplyr::filter

# compare before/after
#pearson correlation

library(ggpubr)
library(ggcorrplot)
library(ggridges)
library(ggh4x)
library(ggplot2)



mean_center <- function(x){
  x-mean(x)
}

#' Regress surrogates of unwanted covariance
#'
#' @param data A tibble.
#' @param markers Vector of marker values to normalise.
#' @param surrogates Vector of surrogates to use for normalisation.
#' @param approach Should be "per_sample" if XXXX 
#' @param out_ruc_obj Some object 
#' @param clust Some clust 
#' 
#' @return A ggplot object that plots 
#' @examples
#'
#' @import ggpubr
#' @import ggcorrplot
#' @import ggridges
#' @import ggh4x
#' @import ggplot2
#' @export
plot_pearson_corr <- function(data, surrogates, markers, approach, out_ruc_obj, clust){

  if (approach == "per_sample") {

    data$sample <- data[[out_ruc_obj$col_name_sample]] #create column "sample" for simplicity
    out_ruc_obj$data_reg$sample <- out_ruc_obj$data_reg[[out_ruc_obj$col_name_sample]] #create column "sample" for simplicity

    p <- list()


    for (l in unique(data$sample)){

      t <- paste0(out_ruc_obj$col_name_sample," ", l)

      a <- data %>%
        filter(sample == l) %>%
        select(surrogates,markers) %>%
        mutate_all(asinh) %>%
        cor(method = "pearson") %>%
        ggcorrplot(type = "upper",outline.color = "white", show.legend = TRUE, title = paste0("Before, ",t), hc.order = clust)

      b <- out_ruc_obj$data_reg %>%
        filter(sample == l) %>%
        select(surrogates,markers) %>%
        mutate_all(asinh) %>%
        cor(method = "pearson") %>%
        ggcorrplot(type = "upper", outline.color = "white", show.legend = TRUE, title = "After", hc.order = clust)


      p[[l]] <- ggarrange(plotlist = list(a,b), ncol = 2)


    }

    } else {

      a <- data %>%
        select(surrogates,markers) %>%
        mutate_all(asinh) %>%
        cor(method = "pearson") %>%
        ggcorrplot(type = "upper",outline.color = "white", show.legend = TRUE, title = "Before", hc.order = clust)


      b <- out_ruc_obj$data_reg %>%
        select(surrogates,markers) %>%
        mutate_all(asinh) %>%
        cor(method = "pearson") %>%
        ggcorrplot(type = "upper",outline.color = "white", show.legend = TRUE, title = "After", hc.order = clust)


      p <- ggarrange(plotlist = list(a,b), ncol = 2)

      }



  return(p)

}

# scatter plots with regression sample

#' @export
plot_scatter_fit <- function(data, m, surrogates, out_ruc_obj, frac){


  data$sample <- data[[out_ruc_obj$col_name_sample]] #create column "sample" for simplicity
  out_ruc_obj$data_reg$sample <- out_ruc_obj$data_reg[[out_ruc_obj$col_name_sample]]

  if (out_ruc_obj$center_surr == "per_sample") {

    dt_before <- data %>%
      mutate_at(vars(m,surrogates),asinh) %>%
      group_by(sample) %>%
      mutate_at(vars(surrogates),mean_center)


  } else if (out_ruc_obj$center_surr == "across_samples") {

    dt_before <- data %>%
      mutate_at(vars(m,surrogates),asinh) %>%
      mutate_at(vars(surrogates),mean_center) %>%
      ungroup()

  }

  dt_after <- out_ruc_obj$data_reg %>%
    mutate_at(vars(m,surrogates),asinh) %>%
    mutate_at(vars(surrogates),mean_center) %>%
    ungroup()

  ids_plot <- dt_before %>% group_by(sample) %>% sample_frac(frac) %>% pull(id)

  if (out_ruc_obj$model == "interaction"){

      fit <- out_ruc_obj$eff_coefficients %>%
        filter(marker == m, surrogate != "FALSE") %>%
        mutate(slope = eff_value) %>%
        select(marker,sample,surrogate,slope) %>%
        left_join(
          out_ruc_obj$eff_coefficients %>%
            filter(marker == m, surrogate == "FALSE") %>%
            mutate(intercept = eff_value) %>%
            select(marker,sample,intercept),
          by = c("marker","sample"))

  } else if (out_ruc_obj$model == "offset") {

    fit <- out_ruc_obj$eff_coefficients %>%
      filter(marker == m, surrogate != "FALSE") %>%
      mutate(slope = eff_value) %>%
      select(marker,surrogate,slope) %>%
      left_join(
        out_ruc_obj$eff_coefficients %>%
          filter(marker == m, surrogate == "FALSE") %>%
          mutate(intercept = eff_value) %>%
          select(marker,sample,intercept),
        by = c("marker"))

  } else {

    fit <- out_ruc_obj$eff_coefficients %>%
      filter(marker == m, surrogate != "FALSE") %>%
      mutate(slope = eff_value) %>%
      select(marker,surrogate,slope) %>%
      left_join(
        out_ruc_obj$eff_coefficients %>%
          filter(marker == m, surrogate == "FALSE") %>%
          mutate(intercept = eff_value) %>%
          select(marker,intercept),
        by = c("marker"))

  }



  if (out_ruc_obj$model %in% c("interaction", "offset")) {

    a <- dt_before %>%
      filter(id %in% ids_plot) %>%
      ungroup() %>%
      select(m,surrogates, sample) %>%
      pivot_longer(names_to = "surrogate",values_to = "surr_value",surrogates) %>%
      ggplot(aes_string(x = "surr_value", y = m, fill = "sample", color = "sample")) +
      geom_point(alpha = 0.1, size = 0.1)  +
      facet_wrap(~surrogate, scale = "free_x", nrow = 1) +
      theme(strip.background = element_rect(fill = "white")) +
      xlab("") +
      guides(x = guide_axis(angle = 45), colour = guide_legend(override.aes = list(size=3, alpha = 1))) +
      ggtitle("Before") +
      theme(axis.line = element_line()) +
      coord_axes_inside(labels_inside = TRUE) +
      scale_color_viridis_d()+
      geom_abline(mapping = aes(intercept = intercept, slope = slope, col = sample), data = fit)

    b <- dt_after %>%
      filter(id %in% ids_plot) %>%
      select(m,surrogates, sample) %>%
      pivot_longer(names_to = "surrogate",values_to = "surr_value",surrogates) %>%
      ggplot(aes_string(x = "surr_value", y = m, fill = "sample", color = "sample")) +
      geom_point(alpha = 0.1, size = 0.1)  +
      facet_wrap(~surrogate, scale = "free_x", nrow = 1) +
      theme(strip.background = element_rect(fill = "white")) +
      xlab("surrogates") +
      geom_vline(xintercept = 0) +
      geom_hline(yintercept = 0) +
      guides(x = guide_axis(angle = 45), colour = guide_legend(override.aes = list(size=3, alpha = 1))) +
      ggtitle("After") +
      theme(axis.line = element_line()) +
      coord_axes_inside(labels_inside = TRUE) +
      scale_color_viridis_d()



  } else {

    a <- dt_before %>%
      filter(id %in% ids_plot) %>%
      ungroup() %>%
      select(m,surrogates, sample) %>%
      pivot_longer(names_to = "surrogate",values_to = "surr_value",surrogates) %>%
      ggplot(aes_string(x = "surr_value", y = m)) +
      geom_point(alpha = 0.1, size = 0.1)  +
      facet_wrap(~surrogate, scale = "free_x", nrow = 1) +
      theme(strip.background = element_rect(fill = "white")) +
      xlab("") +
      guides(x = guide_axis(angle = 45), colour = guide_legend(override.aes = list(size=3, alpha = 1))) +
      ggtitle("Before") +
      theme(axis.line = element_line()) +
      coord_axes_inside(labels_inside = TRUE) +
      scale_color_viridis_d()+
      geom_abline(mapping = aes(intercept = intercept, slope = slope), data = fit)



    b <- dt_after %>%
      filter(id %in% ids_plot) %>%
      select(m,surrogates, sample) %>%
      pivot_longer(names_to = "surrogate",values_to = "surr_value",surrogates) %>%
      ggplot(aes_string(x = "surr_value", y = m)) +
      geom_point(alpha = 0.1, size = 0.1)  +
      facet_wrap(~surrogate, scale = "free_x", nrow = 1) +
      theme(strip.background = element_rect(fill = "white")) +
      xlab("surrogates") +
      geom_vline(xintercept = 0) +
      geom_hline(yintercept = 0) +
      guides(x = guide_axis(angle = 45), colour = guide_legend(override.aes = list(size=3, alpha = 1))) +
      ggtitle("After") +
      theme(axis.line = element_line()) +
      coord_axes_inside(labels_inside = TRUE) +
      scale_color_viridis_d()


  }


    p <-  ggarrange(plotlist = list(a,b), nrow = 2, common.legend = TRUE)

  return(p)

}


# scatter plots with regression sample
#' @export
plot_residuals <- function(data, markers, surrogates, out_ruc_obj, frac){

  data$sample <- data[[out_ruc_obj$col_name_sample]] #create column "sample" for simplicity

  dt <- data %>%
    mutate_at(vars(markers),asinh) %>%
    select(id,markers,sample) %>%
    pivot_longer(names_to = "marker", values_to = "original_values", markers) %>%
    left_join(out_ruc_obj$model_residuals %>% pivot_longer(names_to = "marker", values_to = "residuals", markers), by = c("id","marker"))


  ids_plot <- data %>% group_by(sample) %>%
    sample_frac(frac) %>% pull(id)

  p <-  dt %>% filter(id %in% ids_plot) %>%
    ggplot(aes(x = original_values, y = residuals, color = sample)) +
    geom_point(alpha = 0.5, size = 0.1)  +
    facet_wrap(~marker) +
    geom_hline(yintercept = 0) +
    guides(colour = guide_legend(override.aes = list(size=3, alpha = 1)))

  return(p)

}


# model coefficients


plot_model_coeffs <- function(out_ruc_obj, markers){

  out_ruc_obj$data_reg$sample <- out_ruc_obj$data_reg[[out_ruc_obj$col_name_sample]]

  if (out_ruc_obj$model == "simple") {
    a <- out_ruc_obj$eff_coefficients %>%
      filter(surrogate == "FALSE", marker %in% markers) %>%
      ggplot(aes(y = marker, x = eff_value)) +
      geom_col() +
      ggtitle("Intercept")
  } else {

    a <- out_ruc_obj$eff_coefficients %>%
      filter(surrogate == "FALSE", marker %in% markers) %>%
      ggplot(aes(y = sample, x = eff_value, fill = sample)) +
      geom_col() +
      facet_wrap(~marker, scale = "free") +
      ggtitle("Intercept")

  }



  if(out_ruc_obj$model == "interaction"){

    b <- out_ruc_obj$stand_slopes %>%
      filter(marker %in% markers) %>%
      ggplot(aes(y = sample, fill = surrogate, x = stand_value)) +
      geom_col(position = "dodge") +
      facet_wrap(~marker) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      ggtitle("Standardized slope")


  }else{

    b <- out_ruc_obj$stand_slopes %>%
      filter(marker %in% markers) %>%
      ggplot(aes(y = surrogate, fill = surrogate, x = stand_value)) +
      geom_col(position = "dodge") +
      facet_wrap(~marker, nrow = 3) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      ggtitle("Standardized slope")

  }


  p <-  ggarrange(plotlist = list(a,b), nrow = 2)

  return(p)
}


# adjusted R2: goodness of fit
#' @export
plot_adjR2 <- function(out_ruc_obj, markers){

  p <- out_ruc_obj$adjr2 %>%
    filter(marker %in% markers) %>%
    ggplot(aes(y = reorder(marker,adj_r_squared), x = adj_r_squared)) +
    geom_col()

  return(p)

}



# distributions
#' @export
plot_distributions <- function(data, markers, out_ruc_obj, per_sample){

  data$sample <- data[[out_ruc_obj$col_name_sample]] #create column "sample" for simplicity
  data$data_type<- "before"
  out_ruc_obj$data_reg$sample <- out_ruc_obj$data_reg[[out_ruc_obj$col_name_sample]]
  out_ruc_obj$data_reg$data_type <- "after"

  dt <- data %>% select(markers,out_ruc_obj$surrogates,sample,data_type) %>%
    rbind(out_ruc_obj$data_reg %>% select(markers,out_ruc_obj$surrogates,sample,data_type)) %>%
    mutate(data_type = factor(data_type, levels = c("before","after")))

  if (per_sample == TRUE) {

    p <- dt %>%
      mutate_at(vars(markers),asinh) %>%
      pivot_longer(names_to = "markers", markers) %>%
      ggplot(aes(x = value, y = sample, color = data_type)) +
      geom_density_ridges(alpha = 0) +
      scale_color_manual(values = c("black","red")) +
      facet_wrap(~markers, scale = "free")

  } else {

    p <- dt %>%
      mutate_at(vars(markers),asinh) %>%
      pivot_longer(names_to = "markers", markers) %>%
      ggplot(aes(x = value, color = data_type)) +
      geom_density(alpha = 0) +
      scale_color_manual(values = c("black","red")) +
      facet_wrap(~markers, scale = "free")

  }

  return(p)


}

# mean values
#' @export
plot_mean_marker <- function(data,markers,out_ruc_obj){


  data$sample <- data[[out_ruc_obj$col_name_sample]] #create column "sample" for simplicity
  data$data <- "before"
  out_ruc_obj$data_reg$sample <- out_ruc_obj$data_reg[[out_ruc_obj$col_name_sample]]
  out_ruc_obj$data_reg$data <- "after"

  dt <- data %>% rbind(out_ruc_obj$data_reg) %>% mutate(data = factor(data, levels = c("before","after")))

  p <- dt %>%
    mutate_at(vars(markers),asinh) %>%
    group_by(sample,data) %>%
    summarize_at(vars(markers),mean) %>%
    pivot_longer(names_to = "marker", values_to = "mean_value", markers) %>%
    mutate(marker = factor(marker, levels = markers)) %>%
    ggplot(aes(y = sample, x = mean_value, fill  = data)) +
    geom_col(position = "dodge") +
    facet_wrap(~marker, scale = "free") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    scale_fill_manual(values = c("black","red"))



  return(p)

}


