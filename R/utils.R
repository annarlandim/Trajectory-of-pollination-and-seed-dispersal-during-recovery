hulls_hbt <- function(scores){
  
  hulls <- scores %>%
    group_by(Treatment3) %>%
    summarise(chull_indices = list(chull(RC1, RC2)))
  hull_points <- list()
  for (i in 1:nrow(hulls)) {
    indices <- c(hulls$chull_indices[[i]], hulls$chull_indices[[i]][1]) 
    points <- scores[scores$Treatment3 == hulls$Treatment3[i], ][indices, ]
    points$Treatment3 <- hulls$Treatment3[i] 
    hull_points[[i]] <- points
  }
  hull_points_df <- do.call(rbind, hull_points)
  return(hull_points_df)
}

calc_dist <- function(point1, point2) {
  sqrt(sum((point1 - point2) ^ 2))
}

originality <- function(unique_plot, scores){
  
  orig <- unique_plot %>%
    left_join(distinct(scores, !!sym(names(scores)[1]), !!sym(names(scores)[2]), .keep_all = TRUE), 
              by = c(names(scores)[1], names(scores)[2]), relationship = "many-to-many")
  
  centroids <- scores %>%
    group_by(Treatment3) %>%
    summarise(RC1 = mean(RC1), RC2 = mean(RC2))  
  
  orig$orig <- mapply(function(x, y, trat) {
    center <- centroids[centroids$Treatment3 == trat, c("RC1", "RC2")]
    calc_dist(c(x, y), center)
  }, orig$RC1, orig$RC2, orig$Treatment3)
  
  return(orig %>% mutate(count.t3 = NULL, count.plants.t3 = NULL))
}
