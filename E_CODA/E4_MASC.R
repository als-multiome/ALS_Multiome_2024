MASC <- function(
    dataset, 
    cluster, 
    contrast, 
    random_effects = NULL, 
    fixed_effects = NULL,
    verbose = FALSE, 
    return_models = TRUE
){
 
  # Check inputs
  if (is.factor(dataset[[contrast]]) == FALSE) {
    message("Specified contrast term is not coded as a factor in dataset")
  }
  
  
  # Generate design matrix from cluster assignments
  cluster <- as.character(dataset[[cluster]])
  designmat <- model.matrix(~ cluster + 0, data.frame(cluster = cluster))
  dataset <- cbind(designmat, dataset)
  
  
  # Create output list to hold results
  res <- vector(mode = "list", length = length(unique(cluster)))
  names(res) <- attributes(designmat)$dimnames[[2]]
  
  # Create model formulas
  if (!is.null(fixed_effects) && !is.null(random_effects)) {
    model_rhs <- paste0(
      c(
        paste0(fixed_effects, collapse = " + "),
        paste0("(1|", random_effects, ")", collapse = " + ")
      ),
      collapse = " + "
    )
    if (verbose == TRUE) {
      message(paste("Using null model:", "cluster ~", model_rhs))
    }
  } else {
      stop("No random or fixed effects specified")
    }
  
  # Initialize list to store model objects for each cluster
  cluster_models <- vector(mode = "list",
                           length = length(attributes(designmat)$dimnames[[2]]))
  names(cluster_models) <- attributes(designmat)$dimnames[[2]]
  
  
  
  # Run nested mixed-effects models for each cluster
  for (i in seq_along(attributes(designmat)$dimnames[[2]])) {
    test_cluster <- attributes(designmat)$dimnames[[2]][i]
    if (verbose == TRUE) {
      message(paste("Creating logistic mixed models for", test_cluster))
    }
    null_fm <- as.formula(
      paste0(
        c(
          paste0(test_cluster, " ~ 1 + "),
          model_rhs
        ), 
        collapse = ""
      )
    )
    full_fm <- as.formula(
      paste0(
        c(
          paste0(test_cluster, " ~ ", contrast, " + "),
          model_rhs
        ), 
        collapse = ""
      )
    )
    
    # Run null and full mixed-effects models
    null_model <- lme4::glmer(
      formula = null_fm, 
      data = dataset,
      family = binomial, 
      nAGQ = 1, 
      verbose = 0,
      control = glmerControl(optimizer = "bobyqa")
    )
    
    full_model <- lme4::glmer(
      formula = full_fm, 
      data = dataset,
      family = binomial, 
      nAGQ = 1, 
      verbose = 0,
      control = glmerControl(optimizer = "bobyqa")
    )
    
    model_lrt <- anova(null_model, full_model)
    
    # calculate confidence intervals for contrast term beta
    contrast_lvl2 <- paste0(
      contrast, 
      levels(dataset[[contrast]])[2]
    )
    contrast_ci <- confint.merMod(
      full_model, 
      method = "Wald",
      parm = contrast_lvl2
    )
    
    # Save model objects to list
    cluster_models[[i]]$null_model <- null_model
    cluster_models[[i]]$full_model <- full_model
    cluster_models[[i]]$model_lrt <- model_lrt
    cluster_models[[i]]$confint <- contrast_ci
  }
  
  # Organize results into output dataframe
  output <- data.frame(
    cluster = attributes(designmat)$dimnames[[2]],
    size = colSums(designmat)
  )
  
  output$model.pvalue <- sapply(
    cluster_models, 
    function(x) x$model_lrt[["Pr(>Chisq)"]][2]
  )
  
  output[[paste(contrast_lvl2, "OR", sep = ".")]] <- sapply(
    cluster_models, 
    function(x){
      exp(fixef(x$full)[[contrast_lvl2]])
    }
  )
  
  output[[paste(contrast_lvl2, "OR", "95pct.ci.lower", sep = ".")]] <- sapply(
    cluster_models, 
    function(x){
      exp(x$confint[contrast_lvl2, "2.5 %"])
    } 
  )
  
  output[[paste(contrast_lvl2, "OR", "95pct.ci.upper", sep = ".")]] <- sapply(
    cluster_models, 
    function(x){
      exp(x$confint[contrast_lvl2, "97.5 %"])
    } 
  )
  
  output$fdr <- p.adjust(
    output$model.pvalue, 
    method = "fdr"
  )
  
  if(return_models){
    return(
      list(
        output, 
        cluster_models
      )
    )
  } else {
    return(
      list(
        output
      )
    )
  }
  
  
}
  