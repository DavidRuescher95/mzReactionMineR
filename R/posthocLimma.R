#' posthocLimma
#'
#' A wrapper function that performs an ANOVA test on a summarized experiment using limma.
#'
#' @param object A summarized experiment object.
#' @param assay Character. What is the name of the assay the ANOVA should be executed on.
#' @param blocking_variables Character vector. The names of the columns in colData(object) that should be used as blocking factors
#' @param contrast_variable Character. The name of the column in colData(object) that should be used as test variables.
#' @param controls Character vector. The names of the levels in the contrast_variable that should be used as controls.
#' @param padj_method Character. The method used to adjust the p-values. Default is "fdr".
#' @return Returns a data.frame that contains the ANOVA table.
#' @examples
#' posthocLimma(object = se, assay = "vsn", blocking_variables = c("batch"), test_variable = "condition", controls = c("control1", "control2"))
#' @export

posthocLimma <- function(object = NULL, 
                       assay = NULL, 
                       blocking_variables = NULL, 
                       contrast_variable = NULL,
                       controls = NULL,
                       padj_method = "fdr") {
  
  sample_data <- as.data.frame(colData(object))
  
  formula_string <- paste(
    " ~ 0 +",
    paste(contrast_variable, collapse = " + ")
  )
  
  X <- model.matrix(
    data = sample_data,
    as.formula(formula_string)
  )
  
  for(name in contrast_variable){
    colnames(X) <- gsub(name, "", colnames(X))
  }
  
  if(!is.null(blocking_variables)) {
    
    Z <- list()
    
    for(block_name in blocking_variables) {
      
      block <- as.factor(sample_data[,block_name])
      
      contrasts(block) <- contr.sum(levels(block))
      
      Z[[block_name]] <- model.matrix(~block)[,-1,drop=FALSE]
      
      colnames(Z[[block_name]]) <- levels(block)[-1]
      
    }
    
    do.call(cbind, Z) -> Z
    
    X <- cbind(Z, X)
    
  }
  
  treatments <- unique(sample_data[,contrast_variable][!sample_data[,contrast_variable] %in% controls])
  
  contrasts = expand.grid(
    controls,
    treatments
  )
  contrasts <- contrasts %>%
    mutate(
      contrast = paste0(Var2, "-", Var1)
    )
  contrasts <- contrasts$contrast
  
  # generate contrast matrix. Limma function
  
  contrast_matrix <- limma::makeContrasts(
    contrasts = contrasts,
    levels = colnames(X)
  ) 
  
  fit <- limma::lmFit(
    assays(object)[[assay]] %>% replace(is.na(.),0),
    design = X
  )
  
  fit_contrast <- limma::contrasts.fit(fit, contrast = contrast_matrix)
  
  fit2 <- limma::eBayes(fit_contrast)
  
  pairwise_res <- list()
  
  for(n in seq_along(contrasts)) {
    
    pairwise_res[[contrasts[n]]] <- data.frame(
      as.data.frame(rowData(object)),
      limma::topTable(
        fit2,
        n=Inf, 
        sort.by="none",
        adjust.method = padj_method,,
        coef = n
      )
    ) 
    
  }
  
  pairwise_res <- pairwise_res %>%
    bind_rows(.id = "contrast")
  
  return(pairwise_res)
  
}