predict_lin <- function(object, g.new, e.new = NULL, y.new,
                        robust = TRUE, quant, ...) {
  
  if (robust) {
    
    if (is.null(quant)) {
      stop("quant must be specified when robust = TRUE.")
    }
    
    out <- predict_lin_robust(
      obj   = object,
      g.new = g.new,
      e.new = e.new,
      y.new = y.new,
      quant = quant,
      ...
    )
    
  } else {
    
    if (!is.null(quant)) {
      stop("quant must be NULL when robust = FALSE.")
    }
    out <- predict_lin_nonrobust(
      obj   = object,
      g.new = g.new,
      e.new = e.new,
      y.new = y.new,
      ...
    )
    
  }
  
  return(out)
}