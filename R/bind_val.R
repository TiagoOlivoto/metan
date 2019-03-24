bind_val = function(...){
  dots = list(...)
  data = do.call(rbind, lapply(dots, function(x){
    x$RMSPD
    })
  )
  return(structure(list(RMSPD = data),
                   class = "validation.AMMIF"))
}
