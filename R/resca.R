resca = function(values,
                 new_min = 0,
                 new_max = 100){
new_v = function(v){
(new_max - new_min) / (max(values) - min(values)) * (v - max(values)) + new_max
}
return(sapply(values, new_v))
}
