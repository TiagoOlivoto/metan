mmat = function (data, row, col, value){
data = data.frame(R = eval(substitute(row), eval(data)),
                  C = eval(substitute(col), eval(data)),
                  V = eval(substitute(value), eval(data)))
return(data.frame(tapply(data[, 3], data[, c(1, 2)], mean)))
}
