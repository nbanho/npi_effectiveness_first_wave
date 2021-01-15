require(tikzDevice)

# Settings for tikzdevice
if (Sys.info()['sysname'] == 'Darwin') {
  options("tikzLatex" = Sys.which("pdflatex"))
} else {
  options("tikzLatex" = 'C:/texlive/2018/bin/win32/pdflatex')
}




# Round function
round_k <- function(x, k = 2) {
  trimws(format(round(x, k), nsmall = k))
}

# Check if value is double
is.numeric.double <- function(x) {
  if (is.numeric(x)) {
    if (all(x %% 1 == 0)) {
      return(F)
    }
    else {
      return(T)
    }
  } 
  else {
    return(F)
  }
}

# Remove leading zero
rzero <- function(x) {
  x <- as.character(x)
  fc <- substr(x, 1, 1)
  if (fc == "0")
    return(substr(x, 2, nchar(x)))
  else
    return(paste0("$-$", substr(x, 3, nchar(x))))
}

# Math tex sign
tex_sign <- function(x) {
  x_new <- character(length(x))
  for (i in 1:length(x)) {
    is_negative <- ifelse(substr(x[i], 1, 1) == "-", T, F)
    if (is_negative) {
      x_new[i] <- paste0("$-$", substr(x[i], 2, nchar(x[i])))
    } else {
      x_new[i] <- x[i]
    }
  }
  return(x_new)
}