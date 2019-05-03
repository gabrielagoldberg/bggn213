# Section 1: Improving code by writing functions
  # (A. Can you improve this analysis code?
  df <- data.frame(a=1:10, b=seq(200,400,length=10),c=11:20,d=NA)
  df$a <- (df$a - min(df$a)) / (max(df$a) - min(df$a))
  df$b <- (df$b - min(df$a)) / (max(df$b) - min(df$b))
  df$c <- (df$c - min(df$c)) / (max(df$c) - min(df$c))
  df$d <- (df$d - min(df$d)) / (max(df$a) - min(df$d)) 
  
  # Simplify to work with a generic vector named "x"
  x <- (x - min(x)) / (max(x) - min(x))
  ## Note that we call the min() function twice…
  x <- (x - min(x)) / (max(x) - min(x))
  ## Note that we call the min() function twice…
  xmin <- min(x)
  x <- (x - xmin) / (max(x) - xmin)
  ## Further optimization to use the range() function…
  rng <- range(x)
  x <- (x - rng[1]) / (rng[2] - rng[1])
  ## You need a “name”, “arguments” and “body”…
  rescale <- function(x) {
    rng <-range(x)
    (x - rng[1]) / (rng[2] - rng[1])
  }
  
  # Test on something you know the answer to
  rescale(1:10)
  rescale(c(3, 5))
  rescale

  # Randomness
  rescale <- function(x, na.rm=TRUE, plot=FALSE) {
    if(na.rm) {
      rng <-range(x, na.rm=na.rm)
    } else {
      rng <-range(x)
    }
    print("Hello")
    answer <- (x - rng[1]) / (rng[2] - rng[1])
    print("is it me you are looking for?")
    if(plot) {
      plot(answer, typ="b", lwd=4)
    }
    print("I can see it in ...")
  }  
rescale(1)  
rescale(5)

