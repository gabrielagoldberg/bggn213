Untitled
================
Gabriela Goldberg
4/24/2019

## More on function writing\!

First we will revisit our function from last class.

``` r
source("http://tinyurl.com/rescale-R")
```

Test the **rescale()**
    function

``` r
rescale(1:10)
```

    ##  [1] 0.0000000 0.1111111 0.2222222 0.3333333 0.4444444 0.5555556 0.6666667
    ##  [8] 0.7777778 0.8888889 1.0000000

``` r
#rescale( c(1:10, "string"))
```

Let’s fix this with another line of code.

``` r
#rescale2( c(1:10, "string"))
```

Now the error will tell us what is wrong.

## Function practice\!

Write a function to identify NA elements in two vectors.

Start with a simple example input where I know what the answer should
be.

``` r
x <- c( 1, 2, NA, 3, NA)
y <- c(NA, 3, NA, 3, 4)
```

``` r
is.na(x)
```

    ## [1] FALSE FALSE  TRUE FALSE  TRUE

``` r
is.na(y)
```

    ## [1]  TRUE FALSE  TRUE FALSE FALSE

I am looking for the positions where it is TRUE in both vectors…

``` r
is.na(x) & is.na(y)
```

    ## [1] FALSE FALSE  TRUE FALSE FALSE

Take the sum to find how many TRUEs

``` r
sum(is.na(x) & is.na(y))
```

    ## [1] 1

This is my working snippet of code that I can use as the body of my
first function\!

``` r
both_na <- function(x,y) {
  sum(is.na(x) & is.na(y))
}
```

TEST\!\!

``` r
both_na(x, y)
```

    ## [1] 1

``` r
both_na( c(NA, NA, NA, 5, 6, 7, 9, NA), c(NA, 5, 6, 7, NA, NA, NA))
```

    ## Warning in is.na(x) & is.na(y): longer object length is not a multiple of
    ## shorter object length

    ## [1] 2

Check to see if the lengths of our inputs are equal

``` r
x <- c(NA, NA, NA)
y <- c(NA, NA, NA, NA, NA, NA)
length(x) != length(y)
```

    ## [1] TRUE

Function **both\_na2()** will tell us to make the vectors the same
length.

``` r
#both_na2(x, y)
```

Function **both\_na3()** makes the code better and gives us more
information.

``` r
x <- c(NA, NA, NA, 1, 2, 3)
y <- c(NA, 3, 4, 5, NA, NA)
both_na3(x, y)
```

    ## Found 1 NA's at position(s):1

    ## $number
    ## [1] 1
    ## 
    ## $which
    ## [1] 1

## Your turn\!

Write a function **grade()** to determine an overall grade from a vector
of student homework assignment scores dropping the lowest single
alignment score.

Student 1

``` r
s1 <- c(100, 100, 100, 100, 100, 100, 100, 100, 90)
```

Student 2

``` r
s2 <- c(100, NA, 90, 90, 90, 90, 97, 80)
```

Let’s first find the minimum of each student’s scores.

``` r
min(s1)
```

    ## [1] 90

``` r
min(s2)
```

    ## [1] NA

We can then add up all the scores and subtract them from the total and
divide them by the length of the vector minus 1.

``` r
(sum(s1) - min(s1)) / (length(s1) - 1)
```

    ## [1] 100

``` r
(sum(s2) - min(s2)) / (length(s2) - 1)
```

    ## [1] NA

Let’s write the function\!

``` r
grade <- function(x) {
  (sum(x) - min(x)) / (length(x) - 1)
}
```

Let’s test it on student 1

``` r
grade(s1)
```

    ## [1] 100

Does it work on student 2?

``` r
grade(s2)
```

    ## [1] NA

Let’s try to fix this with a condition statement…

``` r
grade <- function(x) {
  (sum(x, na.rm = TRUE) - min(x, na.rm = TRUE)) / (length(x) - 1)
}
```

Does it work on student 2?

``` r
grade(s2)
```

    ## [1] 79.57143

Great\! Now let’s try to do this on a bigger file\!

``` r
url <- "http://tinyurl.com/gradeinput"
```

``` r
classgrades <- read.csv(url, row.names = 1)
```

This csv file has 20 rows (students) and 5 columns
    (grades).

``` r
grade(classgrades[1,])
```

    ## [1] 91.75

``` r
finalgrades <- apply(classgrades, 1, grade)
```

``` r
sort(finalgrades, decreasing = TRUE)
```

    ##  student-7  student-8 student-13  student-1 student-12 student-16 
    ##      94.00      93.75      92.25      91.75      91.75      89.50 
    ##  student-6  student-5 student-17  student-9 student-14 student-11 
    ##      89.00      88.25      88.00      87.75      87.75      86.00 
    ##  student-3 student-19 student-20  student-2 student-18  student-4 
    ##      84.25      82.75      82.75      82.50      72.75      66.00 
    ## student-15 student-10 
    ##      62.50      61.00

## Now let’s do a real-world example

Find common genes in two data sets and return their associated data
(from each data set)

Start with a simple version of the problem

``` r
df1 <- data.frame(IDs=c("gene1", "gene2", "gene3"),
 exp=c(2,1,1),
 stringsAsFactors=FALSE)
df2 <- data.frame(IDs=c("gene2", "gene4", "gene3", "gene5"),
 exp=c(-2, NA, 1, 2),
 stringsAsFactors=FALSE)
```

simplify further to single vectors

``` r
x <- df1$IDs
y <- df2$IDs
```

find the intersection of both vectors

``` r
intersect(x, y)
```

    ## [1] "gene2" "gene3"

That function doesn’t give us the option to give us all the information
along with the gene.

``` r
x
```

    ## [1] "gene1" "gene2" "gene3"

``` r
y
```

    ## [1] "gene2" "gene4" "gene3" "gene5"

``` r
x %in% y
```

    ## [1] FALSE  TRUE  TRUE

``` r
x[x %in% y]
```

    ## [1] "gene2" "gene3"

Let’s bind these columns together

``` r
gene_intersect <- function(x, y) {
  cbind( x[x %in% y], y[y %in% x])
}
```
