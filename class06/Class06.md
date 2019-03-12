class06
================
Brandon Gonzalez
January 24, 2019

Class 6 work - reading files again
----------------------------------

We will be opening a couple of different files. We will use the fucntion **read.table()** and friends of it to read some example flat files.

``` r
table1 <- read.csv(file="./data_files/test1.txt", header = TRUE)
table1
```

    ##   Col1 Col2 Col3
    ## 1    1    2    3
    ## 2    4    5    6
    ## 3    7    8    9
    ## 4    a    b    c

``` r
table2 <- read.table(file="./data_files/test2.txt",header = TRUE, sep="$")
table2
```

    ##   Col1 Col2 Col3
    ## 1    1    2    3
    ## 2    4    5    6
    ## 3    7    8    9
    ## 4    a    b    c

``` r
table3 <- read.table(file="./data_files/test3.txt", header = FALSE)
```

    ## Warning in read.table(file = "./data_files/test3.txt", header = FALSE):
    ## incomplete final line found by readTableHeader on './data_files/test3.txt'

``` r
table3
```

    ##   V1 V2 V3
    ## 1  1  6  a
    ## 2  2  7  b
    ## 3  3  8  c
    ## 4  4  9  d
    ## 5  5 10  e

``` r
ncol(table3)
```

    ## [1] 3

Section 2. R Functions
----------------------

Basic function structure: name.of.function &lt;- function( arg1, arg2,...){ expr return(something) }

``` r
add <- function(x, y=1){
  #y default value is 1, optional arg
  x + y
}

x <- 1
#using this function:
add(x, y=100)
```

    ## [1] 101

``` r
add(x)
```

    ## [1] 2

``` r
add(c(1,2,3))
```

    ## [1] 2 3 4

``` r
add(c(1,2,3),4)
```

    ## [1] 5 6 7

``` r
rescale <-function(x, na.rm=TRUE){
  rng <- range(x, na.rm = na.rm)
  x <- (x - rng[1]) / (rng[2] - rng[1])
  return(x)
}

x = 1:10
rescale(x)
```

    ##  [1] 0.0000000 0.1111111 0.2222222 0.3333333 0.4444444 0.5555556 0.6666667
    ##  [8] 0.7777778 0.8888889 1.0000000

``` r
rescale3 <- function(x, na.rm=TRUE, plot=FALSE) {
 rng <-range(x, na.rm=na.rm)
 print("Hello")
 answer <- (x - rng[1]) / (rng[2] - rng[1])
 print("is it me you are looking for?")
 if(plot) {
 plot(answer, typ="b", lwd=4)
 }
 print("I can see it in ...")
 return(answer)
}
```

``` r
rescale3(c(1:15,NA,10))
```

    ## [1] "Hello"
    ## [1] "is it me you are looking for?"
    ## [1] "I can see it in ..."

    ##  [1] 0.00000000 0.07142857 0.14285714 0.21428571 0.28571429 0.35714286
    ##  [7] 0.42857143 0.50000000 0.57142857 0.64285714 0.71428571 0.78571429
    ## [13] 0.85714286 0.92857143 1.00000000         NA 0.64285714
