# thrreg

Fork of thrreg package originally by Marcel Kremer, from mlkremer/thrreg, itself a port of Bruce Hansen's code from his website.

This is just aimed at adding some slight modifications so that its interface is more in line with other R packages.

## Overview
thrreg estimates the threshold regression model by Hansen (2000). 
It provides three functions:
* `thr_test_hom()` tests for a threshold in a linear regression under homoskedasticity.
* `thr_test_het()` tests for a threshold in a linear regression under heteroskedasticity.
* `thr_est()` estimates the threshold and the regression parameters of the threshold model.

## Installation

Install development version from GitHub:
```r
# install.packages("devtools")
devtools::install_github("zeyuz35/thrreg")
```

## Usage
This example performs the empirical work reported in Hansen (2000) and prints 
the LaTeX output to `output/output.tex`:
```r
library(thrreg)

## Load data
data <- dur_john

## Define variable names for output (optional)
var_names <- c("GDP\\_Gwth", "Log\\_GDP\\_1960", "Log\\_Inv/GDP",
               "Log\\_Pop\\_Gwth", "Log\\_School", "GDP\\_1960", "Literacy")
xi <- 2:5
h <- 1


## First Level
sink("output/output.tex")
cat("\\section{First Sample Split}", "\n\n")

## Test for a Sample Split
test_1_gdp <- thr_test_het(data, "GDPGwth", xi, "GDP1960", var_names)
test_1_lit <- thr_test_het(data, "GDPGwth", xi, "Literacy", var_names)

## Estimate Sample Split
qhat_1 <- thr_est(data, "GDPGwth", xi, "GDP1960", h, test_1_gdp$p_value,
                  var_names, digits.thr = 0, header = var_names[1])
sink()


## Second Level
dat_2l <- subset(data, GDP1960 <= qhat_1)
dat_2u <- subset(data, GDP1960 > qhat_1)

sink("output/output.tex", append = T)
cat(paste0("\\section{Second Sample Split: Subsample, Incomes below ", qhat_1,
           "}", "\n\n"))

## Test for a Sample Split
test_2l_gdp <- thr_test_het(dat_2l, "GDPGwth", xi, "GDP1960", var_names)
test_2l_lit <- thr_test_het(dat_2l, "GDPGwth", xi, "Literacy", var_names)

cat(paste0("\\section{Second Sample Split: Subsample, Incomes above ", qhat_1,
           "}", "\n\n"))

## Test for a Sample Split
test_2u_gdp <- thr_test_het(dat_2u, "GDPGwth", xi, "GDP1960", var_names)
test_2u_lit <- thr_test_het(dat_2u, "GDPGwth", xi, "Literacy", var_names)

## Estimate Sample Split
qhat_2 <- thr_est(dat_2u, "GDPGwth", xi, "Literacy", h, test_2u_lit$p_value,
                  var_names, digits.thr = 0, header = var_names[1])
sink()


## Third Level
dat_3l <- subset(dat_2u, Literacy <= qhat_2)
dat_3u <- subset(dat_2u, Literacy > qhat_2)

sink("output/output.tex", append = T)
cat(paste0("\\section{Third Sample Split: Subsample, Incomes above ", qhat_1,
           ", Literacy below ", qhat_2, "}", "\n\n"))

## Test for a Sample Split
test_3l_gdp <- thr_test_het(dat_3l, "GDPGwth", xi, "GDP1960")
test_3l_lit <- thr_test_het(dat_3l, "GDPGwth", xi, "Literacy")

cat(paste0("\\section{Third Sample Split: Subsample, Incomes above ", qhat_1,
           ", Literacy above ", qhat_2, "}", "\n\n"))

## Test for a Sample Split
test_3u_gdp <- thr_test_het(dat_3u, "GDPGwth", xi, "GDP1960")
test_3u_lit <- thr_test_het(dat_3u, "GDPGwth", xi, "Literacy")

sink()
```

Compile LaTeX output in directory `output/` via:
```latex
\documentclass[a4paper]{scrartcl} %twocolumn %scrartcl, scrreprt, scrbook
\usepackage[margin=1in]{geometry}
\usepackage{hyperref}
\usepackage{booktabs}
\usepackage[table]{xcolor}
\usepackage{siunitx}

\setlength{\parindent}{0in}

\newcommand{\RedA}{red}
\newcommand{\RedB}{red!60}
\newcommand{\RedC}{red!30}
\newcommand{\BlueA}{blue!75}
\newcommand{\BlueB}{blue!55}
\newcommand{\BlueC}{blue!30}

\begin{document}
\input{output.tex}
\end{document}
```

## References
thrreg is based on the 
[source code](https://www.ssc.wisc.edu/~bhansen/progs/ecnmt_00.html) 
provided by Bruce E. Hansen.

Hansen, B. E. (2000). Sample splitting and threshold estimation.
*Econometrica*, 68(3):575--603. 
https://onlinelibrary.wiley.com/doi/10.1111/1468-0262.00124. 
Open access: https://www.ssc.wisc.edu/~bhansen/papers/ecnmt_00.pdf.
