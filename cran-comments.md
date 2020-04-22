# Release summary
This is a path release (v1.5.1) that fixes the error with `donttest{}` examples in `fai_blup()` function.  I was confusing about the origin of the error because in the last submission (v1.5.0) when I ran `R CMD check --as-cran --run-donttest metan_1.5.0.tar.gz`, no error message was generated. Please, see the results below.

```r{}
* using log directory 'D:/Desktop/metan/metan.Rcheck'
* using R version 3.6.2 (2019-12-12)
* using platform: x86_64-w64-mingw32 (64-bit)
* using session charset: ISO8859-1
* using options '--run-donttest --as-cran'
* checking for file 'metan/DESCRIPTION' ... OK
* checking extension type ... Package
* this is package 'metan' version '1.5.0'
* package encoding: UTF-8
* checking CRAN incoming feasibility ... WARNING
Maintainer: 'Tiago Olivoto <tiagoolivoto@gmail.com>'

Insufficient package version (submitted: 1.5.0, existing: 1.5.0)

Days since last update: 1
* checking package namespace information ... OK
* checking package dependencies ... OK
* checking if this is a source package ... OK
* checking if there is a namespace ... OK
* checking for .dll and .exe files ... OK
* checking for hidden files and directories ... OK
* checking for portable file names ... OK
* checking whether package 'metan' can be installed ... OK
* checking package directory ... OK
* checking for future file timestamps ... OK
* checking 'build' directory ... OK
* checking DESCRIPTION meta-information ... OK
* checking top-level files ... OK
* checking for left-over files ... OK
* checking index information ... OK
* checking package subdirectories ... OK
* checking R files for non-ASCII characters ... OK
* checking R files for syntax errors ... OK
* checking whether the package can be loaded ... OK
* checking whether the package can be loaded with stated dependencies ... OK
* checking whether the package can be unloaded cleanly ... OK
* checking whether the namespace can be loaded with stated dependencies ... OK
* checking whether the namespace can be unloaded cleanly ... OK
* checking loading without being on the library search path ... OK
* checking use of S3 registration ... OK
* checking dependencies in R code ... OK
* checking S3 generic/method consistency ... OK
* checking replacement functions ... OK
* checking foreign function calls ... OK
* checking R code for possible problems ... OK
* checking Rd files ... OK
* checking Rd metadata ... OK
* checking Rd line widths ... OK
* checking Rd cross-references ... OK
* checking for missing documentation entries ... OK
* checking for code/documentation mismatches ... OK
* checking Rd \usage sections ... OK
* checking Rd contents ... OK
* checking for unstated dependencies in examples ... OK
* checking contents of 'data' directory ... OK
* checking data for non-ASCII characters ... OK
* checking data for ASCII and uncompressed saves ... OK
* checking installed files from 'inst/doc' ... OK
* checking files in 'vignettes' ... OK
* checking examples ... NOTE
Examples with CPU or elapsed time > 5s
                user system elapsed
corr_plot      44.14   0.08   45.48
cv_ammi        12.11   0.00   12.35
cv_ammif       11.59   0.01   11.81
get_model_data 10.84   0.05   11.31
waasb          10.02   0.06   10.95
desc_stat       7.82   0.02    7.97
clustering      6.99   0.01    7.25
waas            5.86   0.02    6.34
* checking for unstated dependencies in vignettes ... OK
* checking package vignettes in 'inst/doc' ... OK
* checking re-building of vignette outputs ... OK
* checking PDF version of manual ... OK
* checking for detritus in the temp directory ... OK
* DONE
Status: 1 WARNING, 1 NOTE

```
I just realized that the error in `fai_blup()` example occurred because I used a data set with a non-significant genotype effect to compute an index that requires a significant genotype effect for all traits. Now, the correct data set is used in that example. I'd like to thank Brian D. Ripley for his kind e-mail alerting me about this error.


# Test environments
- local OS X install, R 3.6.3
- Rhub
   - Windows Server 2008 R2 SP1, R-devel, 32/64 bit
   - Ubuntu Linux 16.04 LTS, R-release, GCC
   - Fedora Linux, R-devel, clang, gfortran
- win-builder (devel and release)


# R CMD check results (run with --run-donttest --as-cran)
0 errors | 0 warnings | 0 notes






