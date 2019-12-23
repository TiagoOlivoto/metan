# Resubmission `metan`
This is a resubmission for the new package `metan`. Now in version 1.2.0, I have incorporated the feedback kindly provided by Jelena Saf on December 20, 2019. See responses inline below. 

> 
* Please omit the redundant 'Provides functions for' from your description.

> 
* You write information messages to the console that cannot be easily suppressed. It is more R like to generate objects that can be used to extract the information a user is interested in, and then print() that object. Instead of print()/cat() rather use message()/warning()  or if(verbose)cat(..) if you really have to write text to the console. (except for print() and summary() functions) F.i.: corr_plot.R

In this version I have made the following changes:
* Omitted *"Provides functions for"* from description file.
* In the function `corr_plot()` the message *"The factors ... where excluded to perform the analysis. Only numeric variables were used."* was deleted. Note that the progress bar can be supressed with the argument `progress`.
* Some new functions were also included in this version. Please, see "News" for more details.


# Test environments

- local OS X install, R 3.6.1
- Ubuntu Linux 16.04 LTS, R-release, GCC
- Rhub
   - Windows Server 2008 R2 SP1, R-devel, 32/64 bit
   - Ubuntu Linux 16.04 LTS, R-release, GCC
   - Fedora Linux, R-devel, clang, gfortran
- win-builder (devel and release)
   - One note (Maintainer: 'Tiago Olivoto <tiagoolivoto@gmail.com>')

# R CMD check results
0 errors | 0 warnings | 0 notes


Thank for your time in helping me with this submission.

Best Regards

Tiago Olivoto  

