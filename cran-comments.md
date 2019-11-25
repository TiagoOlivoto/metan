# Resubmission
This is a resubmission for the new package `metan`. Now in version 1.1.1, I have incorporated the feedback kindly provided by Martina Schmirl on November 22th, 2019. See responses inline below. 

> 
A CRAN package check should not take longer than 10 minutes. Please considerably reduce the check time of your package to stay below the threshold of 10 minutes. This can for example be achieved by:
> 
- Omitting the less important and lengthy tests by only running them conditionally if some environment variable is set that you only define on your machine.
> 
- Using precomputed results for the computational intensive parts in vignettes.
> 
Ensure that even after this change the package still includes illustrative examples of the use of the user-facing functions. Please aim at trying to test and also illustrate the use of as much functionality of the package as possible while making considerate use of the computational resources!
> 
Please make sure that you do not change the user's options, par or working directory. If you really have to do so, please ensure with an *immediate* call of on.exit() that the settings are reset when the function is exited. e.g.:
...
oldoptions <- options(SCIPEN = 100)   # code line i
on.exit(options(oldoptions))          # code line i+1
...
e.g.: print.ecovalence.R
If you're not familiar with the function, please check ?on.exit. This function makes it possible to restore options before exiting a function even if the function breaks. Therefore it needs to be called immediately after the option change within a function.
> 
Please fix and resubmit.



* Done. I have reviewed the code and used `\donttest{}` to omitt lengthy tests and used precomputed results in the vignettes. The R CMD check results was now 3m 44.8s and the check time <https://win-builder.r-project.org/2744135XttrV> 587 seconds (~9.8 min)

* I've added `on.exit()` to ensure that the settings are reset when the function is exited.


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

