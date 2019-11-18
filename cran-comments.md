# Resubmission
This is a resubmission for the new package `metan`. Now in version 1.1.0, I have incorporated the feedback kindly provided by Martina Schmirl on Septempber 27th, 2019. See responses inline below. 

- Please ensure that your functions do not write by default or in your examples/vignettes/tests in the user's home filespace. That is not allow by CRAN policies. Please only write/save files if the user has specified a directory. In your examples/vignettes/tests you can write to tempdir().

   * Done! I have reviewed twice and no functions write files by default. In addition, I have updated the vignette "vignettes_stability.Rmd" that was writing a file in the working directory.


After committing the changes I ran:

- `devtools::spell_check()`
- `devtools::check()`
   - 0 errors | 0 warnings | 0 notes
   
- `devtools::check_rhub()`
   - Windows Server 2008 R2 SP1, R-devel, 32/64 bit: CRAN incoming feasibility;
   - Ubuntu Linux 16.04 LTS, R-release, GCC:
   - Fedora Linux, R-devel, clang, gfortran
   - One note (Maintainer: 'Tiago Olivoto <tiagoolivoto@gmail.com>')
   
- `devtools::check_win_devel()`
   - CRAN incoming feasibility: One note (Maintainer: 'Tiago Olivoto <tiagoolivoto@gmail.com>')

- `devtools::check_win_release()`
   - CRAN incoming feasibility: One note (Maintainer: 'Tiago Olivoto <tiagoolivoto@gmail.com>')

# Test environments

- local OS X install, R 3.6.1
- Ubuntu Linux 16.04 LTS, R-release, GCC
- Rhub
   - Windows Server 2008 R2 SP1, R-devel, 32/64 bit
   - Ubuntu Linux 16.04 LTS, R-release, GCC
   - Fedora Linux, R-devel, clang, gfortran
- win-builder (devel and release)

# R CMD check results
0 errors | 0 warnings | 0 notes


Thank for your time in helping me with this submission.

Best Regards

Tiago Olivoto  

