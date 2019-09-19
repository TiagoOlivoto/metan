# Resubmission
This is a resubmission. In this version, I have made the following changes based on the Uwe Ligges' comments. 

- Fixed invalid link in README.md
- Updated DESCRIPTION file
- updated cran-comments.md

After committing the changes I ran:

- `devtools::spell_check()`
- `devtools::check()`
   - 0 errors | 0 warnings | 0 notes
   
- `devtools::check_rhub()`
   - Build ID (metan_1.0.0.tar.gz-844323c730174fbc93ac7742490cd86b): CRAN incoming feasibility: One note (Maintainer: 'Tiago Olivoto <tiagoolivoto@gmail.com>')
   - Build ID (metan_1.0.0.tar.gz-f05d6122c07b42448691b3c4c2d7205a): metan 1.0.0: OK.
   
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


# Other questions

I have a vignette (vignettes/vignettes_metan.Rmd) that I would like to include when installing package but when I run `devtools::check()` the files in inst/doc are deleted. How should I proceed in this case?

Thank for your time in helping me with this submission.

Best Regards

Tiago Olivoto

