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
   - CRAN incoming feasibility: One note (Maintainer: 'Tiago Olivoto <tiagoolivoto@gmail.com>')
   
   - After twenty minutes I've received other email with an error (metan_1.0.0.tar.gz-f4b74c8beceb45098d5ae31e899dbeee) saying: Packages suggested but not available: 'kableExtra', 'roxygen2'. This sounds strange because I have these two packages installed and updated. What would a possible reason for this error?
   - A few minutes latter other email (metan_1.0.0.tar.gz-639d6271c76b43c8ba4ff2b56104ed95): metan 1.0.0: OK


- `devtools::check_win_devel()`
   - CRAN incoming feasibility: One note (Maintainer: 'Tiago Olivoto <tiagoolivoto@gmail.com>')

- `devtools::check_win_release()`

# Test environments

- local OS X install, R 3.6.1
- Ubuntu Linux 16.04 LTS, R-release, GCC
- Rhub
   - Windows Server 2008 R2 SP1, R-devel, 32/64 bit
   - Ubuntu Linux 16.04 LTS, R-release, GCC
   - Fedora Linux, R-devel, clang, gfortran
   - One note
      - Maintainer: 'Tiago Olivoto <tiagoolivoto@gmail.com>'
- win-builder (devel and release)

# R CMD check results



# Other questions

I have a vignette (vignettes/vignettes_metan.Rmd) that I would like to include when installing package but when I run `devtools::check()` the files in inst/doc are deleted. How should I proceed in this case?

Thank for your time in helping me with this submission.

Best Regards

Tiago Olivoto

