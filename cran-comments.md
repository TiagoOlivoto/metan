# Release summary
This is a path release (v1.7.1) that includes minor improvements. Please, see the `news.md` file for more details.

# Test environments
- local OS X install, R 4.0.0
- Rhub
   - Windows Server 2008 R2 SP1, R-devel, 32/64 bit
   - Ubuntu Linux 16.04 LTS, R-release, GCC
   - Fedora Linux, R-devel, clang, gfortran
- win-builder (devel and release)


# R CMD check results (run with --run-donttest --as-cran)
0 errors | 0 warnings | 2 notes

> 
checking CRAN incoming feasibility ... Note_to_CRAN_maintainers
Maintainer: 'Tiago Olivoto <tiagoolivoto@gmail.com>'

**I think the following note doesn't matter.**

> 
checking top-level files ... NOTE
Files 'README.md' or 'NEWS.md' cannot be checked without 'pandoc' being installed.

I can put these in .Rbuildignore but CRAN supports NEWS.md now.





