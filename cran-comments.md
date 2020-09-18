# Release summary
This is a minor release (v1.9.0) that list packages providing the Rd macros in 'Imports' instead of 'Suggests' as suggested by the CRAN team and add new functions for stability analysis. See `News.md` for more information.

# Test environments
- local OS X install, R 4.0.0
- Rhub
   - Windows Server 2008 R2 SP1, R-devel, 32/64 bit
   - Ubuntu Linux 16.04 LTS, R-release, GCC
   - Fedora Linux, R-devel, clang, gfortran
- win-builder (devel and release)


# R CMD check results (run with --run-donttest --as-cran)
0 errors | 0 warnings | 3 notes

> 
checking CRAN incoming feasibility ... Note_to_CRAN_maintainers
Maintainer: 'Tiago Olivoto <tiagoolivoto@gmail.com>'

**I think the following note doesn't matter.**

>
unable to verify current time

**I believe this requires a fix in the check function by the R development team or the web-resource coming online again.**

> 
checking top-level files ... NOTE
Files 'README.md' or 'NEWS.md' cannot be checked without 'pandoc' being installed.

**I can put these in .Rbuildignore but CRAN supports NEWS.md now.**





