# Release summary
This is a minor release (v1.12.0) that add new functions and minor improvements. See `News.md` for more information.

# Test environments
* GitHub Actions (ubuntu-16.04): devel, release, oldrel, 3.5, 3.4, 3.3
* GitHub Actions (windows): release, oldrel
* GitHub Actions (macOS): release
* win-builder: devel


# R CMD check results (run with --run-donttest --as-cran)
0 errors | 0 warnings | 2 notes

* checking CRAN incoming feasibility ... Note_to_CRAN_maintainers Maintainer: 'Tiago Olivoto <tiagoolivoto@gmail.com>'

   - I think the following note doesn't matter.

* checking top-level files ... NOTE Files 'README.md' or 'NEWS.md' cannot be checked without 'pandoc' being installed.
   - I can put these in `.Rbuildignore` but CRAN supports `NEWS.md` now.
