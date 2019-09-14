# Resubmission
This is a resubmission. In this version, I have made the following changes based on the Martina Schmirl's comments. Each comment was labeled and my responses were inserted following each one.

1. Please do not start the description with "This package", package name, title or similar. Just start with "Provides ...".
   * Thanks for the suggestion. In this version I started the description with "Provides...."

2. If there are references describing the methods in your package, please add these in the description field of your DESCRIPTION file in the form authors (year) <doi:...> authors (year) <arXiv:...> authors (year, ISBN:...) or if those are not available: <https:...> with no space after 'doi:', 'arXiv:', 'https:' and angle brackets for auto-linking.
   * There is no specific references to cite in the DESCRIPTION file.
   
3. Please add \\value to .Rd files and explain the functions results in the documentation. (See: Writing R Extensions
<https://cran.r-project.org/doc/manuals/r-release/R-exts.html#Documenting-functions>)
   * Greate comment! Thanks! I have included/updated the \\value section of the following functions:
      * `AMMI_indexes()`
      * `Annichiarico()`
      * `anova_ind()`
      * `as.lpcor()`
      * `as.split_factors()`
      * `bind_cv()`
      * `can_cor()`
      * `comb_vars()`
      * `corr_ci()`
      * `corr_plot()`
      * `covcor_design()`
      * `cv_ammi()`
      * `cv_ammif()`
      * `cv_blup()`
      * `desc_stat()`
      * `ecovalence()`
      * `fai_blup()`
      * `ge_factanal()`
      * `ge_plot()`
      * `ge_reg()`
      * `ge_stats()`
      * `get_model_data()`
      * `is.lpcorr()`
      * `is.split_factors()`
      * `mahala()`
      * `mahala_design()`
      * `make_mat()`
      * `make_sym()`
      * `mtsi()`
      * `pairs_mantel()`
      * `plot.*()` and `plot_*()` functions
      * `rbind_fill()`
      * `resca()`
      * `resp_surf()` 
      * `waas()`
      * `wsmp()`
      * `waasb()`
      
4. \\dontrun{} should be only used if the example really cannot be executed (e.g. because of missing additional software, missing API keys, ...) by the user. That's why wrapping examples in \\dontrun{} adds the comment ("# Not run:") as a warning for the user. Seems not necessary. Please unwrap the examples if they are executable in < 5 sec, or create additionally small toy examples to allow automatic testing (then replace \\dontrun with \\donttest). You could also replace \\dontrun{} with \\donttest, but it would be preferable to have automatic checks for functions.

* In this version I have:
   * Deleted \\dontrun{} of the following functions:
      * `bind_cv()`
      * `clustering()`
      * `comb_vars()`
      * `corr_ci()`
      * `corr_plot()`
      * `covcor_design()`
      * `desc_stat()`
      * `ecovalence()`
      * `path_coefff()`
      * `plot.fai_blup()`
      * `plot.mtsi()`
      * `plot.wsmp()`
      * `plot_ci()`
      * `summary.mtsi()`
      * `summary.can_cor()`
      * `summary.path_coeff()`
      * `wsmp()`
   * Changed \\dontrun{} with \\donttest{} in the following functions:
      * `cv_ammi()`
      * `cv_ammif()`
      * `cv_blup()`
      * `plot.cv_ammif()`
      

5. You write information messages to the console that cannot be easily suppressed. It is more R like to generate objects that can be used to extract the information a user is interested in, and then print() that object. Instead of print()/cat() rather use message()/warning()  or if(verbose)cat(..) if you really have to write text to the console. (except for print() and summary() functions)

* Thank you for the suggestion! In this version I have included "verbose" argument that allows supressing messages to the console in the following functions.
   * `corr_ss()`
   * `find_outliers()`
* For all functions that write information message to the console the argument `verbose = FALSE` may be used to easily suppress those messages. Thank you!

I have a vignette that I want to build it when the package is installed but when I run `devtools::check()` the files in inst/doc are deleted. How should I proceed in this case?

Thank for your time in reviewing my R package
Best Regards
Tiago Olivoto

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
0 errors | 0 warnings | 0 notes
