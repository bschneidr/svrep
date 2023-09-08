## R CMD check results

0 errors | 0 warnings | 0 notes

* Update since last submission: 

  * Fixed issue caused by yesterday's updates on r-devel, which make `crossprod()` now a primitive. As a result, on r-devel the 'Matrix' package no longer exports `crossprod()` and this in turn causes an error on r-devel for the 'survey' package (and 'svrep', by extension). The 'Matrix' package authors are fixing the issue so that `crossprod()` will again be exported. For now, I've conditionally disabled affected tests in 'svrep' so that they will be skipped until the Matrix package again exports `crossprod()`.

* As noted in last submission (rejected because of new issues caused by yesterday's updates on r-devel):

  * Tracked down and removed 12 unmarked UTF-8 strings that caused the NOTE for "checking data for non-ASCII characters".

  * A small handful of tests are now disabled for Atlas/MKL. It is prohibitively difficult to reproduce the test errors from CRAN's Atlas/MKL checks on r-devel and the extensive tests are passing on all other tested platforms architectures. There are known inconsistencies with runs of Atlas/MKL for the same inputs and the same use of `set.seed()`- avoiding these cannot be done from within a given R package's test suite. Please see this recent discussion on r-help: https://www.mail-archive.com/r-help@r-project.org/msg267775.html

  * More updates to resolve the LaTeX issues that only appear on "r-oldrel-macos". I'm not sure why there was a LaTeX issue though for "r-oldrel-macos"; there wasn't any issue for the other flavors. I removed uses of `\begin{cases}...\end{cases}` that the "r-oldrel-macos" manual building couldn't handle.


## Downstream dependencies

* Check OK for downstream dependencies: 'nrba' (direct dependency) and 'idcnrba' (indirect dependency)
