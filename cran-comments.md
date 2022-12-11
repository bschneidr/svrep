## R CMD check results

0 errors | 0 warnings | 1 note

* As in previous accepted submissions, I see a note
about DOI links in the DESCRIPTION. The links are valid,
but it seems the URL checker isn't working correctly for
these DOI links.

* I have reduced the runtimes of unit tests
and examples, since automated checks showed overall checktime was over 10
minutes on "r-devel-windows-x86_64" according to the automated CRAN checks
when I last submitted.

## Downstream dependencies

* This package does not yet have any downstream dependencies.
