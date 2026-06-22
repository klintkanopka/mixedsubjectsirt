## Submission

This is a new submission of the package `mixedsubjectsirt` (version 1.0.0).

## Test environments

* local: macOS (aarch64), R 4.5.2
* win-builder: R-devel and R-release
* R-hub: (Linux, Windows, macOS)
* GitHub Actions: macOS-release, windows-release, ubuntu-devel/release/oldrel-1

## R CMD check results

0 errors | 0 warnings | 1 note

The note is the expected "New submission" note for a first-time submission:

```
Maintainer: 'Klint Kanopka <klint.kanopka@nyu.edu>'
New submission
```

On win-builder the incoming check additionally flags possibly misspelled words
in the Description. These are author surnames cited in the method references
(Angelopoulos, Bates, Broska, Duchi, Fannjiang, Howes, Zrnic) and the standard
psychometric term "unidimensional"; all are spelled correctly.

The method references in the Description are given in the requested
`<doi:...>` form. Publisher landing pages for some references return HTTP 403
to automated agents (bot blocking); the DOIs themselves resolve correctly.

## Downstream dependencies

There are currently no downstream dependencies (this is a new package).
