
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->
<!-- badges: end -->

# Astrocytes_editing

## Background

This repository contains analysis code for the publication found
[here](), which explored the response of astrocytes, neurons and
astro-neuronal co-cultures to alpha-synuclein oligomers. The analysis
includes work shared between two repositories:

- The repository at <https://github.com/karishdsa/ipscAstrNeurCocul>
  contains single cell and bulk RNA sequencing analysis of the cell
  culture models utilised in the study.
- This current repository contains the analysis of RNA editing in the
  cell culture models, as well as in a [published
  dataset](https://link.springer.com/article/10.1007/s00401-021-02343-x)
  of post-mortem brain samples. The results of this repository can be
  viewed interactively at

## Code contents

Within this repository you will find:

| Directory            | Description                                                                                                                                             |
|----------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------|
| [docs](docs)         | Contains all `.Rmd`s and their corresponding `.html`s describing analyses performed for this project. These can be view interactively at: [link](#TODO) |
| [scripts](scripts)   | Contains all non R and R related scripts.                                                                                                               |
| [raw_data](raw_data) | External tables used in analyses.                                                                                                                       |
| [renv](renv)         | `renv`-related scripts                                                                                                                                  |
| [results](results)   | Results from all analyses.                                                                                                                              |
| [rscripts](rscripts) | Contains analysis scripts. Each script contains a one-line description and is also referenced in its corresponding `.Rmd`.                              |

## Reproducibility

<!-- Modify selection below depending on how package dependencies have been managed. -->

### `renv`

<!-- Consider using renv for reproducibility. Delete this section if you will not be doing this. -->

This repository uses [`renv`](https://rstudio.github.io/renv/index.html)
to create a reproducible environment for this R project.

1.  When you first launches this project, `renv` should automatically
    bootstrap itself, thereby downloading and installing the appropriate
    version of `renv` into the project library.
2.  After this has completed, you can use `renv::restore()` to restore
    the project library locally on your machine.

For more information on collaborating with `renv`, please refer to this
[link](https://rstudio.github.io/renv/articles/collaborating.html).

## License

The code in this repository is released under an MIT license. This
repository is distributed in the hope that it will be useful to the
wider community, but without any warranty of any kind. Please see the
[LICENSE](LICENSE) file for more details.

## Citation

The RNA editing analysis was undertaken with
[JACUSA2](https://github.com/dieterich-lab/JACUSA2). Converting the
output to bed format was undertaken using a custom script written by
[George
Young](https://www.linkedin.com/in/george-young-6153a11bb/?originalSubdomain=uk).
