reproduce
=========

[![DOI](https://zenodo.org/badge/140708206.svg)](https://zenodo.org/badge/latestdoi/140708206)

Assessing the quality and reproducibility of absolute proteomic data. Companion repository to the pre-print [_"Benchmarking accuracy and precision of intensity-based absolute quantification of protein abundances in Saccharomyces cerevisiae"_](https://www.biorxiv.org/content/10.1101/2020.03.23.998237v1).

This repository is administered by Benjamín J. Sánchez ([@BenjaSanchez](https://github.com/benjasanchez)), Division of Systems and Synthetic Biology, Department of Biology and Biological Engineering, Chalmers University of Technology.

Requirements
------------

* [R](https://www.r-project.org/) (tested with v3.6)
* [RStudio](https://rstudio.com/) (tested with v1.2)
* [Rtools](https://cran.r-project.org/bin/windows/Rtools/) (in case of a windows setup)
* [Git](https://git-scm.com/) (or any Git client of your choice, e.g. [Github Desktop](https://desktop.github.com/))

Installation
------------

* Clone locally this repository.
* From the command window in RStudio, run each line in `requirements.txt`.

Usage
-----

All analysis can be reproduced by knitting the RMarkdown file `code/iBAQstudy.Rmd`

Repository Structure
--------------------

    |- code/            # all programmatic code relating to the project
    |  +- templates/    # scripts for generating template files
    |
    |- data/            # all data from the study
    |  |- raw_internal/ # raw data generated in-lab or by collaborators, will not be altered
    |  |- raw_external/ # data from third-party sources, databases etc, will not be altered
    |     +- colormaps/ # color palettes used for all figures
    |
    |- doc/             # documentation for the study and other explanatory material
    |  +- paper/        # contains the generated pdf from knitting the markdown file
    |
    |- results          # all output from workflows and analyses
    |  |- figures/      # graphs, designated for manuscript figures
    |  +- pictures/     # diagrams, images, and other non-graph graphics
    |
    |- .gitignore       # files that will not sync to Github
    |- LICENSE          # license
    |- README.md        # the top level description of content
    |- reproduce.Rproj  # contains project information used to customize the behavior of RStudio  
    +- requirements.txt # the requirements file for reproducing the analysis environment

Acknowledgements
----------------

The initial file and directory structure of this project was developed by a group of participants in the Reproducible Science Curriculum Workshop, held at [NESCent] in December 2014. The structure is based on, and heavily follows the one proposed by [Noble 2009], with a few but small modifications. The [original repository] has been modified to the [reproducible-research-init repository] and adapted to this project.

[original repository]: https://github.com/Reproducible-Science-Curriculum/rr-init
[reproducible-research-init repository]: https://github.com/EngqvistLab/reproducible-research-init
[NESCent]: http://nescent.org
[Noble 2009]: http://dx.doi.org/10.1371/journal.pcbi.1000424
