# treeIS
This repository houses code to replicate figures in a paper examining importance sampling for large coalescent trees.

## Requirements

All tested on Ubuntu 24.04.

- the g++ compiler (tested on version 13.2.0),
- the [GNU Scientific Library](https://www.gnu.org/software/gsl/) (tested on version 2.7.1),
- [R](https://www.r-project.org/) (tested on version 4.3.3),
- the [viridis](https://cran.r-project.org/web/packages/viridis/index.html) package (tested on version 0.6.5),
- while not strictly necessary, [RStudio](https://posit.co/) or a similar IDE will also make life easier
  (tested on RStudio version 2023.06.1).

## Compilation

Navigate to `Finite-alleles` and call the following:

- `make`
- `make gt`

Then navigate to `Infinite-sites` and call the following:

- `make`
- `make sd`
- `make gt`
- `make prec`

Each compilation should be near-instantaneous.

## Usage

### Finite alleles

Execute all the shell scripts in `Finite-alleles/Scripts`.
Then execute `plots.R` in `Finite-alleles/Results`.

*Warning*: These shell scripts have considerable runtimes of between 10 minutes and 4 hours each.

### Infinite sites

Execute `precompute.sh` in `Infinite-sites/Scripts`, which can take around 2 hours.
Then execute the remaining shell scripts in the same folder.
The longest script should last around 3 hours.

## Caveat

While this repository contains working implementations of the [Griffiths-Tavare](https://doi.org/10.1006/tpbi.1994.1023), [Stephens-Donnelly](https://doi.org/10.1111/1467-9868.00254), and [Hobolth-Uyenoyama-Wiuf](https://doi.org/10.2202/1544-6115.1400) importance sampling schemes for the coalescent, they are not intended for stand-alone use.
The source files include hard-coded paths to data files or other resources, and will not work correctly unless called via the provided shell scripts at their provided locations in the filesystem.

## Accompanying manuscript

Link to be posted.
