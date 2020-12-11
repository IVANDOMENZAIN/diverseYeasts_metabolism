# diverseYeasts_metabolism

Metabolic modelling simulation and analysis tools for diverse yeast species.

- Brief Description

This repository contains a collection of scripts for simulation, modification, data integration and visualization using Genome-scale mEtabolic Models (GEMs) for probing metabolic capabilities of yeast organisms.

- KeyWords

 **Utilisation:** GEMs reconstruction, FBA, metabolic engineering, data integration; **Field:** Constraint-based methods; **Omic Source:** NA

Last update: 2020-12-11

This repository is administered by [@IVANDOMENZAIN](https://github.com/IVANDOMENZAIN), Division of Systems and Synthetic Biology, Department of Biology and Biological Engineering, Chalmers University of Technology

## Installation
### Required Software
* The scripts in this repository use the [RAVEN toolbox for MATLAB (version 2.0)](https://github.com/SysBioChalmers/RAVEN) as a dependency.
Alternatively, the models can also be directly used with the [COBRA toolbox for MATLAB](https://github.com/opencobra/cobratoolbox).
### Dependencies - Recommended Software
* Please see the [RAVEN toolbox repository](https://github.com/SysBioChalmers/RAVEN) for dependencies regarding RAVEN.
### Installation Instructions
* Clone master branch from this [GitHub repository](https://github.com/IVANDOMENZAIN/diverseYeasts_metabolism).
### Reconstruction
* The models in this repository have been reconstructed based on homology of genes with *S. cerevisiae* and using the model [YeastGEM](https://github.com/SysBioChalmers/yeast-GEM) as a template.

## Development Guidelines

Anybody is welcome to contribute to the development of this modeling and simulation Toolbox, but please abide by the following guidelines.

Each function should start with a commented section describing the function and explaining the parameters. Existing functions can clarify what style should be used. When making *any* changes to an existing function (`*.m`-file), change the date and name of developer near the bottom of this commented section in the *last modified* line.

### Bugfixes, new features and functions
* For any development, whether bugfixes, introducing new functions or new/updated features for existing functions: make a separate branch from `devel` and name the branch for instance after the function/feature you are fixing/developing. If you work on a fix, start the branch name with `fix/`, if you work on a feature, start the branch name with `feat/`. Examples: `fix/format_reactions` or `feat/new_algorithms`.
* Make commits to this branch while developing. Aim for backwards compatibility, and try to avoid very new MATLAB functions when possible, to accommodate users with older MATLAB versions.
* When you are happy with your new function/feature, make a pull request to the `devel` branch. Also, see [Pull request](#pull-request) below.

### Semantic commits
Use semantic commit messages to make it easier to show what you are aiming to do:
* `chore`: updating binaries (model `MATLAB` structures), UniProt databases, physiology and protemics data files, etc.
* `doc`: updating documentation (in `doc` folder) or explanatory comments in functions.
* `feat`: new feature added, e.g. new function introduced / new parameters / new algorithm / etc.
* `fix`: bugfix.
* `refactor`: see [code refactoring](https://en.wikipedia.org/wiki/Code_refactoring).
* `style`: minor format changes of functions (spaces, semi-colons, etc., no code change).

Examples:
```
feat: exportModel additional export to YAML
chore: update UniProt database for CENPK113-7D
fix: optimizeProb parsing results from Gurobi
```
More detailed explanation or comments can be left in the commit description.

### Pull request
* No changes should be directly commited to the `master` or `devel` branches. Commits are made to side-branches, after which pull requests are made for merging with `master` or `devel`.
* The person making the pull request and the one accepting the merge _cannot_ be the same person.
* A merge with the master branch invokes a new release.
