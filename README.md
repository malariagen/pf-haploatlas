# Mutation Discovery App

A project to develop an interactive web app for eventual public use. The overall aim is to give Pf7(+) users the ability to look at genetic, spatial, and temporal trends for any Pf gene of interest in a quick and intuitive way. 

### What will the app do?

We would like the app to be able to acheive the following, for any of the ~5000 genes in the Pf genome:

 - Show all discovered haplotypes for a selected gene
 - Generate a bar plot of sample counts per haplotype
 - Generate a stacked bar plot of major subpopulation proportions per haplotype
 - Generate an UpSet plot to detail the mutations within haplotypes
 - Have an area within the above plots which contains clickable features, which lead to:
   - An abacus plot per haplotype, detailing the haplotype's frequency across administrative divisons and years

### What features are key to the app?
- It should run very quickly 
- It will be interactive

## Resources
This app will be built off of code which exists from generating Figure 3A and Supplementary Figure 9 from the [Pf7 paper](https://wellcomeopenresearch.org/articles/8-22/v1)

The current code for this is located [here](https://github.com/malariagen/PDNA_SP_markers/blob/9b73a05d2c573e66ea70b282adb7055a66e2022e/work/3_haplotype_plots/20221221_haplotype_plots.ipynb)


## Project Management

We will use issues to structure the work, and management of issues can be dealt with at the [issue board](https://gitlab.com/malariagen/gsp/mutation-discovery-app/-/boards).

Labels are used to classify issues:

`status::doing` for things people are currently working on.

`status::currnt_sprint` for things planned for the current work block

`status::backlog` for things that will be tackled in the future

If you want to work on an issue, assign it to yourself and move it to the status::doing board. This just signals
other people you are actively working on that issue, and they should not try to tackle it. You can assign yourself an
issue and leave it in the current-sprint/backlog boards to signal other people you'd like to work on those in the future
(but other people may take it if required for a deadline).

## Working in the Repo

There is currently a folder `work` to contain code for this project. We can add new folders as required. 

We envision that the work will require two steps:
1. Precomputing data summaries for use by the app
2. Building the app 

### Branches and issues

As with most other malariagen/gsp work, we will complete work within the scope of a single issue, and create a branch to work on each. Name the branch using the format `<issue_number>_<issue_name>`, for example: 01_sample_summary_table. Use common sense to compress very long issue names into something meaningful.

Before closing an issue, please report your findings or a summary of the analysis results in a comment, it could also be a good idea to include links to any relevant notebook.

Please write meaningful commit messages (which issue number it addresses, etc.).

