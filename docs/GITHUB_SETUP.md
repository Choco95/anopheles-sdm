# GitHub repository and local setup

The public repository is [Choco95/anopheles-sdm](https://github.com/Choco95/anopheles-sdm).

## Clone the repository

In PowerShell or another terminal, run:

```bash
git clone https://github.com/Choco95/anopheles-sdm.git
cd anopheles-sdm
```

Open `anopheles-sdm.Rproj` in RStudio and follow the [README](../README.md) and [SOP](SOP.md). Input datasets and generated model results are not included. Complete the documented input preparation and local validation before running the analysis.

## Save changes from your computer

Before committing, review the changed files and confirm that only intended scripts, documentation and metadata are staged:

```bash
git status --short
git add README.md docs scripts metadata
git diff --cached
git commit -m "Update Anopheles workflow documentation"
git push origin main
```

Adjust the staged paths to match your changes. The repository is public: do not commit credentials, restricted data or local run outputs. If Git requests author details, configure your preferred Git identity. Complete authentication using GitHub's secure sign-in flow.

For collaboration, create a branch and open a pull request. Do not force-push or overwrite other contributors' changes. After publishing changes, verify that the README renders and relative SOP and data-guide links work.

The repository contains a documented workflow draft; publication does not establish that the models have been reproduced. Record the R environment and complete the QA checklist for each locally validated run.
