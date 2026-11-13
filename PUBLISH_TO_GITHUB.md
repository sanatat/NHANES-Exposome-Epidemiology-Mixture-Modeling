# Publishing to GitHub

The archive is ready for a new repository, but **has not been uploaded to GitHub**. Review code, attributions, and the public README before publishing. Only the repository code and docs (not `data/`, `results/`, credentials, or raw EPA files) should be committed.

## Using the GitHub website

Create an empty repository (for example `nhanes-exposome-epidemiology`), then upload the **contents** of this folder or use the command-line workflow below. Avoid a duplicate README created by GitHub when you initialize the remote.

## Local command line

```bash
cd ~/nhanes-exposome-github
git init
git add .
git status --short
git commit -m "Initial reproducible NHANES exposome analysis"
git branch -M main
# Replace YOUR_USERNAME with your actual GitHub account and use the new empty repo URL:
git remote add origin https://github.com/YOUR_USERNAME/nhanes-exposome-epidemiology.git
git push -u origin main
```

Optional with GitHub CLI installed and authenticated:

```bash
gh auth login
gh repo create nhanes-exposome-epidemiology --public --source=. --remote=origin --push
```

**Before pushing**, run `git status --short` and `git ls-files` and check that private data and local logs are absent. In particular, `data/` and `results/` should contain only `.gitkeep`.
