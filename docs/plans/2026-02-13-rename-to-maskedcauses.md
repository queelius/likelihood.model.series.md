# Rename `likelihood.model.series.md` to `maskedcauses` — Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Rename the R package from `likelihood.model.series.md` to `maskedcauses` across all repos, GitHub, r-universe, and CRAN-facing metadata.

**Architecture:** The rename is mechanical but wide: ~620 references across 4+ repos plus GitHub infrastructure. We work inside-out — package internals first, then downstream dependents, then external references (blog, papers, r-universe). The GitHub repo rename is done via `gh` CLI and GitHub automatically redirects the old URL. The `Remotes` field is removed since all deps are on CRAN.

**Tech Stack:** R (devtools, pkgdown, testthat), git, gh CLI, Hugo

**Naming convention:** The old name `likelihood.model.series.md` becomes `maskedcauses` everywhere. The underscore variant `likelihood_model_series_md` becomes `maskedcauses`. GitHub URLs change from `queelius/likelihood.model.series.md` to `queelius/maskedcauses`.

---

## Phase 1: Package Internals (this repo)

### Task 1: Rename DESCRIPTION

**Files:**
- Modify: `DESCRIPTION`

**Step 1: Update package metadata**

Change these fields in `DESCRIPTION`:
```
Package: maskedcauses
Title: Likelihood Models for Systems with Masked Component Cause of Failure
URL: https://queelius.github.io/maskedcauses/
BugReports: https://github.com/queelius/maskedcauses/issues
```

Remove the entire `Remotes` block (lines 29-31):
```
Remotes:
    github::queelius/algebraic.mle,
    github::queelius/algebraic.dist,
    github::queelius/likelihood.model
```

All three packages are on CRAN — `Remotes` is no longer needed and CRAN rejects it.

**Step 2: Run `Rscript -e "devtools::document()"`**

This regenerates NAMESPACE and Rd files with the new package name.

Expected: clean completion, NAMESPACE header changes from `likelihood.model.series.md` to `maskedcauses`.

**Step 3: Commit**

```bash
git add DESCRIPTION NAMESPACE man/
git commit -m "Rename package to maskedcauses, remove Remotes field"
```

---

### Task 2: Update pkgdown config and rebuild

**Files:**
- Modify: `_pkgdown.yml`

**Step 1: Update URL in `_pkgdown.yml`**

```yaml
url: https://queelius.github.io/maskedcauses/
```

**Step 2: Commit**

```bash
git add _pkgdown.yml
git commit -m "Update pkgdown URL to maskedcauses"
```

---

### Task 3: Update vignettes

**Files:**
- Modify: `vignettes/exponential_series.Rmd`
- Modify: `vignettes/weibull_homogeneous_series.Rmd`
- Modify: `vignettes/weibull_series.Rmd` (if it has library() call)
- Modify: `vignettes/censoring_comparison.Rmd`
- Modify: `vignettes/framework.Rmd`

**Step 1: Replace all `library(likelihood.model.series.md)` with `library(maskedcauses)`**

In each vignette, find and replace:
- `library(likelihood.model.series.md)` -> `library(maskedcauses)`
- Any prose mentioning `likelihood.model.series.md` as the package name -> `maskedcauses`
- `vignette("...", package = "likelihood.model.series.md")` -> `vignette("...", package = "maskedcauses")`

**Step 2: Commit**

```bash
git add vignettes/
git commit -m "Update vignettes to use maskedcauses package name"
```

---

### Task 4: Update test infrastructure

**Files:**
- Modify: `tests/testthat.R`

**Step 1: Replace package name in test harness**

In `tests/testthat.R`:
```r
library(maskedcauses)
testthat::test_check("maskedcauses")
```

**Step 2: Run full test suite**

```bash
Rscript -e "devtools::test()"
```

Expected: 798 tests, all passing. The S3 methods and function names don't change — only the package identity.

**Step 3: Commit**

```bash
git add tests/
git commit -m "Update test harness to maskedcauses"
```

---

### Task 5: Update README

**Files:**
- Modify: `README.Rmd`
- Regenerate: `README.md`

**Step 1: Replace all references in `README.Rmd`**

Find and replace:
- `likelihood.model.series.md` -> `maskedcauses` (package name in prose, code, URLs)
- `queelius/likelihood.model.series.md` -> `queelius/maskedcauses` (GitHub URLs)
- `queelius.github.io/likelihood.model.series.md` -> `queelius.github.io/maskedcauses` (pkgdown URL)
- `queelius.r-universe.dev/likelihood.model.series.md` -> `queelius.r-universe.dev/maskedcauses` (r-universe URL)

**Step 2: Regenerate README.md**

```bash
Rscript -e "devtools::build_readme()"
```

If `build_readme()` is not available, use:
```bash
Rscript -e "rmarkdown::render('README.Rmd', output_format = 'github_document')"
```

**Step 3: Commit**

```bash
git add README.Rmd README.md
git commit -m "Update README for maskedcauses rename"
```

---

### Task 6: Update CLAUDE.md

**Files:**
- Modify: `CLAUDE.md`

**Step 1: Replace all references**

Find and replace throughout:
- `likelihood.model.series.md` -> `maskedcauses`

Keep the architecture documentation accurate — the model class names (`exp_series_md_c1_c2_c3`, etc.) and function names do NOT change. Only the package name changes.

**Step 2: Commit**

```bash
git add CLAUDE.md
git commit -m "Update CLAUDE.md for maskedcauses rename"
```

---

### Task 7: Update paper (internal)

**Files:**
- Modify: `paper/Makefile` (change `MAIN = likelihood_model_series_md` to `MAIN = maskedcauses`)
- Modify: `paper/likelihood_model_series_md.tex` (rename file to `paper/maskedcauses.tex`, update internal references)

**Step 1: Rename and update**

```bash
mv paper/likelihood_model_series_md.tex paper/maskedcauses.tex
```

Update `paper/Makefile`:
```makefile
MAIN = maskedcauses
```

Update references inside `paper/maskedcauses.tex`:
- Package name mentions -> `maskedcauses`
- GitHub URLs -> `queelius/maskedcauses`
- pkgdown URLs -> `queelius.github.io/maskedcauses`

**Step 2: Commit**

```bash
git add paper/
git commit -m "Update internal paper for maskedcauses rename"
```

---

### Task 8: Update simulations

**Files:**
- Modify: `simulations/R/sim_helpers.R`

**Step 1: Replace `library(likelihood.model.series.md)` with `library(maskedcauses)`**

**Step 2: Commit**

```bash
git add simulations/
git commit -m "Update simulation helpers for maskedcauses rename"
```

---

### Task 9: Full verification

**Step 1: Run R CMD check**

```bash
Rscript -e "devtools::check()"
```

Expected: 0 errors, 0 warnings, 0 notes (or only pre-existing notes).

**Step 2: Run test suite with coverage**

```bash
Rscript -e "devtools::test()"
```

Expected: 798 tests passing.

**Step 3: Build pkgdown site**

```bash
Rscript -e "pkgdown::build_site()"
```

Verify the built site in `docs/` uses `maskedcauses` throughout.

**Step 4: Commit docs**

```bash
git add docs/
git commit -m "Rebuild pkgdown site as maskedcauses"
```

---

## Phase 2: GitHub Repo Rename

### Task 10: Rename the GitHub repository

**Step 1: Rename via gh CLI**

```bash
gh repo rename maskedcauses --repo queelius/likelihood.model.series.md --yes
```

GitHub automatically creates a redirect from the old URL. All existing links will continue to work, but new references should use the new URL.

**Step 2: Update local git remote**

```bash
git remote set-url origin https://github.com/queelius/maskedcauses.git
```

**Step 3: Push all commits**

```bash
git push origin master
```

**Step 4: Verify redirect works**

```bash
curl -sI https://github.com/queelius/likelihood.model.series.md | head -5
```

Expected: 301 redirect to `https://github.com/queelius/maskedcauses`.

---

## Phase 3: r-universe

### Task 11: Update r-universe config

**Files:**
- Modify: `/home/spinoza/github/rlang/queelius.r-universe.dev/packages.json`

**Step 1: Update the entry**

Change:
```json
{
  "package": "likelihood.model.series.md",
  "url": "https://github.com/queelius/likelihood.model.series.md"
}
```
To:
```json
{
  "package": "maskedcauses",
  "url": "https://github.com/queelius/maskedcauses"
}
```

**Step 2: Commit and push**

```bash
cd /home/spinoza/github/rlang/queelius.r-universe.dev
git add packages.json
git commit -m "Rename likelihood.model.series.md to maskedcauses"
git push
```

r-universe rebuilds automatically on push.

---

## Phase 4: Downstream Dependent — dfr.lik.series.md

### Task 12: Update dfr.lik.series.md

**Files:**
- Modify: `/home/spinoza/github/rlang/dfr_lik_series_md/DESCRIPTION`
- Modify: `/home/spinoza/github/rlang/dfr_lik_series_md/R/methods.R`
- Modify: `/home/spinoza/github/rlang/dfr_lik_series_md/tests/testthat/test-cross-validate.R`

**Step 1: Update DESCRIPTION**

In `Suggests`, change `likelihood.model.series.md` to `maskedcauses`.
In `Remotes`, change `github::queelius/likelihood.model.series.md` to `github::queelius/maskedcauses`.

**Step 2: Update R/methods.R comment**

Change the comment referencing `likelihood.model.series.md` to `maskedcauses`.

**Step 3: Update test-cross-validate.R**

Replace all occurrences:
- `skip_if_not_installed("likelihood.model.series.md")` -> `skip_if_not_installed("maskedcauses")`
- `likelihood.model.series.md::exp_series_md_c1_c2_c3()` -> `maskedcauses::exp_series_md_c1_c2_c3()`

**Step 4: Run tests**

```bash
cd /home/spinoza/github/rlang/dfr_lik_series_md
Rscript -e "devtools::test()"
```

**Step 5: Commit and push**

```bash
git add DESCRIPTION R/methods.R tests/
git commit -m "Update dependency from likelihood.model.series.md to maskedcauses"
git push
```

---

## Phase 5: External References — Papers

### Task 13: Update masked-causes-in-series-systems paper

**Files:**
- Modify: `/home/spinoza/github/papers/masked-causes-in-series-systems/refs.bib`

**Step 1: Update bibliography entry**

Change the package citation to reference `maskedcauses` and the new GitHub URL.

**Step 2: Commit and push**

```bash
cd /home/spinoza/github/papers/masked-causes-in-series-systems
git add refs.bib
git commit -m "Update package citation to maskedcauses"
git push
```

---

### Task 14: Update masked-series-companions

**Files:**
- Modify: `/home/spinoza/github/papers/masked-series-companions/CLAUDE.md`
- Modify: `/home/spinoza/github/papers/masked-series-companions/README.md`
- Modify: `/home/spinoza/github/papers/masked-series-companions/01-identifiability/brief.md`
- Modify: `/home/spinoza/github/papers/masked-series-companions/03-observation-composition/brief.md`
- Modify: `/home/spinoza/github/papers/masked-series-companions/05-weibull-companion/brief.md`

**Step 1: Find-and-replace `likelihood.model.series.md` with `maskedcauses` in all 5 files**

**Step 2: Commit and push**

```bash
cd /home/spinoza/github/papers/masked-series-companions
git add .
git commit -m "Update package references to maskedcauses"
git push
```

---

## Phase 6: External References — Metafunctor Blog

### Task 15: Update metafunctor project page

**Files:**
- Modify: `/home/spinoza/github/repos/metafunctor/content/projects/likelihood.model.series.md/index.md`

**Step 1: Update frontmatter and body**

Update GitHub URL, pkgdown URL, and r-universe URL references.
Keep the project slug as-is (Hugo uses directory name as slug, and changing the directory would break existing links). Add a note that the package was renamed.

Alternatively, rename the directory:
```bash
mv content/projects/likelihood.model.series.md content/projects/maskedcauses
```
And set up an alias in frontmatter:
```yaml
aliases:
  - /projects/likelihood.model.series.md/
```

**Step 2: Update all internal links in posts that reference `/projects/likelihood.model.series.md/`**

Search for this path in all blog posts and update to `/projects/maskedcauses/`.

---

### Task 16: Update metafunctor blog posts

**Files (all in `/home/spinoza/github/repos/metafunctor/content/post/`):**
- Modify: `2026-02-13-observation-functors/index.md` (~9 references)
- Modify: `2026-02-05-likelihood-model-series-md/index.md` (~13 references)
- Modify: `2026-02-07-ecosystem-overview/index.md` (~4 references)
- Modify: `2025-12-03-mdrelax/index.md` (~2 references)
- Modify: `2025-12-03-weibull-model-selection/index.md` (~2 references)
- Modify: `2024-04-15-expo-masked-fim/index.md` (~1 reference)
- Modify: `2024-05-15-algebraic.mle/index.md` (~2 references)
- Modify: `2022-04-weibull-survival-cancer/index.md` (~1 tag reference)
- Modify: `2019-08-reliability-censored-data/index.md` (~1 tag reference)

**Step 1: In each file, replace:**
- `likelihood.model.series.md` -> `maskedcauses` (package name in prose and code)
- `queelius/likelihood.model.series.md` -> `queelius/maskedcauses` (GitHub URLs)
- `queelius.github.io/likelihood.model.series.md` -> `queelius.github.io/maskedcauses` (pkgdown URLs)

**Context-sensitive replacements:**
- In `linked_project` frontmatter lists: `likelihood.model.series.md` -> `maskedcauses`
- In prose, consider keeping the old name parenthetically on first mention in older posts: "maskedcauses (formerly likelihood.model.series.md)"

**Step 2: Commit**

```bash
cd /home/spinoza/github/repos/metafunctor
git add content/post/
git commit -m "Update blog posts for maskedcauses rename"
```

---

### Task 17: Update metafunctor series and paper pages

**Files:**
- Modify: `/home/spinoza/github/repos/metafunctor/content/series/statistical-reliability/_index.md`
- Modify: `/home/spinoza/github/repos/metafunctor/content/papers/likelihood.model.series.md/index.md`
- Modify: `/home/spinoza/github/repos/metafunctor/content/projects/queelius.r-universe.dev/index.md`

**Step 1: Update all package name references**

In the series `_index.md`, update the `associations.projects` entry and the ecosystem diagram.

For the paper page, consider renaming the directory or adding an alias:
```bash
mv content/papers/likelihood.model.series.md content/papers/maskedcauses
```
With alias:
```yaml
aliases:
  - /papers/likelihood.model.series.md/
```

**Step 2: Build and verify Hugo site**

```bash
cd /home/spinoza/github/repos/metafunctor
hugo
```

Expected: clean build, no broken internal links.

**Step 3: Commit**

```bash
git add content/
git commit -m "Update series, paper, and project pages for maskedcauses rename"
```

---

### Task 18: Rebuild and deploy metafunctor

**Step 1: Full Hugo build**

```bash
cd /home/spinoza/github/repos/metafunctor
hugo
```

**Step 2: Commit and push**

```bash
git add docs/
git commit -m "Rebuild site with maskedcauses references"
git push
```

---

## Phase 7: Final Verification

### Task 19: End-to-end verification

**Step 1: Install the renamed package from local source**

```bash
cd /home/spinoza/github/rlang/likelihood.model.series.md
Rscript -e "devtools::install()"
```

Verify: `library(maskedcauses)` loads without error.

**Step 2: Run R CMD check**

```bash
Rscript -e "devtools::check()"
```

Expected: 0 errors, 0 warnings.

**Step 3: Verify downstream**

```bash
cd /home/spinoza/github/rlang/dfr_lik_series_md
Rscript -e "devtools::test()"
```

**Step 4: Verify r-universe builds**

After pushing the r-universe config, check:
```bash
curl -s https://queelius.r-universe.dev/api/packages/maskedcauses | head -20
```

(May take a few minutes to rebuild.)

---

## Scope Exclusions

Things we are NOT doing in this rename:

- **Function names stay the same**: `exp_series_md_c1_c2_c3()`, `wei_series_md_c1_c2_c3()`, etc. These are the API. Renaming them would break every user's code.
- **S3 class names stay the same**: `"series_md"`, `"likelihood_model"`. These are internal dispatch mechanisms.
- **Column names stay the same**: `t`, `omega`, `x1`, etc.
- **Local directory name**: The directory `/home/spinoza/github/rlang/likelihood.model.series.md/` can optionally be renamed to `maskedcauses/` for clarity, but this is cosmetic and not required for the package to work.
- **CRAN submission**: That's a separate task after the rename is complete and verified.

---

## Risk Notes

1. **GitHub redirect**: After `gh repo rename`, old URLs redirect automatically. No immediate breakage. But the redirect may stop working if a new repo is ever created with the old name.

2. **r-universe build lag**: After pushing the config change, r-universe takes ~30 minutes to rebuild. During this window, `install.packages("maskedcauses", repos = "https://queelius.r-universe.dev")` won't work yet.

3. **Installed package conflict**: If the old `likelihood.model.series.md` is installed alongside the new `maskedcauses`, R won't care (different package names). But users should `remove.packages("likelihood.model.series.md")` to avoid confusion.

4. **CRAN has never seen this package**: The rename happens before the first CRAN submission, so there's no "old CRAN version" to worry about. Clean slate.
