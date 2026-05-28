# Setup Instructions

## 1. Install Quarto

Download from https://quarto.org/docs/get-started/ and install (it's a single binary, ~100MB).

Verify: `quarto --version`

## 2. Replace your sebastianlech.github.io repo

Since you're starting fresh, the cleanest approach:

```bash
# Clone your existing repo
git clone https://github.com/SebastianLech/sebastianlech.github.io
cd sebastianlech.github.io

# Delete all existing Jekyll files
git rm -rf .

# Copy everything from _quarto-site/ into the repo root
cp -r /path/to/probability/_quarto-site/. .

# Preview locally first
quarto preview

# When it looks good, push
git add .
git commit -m "migrate to Quarto"
git push origin main
```

## 3. Enable GitHub Pages

In your repo on GitHub:
- Settings → Pages
- Source: **Deploy from a branch**
- Branch: **gh-pages** / `/(root)`
- Save

The first push will trigger the GitHub Action which builds and deploys.
Your site will be live at https://sebastianlech.github.io in ~2 minutes.

## 4. Add the probability notes as a submodule (optional but clean)

If you want to keep the probability repo separate from the main site:

```bash
cd sebastianlech.github.io
git submodule add https://github.com/SebastianLech/probability notes/probability-src
```

Then symlink or copy `.qmd` files from the submodule into `notes/probability/` as part of the build.
Simpler: just keep `.qmd` files directly in the main site repo and use the probability repo for `.tex` source only.

## 5. Adding new content

**New probability chapter:**
- Create `notes/probability/ch02-combinatorics.qmd`
- Add a row to the table in `notes/probability/index.qmd`
- Push — the site rebuilds automatically

**New article or notebook:**
- Drop a `.qmd` or `.ipynb` file into `articles/`
- It appears automatically in the articles listing

**New Python project:**
- Drop a `.qmd` or `.ipynb` file into `projects/`
- Update `projects/index.qmd` with a card

## 6. Local workflow

```bash
# Live preview with hot reload
quarto preview

# One-off render
quarto render
```

Quarto caches executed notebook output (`freeze: auto` in _quarto.yml),
so re-renders are fast — only changed files re-execute.
