# Publishing CellAgent

This directory is ready to be published as a standalone GitHub repository.

## Recommended Repository

- Owner: `xuanyuelingwu`
- Repository name: `CellAgent`
- Suggested visibility: public, unless the project should remain private

The package metadata in `pyproject.toml` currently points to:

```text
https://github.com/xuanyuelingwu/CellAgent
```

Update those URLs before publishing if you choose a different repository name.

## Publish With Git

Create an empty GitHub repository first, then run:

```bash
cd CellAgent
git init
git branch -M main
git add .
git commit -m "Initial standalone CellAgent release"
git remote add origin https://github.com/xuanyuelingwu/CellAgent.git
git push -u origin main
```

## Publish With GitHub CLI

If `gh` is installed and authenticated:

```bash
cd CellAgent
git init
git branch -M main
git add .
git commit -m "Initial standalone CellAgent release"
gh repo create xuanyuelingwu/CellAgent --public --source . --remote origin --push
```

Use `--private` instead of `--public` if needed.

## Publish Without Git

If Git is unavailable, upload the prepared zip archive through the GitHub web UI:

1. Create an empty repository named `CellAgent`.
2. Open the repository page in GitHub.
3. Choose "uploading an existing file".
4. Upload the contents of the prepared standalone repo directory, not the parent directory itself.

GitHub web upload is fine for the first import, but Git is strongly recommended after that.
