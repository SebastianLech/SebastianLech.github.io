## Workflow

### Preview locally
```bash
export PATH="$PATH:/Applications/quarto/bin"
cd ~/Desktop/projects/sebastianlech.github.io
quarto preview
```
Opens a live preview at `localhost:4848` with hot reload on save.

### Publish changes
```bash
export PATH="$PATH:/Applications/quarto/bin"
cd ~/Desktop/projects/sebastianlech.github.io
quarto publish gh-pages
```