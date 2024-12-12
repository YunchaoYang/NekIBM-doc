# NekIBM-doc
This is the repository of the NekIBM documentation adapted from the original [Nek5000 docs](https://github.com/Nek5000/NekDoc) written using the [Sphinx](http://www.sphinx-doc.org/) documentation framework.


## How to update

### 1. Prepare sphinx environment
```bash
module load conda
conda activate a_python_env # torch-timm
pip install sphinx
pip install recommonmark # add markdown 
pip install sphinx_rtd_theme

```

### 2. modify conf.py
```bash
emacs source/conf.py # language = 'en'
```

### 3. modify rst files

```bash
git add --all
git commit -m "update message"
git push origin master
```

### 4. publish on GitHub Pages

To update the GitHub Page, a contributor must have write permissions to the main NekDoc repository.
The GitHub Page should not contain any edits that are newer than the master branch of the main repository.

Workflow:

Checkout the latest master
run 

```
make gh-pages
```
