# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import sys
import os

sys.path.insert(0, os.path.abspath('.'))
sys.path.insert(0, os.path.abspath('..'))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'planetMagFields'
copyright = '2025, Ankit Barik'
author = 'Ankit Barik'
release = '1.5.6'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.mathjax',
    'sphinx.ext.intersphinx',
    'sphinx.ext.viewcode',
    'sphinx.ext.coverage',
    'sphinx.ext.napoleon',
    'sphinxcontrib.bibtex',
    'sphinx_copybutton'
    ]
bibtex_bibfiles = ['bib.bib']
bibtex_default_style = 'unsrt'
templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'furo'
html_static_path = ['_static']

pygments_style = 'sphinx'
pygments_dark_style = "monokai"
master_doc = 'index'
source_suffix = '.rst'

latex_elements = {}
latex_elements['preamble']=r"""

   \usepackage{amsmath}
   \usepackage{amssymb}
"""

copybutton_exclude = '.linenos, .gp, .go'
copybutton_prompt_text = r">>> |\.\.\. |\$ |In \[\d*\]: | {2,5}\.\.\.: |Out\[\d*\]: | {5,8}: "
copybutton_prompt_is_regexp = True
