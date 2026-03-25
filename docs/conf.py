# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

##---------- This part needed for readthedocs -----------##

import os
import sys

# Add the parent directory of 'docs' to the path
# Assuming structure: project_root/docs/conf.py and project_root/src/my_package/
# sys.path.insert(0, os.path.abspath('../gastronomy')) 
# OR if your code is in the root:
sys.path.insert(0, os.path.abspath('..'))

#autodoc_mock_imports = ['numpy', 'astropy']


# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'Gastronomy'
copyright = '2026, Lia Corrales'
author = 'Lia Corrales'
release = '1.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx.ext.autodoc', 'numpydoc']

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

html_baseurl = '/'

# Make sure theme options are set
html_theme_options = {
    'navigation_depth': 4,
    'collapse_navigation': False,
}
