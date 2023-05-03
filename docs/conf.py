# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'pybalonor'
copyright = '2023, Malte Jörn Ziebarth'
author = 'Malte Jörn Ziebarth'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx.ext.mathjax',
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinxcontrib.cairosvgconverter'
]

mathjax_path = "https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.5/MathJax.js?config=TeX-MML-AM_CHTML"

autodoc_typehints = "none"
add_module_names = False

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_static_path = ['_static']
html_css_files = [
    'css/custom.css',
]
html_theme_options = {
    "sidebar_width" : "7cm",
}

# -- Options for LaTeX output ------------------------------------------------
latex_elements = {
    'preamble' : r'''
    \usepackage{enumitem}
    \setlistdepth{10}
    '''
}
