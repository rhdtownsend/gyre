# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# http://www.sphinx-doc.org/en/master/config

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
sys.path.insert(0, os.path.abspath('exts'))

import sphinx_rtd_theme


# -- Project information -----------------------------------------------------

project = 'GYRE'
author = 'Rich Townsend & The GYRE Team'
version = "5.2"
release = "5.2"
copyright = '2020, Rich Townsend & The GYRE Team'


# -- General configuration ---------------------------------------------------

# Numbered figures
numfig = True

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx_rtd_theme',
    'sphinx.ext.mathjax',
    'sphinx.ext.extlinks',
    'sphinxcontrib.bibtex',
    'sphinx-prompt',
    'sphinx_substitution_extensions',
    'nml_roles'
]


# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']


# -- Additional configuration ------------------------------------------------

# sphinx_rtd options
html_theme_options = {
    'collapse_navigation': True,
    'sticky_navigation': True,
    'navigation_depth': 4,
    'includehidden': True,
    'titles_only': False,
    'logo_only': True
}

# Set master doc
master_doc = 'index'

# Set logo
html_logo = 'gyre-logo-200.png'

# Set up Extlinks
extlinks = {
    'wiki': ('https://en.wikipedia.org/wiki/%s', ''),
    'ads': ('https://ui.adsabs.harvard.edu/abs/%s/abstract', ''),
    'repo': ('https://github.com/rhdtownsend/gyre/blob/release-{:s}/%s'.format(release), '%s')
}

# Set up additional roles
from docutils.parsers.rst import roles, nodes

roles.register_generic_role('nml_n', nodes.literal)
roles.register_generic_role('nml_v', nodes.literal)

# Set substitutions for sphinx_substitution_extensions 
substitutions = [
    ('|release|', release)
]

# Set site-wide targets
targets = {
    'github-tarball': 'https:///github.com/rhdtownsend/gyre/archive/v{0:s}.tar.gz'.format(release),
    'gyre-forums': 'http://www.astro.wisc.edu/~townsend/gyre-forums/',
    'mesa-sdk': 'http://www.astro.wisc.edu/~townsend/static.php?ref=mesasdk',
    'mesa': 'http://mesa.sourceforge.net/',
}

rst_prolog = '\n'.join(['.. _{:s}: {:s}'.format(x, targets[x]) for x in targets])

# Latex macros

mathjax_config = {                  
    'TeX': {                        
        'Macros': {                 
            'Msun': r'{\rm M}_{\odot}',
            'deriv': [r'\frac{{\rm d}^{#3}#1}{{\rm d}#2^{#3}}', 3],
            'pderiv': [r'\frac{\partial^{#3}#1}{\partial#2^{#3}}', 3],
            'ii': r'{\rm i}',
            'disc': r'\mathcal{D}'
        }
    }
}
