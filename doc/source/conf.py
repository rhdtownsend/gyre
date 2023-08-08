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
import re

sys.path.insert(0, os.path.abspath('exts'))

import sphinx_rtd_theme


# -- Project information -----------------------------------------------------

project = 'GYRE'
author = 'Rich Townsend & The GYRE Team'
version = "master"
release = "master"
branch = "master"
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
    'sphinx.ext.intersphinx',
    'sphinxcontrib.spelling',
    'sphinxcontrib.cairosvgconverter',
    'sphinx-prompt',
    'sphinx_substitution_extensions',
    'ads_cite',
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

# CSS
html_css_files = ["_static/theme_overrides.css"] # override wide tables in RTD theme

# Set master doc
master_doc = 'index'

# Set logo
html_logo = 'gyre-logo.png'

# Set up Extlinks
extlinks = {
    'wiki': ('https://en.wikipedia.org/wiki/%s', None),
    'netlib': ('https://www.netlib.org/%s', None),
    'git': ('https://github.com/%s', None),
    'repo': ('https://github.com/rhdtownsend/gyre/blob/{:s}/%s'.format(branch), None)
}

# Set site-wide targets

if release == 'master':
    tarball = 'master'
else:
    tarball = 'v{0:s}'.format(release)

targets = {
    'github-tarball': 'https://codeload.github.com/rhdtownsend/gyre/tar.gz/{0:s}'.format(tarball),
    'gyre-forums': 'http://www.astro.wisc.edu/~townsend/gyre-forums/',
    'mesa-sdk': 'http://www.astro.wisc.edu/~townsend/static.php?ref=mesasdk',
    'mesa': 'http://mesa.sourceforge.net/',
}

rst_prolog = '\n'.join(['.. _{:s}: {:s}'.format(x, targets[x]) for x in targets])

# Add substitutions for sphinx_substitution_extensions

rep_exts = {"release": release,
            "author": author}

for rep_ext_key, rep_ext_val in rep_exts.items():
    rst_prolog += "\n.. |{:s}| replace:: {:s}".format(rep_ext_key, rep_ext_val)

# Mathjax & Latex macros

macros = {}

with open('macros.def', encoding='utf-8') as f:
    line = f.readline()
    while line:
        key, value = line.rstrip().split('\t')
        macros[key] = value
        line = f.readline()

mathjax_macros = {}
latex_preamble = ''

for key, value in macros.items():
    argnums = re.findall('#(\d)', value)
    if argnums:
        n_args = int(max(argnums))
        mathjax_macros[key] = [value, n_args]
        latex_preamble += f'\\newcommand{{\\{key}}}[{n_args}]{{{value}}}\n'
    else:
        mathjax_macros[key] = value
        latex_preamble += f'\\newcommand{{\\{key}}}{{{value}}}\n'

#mathjax_config = {                  
#    'TeX': { 
#        'Macros': mathjax_macros
#    }
#}
mathjax3_config = {                  
    'tex': { 
        'macros': mathjax_macros
    }
}

latex_elements = {
    'preamble': latex_preamble
}

# Enable email obfuscation
email_automode = True

# Intersphinx mapping
intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'numpy': ('https://numpy.org/doc/stable', None),
    'scipy': ('https://docs.scipy.org/doc/scipy', None),
    'matplotlib': ('https://matplotlib.org/stable', None),
    'astropy': ('https://docs.astropy.org/en/latest', None),
    'pygyre': ('https://pygyre.readthedocs.io/en/latest', None)
}

# Equation number formatting
math_eqref_format = '{number}'
                       
# Spelling
spelling_word_list_filename='spelling_wordlist.txt'
