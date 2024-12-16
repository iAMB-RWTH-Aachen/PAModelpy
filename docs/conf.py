
# Set the source directory for your Sphinx documentation files
source_suffix = ['.rst', '.md']  # If using Markdown files with MyST

# -- Project information

project = 'PAModelpy'
copyright = '2024, iAMB, RWTH Aachen University'
author = 'Samira van den Bogaard'

release = '0.4.2'
version = '0.4.2'

# -- General configuration

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
    'myst_parser',  # Add this for Markdown support
]

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
}
intersphinx_disabled_domains = ['std']

templates_path = ['_templates']

# -- Options for HTML output

html_theme = 'sphinx_rtd_theme'

# -- Options for EPUB output
epub_show_urls = 'footnote'

# Enable specific MyST features
myst_enable_extensions = [
    "colon_fence",  # Support for ::: directives
    "deflist",      # Support for definition lists
    "html_admonition",  # HTML-style admonitions
    "html_image",   # Use HTML-style <img> tags
]