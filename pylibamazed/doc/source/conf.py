# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys

sys.path.insert(0, os.path.abspath("../../python/pylibamazed"))


# -- Project information -----------------------------------------------------

project = "PyLibAmazed"
copyright = "2021, CESAM-LAM"
author = "CESAM-LAM"

# The full version, including alpha/beta/rc tags
release = "0.28.0"


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "myst_parser",
    "sphinx.ext.mathjax",
    "sphinx.ext.autodoc",
    "sphinx_autodoc_typehints",
    "sphinx.ext.autosummary",
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_rtd_theme"

html_theme_options = {"titles_only": False, "navigation_depth": 4}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]

html_css_files = [
    "css/custom.css",
]

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

myst_gfm_only = True

# Autodoc and autodoc typehints parameters
add_module_names = False
always_document_param_types = True
typehints_use_signature = True
typehints_use_signature_return = True
typehints_defaults = "comma"


def skip_member(app, what, name, obj, skip, opts):
    if isinstance(obj, property):
        if getattr(obj.fget, "is_doc_method", False):
            return False
        else:
            return True

    elif callable(obj):
        # Explicitly show init methods
        if getattr(obj, "__name__", "") == "__init__":
            return False
        # Show all methods that are marked as API methods
        elif getattr(obj, "is_doc_method", False):
            return False
        # If method comes from c++, keep default behavior
        # Therefore, for c++ methods need to specify all elements to document in autodoc
        else:
            return True
    # Returns none for other elements -> default behavior


def setup(app):
    app.connect("autodoc-skip-member", skip_member)
