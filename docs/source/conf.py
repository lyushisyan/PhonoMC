# Configuration file for the Sphinx documentation builder.

import os
import sys

sys.path.insert(0, os.path.abspath("../.."))

project = "PhonoMC"
copyright = "2026, Shixian Liu"
author = "Shixian Liu"
release = "development"

extensions = [
    "sphinx.ext.mathjax",
]

templates_path = ["_templates"]
exclude_patterns = []

source_suffix = ".rst"
master_doc = "index"
language = "en"

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
html_logo = "_static/PhonoMC-logo-1040x800.png"
html_theme_options = {
    "navigation_depth": 4,
    "collapse_navigation": False,
    "sticky_navigation": True,
    "titles_only": False,
    "logo_only": True,
}

html_title = "PhonoMC Documentation"
html_short_title = "PhonoMC"
