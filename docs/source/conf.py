import os
import sys

sys.path.append(os.path.abspath("./ext"))
sys.path.append(os.path.abspath("./conf_params"))

project = "Simple VOF Solver"
author = "Naoki Hori"
copyright = f"2022-2023, {author}"

extensions = [
        "myliteralinclude",
        "sphinx.ext.mathjax",
]

from alabaster_params import html_theme
from alabaster_params import html_theme_options

from mathjax_params import mathjax_path
from mathjax_params import mathjax3_config

