# stdlib
import os
import sys

# Remove current dir from sys.path so Python doesn't try to import the source dir
sys.path.remove(os.getcwd())

# this package
# from .engines import search  # noqa: E402,F401
from .spectra import spectra  # noqa: E402,F401
