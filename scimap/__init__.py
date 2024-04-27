import importlib.metadata

try:
    __version__ = importlib.metadata.version('scimap')
except importlib.metadata.PackageNotFoundError:
    __version__ = '(local)'

import scimap

print(f"Running SCIMAP  {scimap.__version__}")

from . import preprocessing as pp
from . import plotting as pl
from . import tools as tl
from . import helpers as hl
from . import external as ex
