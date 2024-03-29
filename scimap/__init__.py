from pkg_resources import get_distribution, DistributionNotFound

try:
    __version__ = get_distribution('scimap').version
except DistributionNotFound:
    __version__ = '(local)'

import scimap
print(f"Running SCIMAP  {scimap.__version__}")

from . import preprocessing as pp
from . import plotting as pl
from . import tools as tl
from . import helpers as hl