from pkg_resources import get_distribution, DistributionNotFound

try:
    __version__ = get_distribution('scimap').version
except DistributionNotFound:
    __version__ = '(local)'


from . import preprocessing as pp