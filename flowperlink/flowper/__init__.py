"""Sub-package for linking FLOwPER datasets with hydrography flowlines"""

from ..__about__ import __author__, __email__, __version__

# Import classes for easier access
from .link import FlowperLink

def describe():
    return "This is the flowper sub-package of flowperlink for linking FLOwPER datasets with hydrography flowlines."
