"""Sub-package for linking FLOwPER datasets with hydrography flowlines"""

# Contact Information
__author__ = "Steven Pestana"
__email__ = "spestana@usgs.gov"

# provide version, PEP - three components ("major.minor.micro")
__version__ = '0.0.1'

# Import classes for easier access
from .link import FlowperLink

def describe():
    return "This is the flowper sub-package of flowperlink for linking FLOwPER datasets with hydrography flowlines."
