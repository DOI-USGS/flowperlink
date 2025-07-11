"""Sub-package for linking point datasets with TerrainWorks hydrography points"""

# Contact Information
__author__ = "Steven Pestana"
__email__ = "spestana@usgs.gov"

# provide version, PEP - three components ("major.minor.micro")
__version__ = '0.0.1'

# Import classes for easier access
from .link import TerrainWorksLink

def describe():
    return "This is the terrainworks sub-package of flowperlink for linking point datasets with TerrainWorks hydrography points."
