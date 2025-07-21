"""Sub-package for linking point datasets with TerrainWorks hydrography points"""

from ..__about__ import __author__, __email__, __version__

# Import classes for easier access
from .link import TerrainWorksLink

def describe():
    return "This is the terrainworks sub-package of flowperlink for linking point datasets with TerrainWorks hydrography points."
