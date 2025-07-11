# flowperlink

![Generic badge](https://img.shields.io/badge/Version-0.0.1-<COLOR>.svg)

The flowperlink package extends the functionality of hydrolink by adding new features specific for linking FLOwPER (streamflow permanence) observation data to various hydrography datasets (e.g. NHDPlusHR, TerrainWorks), for ready ingestion of the observation data into the Probability of Streamflow Permanence (PROSPER) model ([Jaeger et al., 2019](https://doi.org/10.1016/j.hydroa.2018.100005)).

## Features

* The `FlowperLink` class, derived from hydrolink's `CustomHydrography`, allows the snapping of FLOwPER point observations using tributary junction information and adjustable buffer distances to account for GPS accuracy and hydrography linework uncertainty.
* The `TerrainWorksLink` class, for snapping FLOwPER and other streamflow permanence observations to high-resolution terrain-derived hydrography from TerrainWorks.
* A preprocessing script, `preprocess_flowlines.py`, to derive tributary junction information from NHD hydrography datasets is included, to allow snapping FLOwPER points using the tributary junction fields in the observation data.
* Other helper scripts: `download_flowper.py` for downloading FLOwPER data from ScienceBase via command line, and `merge_flowper.py` to merge multiple FLOwPER shapefiles into a single GeoDataFrame after downloading.

## Contact Information

Steven Pestana - spestana@usgs.gov

## USGS Software Release Information

## Quick Start

### Installation

It is recommended to use [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file) to create a new python environment for using flowperlink. Activate your new environment before installing the flowperlink package.

To install flowperlink from the main branch of this repository using pip:

```bash
pip install git+https://code.usgs.gov/streamflow-permanence/flowperlink.git
```

To install a different branch of flowperlink, such as a development branch:

```bash
pip install git+https://code.usgs.gov/streamflow-permanence/flowperlink.git@branch_name
```

The flowperlink package can also be installed from the source code directly. This may be of interest if you are developing new features for the package.

```bash
# clone the repository to get the source code
git clone https://code.usgs.gov/streamflow-permanence/flowperlink.git
# navigate to the flowperlink directory
cd flowperlink
# create the conda environment
conda env create -f environment.yml
# activate the new environment
conda activate flowperenv
# install flowperlink
pip install -e .
```
### Command line tools

#### Downloading FLOwPER data

The commands below illustrate usage of the `download_flowper.py` and `merge_flowper.py` scripts to access FLOwPER observations from ScienceBase and prepare them for linking to hydrography datasets.

```bash
# Download FLOwPER shapefiles
python download_flowper.py --download_directory D:/FLOwPER_data/downloads/
```

```bash
# Merge FLOwPER shapefiles that have been downloaded
python merge_flowper.py --input_directory D:/FLOwPER_data/downloads/ --output_filepath D:/FLOwPER_data/FLOwPER_merged.shp --crs 4326
```


#### Preprocessing hydrography flowlines

The command below illustrates usage of the `preprocess_flowlines.py` script with an NHDPlusHR hydrography dataset, where tributary junction in formation is retrieved from the NHDFlowline, NHDPlusFlow, and NHDPlusFlowlineVAA layers of the geopackage file. In this example, a distance of 30 m is used to define the portion of flowlines near tributary junctions, and the output is written to a new geopackage file. 

```bash
python preprocess_flowlines.py 
    --flowlines_filepath NHDPLUS_H_1701_HU4_GPKG.gpkg 
    --flowline_layer NHDFlowline 
    --flowlines_identifier Permanent_Identifier 
    --flow_layer NHDPlusFlow 
    --from_identifier FromPermID 
    --to_identifier ToPermID  

    --flowlineVAA_layer NHDPlusFlowlineVAA 
    --flowlinesVAA_identifier NHDPlusID 
    --divergence_field Divergence 
    --trib_jcn_dist 30 
    --output_filepath NHDPLUS_H_1701_HU4_GPKG_preprocessed.gpkg 
```



#### Snapping with tributary junction informaiton

```python
from flowperlink.flowper import FlowperLink

nhdplushr = FlowperLink(points = 'FLOwPER_points.shp',   
                  flowlines = 'NHDPLUS_H_1701_HU4_GPKG_preprocessed.gpkg', 
                  source_identifier = 'GlobalID',  
                  flowlines_identifier = ' Permanent_Identifier', 
                  water_name = 'Strm_Nm_Sp',  
                  flowline_name = 'GNIS_name', 
                  buffer_m = 'AccuracyH', 
                  buffer_multiplier = 10, 
                  default_buffer = 100, 
                  no_stream_name_min_buffer = 10, 
                  yes_stream_name_min_buffer = 15, 
                  max_buffer_distance = 100, 
                  keep_points_attributes = ['TribJncTyp'], 
                  keep_flowlines_attributes = ['mainstem_flag', 'trib_jcn']) 

nhdplushr.hydrolink_method(method = 'name_match', 
                     hydro_type = 'flowline',  
                     trib_jcn = 'TribJncTyp',  
                     outfile_name = 'nhdplushr_snapped_output.gpkg',  
                     similarity_cutoff = 0.6) 

nhdplushr.write_connecting_lines(outfile_name='nhdplushr_snapped_connectors.gpkg') 

nhdplushr.buffered_points_gdf.to_file('nhdplushr_snapped_buffer_pts.gpkg') 
```

#### Snapping to TerrainWorks hydrography

```python
from flowperlink.terrainworks import TerrainWorksLink

tw = TerrainWorksLink(points = 'FLOwPER_points.shp',   
                      flowlines = 'Nodes_UpperDeschutes.gdb',  
                      source_identifier='OBJECTID', 
                      flowlines_identifier='NODE_ID',  
                      water_name = None,  
                      flowline_name = None, 
                      buffer_m = 100, 
                      buffer_multiplier = 1, 
                      default_buffer = 100, 
                      no_stream_name_min_buffer = 10, 
                      yes_stream_name_min_buffer = 15, 
                      max_buffer_distance = 100, 
                      flowline_grid_offsets = (1, 1)) 

tw.hydrolink_method(method = 'closest', 
                    trib_jcn = None, 
                    hydro_type = 'flowline',  
                    outfile_name = 'tw_snapped_output.gpkg', 
                    similarity_cutoff = 0.6) 

tw.write_connecting_lines(outfile_name='tw_snapped_connectors.gpkg') 

tw.buffered_points_gdf.to_file('tw_snapped_buffer_pts.gpkg') 
```
## Testing

Test scripts can be ran using [pytest](https://docs.pytest.org/en/stable/). To run these tests, clone the repository (see instructions above in the Installation section). Pytest will execute test scripts located within the `test` directory where the filename has the prefix `test_`.

```bash
pytest
```

# Hydrography datasets

Table based on [this](https://www.usgs.gov/news/which-nhd-product-do-you-need-and-which-do-you-have)

| | [NHD](https://www.usgs.gov/national-hydrography/national-hydrography-dataset?qt-science_support_page_related_con=0#qt-science_support_page_related_con) | [NHDPlus HR](https://www.usgs.gov/national-hydrography/nhdplus-high-resolution) | [NHDPlus V2](https://www.epa.gov/waterdata/nhdplus-national-hydrography-dataset-plus) | [NHDPlus V1](https://nhdplus.com/NHDPlus/NHDPlusV1_home.php) | [TerrainWorks](https://terrainworks.com/)
| ------ | ------ | ------ | ------ | ------ | ------ |
|    Currency    |     As of October 1, 2023, the NHD was retired. NHD data will continue to be available, but no longer maintained.    |    Based on snapshot of 1:24K NHD from 2015-2019 timeframe (differs by area)    |   Based on snapshot of 1:100K NHD from 2011 timeframe     |   NHDPlusV1 is replaced in its entirety by NHDPlusV2     |   |
|   Level of Detail     |    1:24K or more detailed (1:63,360 or more detailed in AK)    |    1:24K or more detailed    |    1:100K (more detailed for HI, U.S. Virgin Islands, Puerto Rico, and Pacific Island territories)    |    1:100K  | 2 m |
|   Geography     |    CONUS, AK, HI, Puerto Rico, U.S. Virgin Islands, and several Pacific Island territories    |   CONUS, HI, Puerto Rico, U.S. Virgin Islands, and several Pacific Island territories; in development in AK     |    CONUS, HI, Puerto Rico, U.S. Virgin Islands, and several Pacific Island territories; not available in AK    |        | parts of Oregon |
|   Additional Attributes     |    None    |   Value-Added Attributes, Mean Annual Flow Estimates     |    Value-Added Attributes, Mean Annual and Mean Monthly Flow Estimates    |        |     | 
|    Additional Features    |   None     |   Elevation-derived catchment areas based on 1:24K NHD and 10 meter DEM      |  Elevation-derived catchment areas based on 1:100K NHD and 30 meter DEM      |        |    |
|    Link    |    [download by HUC](https://prd-tnm.s3.amazonaws.com/index.html?prefix=StagedProducts/Hydrography/NHD/)    |   [download by HUC](https://prd-tnm.s3.amazonaws.com/index.html?prefix=StagedProducts/Hydrography/NHDPlusHR/)     |    [download by HUC](https://www.epa.gov/waterdata/get-nhdplus-national-hydrography-dataset-plus-data#Download)    |        |    |
|    Formats    |    GDB, GPKG, SHP    |    GDB, GPKG, Raster    |    GDB, SHP, Raster     |        |    |
|    [Hydrolink 1.0](https://code.usgs.gov/sas/bioscience/hlt/hydrolink)    |   "NHD HR" [MapServer](https://hydromaintenance.nationalmap.gov/arcgis/rest/services/HEM/NHDHigh/MapServer)     |   -    |   "NHD MR" [MapServer](https://watersgeo.epa.gov/arcgis/rest/services/NHDPlus/NHDPlus/MapServer)      |        |  n/a  |
|    [Hydrolink 2.0 in dev](https://code.usgs.gov/spestana/hydrolink/-/tree/custom-hydrography?ref_type=heads)  |   [download](https://prd-tnm.s3.amazonaws.com/index.html?prefix=StagedProducts/Hydrography/NHD/) to use CustomHydrography module  |   [download](https://prd-tnm.s3.amazonaws.com/index.html?prefix=StagedProducts/Hydrography/NHDPlusHR/) to use CustomHydrography module     |    [download](https://www.epa.gov/waterdata/get-nhdplus-national-hydrography-dataset-plus-data#Download) to use CustomHydrography module    |        |  n/a  |
|    ["FLOwPERlink"](#)    |        |        |        |        |  use ["TerrainWorksLink"](#)  |
|       |        |        |        |        |