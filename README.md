# flowperlink

![Generic badge](https://img.shields.io/badge/Version-1.0.0-<COLOR>.svg)

The flowperlink package extends the functionality of [hydrolink](https://code.usgs.gov/sas/bioscience/hlt/hydrolink) by adding new features specific for linking [FLOwPER (stream FLOw PERmanence) observation data](https://www.usgs.gov/national-hydrography/national-hydrography-dataset) to various hydrography datasets (e.g. the [National Hydrography Dataset (NHD)](https://www.usgs.gov/national-hydrography/national-hydrography-dataset), TerrainWorks), for ready ingestion of the observation data into the PRObability of Streamflow PERmanence (PROSPER) model ([Jaeger et al., 2019](https://doi.org/10.1016/j.hydroa.2018.100005)).

## Features

* The `FlowperLink` class, derived from hydrolink's `CustomHydrography`, allows the snapping of FLOwPER point observations using tributary junction information and adjustable buffer distances to account for uncertainty in GPS accuracy and hydrography linework.
* The `TerrainWorksLink` class, for snapping FLOwPER and other streamflow permanence observations to high-resolution terrain-derived hydrography from TerrainWorks.
* A preprocessing script, `preprocess_flowlines.py`, to derive tributary junction information from NHD hydrography datasets is included, to allow snapping FLOwPER points using the tributary junction fields in the observation data.
* Other helper scripts: `download_flowper.py` for downloading FLOwPER data from ScienceBase via command line, and `merge_flowper.py` to merge multiple FLOwPER shapefiles into a single GeoDataFrame after downloading.

## Contact Information

Steven Pestana, Physical Scientist, Washington Water Science Center, Tacoma, Washington, [spestana@usgs.gov](mailto:spestana@usgs.gov), https://orcid.org/0000-0003-3360-0996

## USGS Software Release Information

Pestana, S., Chelgren, N., Heaston, E., Jaeger, K., Labriola, L., Sando, R., Wurster., P. (2025). flowperlink Version-1.0.0. U.S. Geological Survey software release. https://doi.org/10.5066/P1DN3GHB

IPDS-`181392`

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

The flowperlink package can also be installed from the source code directly. This may be of interest if you are developing new features for the package. Hydrolink must be installed separately in this case.

```bash
# clone the repository to get the source code
git clone https://code.usgs.gov/streamflow-permanence/flowperlink.git
# navigate to the flowperlink directory
cd flowperlink
# create the conda environment
conda env create -f environment.yml
# activate the new environment
conda activate flowperenv
# clone and install hydrolink (must be version 2.0 or above)
cd ..
git clone https://code.usgs.gov/sas/bioscience/hlt/hydrolink.git
cd hydrolink
pip install -e .
# navigate back to flowperlink
cd ../flowperlink
# install flowperlink
pip install -e .
```
### Command line tools

#### Downloading FLOwPER data

The commands below illustrate usage of the `download_flowper.py` and `merge_flowper.py` scripts to access FLOwPER observations from ScienceBase and prepare them for linking to hydrography datasets.

Download a default list of FLOwPER observations from ScienceBase:

```bash
python -m flowperlink.download_flowper --download_directory ./FLOwPER_data/downloads/
```

Or, download from a custom list of ScienceBase datasets in a csv file:

```bash
python -m flowperlink.download_flowper
    --download_directory ./FLOwPER_data/downloads/ 
    --input_filepath flowper_records.csv 
    --name_col dataset_name 
    --id_col sciencebase_ID
```

Merge the downloaded datasets into one shapefile:

```bash
# Merge FLOwPER shapefiles that have been downloaded
python -m flowperlink.merge_flowper
    --input_directory ./FLOwPER_data/downloads/ 
    --output_filepath ./FLOwPER_data/FLOwPER_merged.shp 
    --crs 4326
```


#### Preprocessing hydrography flowlines

The command below illustrates usage of the `preprocess_flowlines.py` script with the National Hydrography Dataset Plus, High Resolution (NHDPlusHR), where tributary junction information is retrieved from the NHDFlowline, NHDPlusFlow, and NHDPlusFlowlineVAA layers of the geopackage file. In this example, a distance of 30 m is used to define the portion of flowlines near tributary junctions, and the output is written to a new geopackage file. The NHDPlusHR file geopackage file used in this example can be downloaded [here](https://prd-tnm.s3.amazonaws.com/StagedProducts/Hydrography/NHDPlusHR/VPU/Current/GPKG/NHDPLUS_H_1701_HU4_GPKG.zip).

```bash
python -m flowperlink.preprocess_flowlines
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

Start an interactive python shell to run the following code (or run the code from your own python script or jupyter notebook):

```bash
conda activate flowperenv
python
```

```python
from flowperlink.flowper import FlowperLink

nhdplushr = FlowperLink(points = './FLOwPER_data/FLOwPER_merged.shp',   
                  flowlines = 'NHDPLUS_H_1701_HU4_GPKG_preprocessed.gpkg', 
                  points_identifier = 'GlobalID',  
                  flowlines_identifier = 'Permanent_Identifier', 
                  water_name = 'Strm_Nm_Sp',  
                  flowline_name = 'GNIS_Name', 
                  buffer_m = 'AccuracyH', 
                  buffer_multiplier = 10, 
                  replace_nodata_buffer_with = 100, 
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

tw = TerrainWorksLink(points = './FLOwPER_data/FLOwPER_merged.shp',   
                      flowlines = 'tw_example_flowlines.gpkg',  
                      points_identifier='GlobalID', 
                      flowlines_identifier='NODE_ID',  
                      water_name = None,  
                      flowline_name = None, 
                      buffer_m = 100, 
                      buffer_multiplier = 1, 
                      replace_nodata_buffer_with = 100, 
                      no_stream_name_min_buffer = 10, 
                      yes_stream_name_min_buffer = 15, 
                      max_buffer_distance = 100, 
                      flowline_grid_offsets = (1, 1)) 

tw.hydrolink_method(outfile_name = 'tw_snapped_output.gpkg') 

tw.write_connecting_lines(outfile_name='tw_snapped_connectors.gpkg') 

tw.buffered_points_gdf.to_file('tw_snapped_buffer_pts.gpkg') 
```
## Testing

Test scripts can be ran using [pytest](https://docs.pytest.org/en/stable/). To run these tests, clone the repository (see instructions above in the Installation section). Pytest will execute test scripts located within the `test` directory where the filename has the prefix `test_`.

```bash
pytest
```
