---
title: 'FLOwPERlink: Improved automated snapping of streamflow permanence point observations to hydrography datasets using extensions to the hydrolink package'
tags:
  - python
  - hydrology
  - hydrography
  - streamflow permanence
  - surfacewater
authors:
  - name: Steven Pestana
    orcid: 0000-0003-3360-0996
    corresponding: true
    affiliation: 1
  - name: Laura Labriola
    orcid: 0000-0002-5096-2940
    affiliation: 2
  - name: Author 3
    orcid: 0000-0000-0000-0000
    affiliation: 3
  - name: Author 4
    orcid: 0000-0000-0000-0000
    affiliation: 4
affiliations:
  - name: U.S. Geological Survey, Washington Water Science Center, Tacoma, WA
    index : 1
  - name: U.S. Geological Survey, Oklahoma-Texas Water Science Center, Tacoma, WA
    index: 2
  - name: Affiliation 3
    index: 3
  - name: Affiliation 4
    index: 2

date: 10 July 2025
bibliography: paper.bib
---

# Summary

Simple visual field observations on the presence or absence of surface flow in streams and other waterbodies is a straightforward and cost-effective approach for data collection that is useful to a wide range of applications related to land and water resource management. FLOwPER [@Jaeger:2020], available in the Survey123R application [@ESRI:2015], is one of several feature mapping applications (Stream Tracker [@Kampf:2018], Crowd Water [@Seibert:2019], DRYvER [@Truchy:2023]) that provides a standard protocol for the collection of surface flow presence/absence anywhere in the world. Applying these data to be used in resource management decisions and for inclusion in modeling exercises require accurate georeferencing, or snapping, to a given hydrography. In many cases, an observation taken in the field of surface water presence or absence for a given stream reach is not co-located with the flowline that represents that stream reach in a given hydrography. The recently revised hydrolink tool [@Wieferich:2025] provides automated snapping of point features to any hydrography. In this particular application, the FLOwPERlink package leverages and extends the functionality of hydrolink to preprocess FLOwPER observation data for ready ingestion into the Probability of Streamflow Permanence (PROSPER) model [@Jaeger:2019] that provides predictions on the probability of a stream reach having year-round flow. 

# Statement of need

Georeferencing field observations of surface water presence/absence data to a given hydrography typically occurs manually, is not readily reproducible, can be time consuming, or requires access to proprietary software (e.g. HydroAdd). In addition, the georeferencing procedure can be particularly uncertain or error prone in dense and complex hydrographies, such as forested, rugged landscapes or if geographic information associated with the field observation is poor or unknown.  The FLOwPERlink package presented here leverages the hydrolink tool for automated snapping of FLOwPER observation points but also includes quality control components based on ancillary information available in the FLOwPER survey for improved snapping accuracy to a given hydrography. Specifically, FLOwPERlink uses the geographic horizontal accuracy that is part of the FLOwPER observation record to supplement the snapping based on stream same that is part of the hydrolink package. Where present in the FLOwPER observation record, information about the observation’s location relative to tributary junctions is also used to inform the linking to hydrography flowlines. In addition, information about the snapping procedure is preserved in the output file for evaluation of these data in end-user applications. The FLOwPERlink package can be used to georeference similar data from other crowd-sourced feature mapping applications where geographic horizontal accuracy and/or tributary junction information are preserved as part of the observation record. 

# Software description

FLOwPERlink requires python version 3.9 or newer, and the following packages: 
* hydrolink (version 2.0) 
* geopandas (version 1.0 or newer) 
* sciencebasepy 
* requests 

To leverage the tributary junction information of FLOwPER records, the hydrography dataset must be preprocessed to include mainstem_flag and trib_jcn fields (described below).

## Preprocessing flowlines

To snap FLOwPER observations to flowlines using their tributary junction information (described below), the flowlines within the hydrography dataset must have associated `mainstem_flag` and `trib_jcn fields`. A preprocessing script was developed to derive this tributary junction information from NHD, NHDPlusHR, and NHDPlusV2 datasets. The script, executable from a command line or as part of a separate bash script, preprocesses flowline data by eliminating closed loops, splitting the flowlines at specific distances from junctions, and categorizing the flowlines as either upstream mainstem, downstream mainstem, or tributary based on the network's direction information. 

The required information is pulled from the Flowline, Flow, and Value Added Attributes layers of each respective hydrography dataset. The script begins by removing divergences, and merging relevant attributes from the Flowline, Flow, and Value Added Attributes layers with the `clean_flowlines` function. The NetworkX package (https://networkx.org/) is then used to convert the resulting GeoPandas GeoDataFrame into a Directed Graph network. As a directed graph, the code can navigate the network upstream and downstream and understands how each flowline connects to the rest of the network. 

To classify a flowline as a mainstem or a tributary, at each junction the total length of all lines in the upstream network is computed within the `compute_upstream_network_length` function. This method of identifying mainstems and tributaries assumes that total line length is proportional to drainage area, and that mainstems will have larger drainage areas than tributaries. Though not all hydrography datasets will have the same detail or resolution everywhere (e.g. NHD has 1:24k lines in some locations, and high-resolution terrain derived 3DHP lines elsewhere), this method relies on adjacent basins lines having comparable detail or resolution. 

To label flowlines as near a junction, and upstream or downstream of that nearby junction, the `get_junction_flowlines` and `split_flowline` functions use a `trib_jcn_dist` parameter to split all flowlines a specified distance from each junction. This creates two or three sub-flowlines from every one original flowline. If the original flowline length is greater than twice the `trib_jcn_dist` value, then it will be split into three sub-flowlines, two having lengths of `trib_jcn_dist`, and the third consisting of the remainder of the flowline The sub-flowlines are then labeled as upstream or downstream of a junction, or not near a junction, based on the flow information in the directed graph. Finally, the processed flowlines with these new attributes are saved to a new file ready for hydrolinking. 

However, there still remains the problem of unmapped tributary junctions, where a main stem might be mapped in a hydrography dataset, but a small tributary is not. In these cases, the buffer distance (described below) will control whether or not the point is snapped to a flowline or not snapped at all (if far enough from other flowlines). Observations that are recorded as on a tributary, but where there is no nearby tributary to snap to (or recorded as mainstem, but where there is no nearby mainstem to snap to) are flagged with a processing message in their `trib_jnc_processing_message` field in the process described below. 

### Example usage

The command below illustrates usage of the python script with an NHDPlusHR hydrography dataset, where tributary junction in formation is retrieved from the NHDFlowline, NHDPlusFlow, and NHDPlusFlowlineVAA layers of the geopackage file. In this example, a distance of 30 m is used to define the portion of flowlines near tributary junctions, and the output is written to a new geopackage file. 

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
## Snapping with tributary junction informaiton

A new class, FlowperLink, was developed as a child class to CustomHydrography from hydrolink. This class for snapping FLOwPER, or similar, point observations extends the functionality of CustomHydrography to consider tributary junction information in the snapping process. The method `trib_jcn_match` compares tributary junction information within the point dataset against tributary junction information in the flowlines dataset. 

In FLOwPER records, the tributary junction type field, `TribJncTyp`, provides context about the location of the observation relative to the stream network. The possible values for this field are `On tributary`, `On mainstem upstream`, `On mainstem downstream`, and `No Data` (see the FLOwPER User's Guide or FLOwPER Quick Guide 2.0 [@Jaeger:2020] for more information about FLOwPER data). Based on these values, the `trib_jcn_match` method of the FlowperLink class applies filtering logic to possible flowline matches within the search region to narrow down the number of possible matches. If the value of `TribJncTyp` for a FLOwPER point is `No Data`, all possible flowline matches remain eligible for snapping, and a processing message is recorded to indicate the lack of tributary junction information at this point. If an observation is `On tributary` then only flowlines with a corresponding tributary value in their `mainstem_flag` attribute are retained for the next steps of the snapping process. If an observation is `On mainstem downstream` or `On mainstem upstream` then only flowlines with `mainstem_flag` set as `mainstem` and `trib_jcn` set as `downstream of junction` or `upstream of junction`, respectively, are retained. If none of the previous conditions are met, all possible flowlines are retained, and a processing message is recorded in the `trib_jnc_processing_message` field to indicate that matching using tributary junction information was attempted but was not successful.

Other parameters allow for fine-tuning this buffer distance. The `buffer_multiplier` parameter can be used to multiply the values in column `buffer_m` by a constant value. The horizontal accuracy of GPS point observations can be impacted by the canopy cover of dense forests or surrounding mountain terrain. Therefore, using a buffer multiplier of 10 (e.g. 10 $/sigma$ if we take the RMSE to be equal to the standard error, $/sigma$) is an approach that will help reduce false negatives (a point improperly identified as not belonging to any line) at the expense of increasing false positives (a point improperly identified as belonging to a line). 

The RMSE, or even 10$/sigma$, may be too restrictive in cases where the RMSE is very small (e.g. << 1 m). There is not only uncertainty in the point location but also uncertainty in the hydrography linework (which varies by which hydrography dataset is used). Additionally, the PROSPER models using these observations are resolving a 10-30 m spatial resolution (reach length, or grid cell size). Therefore, the linework uncertainty can be lumped with the point uncertainty by setting a lower limit on the search radius to about 10 m. 

The minimum allowable buffer distance can be set based on whether stream name information (for the `name_match` method) is present (`yes_stream_name_min_buffer`) or not (`no_stream_name_min_buffer`). The minimum buffer distance's dependence on stream name information allows for prioritizing the name matching method when a stream name has been recorded, over the horizontal accuracy of the point's GPS location, by setting `yes_stream_name_min_buffer` > `no_stream_name_min_buffer`. Finally, the `max_buffer_distance` parameter sets the maximum allowable buffer distance. 

### Example usage

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

## Snapping to TerrainWorks hydrography

A separate module was developed to snap FLOwPER and other streamflow permanence observations [@McShane:2017] to a high spatial resolution terrain derived hydrography from TerrainWorks (https://terrainworks.com/). This hydrography necessitated its own class, `TerrainWorksLink`, due to the representation of flow paths not as vector lines, but as points derived from gridded data. Each point in this hydrography dataset represents the lower left coordinate of a 2 m grid cell identified as a flow path from a flow accumulation raster. The TerrainWorks hydrography dataset does not include stream names, nor the attributes required to associate tributary junction information from FLOwPER observations, and therefore only snaps points informed by their geographic horizontal accuracy. Methods within the TerrainWorksLink class follow the same steps as in the FlowperLink class for hydrolinking using horizontal accuracy information, albeit for the point data of the TerrainWorks hydrography. There is an additional option to offset the snapping locations from the input hydrography points. A tuple of x (Eastings) and y (Northings) values in meters (positive or negative) can be passed to the `flowline_grid_offsets` parameter to offset the snapping location. Options for `name_match` and `trib_jcn` methods are included for future development but will currently raise a NotImplementedError if selected. 

### Example usage

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

## Handling FLOwPER observation data

The `download_flowper.py` script, designed to be run in the command line or through a wrapper bash script, is provided to download FLOwPER data from ScienceBase. The script parses user input specifying the download directory and optional input CSV file path, then executes the download process based on the provided parameters. The optional CSV file input argument allows users to specify which datasets to download through a CSV file containing ScienceBase data release IDs, otherwise a default list of FLOwPER data releases from 2019-2023 is downloaded [@Heaston:2025a; @Heaston:2025b; @Heaston:2024; @Dunn:2023; @Heaston:2022; @Heaston:2020]. The `load_flowper_data` function reads this CSV file if provided or returns the default list. The `download_flowper` function then establishes a session with ScienceBase using the `sciencebasepy` package [@Long:2023], checks for the existence of a specified download directory (creating it if necessary), and iterates through the ScienceBase IDs to download relevant files (while filtering out photos and metadata). It also handles the extraction of the downloaded ZIP files. 

The handling of FLOwPER datasets is continued with the `merge_flowper.py` script. Like the download script, it is designed to be run in the command line or through a wrapper bash script. The script merges multiple FLOwPER shapefiles into a single geopandas GeoDataFrame, and saved out to a single file. All input shapefiles are reprojected into a common coordinate reference system specified by the user. 

# Further development directions

# Acknowledgements

Any use of trade, product, or firm names is for descriptive purposes only and does not imply endorsement by the U.S. Government.

# References