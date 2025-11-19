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
  - name: Nathan Chelgren
    orcid: 0000-0003-0944-9165
    affiliation: 3
  - name: Emily Heaston
    orcid: 0000-0002-3949-391X
    affiliation: 3
  - name: Kristin Jaeger
    orcid: 0000-0002-1209-8506
    affiliation: 1
  - name: Roy Sando
    orcid: 0000-0003-0704-6258
    affiliation: 4
  - name: Patrick Wurster
    orcid: 0000-0003-2668-2014
    affiliation: 4
affiliations:
  - index: 1
    name: U.S. Geological Survey, Washington Water Science Center, Tacoma, WA
  - index: 2
    name: U.S. Geological Survey, Oklahoma-Texas Water Science Center, Oklahoma City, OK
  - index: 3
    name: U.S. Geological Survey, Forest and Rangeland Ecosystem Science Center, Corvallis, OR
  - index: 4
    name: U.S. Geological Survey, Wyoming-Montana Water Science Center, Helena, MT
    

date: 21 July 2025
bibliography: paper.bib
---

# Summary

Simple visual field observations on the presence or absence of surface flow in streams and other waterbodies is a straightforward and cost-effective approach for data collection that can be used for a wide range of applications related to land and water resource management. FLOw PERmanence (FLOwPER) [@Jaeger:2020], available in the Survey123® application [@ESRI:2015], is one of several feature mapping applications such as Stream Tracker [@Kampf:2018], Crowd Water [@Seibert:2019], or DRYing riVER networks (DRYvER; [@Truchy:2023]) that provides a standard protocol for the collection of surface flow presence/absence anywhere in the world. Using these data to inform resource management decisions and for inclusion in modeling exercises requires accurate linking of point observations, or snapping, to a given hydrography. In many cases, an observation taken in the field of surface water presence or absence for a given stream reach is not co-located with the line that represents that stream reach in a given hydrography. The recently revised hydrolink tool [@Wieferich:2025] provides automated snapping of point features to any hydrography. The FLOwPERlink package [@Pestana:2025] described here leverages and extends the functionality of hydrolink to prepare FLOwPER observation data for ingestion into the Probability of Streamflow Permanence (PROSPER) model [@Sando:2020; @Jaeger:2019] that provides predictions on the probability of a stream reach having year-round flow.

# Statement of need

Associating field observations of surface water presence/absence data to a given hydrography typically occurs manually, is not readily reproducible, can be time consuming, or requires access to proprietary software (e.g., HydroAdd; [@Demas:2022]). In addition, the georeferencing procedure can be particularly uncertain or error prone in dense and complex hydrographies, such as forested or rugged landscapes, or if the horizontal accuracy of the observation’s GPS coordinates was poor.  The FLOwPERlink package [@Pestana:2025] described here leverages the hydrolink tool [@Wieferich:2025] for automated snapping of FLOwPER observation points and includes quality control components based on ancillary information available in the FLOwPER survey for improved snapping accuracy to a given hydrography. Specifically, FLOwPERlink uses the geographic horizontal accuracy that is part of the FLOwPER observation record to better inform the snapping based on the stream name that is part of the hydrolink package. Where present in the FLOwPER observation record, information about the observation’s location relative to tributary junctions is also used to inform the linking to hydrography flowlines. In addition, information about the snapping procedure is preserved in the output file for evaluation of these data in end-user applications. The FLOwPERlink package can be used to georeference similar data from other crowd-sourced feature mapping applications where geographic horizontal accuracy and/or tributary junction information are preserved as part of the observation record.

# Software description

## Dependencies

FLOwPERlink [@Pestana:2025] requires Python version 3.9 or newer, and the following packages:
  - Hydrolink>=2.0.0 [@Wieferich:2025]
  - geopandas>=1.0.1 [@Bossche:2024]
  - sciencebasepy >= 2.0.18 [@Long:2023]
  - networkx [@Hagberg:2008]
  - setuptools [@PythonPackagingAuthority:2025]
  - requests [@Reitz:2024]

To leverage the tributary junction information of FLOwPER records, any hydrography dataset used must include the `mainstem_flag` and `trib_jcn` fields described below. FLOwPERlink includes a script to do this preprocessing with National Hydrography Dataset hydrographies. 

## Preprocessing flowlines

To snap FLOwPER observations to flowlines using their tributary junction information (described below), the flowlines within the hydrography dataset must have the associated attributes: `mainstem_flag` and `trib_jcn`. A preprocessing script was developed to derive this tributary junction information from the National Hydrography Dataset (NHD; [@USGS:2023]), NHD Plus High Resolution (NHDPlusHR; [@USGS:2025b]), and NHD Plus Version 2 (NHDPlusV2; [@USEPA:2011]) datasets. The script is executable from a command line or as part of a separate Bash script and preprocesses flowline data by eliminating closed loops (i.e., where flowlines representing a stream split and then re-join, creating parallel flow paths), splitting the flowlines at specific distances from junctions, and categorizing the flowlines as either upstream mainstem, downstream mainstem, or tributary based on the network's direction information.

The required information is pulled from the “Flowline”, “Flow”, and “Value Added Attributes” layers of each respective hydrography dataset. The script begins by removing divergences and merging relevant fields from the “Flowline”, “Flow”, and “Value Added Attributes” layers with the `clean_flowlines` function. The NetworkX package [@Hagberg:2008] is then used to convert the resulting GeoPandas GeoDataFrame [@Bossche:2024] into a Directed Graph network. As a directed graph, the code can navigate the network upstream and downstream and understands how each flowline connects to the rest of the network.

To classify a flowline as a mainstem or a tributary, the `compute_upstream_network_length` function computes the total length of all lines in the upstream network at each junction. The line with the longest total upstream length is classified at this junction as the mainstem, and all other lines at this junction are classified as tributaries. This allows a single stream to be classified as a tributary at a downstream node where it joins a larger river but be classified as a mainstem further upstream where it is joined by its own smaller tributaries. This method of identifying mainstems and tributaries assumes that total line length is proportional to drainage area, and that mainstems will have larger drainage areas than tributaries. This method relies on adjacent basin lines having comparable detail or resolution; however, not all hydrography datasets will have the same detail or resolution everywhere. For example, NHD has 1:24k lines in some locations and high-resolution terrain derived 3D Hydrography Program (3DHP) lines elsewhere.

To label flowlines as near a junction and as upstream or downstream of that nearby junction, the `get_junction_flowlines` and `split_flowline` functions use a `trib_jcn_dist` parameter to split all flowlines a specified distance from each junction. This creates two or three sub-flowlines from each original flowline. If the original flowline length is greater than twice the `trib_jcn_dist` value, then it will be split into three sub-flowlines, two having lengths of `trib_jcn_dist`, and the third consisting of the remainder of the flowline. The sub-flowlines are then labeled as upstream or downstream of a junction, or not near a junction, based on the flow information in the directed graph. Finally, the processed flowlines with these new fields (`mainstem_flag` and `trib_jcn`) are saved to a new file ready for hydrolinking.

However, unmapped tributary junctions, where a mainstem might be mapped in a hydrography dataset, still remains a problem. In these cases, the buffer distance parameters (described below) will control whether the point is snapped to a flowline or not snapped at all. Observations that are recorded as on a tributary, but where there is no nearby tributary to snap to (or recorded as mainstem, but where there is no nearby mainstem to snap to) are flagged with a processing message in their `trib_jnc_processing_message` field in the process described below.

### Example usage

The command below illustrates usage of the Python script with an NHDPlusHR hydrography dataset, where tributary junction information is retrieved from the NHDFlowline, NHDPlusFlow, and NHDPlusFlowlineVAA layers of the geopackage file. In this example, a 30 m distance is used to define the portion of flowlines near tributary junctions (`trib_jcn_dist`), and the output is written to a new geopackage file.

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

A new class, FlowperLink, was developed as a child class to CustomHydrography from hydrolink. This class for snapping FLOwPER, or similar point observations extends the functionality of CustomHydrography to consider tributary junction information in the snapping process \autoref{fig:trib_jcn_fig}. The method `trib_jcn_match` compares tributary junction information within the point dataset against tributary junction information in the flowlines dataset.

In FLOwPER records, the tributary junction type field, `TribJncTyp`, provides context about the location of the observation relative to the stream network. The possible values for this field are `On tributary`, `On mainstem upstream`, `On mainstem downstream`, and `No Data` (refer to the FLOwPER User's Guide or FLOwPER Quick Guide 2.0 for more information about FLOwPER data; [@Jaeger:2020]). Based on these values, the `trib_jcn_match` method of the FlowperLink class applies filtering logic to possible flowline matches within the search region to narrow down the number of possible matches. If the value of `TribJncTyp` for a FLOwPER point is `No Data`, all possible flowline matches remain eligible for snapping, and a processing message is recorded to indicate the lack of tributary junction information at this point. If an observation is `On tributary` then only flowlines with a corresponding tributary value in their `mainstem_flag` attribute are retained for the next steps of the snapping process. If an observation is `On mainstem downstream` or `On mainstem upstream` then only flowlines with `mainstem_flag` set as `mainstem` and `trib_jcn` set as `downstream of junction` or `upstream of junction`, respectively, are retained. If none of the previous conditions are met, all possible flowlines are retained, and a processing message is recorded in the `trib_jnc_processing_message` field to indicate that matching using tributary junction information was attempted but was not successful.

![Map of FLOwPER observation points [@Heaston:2025a; @Heaston:2025b; @Heaston:2024; @Dunn:2023; @Heaston:2022; @Heaston:2020] snapped to the NHDPlusHR hydrography dataset [@USGS:2025b] using tributary junction information. Though two of these FLOwPER point locations are closer to the mainstem flowline (solid line), they are snapped to the tributaries (dashed lines) based on the tributary junction information in the FLOwPER records.\label{fig:trib_jcn_fig}]{trib_jcn_fig.jpeg}

## Snapping with horizontal accuracy information

FLOwPER observations include a measure of the horizontal accuracy (root mean square error, RMSE) of the GPS coordinates of the observation point. That information can be used to help guide where on a given hydrography dataset's linework to snap the point observation \autoref{fig:gps_rmse_fig}. The original hydrolink tool (Version 1.0.0; [@Wieferich:2023]) includes a buffer distance parameter (`buffer_m`) to define the radius around a point of interest to search for candidate lines for snapping. The revised hydrolink tool (Version 2.0.0; [@Wieferich:2025]) within the CustomHydrography module can use a dynamic buffer, a different value for each point, and read from a field in the point dataset (in which case `buffer_m` is set to the field name rather than a static buffer distance value). Where a point might be missing a value in this field, the code relies on a `default_buffer` value which can be set to a static value or set to the mean value of the dataset.

Other parameters allow for fine-tuning this buffer distance. The `buffer_multiplier` parameter can be used to multiply the values in field `buffer_m` by a constant value. The horizontal accuracy of GPS coordinates can be impacted by the canopy cover of dense forests or surrounding mountain terrain. Therefore, using a buffer multiplier of 10 (e.g., 10$/sigma$ if the RMSE is taken to be equal to $/sigma$) is an approach that can help reduce false negatives (a point improperly identified as not belonging to any line) at the expense of increasing false positives (a point improperly identified as belonging to a line).

The RMSE, or even 10$/sigma$, may be too restrictive in cases where the RMSE is very small (e.g., << 1 m). In addition to the uncertainty in the point location, there is also uncertainty in the hydrography linework. The linework uncertainty varies based on which hydrography dataset is used due to the different spatial resolutions used to map the hydrography datasets. For example, NHD was mapped at 1:24,000 or larger scale, whereas NHDPlusV2 was mapped at 1:100,000 scale. Additionally, the PROSPER models using these observations resolve a 10-30 m spatial resolution (reach length, or grid cell size). Therefore, the linework uncertainty can be lumped with the point uncertainty by setting a lower limit on the search radius to about 10 m.

The minimum allowable buffer distance can be set based on whether stream name information (for the "name_match" method) is either present (`yes_stream_name_min_buffer`) or not present (`no_stream_name_min_buffer`). The minimum buffer distance's dependence on stream name information allows for prioritizing the name matching method when a stream name has been recorded over the horizontal accuracy of the point's GPS location by setting `yes_stream_name_min_buffer` > `no_stream_name_min_buffer`. Finally, the `max_buffer_distance` parameter sets the maximum allowable buffer distance.

![Map of FLOwPER observation points [@Heaston:2025a; @Heaston:2025b; @Heaston:2024; @Dunn:2023; @Heaston:2022; @Heaston:2020] snapped to the NHDPlusHR hydrography dataset [@USGS:2025b]. The uppermost FLOwPER point snapped to a flowline based on a successful name match. The lowermost FLOwPER point did not snap to any flowline since its search buffer region did not contain any flowlines.\label{fig:gps_rmse_fig}]{gps_rmse_fig.jpeg}

### Example usage

The Python code below shows the use of the FlowperLink class and methods to snap FLOwPER observations to an NHDPlusHR hydrography dataset using a combination of name matching, horizontal accuracy, and tributary junction information present in the FLOwPER observation record. An instance of the FlowperLink class is initiated by passing the input datasets, specifying field names within those datasets, and setting up distance thresholds for the various snapping rules. Calling `hydrolink_method` on this object executes the snapping functions and outputs the result to a file. The last two lines in this example write out connecting lines to help visualize the snap distance and path for each input point, and the buffer region around each input point.

```python
from flowperlink import FlowperLink

nhdplushr = FlowperLink(points = 'FLOwPER_points.shp',  
                        flowlines = 'NHDPLUS_H_1701_HU4_GPKG_preprocessed.gpkg',
                        points_identifier = 'GlobalID', 
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

A separate module was developed to snap FLOwPER and other streamflow permanence observations [@McShane:2017] to a high spatial resolution terrain derived hydrography from TerrainWorks (https://terrainworks.com/; [@ODFW:2024; @TerrainWorks:2014]). This hydrography necessitated its own class, TerrainWorksLink, due to the representation of flow paths not as vector lines, but as points derived from gridded data. Each point in this hydrography dataset represents the lower left coordinate of a 2-m grid cell identified as a flow path from a flow accumulation raster \autoref{fig:tw_fig}. The TerrainWorks hydrography dataset does not include stream names nor the fields required to associate tributary junction information from FLOwPER observations and therefore only snaps points informed by their geographic horizontal accuracy. Methods within the TerrainWorksLink class follow the same steps as in the FlowperLink class for hydrolinking using horizontal accuracy information, albeit for the point data of the TerrainWorks hydrography. There is an additional option to offset the snapping locations from the input hydrography points. A tuple of x (Eastings) and y (Northings) values in meters (positive or negative) can be passed to the `flowline_grid_offsets` parameter to offset the snapping location. Options for "name_match" and "trib_jcn" methods are included for future development but will currently raise a NotImplementedError if selected.

![Map of FLOwPER observation points [@Heaston:2025a; @Heaston:2025b; @Heaston:2024; @Dunn:2023; @Heaston:2022; @Heaston:2020] snapped to a 2-m gridded TerrainWorks hydrography dataset [@ODFW:2024; @TerrainWorks:2014], demonstrating the dynamic search buffer distances used based on the observations’ GPS horizontal accuracy.\label{fig:tw_fig}]{tw_fig.jpeg}

### Example usage

The Python code below shows the use of the TerrainWorksLink class and methods to snap FLOwPER observations to a TerrainWorks hydrography dataset using the closest distance method. Note the additional `flowline_grid_offsets` parameter. In this example, the TerrainWorks hydrography points represent the lower left corner of 2-m grid cells. To snap the FLOwPER point to grid cell centers rather than the corners, a tuple of (1, 1) is passed to this parameter to shift the snap location by 1-m east and 1-m north.

```python
tw = TerrainWorksLink(points = 'FLOwPER_points.shp',  
                      flowlines = 'Nodes_UpperDeschutes.gdb', 
                      points_identifier='OBJECTID',
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

## Handling FLOwPER data

The `download_flowper.py` Python script is used to download FLOwPER data from USGS ScienceBase [@USGS:2025a]. The script is run from the command line or through a wrapper Bash script. The script parses user input specifying the download directory and optional input comma-separated value (CSV) file path, then executes the download process based on the provided parameters. The optional CSV file input argument allows users to specify which datasets to download through a CSV file containing ScienceBase data release IDs, the unique identifiers used within ScienceBase (refer to the sciencebasepy documentation for details about interacting with ScienceBase; [@Long:2023]). If no datasets are specified by the user, a default list of FLOwPER data releases (presented in CSV form below) from 2019-2023 is downloaded instead [@Heaston:2025a; @Heaston:2025b; @Heaston:2024; @Dunn:2023; @Heaston:2022; @Heaston:2020].

```
name        , sb_id
flowper_2019, 5edd08c982ce7e579c6e48db
flowper_2020, 61b79dddd34e78124560f8d5
flowper_2021, 6627e820d34ea70bd5f033f9
flowper_2022, 678fd928d34e28977994d0aa
flowper_2023, 67e1d064d34ee7f142216699
rainier     , 62674439d34e76103cce59d4
```

The `load_flowper_data` function reads this CSV file if provided or returns the default list. The `download_flowper` function then establishes a session with ScienceBase using the sciencebasepy package [@Long:2023], checks for the existence of a specified download directory (creating it if necessary), and iterates through the ScienceBase IDs to download relevant files (while filtering out photographs and metadata). The `download_flowper` function also handles the extraction of the downloaded ZIP files.

The handling of FLOwPER datasets is continued with the `merge_flowper.py` Python script. Similar to `download_flowper.py`, the merging script is designed to be run in the command line or through a wrapper Bash script. The script merges multiple FLOwPER shapefiles into a single Geopandas GeoDataFrame [@Bossche:2024] and exports the merged product to a single file (any file format supported by `geopandas.to_file`, e.g. shapefile, geodatabase, geopackage). All input files are reprojected into a common coordinate reference system specified by the user, or EPSG 4326 if none is provided.

# Acknowledgements

The developers of FLOwPERlink would like to thank Nathan Pasley and Matt Barker of the U.S. Geological Survey for their helpful comments and initial reviews of the software and this paper. The developers also thank Susan Wherry of the U.S. Geological Survey for informing the development of the download script. The development of FLOwPERlink was funded by the U.S. Geological Survey Water Availability and Use Program.

# References