"""Script for preprocessing flowlines to derive tributary junction information."""

"""
FLOwPER observations include a tributary junction field, which could be used to help narrow down
the possible candidate flowlines (and portions of flowlines) that a point should be snapped to.
The flowline data will need to be preprocessed to add tributary junction information that can then
be compared with the FLOwPER data. First, closed loops will need to be removed. Then, at each node,
the flowlines are cut "trib_jcn_dist" from that flowline (if a line is < 2*trib_jcn_dist, then the 
line is split in half). The downstream flowline at this node is classified as the downstream 
mainstem. For the upstream lines, their total upstream length from this node is calculated and 
compared. The line with the longest upstream length is classified as the upstream mainstem, and all 
others as tributaries. The ratio between the upstream lengths will also be computed, if this ratio 
is close to 1 then there is more ambiguity in identifying which upstream line is the mainstem or 
the tributary.
"""

import argparse
import geopandas as gpd
import pandas as pd
import numpy as np
import networkx as nx
from shapely.ops import split, snap
from functools import partial
from typing import Union
from pathlib import Path

def clean_flowlines(flowlines_filepath: Union[str, Path], 
                    flowline_layer: str, 
                    flow_layer: str, 
                    flowlineVAA_layer: str, 
                    divergence_field: str, 
                    flowlines_identifier: str, 
                    flowlinesVAA_identifier: str, 
                    from_identifier: str, 
                    to_identifier: str):
    
    """Reads in NHD or NHDPlus flowlines and associated tables, merges and filters data.
    
    Parameters
    ----------
    flowlines_filepath: str, Path
        Path to a geopackage file containing NHD hydrography
    flowline_layer: str
        Name of layer within geopackage containing the NHD flowlines
    flow_layer: str
        Name of layer within geopackage containing the NHD flow information
    flowlineVAA_layer: str
        Name of layer within geopackage containing the flowline value added attributes information
    divergence_field: str
        Field or column name within the flowlineVAA_layer containing divergence info, for removing closed loops in flowlines
    flowlines_identifier: str
        Field or column name within flowlines layer containing a unique identifier for each flowline
    flowlinesVAA_identifier: str
        Field or column name within flowlineVAA layer containing a unique identifier for each flowline
    from_identifier: str
        Field or column name within flow layer containing flow from information
    to_identifier: str
        Field or column name within flow layer containing flow to informaiton
    
    Returns
    -------
    gdf: gpd.GeoDataFrame
        A merged geopandas GeoDataFrame containing flowlines, flow to and from fields, and reprojected to local UTM coordinates

    """
    
    flowlines_gdf = gpd.read_file(flowlines_filepath, layer=flowline_layer)
    # set datatype for the flowlines_identifier column to string
    flowlines_gdf[flowlines_identifier] = flowlines_gdf[flowlines_identifier].astype('str')

    flow_gdf = gpd.read_file(flowlines_filepath, layer=flow_layer)
    
    if flowlineVAA_layer != None:
        flowlinesVAA_gdf = gpd.read_file(flowlines_filepath, layer=flowlineVAA_layer)
        #if divergence_field != None:
        #    # remove where Divergence = 2 (minor path divergence)
        #    flowlinesVAA_gdf = flowlinesVAA_gdf.loc[flowlinesVAA_gdf[divergence_field] != 2]
        # remove flowlines that were a minor path divergence
        #gdf = flowlines_gdf[flowlines_gdf[flowlinesVAA_identifier].isin(set(flowlinesVAA_gdf[flowlinesVAA_identifier].unique()))]
        # some lines don't have any information in the flowlinesVAA table. Therefore the above line was removing them from the gdf, which we don't want
        # set datatype for the flowlines_identifier column to string
        flowlinesVAA_gdf[flowlinesVAA_identifier] = flowlinesVAA_gdf[flowlinesVAA_identifier].astype('str')
        # join flowlinesVAA to gdf on perm IDs
        gdf = flowlines_gdf.merge(flowlinesVAA_gdf, left_on=flowlines_identifier, right_on=flowlinesVAA_identifier, how='left')
    else:
        gdf = flowlines_gdf

    # join flow to new gdf on perm IDs
    gdf = gdf.merge(flow_gdf[[from_identifier, to_identifier]], left_on=flowlines_identifier, right_on=from_identifier, how='left')#[[flowlines_identifier, from_identifier, to_identifier, 'geometry']]

    # reproject to local UTM
    gdf = gdf.to_crs(gdf.estimate_utm_crs())
    # compute line lengths
    gdf['length'] = gdf.length

    return gdf


def compute_upstream_network_length(gdf: gpd.GeoDataFrame, 
                                    flowlines_identifier: str, 
                                    from_identifier: str, 
                                    to_identifier: str, 
                                    measure: str  = 'length'):
    
    """For each flowline, compute the sum of the lengths of all upstream contributing flowlines.

    This function constructs a directed graph from the flowlines in the provided GeoDataFrame 
    and calculates the total length of all upstream flowlines contributing to each flowline.
    
    Parameters
    ----------
    gdf: gpd.GeoDataFrame
        A merged geopandas GeoDataFrame containing flowlines, flow to and from fields, and reprojected to local UTM coordinates
    flowlines_identifier: str
        Field or column name within flowlines layer containing a unique identifier for each flowline
    from_identifier: str
        Field or column name within flow layer containing flow from information
    to_identifier: str
        Field or column name within flow layer containing flow to informaiton
    measure: str
        Field or column name with flowline length, defaults to 'length'
    
    Returns
    -------
    gdf: gpd.GeoDataFram
        Geopandas GeoDataFrame with origial data, now also containing new columns with upstream 
        lengths, number of upstream and downstream nodes computed for each flowline

    Notes
    ------------
    The computation may take significant time depending on the size of the input GeoDataFrame. 
    For example, it may take about 5 minutes with 200,000 flowlines across a HUC4.
    """

    # Create directed graph
    graph = nx.from_pandas_edgelist(gdf, source=from_identifier, target=to_identifier, edge_attr=True, create_using=nx.DiGraph())

    # Initialize dict for upstream lengths
    upstream_lengths = {}
    # Initialize dict for max_upstream_lengths, lengths, mainstems
    max_upstream_lengths = {}
    lengths = {}
    upstream_mainstems = {}
    # Initialize dicts for number of nodes to detect ends
    n_upstream_nodes = {}
    n_downstream_nodes = {}

    # Here is an opportunity to make the code more efficient, avoid looping through
    # all nodes or at least distribute across CPUs in parallel since each
    # node is computed independently. Alternatively, start at the 1st order
    # streams and traverse the network downstream summing up lengths could cut down
    # on total compute time.
    for node in graph.nodes():
        # initialize total upstream length to be 0
        total_upstream_length = 0
        # get all upstream and downstream nodes
        upstream_nodes = nx.ancestors(graph, node)
        downstream_nodes = nx.descendants(graph, node)
        # sum up the lengths for all the edges between upstream nodes
        for upstream_node in upstream_nodes:
            for _, _, data in graph.out_edges(upstream_node, data=True):
                lengths[upstream_node] = data[measure] # get the lengths in a dict to use later
                total_upstream_length += data[measure] # add this length to running total of upstream lengths
        # print(node)
        # print(list(graph.successors(node)))
        # (_, _, data) = graph.out_edges(node, data=True)
        
        upstream_lengths[node] = total_upstream_length #+ data[measure] # add the length of this line to the total
        n_upstream_nodes[node] = len(upstream_nodes)
        n_downstream_nodes[node] = len(downstream_nodes)
        
    # compare upstream line lengths at each node to find mainstem upstream
    for node in graph.nodes():
        upstream_nodes = graph.predecessors(node) # predecessor node just one step upstream
        # initialize max_upstream_length, and max_upstream_node
        this_max_upstream_length = 0
        this_max_upstream_node = None
        for upstream_node in upstream_nodes:
            if upstream_lengths[upstream_node] > this_max_upstream_length:
                this_max_upstream_length = upstream_lengths[upstream_node]
                this_max_upstream_node = upstream_node
        max_upstream_lengths[node] = this_max_upstream_length
        # if n_upstream_nodes == 0, we are at the upstream-most junction
        # therefore the line with the longest length at each junction we'll call the "mainstem"
        # This needs some work still, right now none of the upstream-most lines are labeled "mainstem"
        # when we would want the longest line at the last node labeled "mainstem"
        if n_upstream_nodes[node] == 0:
            this_max_length = 0
            for upstream_node in upstream_nodes:
                if lengths[upstream_node] > this_max_length:
                    this_max_upstream_node = upstream_node
        # label the mainstem here
        upstream_mainstems[this_max_upstream_node] = 'mainstem'


    # Map dict of computed lengths to the gdf, producing a DataSeries for a new column
    upstream_lengths_ds = gdf[flowlines_identifier].map(upstream_lengths) + gdf[measure]
    gdf['total_upstream_length'] = upstream_lengths_ds # add new column

    # Map the dict of max upstream lengths (minus the length of adjacent lines)
    max_upstream_lengths_ds = gdf[flowlines_identifier].map(max_upstream_lengths)
    gdf['max_upstream_length'] = max_upstream_lengths_ds
    upstream_mainstems_ds = gdf[flowlines_identifier].map(upstream_mainstems)
    gdf['mainstem_flag'] = upstream_mainstems_ds
    gdf['mainstem_flag'] = gdf['mainstem_flag'].fillna('tributary') # where not mainstem, label tributary
    
    # Map the dicts of the number of upstream and downstream nodes and add new columns for each
    n_upstream_nodes_ds = gdf[flowlines_identifier].map(n_upstream_nodes)
    gdf['n_upstream_nodes'] = n_upstream_nodes_ds
    n_downstream_nodes_ds = gdf[flowlines_identifier].map(n_downstream_nodes)
    gdf['n_downstream_nodes'] = n_downstream_nodes_ds
    
    return gdf


def split_flowline(row, 
                   trib_jcn_dist: float = 100, 
                   tolerance: float = 1):
    """Split flowlines to identify segments upstream or downstream of junctions.
    
    This function takes a flowline geometry and splits it into segments based on specified 
    distances from junctions. It returns the geometries of the upstream and downstream 
    segments, as well as any middle segment if applicable.

    Parameters
    ----------
    row : GeoDataFrame row
        A row from a GeoDataFrame containing the flowline geometry to be split. 
        It must have a 'geometry' key with a LineString or MultiLineString geometry.
    trib_jcn_dist : float, optional
        The distance from the junction at which to split the flowline. 
        Default is 100 units (e.g., meters).
    tolerance : float, optional
        The tolerance used when snapping points to the flowline. 
        This allows for some wiggle room when interpolating points on the line. 
        Default is 1 unit (e.g., meter).
    
    Returns
    -------
    tuple
        A tuple containing three elements:
        - line_dwnstrm_of_jcn_geom : geometry or None
            The geometry of the segment downstream of a junction, or None if at upstream-most line or the split was not valid.
        - line_upstrm_of_jcn_geom : geometry or None
            The geometry of the segment upstream of a junction, or None if at downstream-most line or the split was not valid.
        - line_not_near_jcn_geom : geometry or None
            The geometry of the middle segment, or None if the line was not long enough 
            to have a middle segment outside of `trib_jcn_dist` from each end.

    Notes
    -----
    If the flowline is shorter than twice the `trib_jcn_dist`, the function will split 
    the line in half. If the line is long enough, it will create segments of length 
    `trib_jcn_dist` from both ends and attempt to return the middle segment if applicable.

    """
    #tolerance because interpolate doesn't fall exactly on the line, we need some wiggle room to snap

    line = row['geometry']
    # if the line isn't long enough, split in half
    if line.length <= 2*trib_jcn_dist:
        pt_dwnstrm_of_jcn = line.interpolate(line.length / 2, normalized=False) # positive value, measuring from upstream end of line
        pt_upstrm_of_jcn = line.interpolate(-line.length / 2, normalized=False) # negative value, measuring from downstream end of line
    # otherwise, grab upstream and downstream reaches of length trib_jcn_dist
    else:
        pt_dwnstrm_of_jcn = line.interpolate(trib_jcn_dist, normalized=False) # positive value, measuring from upstream end of line
        pt_upstrm_of_jcn = line.interpolate(-trib_jcn_dist, normalized=False) # negative value, measuring from downstream end of line

    # split the original line at these points upstream and downstream of junctions
    line_dwnstrm_of_jcn = split(snap(line, pt_dwnstrm_of_jcn, tolerance), pt_dwnstrm_of_jcn)
    line_upstrm_of_jcn = split(snap(line, pt_upstrm_of_jcn, tolerance), pt_upstrm_of_jcn)

    # check for each if a valid geometry split the flowline in two, otherwise return None
    if len(line_upstrm_of_jcn.geoms) == 2:
        line_dwnstrm_of_jcn_geom = line_dwnstrm_of_jcn.geoms[0] # grab the first, upstream-most, of the two lines
    else:
        line_dwnstrm_of_jcn_geom = None
    if len(line_upstrm_of_jcn.geoms) == 2:
        line_upstrm_of_jcn_geom = line_upstrm_of_jcn.geoms[1] # grab the second, downstream-most, of the two lines
    else:
        line_upstrm_of_jcn_geom = None

    # if the line was long enough to have a middle segment outside of trib_jcn_dist from each end
    # we want to return that segment too, otherwise return None
    if line.length > 2*trib_jcn_dist:
        # can either split remainder from line_dwnstrm_of_jcn at point pt_upstrm_of_jcn or remainder from line_upstrm_of_jcn at point pt_dwnstrm_of_jcn
        try:
            line_not_near_jcn_geom = split(snap(line_dwnstrm_of_jcn.geoms[1], pt_upstrm_of_jcn, tolerance), pt_upstrm_of_jcn).geoms[0]
        except:
            line_not_near_jcn_geom = split(snap(line_upstrm_of_jcn.geoms[0], pt_dwnstrm_of_jcn, tolerance), pt_dwnstrm_of_jcn).geoms[1]
    else:
        line_not_near_jcn_geom = None


    if row['n_upstream_nodes'] == 0 and row['n_downstream_nodes'] > 1:
        # only split into segments upstrm_of_jcn and not_near_jcn (we are at an upstream-most flowline, i.e. 1st order)
        # remove line_dwnstrm_of_jcn_geom if it was created
        line_dwnstrm_of_jcn_geom = None
        # set line_not_near_jcn_geom to be the upstream segment from the split with line_upstrm_of_jcn
        line_not_near_jcn_geom = line_upstrm_of_jcn.geoms[0]
    elif row['n_upstream_nodes'] > 0 and row['n_downstream_nodes'] == 1:
        # only split into segments dwnstrm_of_jcn and not_near_jcn (we are at a downstream-most flowline, i.e. highest order)
        # remove line_upstrm_of_jcn_geom if it was created
        line_upstrm_of_jcn_geom = None
        # set line_not_near_jcn_geom to be the downstream segment from the split with line_dwnstrm_of_jcn
        line_not_near_jcn_geom = line_dwnstrm_of_jcn.geoms[1]
    elif row['n_upstream_nodes'] == 0 and row['n_downstream_nodes'] == 1:
        # this flowline is isolated, skip and flag
        # remove line_dwnstrm_of_jcn_geom and line_upstrm_of_jcn_geom
        line_dwnstrm_of_jcn_geom = None
        line_upstrm_of_jcn_geom = None
        # set line_not_near_jcn_geom to be the whole flowline
        line_not_near_jcn_geom = line
    elif row['n_upstream_nodes'] > 0 and row['n_downstream_nodes'] > 1:
        # flowline is within the stream network, not at one end or the other, split into dwnstrm_of_jcn, upstrm_of_jcn, and not_near_jcn
        # change nothing, we've already split this correctly above
        pass

    return line_dwnstrm_of_jcn_geom, line_upstrm_of_jcn_geom, line_not_near_jcn_geom


def get_junction_flowlines(gdf, 
                           flowlines_identifier: str, 
                           trib_jcn_dist: float = 150, 
                           output_separate_files: bool = False):
    
    """Identify and split flowlines near tributary junctions.

    This function uses NetworkX to find and split flowlines that are near tributary junctions 
    based on a specified distance. It creates separate geometries for flowlines that are 
    downstream of junctions, upstream of junctions, and not near junctions. Optionally, it can 
    output separate files for each type of flowline segment for debugging purposes.

    Parameters
    ----------
    gdf : gpd.GeoDataFrame
        A GeoDataFrame containing flowlines with geometries and identifiers. It must include 
        the geometry column and any necessary identifiers.
    flowlines_identifier : str
        The name of the column in the GeoDataFrame that contains unique identifiers for each flowline.
    trib_jcn_dist : float, optional
        The distance from the junction at which to consider flowlines as being near a junction. 
        Default is 150 units (e.g., meters).
    output_separate_files : bool, optional
        If True, the function will output separate shapefiles for each type of flowline segment 
        (downstream, upstream, and not near junctions) for debugging purposes. Default is False.

    Returns
    -------
    gpd.GeoDataFrame
        A GeoDataFrame containing the split flowlines with updated identifiers and geometries 
        for each segment type, including a new column indicating the type of junction proximity.

    Notes
    -----
    The function applies the `split_flowline` function to each row of the input GeoDataFrame 
    to generate the three types of flowline segments. The resulting segments are then merged 
    into a single GeoDataFrame with appropriate labels and identifiers.

    Examples
    --------
    >>> junctions_gdf = get_junction_flowlines(flowlines_gdf, 'flowline_id', trib_jcn_dist=100)
    
    """

    # use networkx to find and split flowlines near tributary junctions
    gdf[['line_dwnstrm_of_jcn', 'line_upstrm_of_jcn', 'line_not_near_jcn']] = gdf.apply(partial(split_flowline, trib_jcn_dist=trib_jcn_dist), axis=1, result_type='expand')

    if output_separate_files == True:
        # optional, for debugging, output separate files for flowline sub-reaches
        gdf.set_geometry('line_dwnstrm_of_jcn').drop(columns=['geometry', 'line_upstrm_of_jcn', 'line_not_near_jcn']).to_file('lines_dwnstrm_of_jcn.gpkg')
        gdf.set_geometry('line_upstrm_of_jcn').drop(columns=['geometry', 'line_dwnstrm_of_jcn', 'line_not_near_jcn']).to_file('lines_upstrm_of_jcn.gpkg')
        gdf.set_geometry('line_not_near_jcn').drop(columns=['geometry', 'line_dwnstrm_of_jcn', 'line_upstrm_of_jcn']).to_file('lines_not_near_jcn.gpkg')

    # for the three line segments, create three new geodataframes and assign geometry column for each line segment
    gdf_a = gdf.drop(columns=['geometry']).rename(columns={'line_dwnstrm_of_jcn': 'geometry'})
    gdf_b = gdf.drop(columns=['geometry']).rename(columns={'line_upstrm_of_jcn': 'geometry'})
    gdf_c = gdf.drop(columns=['geometry']).rename(columns={'line_not_near_jcn': 'geometry'})

    # Add trib junction labels
    gdf_a['trib_jcn'] = 'downstream of junction'
    gdf_b['trib_jcn'] = 'upstream of junction'
    gdf_c['trib_jcn'] = 'not near junction'

    # Append suffixes to flowlines_identifier for each
    gdf_a[flowlines_identifier] = gdf_a[flowlines_identifier].astype(str) + '_dwnstrm_of_jcn'
    gdf_b[flowlines_identifier] = gdf_b[flowlines_identifier].astype(str) + '_upstrm_of_jcn'
    gdf_c[flowlines_identifier] = gdf_c[flowlines_identifier].astype(str) + '_not_near_jcn'

    # Merge these three back together
    junctions_gdf = pd.concat([gdf_a, gdf_b, gdf_c], ignore_index=True).drop(columns=['line_dwnstrm_of_jcn', 'line_upstrm_of_jcn', 'line_not_near_jcn'])
    junctions_gdf = gpd.GeoDataFrame(junctions_gdf, geometry='geometry')

    return junctions_gdf

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Preprocess flowlines to derive tributary junction information.')
    
    # define input arguments
    parser.add_argument('--flowlines_filepath', type=str, required=True, 
                        help='Path to a geopackage file containing NHD hydrography')
    parser.add_argument('--flowline_layer', type=str, required=True, 
                        help='Name of layer within geopackage containing the NHD flowlines')
    parser.add_argument('--flowlines_identifier', type=str, required=True, 
                        help='Field or column name within flowlines layer containing a unique identifier for each flowline')
    parser.add_argument('--flow_layer', type=str, required=True, 
                        help='Name of layer within geopackage containing the NHD flow information')    
    parser.add_argument('--from_identifier', type=str, required=True, 
                        help='Field or column name within flow layer containing flow from information')
    parser.add_argument('--to_identifier', type=str, required=True, 
                        help='Field or column name within flow layer containing flow to informaiton')
    parser.add_argument('--flowlineVAA_layer', type=str, required=False, default=None, 
                    help='Name of layer within geopackage containing the flowline value added attributes information')
    parser.add_argument('--flowlinesVAA_identifier', type=str, required=False, default=None, 
                    help='Field or column name within flowlineVAA layer containing a unique identifier for each flowline')
    parser.add_argument('--divergence_field', type=str, required=False, default = None,
                    help='Field or column name within the flowlineVAA_layer containing divergence info, for removing closed loops in flowlines')
    parser.add_argument('--trib_jcn_dist', type=int, required=False, default=150,
                    help='The distance from the junction at which to consider flowlines as being near a junction. Default is 150 units (e.g., meters).')
    parser.add_argument('--output_filepath', type=str, required=False, default='segmented_flowlines.gpkg',
                    help='Path and file name to save the segmented flowlines. Default is segmented_flowlines.gpkg')
    args = parser.parse_args()
    
    # Read flowlines data and perform cleanup step
    print(f'Reading flowlines data from {args.flowlines_filepath} and performing cleanup step...')
    cleaned_gdf = clean_flowlines(args.flowlines_filepath, 
                                  args.flowline_layer, 
                                  args.flow_layer, 
                                  args.flowlineVAA_layer, 
                                  args.divergence_field, 
                                  args.flowlines_identifier, 
                                  args.flowlinesVAA_identifier, 
                                  args.from_identifier, 
                                  args.to_identifier)
    
    # compute upstream network lengths prior to segmenting flowlines 
    # (otherwise we have ~3x as many flowlines to compute over)
    print('Computing upstream network lengths...')
    cleaned_gdf_w_lengths = compute_upstream_network_length(cleaned_gdf, 
                                                            args.flowlines_identifier, 
                                                            args.from_identifier, 
                                                            args.to_identifier)
    
    # create final geodataframe with all lines segmented
    print('Creating final geodataframe...')
    junctions_gdf = get_junction_flowlines(cleaned_gdf_w_lengths, 
                                           args.flowlines_identifier, 
                                           args.trib_jcn_dist, 
                                           output_separate_files=False)
    
    # save to file
    print(f'Saving to file: {args.output_filepath}...')
    junctions_gdf.set_crs(cleaned_gdf_w_lengths.crs).to_file(args.output_filepath)

    print('Done.')