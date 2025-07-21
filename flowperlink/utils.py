from pathlib import Path
import pandas as pd
import geopandas as gpd
import fiona
from typing import Union, List

def merge_flowlines(DIR: str, layers_of_interest: list = [None], output_filepath: str = None, merge_layers: bool = False):
    '''Merges separate NHD or NHDPLUS geopackage or geodatabase files for a given layer into a new geopackage file.
    
    Parameters
    ----------
    DIR: str
        Path to top level directory of NHD or NHDPLUS data
    layers_of_interest: list
        List of layer names within NHD or NHDPLUS geopackage or geodatabase files to merge
        if None, then all layers available will be merged, defaults to None
    output_filepath: str
        Filepath to save merged geopackage to, defaults to {DIR}/merged.gpkg
    merge_layers: bool
        Set to True to ignore separate layer names in each file and merge all layers together from all files
        Set to False to merge layer by layer for all files
    
    Example
    -------
    # Directory containing NHDPlus geopackage or geodatabase files
    NHDPLUS_DIR = "D:/NHDPLUS/"
    # Layers to load from the NHD geopackage or geodatabase files
    layers_of_interest = ["NHDFlowline", "NHDPlusFlow", "NHDPlusFlowlineVAA"]
    # Merge files
    merge_flowlines(NHDPLUS_DIR, layers_of_interest, "D:/NHDPLUS/NHDPLUS_merged.gpkg")
    
    '''

    # if user set layers_of_interest = None, make this an iterable
    if layers_of_interest == None:
        layers_of_interest = [None]

    # if no output_filepath specified, set to default
    if not Path(output_filepath):
        output_filepath = Path(DIR, f'merged.gpkg')

    # check if a merged NHD or NHDPLUS layer geopackage or geodatabase file exists in this directory
    if not Path(output_filepath).is_file():

        # get a list of gpkgs from our DIR (searches recursively into DIR)
        files = [file for file in Path(DIR).rglob('*.gpkg')]

        # initialize an empty list to hold GeoDataFrames
        gdfs = []

        # initialize a dictionary to hold GeoDataFrames for each layer
        layer_dict = {}
        # initialize a dictionary to hold GeoDataFrames regardless of layer
        all_gdfs = []

        # loop through each file and each layer within the file
        for i, file_path in enumerate(files):
            
            print(f'reading file {i+1} of {len(files)}: {file_path}')

            # get all layers in the geopackage or geodatabase
            layers = fiona.listlayers(file_path)
            print(layers)

            for j, layer in enumerate(layers):

                # only read in layers_of_interest if specified above
                if (layer in layers_of_interest) | (layers_of_interest == [None]):

                    print(f'\treading layer {j+1} of {len(layers)}: {layer}')

                    # fead each layer into a GeoDataFrame
                    gdf = gpd.read_file(file_path, layer=layer)
                    
                    if merge_layers:
                        # if we're merging all layers, just add all gdfs to a list
                        all_gdfs.append(gdf)
                    else:
                        # otherwise, create layer item in dict if it's not yet there
                        if layer not in layer_dict:
                            layer_dict[layer] = [] 
                        # append geodataframe to layer list for merging layer by layer 
                        layer_dict[layer].append(gdf)

        if merge_layers:
            print(f'')

            print(f'merging all {len(all_gdfs)} dataframes...')

            # Set all to same CRS
            all_gdfs = align_crs(all_gdfs, 'EPSG:4269')
            # Concatenate all GeoDataFrames for the current layer
            merged_gdf = pd.concat(all_gdfs, ignore_index=True, sort=False)
            
            # save to the new geopackage
            print(f'saving merged layers to {output_filepath}')
            merged_gdf.to_file(output_filepath, layer='merged_flowlines')

        else:
            # merge and save each layer
            for layer, gdfs in layer_dict.items():

                print(f'merging all {len(gdfs)} {layer} dataframes...')

                # Set all to same CRS
                gdfs = align_crs(gdfs, 'EPSG:4269')
                # Concatenate all GeoDataFrames for the current layer
                merged_gdf = pd.concat(gdfs, ignore_index=True, sort=False)

                print(f'saving merged {layer} to {output_filepath}')
                
                if type(merged_gdf) == pd.DataFrame:
                    # if we have a pandas dataframe (layers without geometry)
                    # convert to geodataframe with geometry=None
                    merged_gdf = gpd.GeoDataFrame(merged_gdf, geometry=None)
                
                # save to the new geopackage
                merged_gdf.to_file(output_filepath, layer=layer)

        print(f"Created {output_filepath}.")

    else:
        
        print(f"File {output_filepath} already exists.")


def align_crs(gdf_list: List[gpd.GeoDataFrame], to_crs: Union[bool, str]=False):
    """
    Align the Coordinate Reference Systems (CRS) of a list of GeoDataFrames.

    This function checks if all GeoDataFrames in the provided list have the same CRS.
    If they do, it prints a message indicating that. If they do not, it can transform
    all GeoDataFrames to a specified CRS or to the CRS of the first GeoDataFrame in the list.

    Parameters
    ----------
    gdf_list : List[gpd.GeoDataFrame]
        A list of GeoDataFrames whose CRS needs to be aligned.

    to_crs : Union[bool, str], optional
        The target CRS to transform the GeoDataFrames to. 
        If set to `False`, all GeoDataFrames will be transformed to the CRS of the first GeoDataFrame.
        If a string representing a valid CRS (e.g., 'EPSG:4326') is provided, 
        all GeoDataFrames will be transformed to that CRS. The default is `False`.

    Returns
    -------
    List[gpd.GeoDataFrame]
        The list of GeoDataFrames with aligned CRS.

    Notes
    -----
    - If all GeoDataFrames have the same CRS, no transformation is performed unless `to_crs` is specified.
    - Ensure that the input GeoDataFrames are valid and have a defined CRS before calling this function.
    """
    crs_list = [gdf.crs for gdf in gdf_list]
    if all(crs == crs_list[0] for crs in crs_list):
        print(f'All geodataframes have the same crs: {crs_list[0]}')
        if to_crs != False:
            print(f'Transforming all geodataframes CRS to: {to_crs}')
            for i, gdf in enumerate(gdf_list):
                gdf_list[i] = gdf.to_crs(to_crs)
    else:
        print(f'Not all geodataframes have the same crs')
        if to_crs != False:
            print(f'Transforming all geodataframes CRS to: {to_crs}')
            for i, gdf in enumerate(gdf_list):
                gdf_list[i] = gdf.to_crs(to_crs)
        else:
            print(f'Transforming all geodataframe CRS to: {crs_list[0]}')
            for i, gdf in enumerate(gdf_list):
                gdf_list[i] = gdf.to_crs(crs_list[0])
    return gdf_list
            