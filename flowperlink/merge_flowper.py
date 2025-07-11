"""Script for merging FLOwPER shapefiles."""

import argparse
from pathlib import Path
import geopandas as gpd
import pandas as pd
from typing import Union

def merge_flowper(input_directory: str, output_filepath: str, crs: Union[int, str]):
    """
    This function merges multiple FLOwPER shapefiles together.

    Parameters
    ----------
    input_directory : str
        The path to a directory containing all input shapefiles. Shapefiles located in subdirectories will also be included.
    output_filepath : str
        The path and filename to save the merged output to.
    crs : int, str
        A common CRS to project all input datasets into. This will be the CRS of the output file.
    """

    # get a list of shapefiles from our DATA_DIR (searches recursively into DATA_DIR)
    files = [file for file in Path(input_directory).rglob('*.shp')]
    print(f'Merging {len(files)} files found in {input_directory}:')
    _ = [print(f'\t{str(file)}') for file in files]

    # read all files into a list of geodataframes, project to common CRS
    gdfs = [gpd.read_file(file_path).to_crs(crs) for file_path in files]

    # merge all together
    merged_gdf = pd.concat(gdfs, ignore_index=True, sort=False)

    # save to a new shapefile
    print(f'Saving to {output_filepath}')
    merged_gdf.to_file(Path(output_filepath), driver='ESRI Shapefile')
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Merge FLOwPER shapefiles.')
    
    # define input arguments
    parser.add_argument('--input_directory', type=str, required=True, 
                        help='The directory containing FLOwPER shapefiles to process')
    parser.add_argument('--output_filepath', type=str, required=False, default='FLOwPER_merged.shp', 
                        help='Path and file name for new shapefile to save the merged FLOwPER data to (default saves to current working directory as FLOwPER_merged.shp)')
    parser.add_argument('--crs', type=str, required=False, default=4326, 
                        help='Coordinate Reference System to use for final output (defaults to EPSG:4326)')
    args = parser.parse_args()
    
    # run the merge_flowper function with user-specified inputs
    merge_flowper(args.input_directory, 
                  args.output_filepath, 
                  args.crs)