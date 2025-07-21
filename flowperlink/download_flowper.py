"""Script for downloading FLOwPER data from ScienceBase."""

import argparse
from pathlib import Path
import pandas as pd
import sciencebasepy as sb
import zipfile
from typing import Union

def load_flowper_data(input_filepath: Union[str, Path] = None, name_col: str = 'name', id_col: str = 'sb_id'):

    """
    Specify the FLOwPER datasets to download from ScienceBase from a csv file, or default list of datasets.
    
    Parameters
    ----------
    input_filepath : Union[str, path], optional
        Path to a csv file containing the names and ScienceBase IDs of FLOwPER data. If None, default data will be used
    name_col : str, optional
        The name of the column to use for the ScienceBase item names. Default is 'name'
    id_col : str, optional
        The name of the column to use for the ScienceBase item IDs. Default is 'sb_id'

    Returns
    -------
    flowper_ids: pd.DataFrame
        A DataFrame containing the names and ScienceBase IDs for FLOwPER datasets
    
    """

    # Define ScienceBase IDs for pages with Flowper data
    flowper_data = [
        {'name': "flowper_2019", 'sb_id': "5edd08c982ce7e579c6e48db"},
        {'name': "flowper_2020", 'sb_id': "61b79dddd34e78124560f8d5"},
        {'name': "flowper_2021", 'sb_id': "6627e820d34ea70bd5f033f9"},
        {'name': "flowper_2022", 'sb_id': "678fd928d34e28977994d0aa"},
        {'name': "flowper_2023", 'sb_id': "67e1d064d34ee7f142216699"},
        {'name': "rainier", 'sb_id': "62674439d34e76103cce59d4"}
    ]
    
    # If a filepath is provided, read the CSV file
    if input_filepath:
        flowper_ids = pd.read_csv(input_filepath, usecols=[name_col, id_col])
    else:
        # Create a DataFrame from the predefined data
        flowper_ids = pd.DataFrame(flowper_data)
    
    # Rename columns to the common format for download_flowper function
    flowper_ids = flowper_ids.rename(columns={'name': name_col, 'sb_id': id_col})
    
    return flowper_ids


def download_flowper(flowper_ids: pd.DataFrame, download_directory: Union[str, Path] = "./sb_data"):
    """
    Download FLOwPER data from ScienceBase

    Parameters
    ----------
    flowper_ids: pd.DataFrame
        A DataFrame containing the names and ScienceBase IDs for FLOwPER datasets
    download_directory : Union[str, Path]
        The path to a directory where all FLOwPER data will be downloaded
    """
    # Initialize a session (replace with your credentials if needed)
    session = sb.SbSession()

    # Check for download directory and create, if missing
    sb_data_dir = Path(download_directory)
    if not sb_data_dir.exists():
        sb_data_dir.mkdir(parents=True)

    # Download and unzip data for each ScienceBase ID
    for index, row in flowper_ids.iterrows():
        # ScienceBase item name and ID for this row
        sb_id = row['sb_id']
        sb_name = row['name']
        
        # create a local download directory for this item
        sb_dir = Path(sb_data_dir / sb_name)
        if not sb_dir.exists():
            sb_dir.mkdir(parents=True)
        
        # get this SB item and it's associated files
        sb_item = session.get_item(sb_id)
        sb_item_files = session.get_item_file_info(sb_item)
        sb_item_files_name = [file['name'] for file in sb_item_files]
        
        # filter out files with "photo" or "metadata" in their filenames
        is_photo = ["photo" in name.lower() for name in sb_item_files_name]
        is_metadata = ["metadata" in name.lower() for name in sb_item_files_name]
        is_data = [not (photo or metadata) for photo, metadata in zip(is_photo, is_metadata)]

        data_files = [sb_item_files[i] for i in range(len(is_data)) if is_data[i]]
        for file in data_files:
            # Download the file
            session.download_file(file['url'], file['name'], sb_dir)
            # Unzip the file if it's a zip file
            if file['name'].endswith('.zip'):
                print(f'extracting {file["name"]}')
                with zipfile.ZipFile(Path(sb_dir / file['name']), 'r') as zip_ref:
                    zip_ref.extractall(sb_dir)

    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Download FLOwPER data.')
    
    # define input arguments
    parser.add_argument('--download_directory', type=str, required=True, 
                        help='The directory where all FLOwPER data will be downloaded')
    parser.add_argument('--input_filepath', type=str, required=False, default=None, 
                        help='Path to a csv file containing the names and ScienceBase IDs of FLOwPER data. If None, default data will be used')
    parser.add_argument('--name_col', type=str, required=False, default='name', 
                        help='Column name within input_filepath with the ScienceBase item name')
    parser.add_argument('--id_col', type=str, required=False, default='sb_id', 
                        help='Column name within input_filepath with the ScienceBase item ID')
    args = parser.parse_args()
    
    # get the ScienceBase IDs for FLOwPER
    if args.input_filepath is not None:
        flowper_ids = load_flowper_data(args.input_filepath, args.name_col, args.id_col)
    else:
        flowper_ids = load_flowper_data()

    # download these FLOwPER data
    download_flowper(flowper_ids, args.download_directory)
 