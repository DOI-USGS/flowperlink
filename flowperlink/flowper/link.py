"""Allows for hydrolinking of FLOwPER point observations to a custom hydrography dataset.

Module that hydrolinks FLOwPER data to a custom hydrography dataset. Classes and methods are
designed to handle all input points together as a geopandas geodataframe. This is a child
class of the CustomHydrography class. Errors for individual points should be captured instead
of causing the process to fail. Currently only point data (not lines or polygons) can be
hydrolinked to flowlines (hydrolinking to waterbodies not yet implemented). The terms
hydrolink(ing), addressing, and snapping are used synonymously throughout this code and both
refer to making a relationship between spatial data and a stream network. This is similar to
the analogy of providing an address on a road network and provides locational context and
position within a stream network.

"""

# Import packages
from hydrolink.custom import CustomHydrography
from pathlib import Path
import geopandas as gpd
import pandas as pd
from typing import Union, Literal
############################################################################################
############################################################################################



class FlowperLink(CustomHydrography):
    """CustomHydrography child class for hydrolinking FLOwPER point data"""

    def hydrolink_method(self,
                         method: Literal['name_match', 'closest'] = 'name_match',
                         trib_jcn: str = 'TribJncTyp',
                         hydro_type: Literal['flowline', 'waterbody'] = 'flowline',
                         outfile_name: Union[str, Path] = 'custom_hydrolink_output.gpkg',
                         similarity_cutoff: float = 0.6):
        """Build hydrolinking pipeline based on specified method and hydro_type.

        Description
        ----------
        Strings together CustomHydrography methods into a pipeline that handles FLOwPER specific
        hydrolink use cases. This pipeline allows for a user to specify method of hydrolink and
        the target hydro_type. At completion data are exported to the outfile specified in "outfile_name".

        Parameters
        ----------
        method: {'name_match', 'closest'}, default 'name_match'
            Method for hydrolinking data. Supported methods are

            - ``'name_match'``: This default method hydrolinks data to the closest flowline feature with a name similarity
            that meets the specified similarity_cutoff. If no flowlines meet similarity cutoff the method hydrolinks
            data to the closest flowline feature.
            - ``'closest'``: This method hydrolinks data to the closest flowline feature.

        trib_jcn: str, default "TribJncTyp"
            Field or column name within FLOwPER data containing tributary junction information for each point. Defaults to
            "TribJncTyp". If trib_jcn is None, the tributary junction matching routine will be skipped. The values within the
            field specified here should be one of the following for FLOwPER V2 records: ['On tributary', 'On mainstem upstream',
            'On mainstem downstream', 'No Data']

        hydro_type: {'waterbody', 'flowline'}, default 'flowline'
            Type of features to hydrolink.

            - ``'flowline'``: This default feature type specifies flowline feature type of flowline.
            Flowline features represent water types such as streams, rivers, canals/ditches.
            Waterbodies also have line representations as flowline type.
            - ``'waterbody'``: Not implemented

        outfile_name: str
            Name and directory of gpkg output file.  default is 'custom_hydrolink_output.gpkg'.

        similarity_cutoff: float
            Values between 0 and 1.0, range of similarity between 0 representing no match to 1.0 being perfect match.
            Defaults to 0.6.

        """
        similarity_cutoff_max = 1.0
        similarity_cutoff_min = 0.6

        if 0.6 > similarity_cutoff or 1.0 < similarity_cutoff:
            self.message = f'Parameter similarity_cutoff should be a float between min: {similarity_cutoff_min} and max: {similarity_cutoff_max}.'
            self.error_handling()

        if self.status == 1:
            if hydro_type == 'waterbody':
                raise NotImplementedError
            if hydro_type == 'flowline':
                # steps to hydrolink points to flowlines:
                print('Create search region around input points...')
                self.buffer_points()
                print('Flowlines to polygons...')
                self.buffer_flowlines()
                print('Find flowlines within search regions...')
                self.intersect_points_flowlines()
                print('Hydrolink points to matching flowlines...')
                self.find_nearest_point()
                if trib_jcn is not None:
                    print('...use tributary junction information...')
                    self.trib_jcn_match()
                if method == 'name_match':
                    print('...use name match method...')
                    self.get_name_similarity()
                    self.select_closest_flowline_w_name_match(similarity_cutoff=similarity_cutoff)
                elif method == 'closest':
                    print('...use closest method...')
                    self.select_closest_flowline()
                print('Handling points that did not link to hydrography...')
                self.unlinked_points_handling()
                print('Writing output file...')
                self.write_hydrolink(outfile_name=outfile_name)
        else:
            print('Status = 0, writing output...')
            self.write_hydrolink(outfile_name=outfile_name)
        print('Done.')

    def trib_jcn_match(self):
        """Uses the FLOwPER tributary junction information to narrow down possible flowlines for snapping

        Description
        -----------
        This method processes the tributary junction information from FLOwPER records contained in the
        `hydrolinked_gdf` GeoDataFrame. It filters flowlines based on the type of tributary junction each
        flowline is associated with. The function uses the `TribJncTyp` field to determine how to filter
        the flowlines and updates the `trib_jnc_processing_message` field with relevant processing messages.

        The filtering logic is as follows:
        - If `TribJncTyp` is 'No Data', no rows are dropped.
        - If `TribJncTyp` is 'On tributary', only rows where `mainstem_flag` is 'tributary' are retained.
        - If `TribJncTyp` is 'On mainstem downstream', only rows where `mainstem_flag` is 'mainstem'
          and `trib_jcn` is 'downstream of junction' are retained.
        - If `TribJncTyp` is 'On mainstem upstream', only rows where `mainstem_flag` is 'mainstem'
          and `trib_jcn` is 'upstream of junction' are retained.
        - If none of the conditions are met, no rows are dropped, but a message is logged.

        Parameters
        ----------
        self : object
            The instance of the class containing the `hydrolinked_gdf` GeoDataFrame and the `source_id`
            attribute.

        Returns
        -------
        None
            This method modifies the `hydrolinked_gdf` GeoDataFrame in place and does not return any value.

        Notes
        -----
        - The method assumes that `hydrolinked_gdf` contains the necessary columns:
          `TribJncTyp`, `mainstem_flag`, and `trib_jcn`.
        - The `trib_jnc_processing_message` field is updated with messages indicating the processing
          status for each unique `global_id`.
        - If all rows are dropped for a given `global_id`, a message will be logged indicating that
          tributary junction matching failed, and no flowlines will be dropped
        """

        # create a field to store trib junction info processing messages per point
        self.hydrolinked_gdf['trib_jnc_processing_message'] = ''

        for global_id in self.hydrolinked_gdf[self.source_id].unique():
            # gdf points directly to this filtered self.hydrolinked_gdf, it is not a copy, this is just for ease of typing
            gdf = self.hydrolinked_gdf[self.hydrolinked_gdf[self.source_id] == global_id]
            # count number of flowlines we could potentially snap to (this is just the number of rows for each global_id)
            num_flowlines = len(gdf)

            # get the TribJncType for this flowper record (this global id)
            flowper_TribJncType = gdf.iloc[0]["TribJncTyp"]

            if (flowper_TribJncType == 'No Data') | (flowper_TribJncType == None):
                # When "TribJncType" == "NoData", no rows should be dropped.
                drop_conditions = None
                trib_jnc_processing_message = 'No tributary junction information'
            elif flowper_TribJncType == 'On tributary':
                # When "TribJncType" == "On tributary", all rows except where "mainstem_flag" == "tributary" should be dropped.
                drop_conditions = gdf['mainstem_flag'] != 'tributary'
                trib_jnc_processing_message = 'Matching with tributary'
            elif flowper_TribJncType == 'On mainstem downstream':
                # When "TribJncType" == "On mainstem downstream", all rows except where "mainstem_flag" == "mainstem" AND "trib_jnc" == "downstream of junction" should be dropped.
                drop_conditions = (gdf['mainstem_flag'] != 'mainstem') | (gdf['trib_jcn'] != 'downstream of junction')
                trib_jnc_processing_message = 'Matching with downstream mainstem'
            elif flowper_TribJncType == 'On mainstem upstream':
                # When "TribJncType" == "On mainstem upstream", all rows except where "mainstem_flag" == "mainstem" AND "trib_jnc" == "upstream of junction" should be dropped.
                drop_conditions = (gdf['mainstem_flag'] != 'mainstem') | (gdf['trib_jcn'] != 'upstream of junction')
                trib_jnc_processing_message = 'Matching with upstream mainstem'
            else:
                # If none of the conditions are met, no rows should be dropped, but we should note this
                drop_conditions = None
                print(f'Tributary junction match attempted for point {global_id}, but an issue was encountered. Field values may have been inconsistent: {flowper_TribJncType}.')
                trib_jnc_processing_message = f'Tributary junction match attempted for point {global_id}, but an issue was encountered. Field values may have been inconsistent. Flowline matching will be skipped for this point.'

            if type(drop_conditions) != type(None):
                # if we have rows to drop
                # find the indices meeting the drop_conditions above
                drop_indices = gdf.loc[drop_conditions].index
                num_drops = len(drop_indices)
                if num_flowlines - num_drops > 0:
                    # only if this won't drop all rows here, proceed with tributary junction matching
                    # and drop rows for this unique self.source_id
                    #print(f'dropping {num_drops} of {num_flowlines} rows')
                    self.hydrolinked_gdf.drop(index=drop_indices, inplace=True)
                if (num_flowlines - num_drops <= 0):
                    # otherwise, note this in the processing message field
                    # In other cases the trib jcn info might have us drop all possible lines
                    # Therefore, don't drop any flowlines and log a message
                    #print(f'cannot drop {num_drops} of {num_flowlines} rows')
                    trib_jnc_processing_message = 'Tributary junction matching did not find any matching flowlines. Flowline matching will be skipped for this point.'

            # apply the processing message
            self.hydrolinked_gdf.loc[self.hydrolinked_gdf[self.source_id] == global_id, 'trib_jnc_processing_message'] = trib_jnc_processing_message

                # In some cases there were no flowlines in the buffer region to begin with
                # Therefore, don't drop any flowlines and log a message
