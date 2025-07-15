"""Allows for hydrolinking of point observations to TerrainWorks hydrography points.

This requires a hydrography represented by point geometries, where each point corresponds to
a pixel in a stream grid raster.

"""

from hydrolink import utils
import geopandas as gpd
import pandas as pd
from typing import Union, Literal
from pathlib import Path
from shapely.geometry import Point, LineString

class TerrainWorksLink():
    """Class for hydrolinking observation point data to TerrainWorks point data"""

    # initialize parent class with additional kwargs
    def __init__(self,
                 points: Union[str, Path, gpd.GeoDataFrame],
                 flowlines: Union[str, Path, gpd.GeoDataFrame],
                 source_identifier: str,
                 flowlines_identifier: str,
                 water_name: str,
                 flowline_name: str,
                 input_crs: Union[int, str] = 4269,
                 buffer_m: Union[str, float] = 1000,
                 buffer_multiplier: float = 1,
                 default_buffer: Union[Literal["mean"], int, float] = 100,
                 no_stream_name_min_buffer: float = 10,
                 yes_stream_name_min_buffer: float = 15,
                 max_buffer_distance: float = 100,
                 flowline_grid_offsets: tuple = None):

        """Initiate attributes for hydrolinking point data to a custom hydrography dataset.

        During initiation of an object the buffer is verified to be less than 2000 meters, though this
        limitation may not be neccessary for custom, local hydrography datasets. Initiation also converts
        supplied point and flowline datasets to a common CRS (defaults to crs=3857) (area for future
        improvement: validate that points are within the bounds of the United States bounding box, or
        that points are within bounding box of the custom hydrography). Initiation should always occur if
        the required parameters are passed. If a process fails during initiation or while updating an
        object with class methods, the object remains available and the message and status of the object
        will be updated accordingly to capture and communicate errors.

        Parameters
        ----------
        points: str, Path, gpd.GeoDataFrame
            Path to a file, or a geopandas geodataframe containing points to hydrolink
        flowlines: str, Path, gpd.GeoDataFrame
            Path to a file, or a geopandas geodataframe containing the TerrainWorks flowline points
        source_identifier: str
            Field or column name within points dataset containing a unique identifier for each point of
            interest
        flowlines_identifier: str
            Field or column name within flowlines dataset containing a unique identifier for each flowline
        water_name: str
            Field or column name within points dataset containing a waterbody name for each point of
            interest, required for name_match method. Optional
        flowline_name: str
            Field or column name within hydrography dataset containing a waterbody name for each flowline,
            required for name_match method. Optional
        input_crs: str, int
            Currently not used (see Notes). Specify coordinate reference system, default is 4269
        buffer_m: str, int, float
            If a number is provided, distance in meters used as buffer to search for candidate flowline
            features to use for hydrolinking, default is 1000, max is 2000. If a string is provided, field
            or column name within points dataset containing a buffer distance unique to each point.
        buffer_multiplier: int, float
            Constant value to multiply the value in buffer_m by, default is 1.
        flowline_grid_offsets: tuple
            Offsets (in meters) to apply to the flowline grid points before snapping. This allows observations
            to be snapped to pixel centers if the flowline points are on pixel corners. Example, if flowline
            points are on the lower left corner of the corresponding stream grid pixel, and the pixels are
            2x2 meters, we can add a 1/2 pixel offset to snap to the grid center with
            flowline_grid_offsets = (x_offset, y_offset) = (1, 1) to shift one meter north and one meter east
            defaults to None and no offets are applied

        Notes
        ----------
        Coordinate reference system recommendations
            Ideally we would use NAD83 (CRS 4269), WGS84 (CRS 4326), Albers (CRS 5070) like the rest of
            hydrolink. Currently using EPSG 4269. Though there is an input_crs parameter, this parameter
            is not being used. Internally, whatever CRS the input points and input flowlines are in,
            they are both reprojected to a common UTM CRS based on the location of the points dataset.
            The hydrolinked results are then reprojected into EPSG 4269 before saving the output file.
        Source name recommendations
            To be most effective the names under water_name should contain no abbreviations and only
            contain official names from USGS Geospatial Names Information System (GNIS).
        Buffer recommendations
            Larger buffers may result in longer run times.
        Expected speed note
            TBD after further testing.

        """

        # column name in the points dataset with a unique ID for each point
        self.source_id = str(source_identifier)
        # column name in the points dataset with a unique ID for each flowline
        if flowlines_identifier is not None:
            self.flowlines_identifier = str(flowlines_identifier)

        # read in points and flowlines to geodataframes
        print('Reading input points...')
        self.input_crs = input_crs
        self.points = points # used later to write a new points file
        self.points_gdf = utils.input_to_gdf(self.points).drop_duplicates()
        print('Reading input flowlines...')
        self.flowlines = flowlines # used later to write a new flowlines file
        self.flowlines_gdf = utils.input_to_gdf(self.flowlines).drop_duplicates()
        if flowlines_identifier is None:
            self.flowlines_gdf['NODE_ID'] = self.flowlines_gdf.index
            self.flowlines_identifier = 'NODE_ID'

        # water name should be the name of a column name in points that contains the stream name
        self.water_name = str(water_name)
        self.flowline_name = str(flowline_name)

        # get options for setting up buffer (search radius) around each point
        self.buffer_m = buffer_m
        self.buffer_multiplier = buffer_multiplier
        self.default_buffer = default_buffer
        self.no_stream_name_min_buffer = no_stream_name_min_buffer
        self.yes_stream_name_min_buffer = yes_stream_name_min_buffer
        self.max_buffer_distance = max_buffer_distance

        # get option for flowline grid offsets (x, y) in meters
        self.flowlines_grid_offsets = flowline_grid_offsets

        self.status = 1  # where 0 is failed, 1 is worked properly
        self.message = ''

        # If buffer is greater than 2000 do not run and set error message
        # we don't really have a limit, but may want to keep it here for keeping computation time down
        if not isinstance(self.buffer_m, str) and self.buffer_m > 2000:
            self.message = 'Maximum buffer is 2000 meters, reduce buffer.'
            self.error_handling()
        # if buffer is a string, this is the column that contains a dynamic buffer value (to also be multiplied by buffer_multiplier)


        # Try converting to NAD83 (crs==4269) ? 3857 coordinate system if different coordinate system provided
        # If fails do not run and set error message
        self.utm_crs = self.points_gdf.estimate_utm_crs()
        print(f'Reproject points and flowlines to common UTM CRS: {self.utm_crs} ...')
        try:
            self.points_gdf.to_crs(self.utm_crs, inplace=True)
            self.flowlines_gdf.to_crs(self.utm_crs, inplace=True)
            # if self.points_gdf.crs != self.input_crs:
            #     self.points_gdf.to_crs(self.input_crs, inplace=True)
            # if self.flowlines_gdf.crs != self.input_crs:
            #     self.flowlines_gdf.to_crs(self.input_crs, inplace=True)


            # Test to make sure coordinates are within U.S. including Puerto Rico and Virgin Islands.
            # This is based on a general bounding box and intended to pick up common issues like missing values, 0 values and positive lon values
            # if (float(self.init_lat) > 17.5 and float(self.init_lat) < 71.5) and (float(self.init_lon) < -64.0 and float(self.init_lon) > -178.5):
            #     pass
            # else:
            #     self.message = f'Coordinates for id: {self.source_id} are outside of the bounding box of the United States.'
            #     self.error_handling()

            # check if any points lie outside the bounds of the input hydrography and alert the user?

        except:
            self.message = f'Issues handling provided coordinate system or coordinates for {self.source_id}. Consider using a common crs like 4269 (NAD83) or 4326 (WGS84).'
            self.error_handling()
        print('Ready to hydrolink.')

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
                print('Find flowlines within search regions...')
                self.intersect_points_flowlines()
                print('Hydrolink points to matching flowlines...')
                if trib_jcn is not None:
                    print('...use tributary junction information...')
                    #self.trib_jcn_match()
                    raise NotImplementedError
                if method == 'name_match':
                    print('...use name match method...')
                    #self.get_name_similarity()
                    #self.select_closest_flowline_w_name_match(similarity_cutoff=similarity_cutoff)
                    raise NotImplementedError
                elif method == 'closest':
                    print('...use closest method...')
                    self.select_closest_flowline()
                #print('Handling points that did not link to hydrography...')
                self.unlinked_points_handling()
                print('Writing output file...')
                self.write_hydrolink(outfile_name=outfile_name)
        else:
            print('Status = 0, writing output...')
            self.write_hydrolink(outfile_name=outfile_name)
        print('Done.')

    def buffer_points(self):

        """Create a buffer around the input points to define search regions for each point.

        Description
        ----------
        Creates buffers around the input points to define search regions. If user specified a fixed
        buffer value, all points will have the same search distance for finding nearby flowlines.
        If user provided a column name in the points dataset, the values in that column will be
        used for each point, such that each point can have a different search distance. The
        buffer values are multiplied by buffer_multiplier (defaults to 1) to allow easy scaling
        of the search distance. Additional checks are performed specific to FLOwPER observations.
        If a point includes stream name information, the buffer is not allowed to be less than
        yes_stream_name_min_buffer (defaults to 15 m) m. If a point does not contain stream name
        information, the buffer is not allowed to be less than no_stream_name_min_buffer (defaults
        to 10 m). A user can also specify how nodata values should be filled in, whether to use a
        mean value or a fixed default buffer value.

        Parameters
        ----------
        buffer_m: str, int, float
            If a number is provided, distance in meters used as buffer to search for candidate flowline
            features to use for hydrolinking. If a string is provided, field or column name within points
            dataset containing a buffer distance unique to each point.
        buffer_multiplier: float
            Constant value to multiply the value in buffer_m by, default is 1.
        default_buffer: "mean", float
            Choose what to do with points with nodata in column buffer_m. Either set to "mean" to use the
            mean value in column buffer_m * buffer_multiplier, or set to a fixed value (float),
            defaults to 1000 m
        water_name: str
            Field or column name within points dataset containing a waterbody name for each point of
            interest, required for name_match method. Optional
        no_stream_name_min_buffer: float
            Minimum allowable buffer size for points without stream name information, defaults to 10 m
        yes_stream_name_min_buffer: float
            Minimum allowable buffer size for points with stream name information, defaults to 15 m
        max_buffer_distance: float
            Maximum allowable buffer size, defaults to 1000 m

        Returns
        -------
        buffered_points_gdf: gpd.GeoDataFrame
            A geopandas GeoDataFrame containing all the same fields as points_gdf, except the point
            geometries have been replaced with the polygon geometries of the buffered search regions
            around each point.

        """

        # if buffer_m is not a string, buffer the points we want to search by some distance within which to search for flowlines
        if not isinstance(self.buffer_m, str):
            self.buffered_points_gdf = self.points_gdf.copy()
            # save the buffer distance used in its own column
            self.buffered_points_gdf['buffer_distance'] = self.buffer_m * self.buffer_multiplier
        # otherwise, buffer by the values in column 'buffer_m'
        else:
            self.buffered_points_gdf = self.points_gdf.copy()
            # make sure this column is interpreted as numeric
            self.buffered_points_gdf[self.buffer_m] = pd.to_numeric(self.buffered_points_gdf[self.buffer_m], errors='coerce')
            # start off with a copy of the values in buffer_m for buffer distances
            self.buffered_points_gdf['buffer_distance'] = self.buffered_points_gdf[self.buffer_m].copy() * self.buffer_multiplier
            # replace nans in column buffer_distance
            if self.default_buffer == "mean":
                # replace any nans with the mean value in column buffer_distance, which is = buffer_m * buffer_multiplier
                self.buffered_points_gdf.fillna({'buffer_distance': self.buffered_points_gdf['buffer_distance'].mean()}, inplace=True)
            else:
                # replace any nans with a default value
                self.buffered_points_gdf.fillna({'buffer_distance': self.default_buffer}, inplace=True)

            if self.water_name == None:
                # if no stream name information is to be used, use the no_stream_name_min_buffer for all points
                self.buffered_points_gdf['buffer_distance'] = self.buffered_points_gdf.apply(
                    lambda row: (
                        max(row['buffer_distance'], self.no_stream_name_min_buffer)
                    ),
                    axis=1
                )
            else:
                # calculate the buffer size for each point based on whether or not there is stream name information for the point
                self.buffered_points_gdf['buffer_distance'] = self.buffered_points_gdf.apply(
                    lambda row: (
                        max(row['buffer_distance'], self.no_stream_name_min_buffer)
                    ) if row[self.water_name] == 'na' else (
                        max(row['buffer_distance'], self.yes_stream_name_min_buffer)
                    ),
                    axis=1
                )
        # enforce a max buffer distance
        self.buffered_points_gdf.loc[self.buffered_points_gdf['buffer_distance'] > self.max_buffer_distance, 'buffer_distance'] = self.max_buffer_distance
        # finally, apply the buffer to each point geometry
        self.buffered_points_gdf['geometry'] = self.buffered_points_gdf['geometry'].buffer(self.buffered_points_gdf['buffer_distance'])



    def intersect_points_flowlines(self): #buffered_points_gdf, flowlines_gdf, flowlines_identifier, points_gdf, source_id):
        """Finds all flowlines that intersect with the search (buffered) region around each point.

        Description
        ----------
        Uses the geopandas overlay method to find all flowlines that are within the search (buffered)
        region around each input point. Merges information about these flowlines with the points
        GeoDataFrame to create a new GeoDataFrame. Later operations on this new GeoDataFrame will
        eliminate flowlines until we find the spot on a flowline we want to link the input point to.

        Parameters
        ----------
        self.buffered_points_gdf: gpd.GeoDataFrame
            The geopandas geodataframe containing search region geometries around each point of interest.
        self.flowlines_gdf: gpd.GeoDataFrame
            The geopandas geodataframe containing the original flowline geometries
        self.points_gdf: gpd.GeoDataFrame
            The geopandas geodataframe containing the original point geometries


        Returns
        -------
        self.hydrolinked_gdf: gpd.GeoDataFrame
            A geopandas GeoDataFrame containing the flowlines within the search region of each point of interest.

        """
        # apply an offset to the flowline points if self.flowline_grid_offsets is set
        if self.flowlines_grid_offsets != None:
            x_offset = self.flowlines_grid_offsets[0]
            y_offset = self.flowlines_grid_offsets[0]
            print(f'Applying an offset of ({x_offset}, {y_offset}) for each flowline (x, y) point geometry')
            self.flowlines_gdf['geometry'] = self.flowlines_gdf['geometry'].apply(lambda geom: Point(geom.x + x_offset, geom.y + y_offset))


        # Find the buffered flowlines that intersect with the buffered points
        closest_flowlines = gpd.sjoin(self.flowlines_gdf,
                                    self.buffered_points_gdf,
                                    predicate='intersects',
                                    how='inner',
                                    lsuffix='line',
                                    rsuffix='point')


        # now replace the buffered flowline geometries with the original line geometries
        self.hydrolinked_gdf = closest_flowlines.merge(self.flowlines_gdf.copy()[[self.flowlines_identifier, 'geometry']],
                                                on=self.flowlines_identifier,
                                                suffixes=('_buffered', ''),
                                                how='left')
        self.hydrolinked_gdf.drop(columns={'geometry_buffered'}, inplace=True)

        # add the point geometry for reference
        geometry_point_gdf = self.points_gdf.rename(columns={'geometry': 'geometry_point'})
        self.hydrolinked_gdf = self.hydrolinked_gdf.merge(geometry_point_gdf[[self.source_id, 'geometry_point']],
                                                            on=self.source_id,
                                                            suffixes=('', '_point'),
                                                            how='left')

        # calculate the distance to this nearest point
        self.hydrolinked_gdf['distance_to_point_on_line'] = self.hydrolinked_gdf['geometry_point'].distance(self.hydrolinked_gdf['geometry'])

    def select_closest_flowline(self):
        """Select closest flowline.

        Description
        -----------
        Selects closest flowline to each point of interest. Requires an initial hydrolinked_gdf
        to operate on, created from running intersect_points_flowlines, and find_nearest_point.

        Parameters
        ----------
        self.hydrolinked_gdf: gpd.GeoDataFrame
            A geopandas GeoDataFrame containing the flowlines within the search region of each point of interest.
            Now also containing the nearest points, distances, and name similarity scores.
        self.source_id: str
            Field or column name within points dataset containing a unique identifier for each point of interest.

        Returns
        -------
        self.hydrolinked_gdf: gpd.GeoDataFrame
            A geopandas GeoDataFrame containing the flowlines and points along the flowlines to snap
            each point of interest to.

        """
        self.hydrolinked_gdf.drop_duplicates(inplace=True)
        # count number of flowlines in search buffer region (this is just the number of rows for each global_id)
        self.hydrolinked_gdf['num_flowlines'] = self.hydrolinked_gdf.groupby(self.source_id)[self.source_id].transform('count')
        # create a field to store processing messages per point
        self.hydrolinked_gdf['processing_message'] = ''

        for global_id in self.hydrolinked_gdf[self.source_id].unique():
            # gdf points directly to this filtered self.hydrolinked_gdf, it is not a copy, this is just for ease of typing
            gdf = self.hydrolinked_gdf[self.hydrolinked_gdf[self.source_id] == global_id]
            # count number of flowlines in search buffer region (this is just the number of rows for each global_id)
            num_flowlines = len(gdf)
            # drop rows that are not the minimum distance
            #print("drow rows that are not the minimum distance")
            drop_conditions = (gdf['distance_to_point_on_line'] != gdf['distance_to_point_on_line'].min())

            # find the indices meeting the drop_conditions above
            drop_indices = gdf.loc[drop_conditions].index
            num_drops = len(drop_indices)

            # check that we are dropping all but 1 row for this unique self.source_id
            if ( num_flowlines - num_drops ) > 1 :
                # proceed for now, but provide a warning
                print(f"WARNING: Multiple flowlines met hydrolinking conditions for id: {global_id}")
                self.hydrolinked_gdf.loc[self.hydrolinked_gdf[self.source_id] == global_id, 'processing_message'] = [f"WARNING: Multiple flowlines met hydrolinking conditions for id: {global_id}"] * num_flowlines

            # drop rows for this unique self.source_id
            #print(f'dropping {num_drops} of {num_flowlines} rows')
            self.hydrolinked_gdf.drop(index=drop_indices, inplace=True)



    def unlinked_points_handling(self):
        """Handles what to do with points that did not link to any hydrography."""
        # Join the original point data back with the linked data to show which points didn't link anywhere
        # Find points (using source_id) that are not in the hydrolinked_gdf (this means no link was found for those points)
        unlinked_points_gdf = self.points_gdf[~self.points_gdf[self.source_id].isin(self.hydrolinked_gdf[self.source_id])].copy()
        # Rename geometry field to match hydrolinked_gdf
        unlinked_points_gdf.rename(columns={'geometry': 'geometry_point'}, inplace=True)
        # Concatenate the un-linked points back into the hydrolinked_gdf
        self.hydrolinked_gdf = pd.concat([self.hydrolinked_gdf, unlinked_points_gdf], ignore_index=True)
        # For the un-linked points, add a warning message and set snap point to source point by default
        self.hydrolinked_gdf.loc[self.hydrolinked_gdf['geometry'].isnull(), 'processing_message'] = f"WARNING: No flowlines found within buffer distance from point. Try increasing buffer. Snap point will be set to source point by default."
        self.hydrolinked_gdf.loc[self.hydrolinked_gdf['geometry'].isnull(), 'geometry'] = self.hydrolinked_gdf.loc[self.hydrolinked_gdf.geometry.isnull(), 'geometry_point']
        # Finally, remove observation points that fall outside a convex hull around the flowline points
        convex_hull = self.flowlines_gdf.unary_union.convex_hull
        convex_hull_buffer = self.max_buffer_distance # add a buffer to the convex hull equal to the maximum allowed point buffer distance (search radius)
        self.hydrolinked_gdf = self.hydrolinked_gdf[self.hydrolinked_gdf['geometry'].within(convex_hull.buffer(convex_hull_buffer))] # this check could be performed at the start of the processing, to avoid finding out only at the end that there are no overlaping points to output



    def write_hydrolink(self, outfile_name='custom_hydrolink_output.gpkg'):
        """Write hydrolink data output to file.

        Parameters
        ----------
        self.hydrolinked_gdf: gpd.GeoDataFrame
            A geopandas GeoDataFrame containing the flowlines and points along the flowlines to snap
            each point of interest to.
        outfile_name: str, Path
            File name and path to write the hydrolinked data to.
        """


        # print(self.hydrolinked_gdf)
        # print(len(self.hydrolinked_gdf))

        # if we have no points snapped here, skip writing a file
        if len(self.hydrolinked_gdf) == 0:
            print("No hydrolinked points to output. Check that there are points actually near this hydrography data.")
        else:
            # make sure that hydrolinked_gdf has a crs set, if not, set it here:
            if self.hydrolinked_gdf.crs == None:
                self.hydrolinked_gdf.set_crs(self.utm_crs, inplace=True)
            # format gdf into standard hydrolink output csv file
            # convert source and snap points to EPSG:4269, extract lats and lons to columns
            source_points_nad83 = self.hydrolinked_gdf.geometry_point.to_crs(4269)
            self.hydrolinked_gdf['source lat nad83'] = source_points_nad83.y
            self.hydrolinked_gdf['source lon nad83'] = source_points_nad83.x
            snap_points_nad83 = self.hydrolinked_gdf.geometry.to_crs(4269)
            self.hydrolinked_gdf['snap lat nad83'] = snap_points_nad83.y
            self.hydrolinked_gdf['snap lon nad83'] = snap_points_nad83.x

            # before saving, drop the geometries of the original input point, and the matched flowline
            self.hydrolinked_gdf.drop(columns=['geometry_point']).set_geometry('geometry').to_crs(4269).to_file(outfile_name)

    def write_connecting_lines(self, outfile_name='custom_hydrolink_connectors.gpkg'):
        """Create connecting lines and write to file.

        Description
        -----------
        Creates line geometries connecting each original input point of interest to the final
        hydrolinked point along the flowline geometries.

        Parameters
        ----------
        self.hydrolinked_gdf: gpd.GeoDataFrame
            A geopandas GeoDataFrame containing the flowlines and points along the flowlines to snap
            each point of interest to.
        outfile_name: str, Path
            File name and path to write the connecting lines data to.
        """
        # if we have no points snapped here, skip writing a file
        if len(self.hydrolinked_gdf) == 0:
            print("No connecting lines to output. Check that there are points actually near this hydrography data.")
        else:
            line_gdf = self.hydrolinked_gdf.apply(lambda row: LineString([Point(row['source lon nad83'], row['source lat nad83']), Point(row['snap lon nad83'], row['snap lat nad83'])]), axis=1)
            line_gdf.set_crs(4269).to_file(outfile_name)


    def error_handling(self):
        """Handle errors throughout hydrolink."""
        self.status = 0
        print(self.message)
