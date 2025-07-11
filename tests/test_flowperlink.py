import pytest
import geopandas as gpd
from shapely.geometry import Point, LineString
from flowperlink.flowper import FlowperLink

@pytest.fixture
def sample_data():
    # Create sample GeoDataFrames for points and flowlines
    points_data = gpd.GeoDataFrame({
        'id': [1, 2, 3, 4, 5, 6, 7],
        'geometry': [Point(-122.4426, 47.60602),
                     Point(-122.4425, 47.60618),
                     Point(-122.4423, 47.60623),
                     Point(-122.4428, 47.60596),
                     Point(-122.4423, 47.60632),
                     Point(-122.4423, 47.60623),
                     Point(-122.4168, 47.49301)],
        'dynamic_buffer': [5, 10, 10, 2, 4, 10, 5],
        'stream_name': ['anduin', 'brandywine', None, None, None, 'Baranduin', 'Isen'],
        'TribJncTyp': ['On mainstem downstream', 'On tributary', 'On mainstem upstream', 'No Data', 'On tributary', 'On tributary', 'No Data'],
    }, crs="EPSG:4326")

    flowlines_data = gpd.GeoDataFrame({
        'reachcode': [1, 2, 3, 4, 5],
        'geometry': [LineString([(-122.4427, 47.60596), (-122.4425, 47.60613)]),
                     LineString([(-122.4425, 47.60613), (-122.4425, 47.60622)]),
                     LineString([(-122.4425, 47.60613), (-122.4423, 47.60622)]),
                     LineString([(-122.4430, 47.60592), (-122.4427, 47.60596)]),
                     LineString([(-122.4425, 47.60622), (-122.4423, 47.60631)])],
        'flowline_name': ['Anduin', 'Baranduin', 'Anduin', 'Anduin', 'Baranduin'],
        'mainstem_flag': ['mainstem', 'tributary', 'mainstem', 'mainstem', 'tributary'],
        'trib_jcn': ['downstream of junction', 'upstream of junction', 'upstream of junction', 'not near junction', 'not near junction']
    }, crs="EPSG:4326")

    return points_data, flowlines_data


def test_initialization_valid(sample_data):
    points_gdf, flowlines_gdf = sample_data
    test_flowperlink = FlowperLink(points = points_gdf,
                                   flowlines = flowlines_gdf,
                                   source_identifier='id',
                                   flowlines_identifier='reachcode',
                                   water_name='stream_name',
                                   flowline_name='flowline_name',
                                   buffer_m = "dynamic_buffer",
                                   buffer_multiplier=1,
                                   default_buffer = 2,
                                   no_stream_name_min_buffer=1,
                                   yes_stream_name_min_buffer=2,
                                   max_buffer_distance=10,
                                   keep_points_attributes=['TribJncTyp'],
                                   keep_flowlines_attributes=['mainstem_flag', 'trib_jcn'])
    assert test_flowperlink.status == 1

def test_buffer_points_fixed(sample_data):
    points_gdf, flowlines_gdf = sample_data
    test_flowperlink = FlowperLink(points = points_gdf,
                                flowlines = flowlines_gdf,
                                source_identifier='id',
                                flowlines_identifier='reachcode',
                                water_name='stream_name',
                                flowline_name='flowline_name',
                                buffer_m = 2,
                                buffer_multiplier=1,
                                default_buffer = 2,
                                no_stream_name_min_buffer=1,
                                yes_stream_name_min_buffer=2,
                                max_buffer_distance=10,
                                keep_points_attributes=['TribJncTyp'],
                                keep_flowlines_attributes=['mainstem_flag', 'trib_jcn'])
    test_flowperlink.buffer_points()
    assert test_flowperlink.buffered_points_gdf.geometry.geom_type[0] == 'Polygon'

def test_buffer_points_dynamic(sample_data):
    points_gdf, flowlines_gdf = sample_data
    test_flowperlink = FlowperLink(points = points_gdf,
                                flowlines = flowlines_gdf,
                                source_identifier='id',
                                flowlines_identifier='reachcode',
                                water_name='stream_name',
                                flowline_name='flowline_name',
                                buffer_m = "dynamic_buffer",
                                buffer_multiplier=1,
                                default_buffer = 2,
                                no_stream_name_min_buffer=1,
                                yes_stream_name_min_buffer=2,
                                max_buffer_distance=10,
                                keep_points_attributes=['TribJncTyp'],
                                keep_flowlines_attributes=['mainstem_flag', 'trib_jcn'])
    test_flowperlink.buffer_points()
    assert test_flowperlink.buffered_points_gdf.geometry.geom_type[0] == 'Polygon'

def test_buffer_flowlines(sample_data):
    points_gdf, flowlines_gdf = sample_data
    test_flowperlink = FlowperLink(points = points_gdf,
                                   flowlines = flowlines_gdf,
                                   source_identifier='id',
                                   flowlines_identifier='reachcode',
                                   water_name='stream_name',
                                   flowline_name='flowline_name',
                                   buffer_m = "dynamic_buffer",
                                   buffer_multiplier=1,
                                   default_buffer = 2,
                                   no_stream_name_min_buffer=1,
                                   yes_stream_name_min_buffer=2,
                                   max_buffer_distance=10,
                                   keep_points_attributes=['TribJncTyp'],
                                   keep_flowlines_attributes=['mainstem_flag', 'trib_jcn'])
    test_flowperlink.buffer_flowlines()
    assert test_flowperlink.buffered_flowlines_gdf.geometry.geom_type[0] == 'Polygon'


def test_intersect_points_flowlines(sample_data):
    """Test that we can intersect the buffered points and lines"""
    points_gdf, flowlines_gdf = sample_data
    test_flowperlink = FlowperLink(points = points_gdf,
                                   flowlines = flowlines_gdf,
                                   source_identifier='id',
                                   flowlines_identifier='reachcode',
                                   water_name='stream_name',
                                   flowline_name='flowline_name',
                                   buffer_m = "dynamic_buffer",
                                   buffer_multiplier=1,
                                   default_buffer = 2,
                                   no_stream_name_min_buffer=1,
                                   yes_stream_name_min_buffer=2,
                                   max_buffer_distance=10,
                                   keep_points_attributes=['TribJncTyp'],
                                   keep_flowlines_attributes=['mainstem_flag', 'trib_jcn'])
    test_flowperlink.buffer_points()
    test_flowperlink.buffer_flowlines()
    test_flowperlink.intersect_points_flowlines()
    columns_match = set(test_flowperlink.buffered_points_gdf.columns).issubset(test_flowperlink.hydrolinked_gdf.columns) and set(test_flowperlink.buffered_flowlines_gdf.columns).issubset(test_flowperlink.hydrolinked_gdf.columns)
    assert not test_flowperlink.hydrolinked_gdf.empty
    assert columns_match

def test_find_nearest_point(sample_data):
    """Test finding the nearest point along flowline to input point"""
    points_gdf, flowlines_gdf = sample_data
    test_flowperlink = FlowperLink(points = points_gdf,
                                   flowlines = flowlines_gdf,
                                   source_identifier='id',
                                   flowlines_identifier='reachcode',
                                   water_name='stream_name',
                                   flowline_name='flowline_name',
                                   buffer_m = "dynamic_buffer",
                                   buffer_multiplier=1,
                                   default_buffer = 2,
                                   no_stream_name_min_buffer=1,
                                   yes_stream_name_min_buffer=2,
                                   max_buffer_distance=10,
                                   keep_points_attributes=['TribJncTyp'],
                                   keep_flowlines_attributes=['mainstem_flag', 'trib_jcn'])
    test_flowperlink.buffer_points()
    test_flowperlink.buffer_flowlines()
    test_flowperlink.intersect_points_flowlines()
    test_flowperlink.find_nearest_point()
    assert 'geometry_point_on_line' in test_flowperlink.hydrolinked_gdf.columns
    assert all(test_flowperlink.hydrolinked_gdf['geometry_point_on_line'].geom_type == 'Point')

def test_get_name_similarity(sample_data):
    """Test computing name similarity between flowline and point"""
    points_gdf, flowlines_gdf = sample_data
    test_flowperlink = FlowperLink(points = points_gdf,
                                   flowlines = flowlines_gdf,
                                   source_identifier='id',
                                   flowlines_identifier='reachcode',
                                   water_name='stream_name',
                                   flowline_name='flowline_name',
                                   buffer_m = "dynamic_buffer",
                                   buffer_multiplier=1,
                                   default_buffer = 2,
                                   no_stream_name_min_buffer=1,
                                   yes_stream_name_min_buffer=2,
                                   max_buffer_distance=10,
                                   keep_points_attributes=['TribJncTyp'],
                                   keep_flowlines_attributes=['mainstem_flag', 'trib_jcn'])
    test_flowperlink.buffer_points()
    test_flowperlink.buffer_flowlines()
    test_flowperlink.intersect_points_flowlines()
    test_flowperlink.find_nearest_point()
    test_flowperlink.get_name_similarity()
    assert 'flowline name similarity' in test_flowperlink.hydrolinked_gdf.columns

def test_trib_jcn_match(sample_data):
    """Test selecting flowline with tributary junction info"""
    points_gdf, flowlines_gdf = sample_data
    test_flowperlink = FlowperLink(points = points_gdf,
                                   flowlines = flowlines_gdf,
                                   source_identifier='id',
                                   flowlines_identifier='reachcode',
                                   water_name='stream_name',
                                   flowline_name='flowline_name',
                                   buffer_m = "dynamic_buffer",
                                   buffer_multiplier=1,
                                   default_buffer = 2,
                                   no_stream_name_min_buffer=1,
                                   yes_stream_name_min_buffer=2,
                                   max_buffer_distance=10,
                                   keep_points_attributes=['TribJncTyp'],
                                   keep_flowlines_attributes=['mainstem_flag', 'trib_jcn'])
    test_flowperlink.buffer_points()
    test_flowperlink.buffer_flowlines()
    test_flowperlink.intersect_points_flowlines()
    test_flowperlink.find_nearest_point()
    test_flowperlink.get_name_similarity()
    test_flowperlink.trib_jcn_match()
    assert test_flowperlink.hydrolinked_gdf.loc[test_flowperlink.hydrolinked_gdf[test_flowperlink.source_id] == 1, 'trib_jnc_processing_message'].iloc[0] == 'Matching with downstream mainstem'
    assert test_flowperlink.hydrolinked_gdf.loc[test_flowperlink.hydrolinked_gdf[test_flowperlink.source_id] == 2, 'trib_jnc_processing_message'].iloc[0] == 'Matching with tributary'
    assert test_flowperlink.hydrolinked_gdf.loc[test_flowperlink.hydrolinked_gdf[test_flowperlink.source_id] == 3, 'trib_jnc_processing_message'].iloc[0] == 'Matching with upstream mainstem'
    assert test_flowperlink.hydrolinked_gdf.loc[test_flowperlink.hydrolinked_gdf[test_flowperlink.source_id] == 4, 'trib_jnc_processing_message'].iloc[0] == 'No tributary junction information'
    assert test_flowperlink.hydrolinked_gdf.loc[test_flowperlink.hydrolinked_gdf[test_flowperlink.source_id] == 5, 'trib_jnc_processing_message'].iloc[0] == 'Matching with tributary'
    assert test_flowperlink.hydrolinked_gdf.loc[test_flowperlink.hydrolinked_gdf[test_flowperlink.source_id] == 6, 'trib_jnc_processing_message'].iloc[0] == 'Matching with tributary'


def test_select_closest_flowline_w_name_match(sample_data):
    """Test selecting flowline with name match from point"""
    points_gdf, flowlines_gdf = sample_data
    test_flowperlink = FlowperLink(points = points_gdf,
                                   flowlines = flowlines_gdf,
                                   source_identifier='id',
                                   flowlines_identifier='reachcode',
                                   water_name='stream_name',
                                   flowline_name='flowline_name',
                                   buffer_m = "dynamic_buffer",
                                   buffer_multiplier=1,
                                   default_buffer = 2,
                                   no_stream_name_min_buffer=1,
                                   yes_stream_name_min_buffer=2,
                                   max_buffer_distance=10,
                                   keep_points_attributes=['TribJncTyp'],
                                   keep_flowlines_attributes=['mainstem_flag', 'trib_jcn'])
    test_flowperlink.buffer_points()
    test_flowperlink.buffer_flowlines()
    test_flowperlink.intersect_points_flowlines()
    test_flowperlink.find_nearest_point()
    test_flowperlink.get_name_similarity()
    test_flowperlink.select_closest_flowline_w_name_match()
    assert 'geometry_point_on_line' in test_flowperlink.hydrolinked_gdf.columns
    assert all(test_flowperlink.hydrolinked_gdf['geometry_point_on_line'].geom_type == 'Point')
    assert test_flowperlink.hydrolinked_gdf.loc[test_flowperlink.hydrolinked_gdf[test_flowperlink.source_id] == 1, 'reachcode'].iloc[0] == 1
    assert test_flowperlink.hydrolinked_gdf.loc[test_flowperlink.hydrolinked_gdf[test_flowperlink.source_id] == 2, 'reachcode'].iloc[0] == 2
    assert test_flowperlink.hydrolinked_gdf.loc[test_flowperlink.hydrolinked_gdf[test_flowperlink.source_id] == 3, 'reachcode'].iloc[0] == 3
    assert test_flowperlink.hydrolinked_gdf.loc[test_flowperlink.hydrolinked_gdf[test_flowperlink.source_id] == 4, 'reachcode'].iloc[0] == 4
    assert test_flowperlink.hydrolinked_gdf.loc[test_flowperlink.hydrolinked_gdf[test_flowperlink.source_id] == 5, 'reachcode'].iloc[0] == 5
    assert test_flowperlink.hydrolinked_gdf.loc[test_flowperlink.hydrolinked_gdf[test_flowperlink.source_id] == 6, 'reachcode'].iloc[0] == 5


def test_select_closest_flowline(sample_data):
    """Test selecting closest flowline to point"""
    points_gdf, flowlines_gdf = sample_data
    test_flowperlink = FlowperLink(points = points_gdf,
                                   flowlines = flowlines_gdf,
                                   source_identifier='id',
                                   flowlines_identifier='reachcode',
                                   water_name='stream_name',
                                   flowline_name='flowline_name',
                                   buffer_m = "dynamic_buffer",
                                   buffer_multiplier=1,
                                   default_buffer = 2,
                                   no_stream_name_min_buffer=1,
                                   yes_stream_name_min_buffer=2,
                                   max_buffer_distance=10,
                                   keep_points_attributes=['TribJncTyp'],
                                   keep_flowlines_attributes=['mainstem_flag', 'trib_jcn'])
    test_flowperlink.buffer_points()
    test_flowperlink.buffer_flowlines()
    test_flowperlink.intersect_points_flowlines()
    test_flowperlink.find_nearest_point()
    test_flowperlink.select_closest_flowline()
    assert 'geometry_point_on_line' in test_flowperlink.hydrolinked_gdf.columns
    assert all(test_flowperlink.hydrolinked_gdf['geometry_point_on_line'].geom_type == 'Point')
    assert test_flowperlink.hydrolinked_gdf.loc[test_flowperlink.hydrolinked_gdf[test_flowperlink.source_id] == 1, 'reachcode'].iloc[0] == 1
    assert test_flowperlink.hydrolinked_gdf.loc[test_flowperlink.hydrolinked_gdf[test_flowperlink.source_id] == 2, 'reachcode'].iloc[0] == 2
    assert test_flowperlink.hydrolinked_gdf.loc[test_flowperlink.hydrolinked_gdf[test_flowperlink.source_id] == 3, 'reachcode'].iloc[0] == 3
    assert test_flowperlink.hydrolinked_gdf.loc[test_flowperlink.hydrolinked_gdf[test_flowperlink.source_id] == 4, 'reachcode'].iloc[0] == 4
    assert test_flowperlink.hydrolinked_gdf.loc[test_flowperlink.hydrolinked_gdf[test_flowperlink.source_id] == 5, 'reachcode'].iloc[0] == 5
    assert test_flowperlink.hydrolinked_gdf.loc[test_flowperlink.hydrolinked_gdf[test_flowperlink.source_id] == 6, 'reachcode'].iloc[0] == 3


def test_write_hydrolink(sample_data, tmp_path):
    """Test writing hydrolink file out"""
    points_gdf, flowlines_gdf = sample_data
    test_flowperlink = FlowperLink(points = points_gdf,
                                   flowlines = flowlines_gdf,
                                   source_identifier='id',
                                   flowlines_identifier='reachcode',
                                   water_name='stream_name',
                                   flowline_name='flowline_name',
                                   buffer_m = "dynamic_buffer",
                                   buffer_multiplier=1,
                                   default_buffer = 2,
                                   no_stream_name_min_buffer=1,
                                   yes_stream_name_min_buffer=2,
                                   max_buffer_distance=10,
                                   keep_points_attributes=['TribJncTyp'],
                                   keep_flowlines_attributes=['mainstem_flag', 'trib_jcn'])
    test_flowperlink.buffer_points()
    test_flowperlink.buffer_flowlines()
    test_flowperlink.intersect_points_flowlines()
    test_flowperlink.find_nearest_point()
    test_flowperlink.select_closest_flowline()
    print(test_flowperlink.hydrolinked_gdf.columns)
    test_flowperlink.write_hydrolink(outfile_name=tmp_path / "output.gpkg")
    assert (tmp_path / "output.gpkg").exists()

def test_unlinked_points_handling(sample_data):
    """Test the handling of unlinked points, should end up with gdf of same length as input points_gdf"""
    points_gdf, flowlines_gdf = sample_data
    test_flowperlink = FlowperLink(points = points_gdf,
                                   flowlines = flowlines_gdf,
                                   source_identifier='id',
                                   flowlines_identifier='reachcode',
                                   water_name='stream_name',
                                   flowline_name='flowline_name',
                                   buffer_m = "dynamic_buffer",
                                   buffer_multiplier=1,
                                   default_buffer = 2,
                                   no_stream_name_min_buffer=1,
                                   yes_stream_name_min_buffer=2,
                                   max_buffer_distance=10,
                                   keep_points_attributes=['TribJncTyp'],
                                   keep_flowlines_attributes=['mainstem_flag', 'trib_jcn'])
    test_flowperlink.buffer_points()
    test_flowperlink.buffer_flowlines()
    test_flowperlink.intersect_points_flowlines()
    test_flowperlink.find_nearest_point()
    test_flowperlink.select_closest_flowline()
    test_flowperlink.unlinked_points_handling()
    assert 'processing_message' in test_flowperlink.hydrolinked_gdf.columns
    assert test_flowperlink.hydrolinked_gdf.loc[test_flowperlink.hydrolinked_gdf[test_flowperlink.source_id] == 7, 'processing_message'].iloc[0] == 'WARNING: No flowlines found within buffer distance from point. Try increasing buffer. Snap point will be set to source point by default.'

def test_error_handling(sample_data):
    """Test error_handling function"""
    points_gdf, flowlines_gdf = sample_data
    test_flowperlink = FlowperLink(points = points_gdf,
                                   flowlines = flowlines_gdf,
                                   source_identifier='id',
                                   flowlines_identifier='reachcode',
                                   water_name='stream_name',
                                   flowline_name='flowline_name',
                                   buffer_m = "dynamic_buffer",
                                   buffer_multiplier=1,
                                   default_buffer = 2,
                                   no_stream_name_min_buffer=1,
                                   yes_stream_name_min_buffer=2,
                                   max_buffer_distance=10,
                                   keep_points_attributes=['TribJncTyp'],
                                   keep_flowlines_attributes=['mainstem_flag', 'trib_jcn'])
    test_flowperlink.message = "Test error"
    test_flowperlink.error_handling()
    assert test_flowperlink.status == 0
    assert test_flowperlink.message == "Test error"
