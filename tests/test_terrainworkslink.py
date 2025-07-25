import pytest
import geopandas as gpd
from shapely.geometry import Point
from flowperlink.terrainworks import TerrainWorksLink


@pytest.fixture
def sample_data():
    # Create sample GeoDataFrames for points and flowlines
    points_data = gpd.GeoDataFrame({
        'id': [1, 2],
        'geometry': [Point(1, 1), Point(2, 2)],
        'dynamic_buffer': [500, 1000],
        'stream_name': ['anduin', 'brandywine']
    }, crs="EPSG:4326")

    flowlines_data = gpd.GeoDataFrame({
        'flow_id': [1, 2],
        'geometry': [Point(0, 0), Point(3, 3)],
        'flowline_name': ['Anduin', 'Baranduin']
    }, crs="EPSG:4326")

    return points_data, flowlines_data



def test_initialization(sample_data):
    points_gdf, flowlines_gdf = sample_data
    test_tw = TerrainWorksLink(
        points=points_gdf,
        flowlines=flowlines_gdf,
        points_identifier='id',
        flowlines_identifier='flow_id',
        water_name='stream_name',
        flowline_name='flowline_name'
    )
    assert test_tw.status == 1

def test_buffer_points_fixed(sample_data):
    points_gdf, flowlines_gdf = sample_data
    test_tw = TerrainWorksLink(
        points=points_gdf,
        flowlines=flowlines_gdf,
        points_identifier='id',
        flowlines_identifier='flow_id',
        water_name='stream_name',
        flowline_name='flowline_name',
        buffer_m=0.5
    )

    test_tw.buffer_points()
    assert 'buffer_distance' in test_tw.buffered_points_gdf.columns
    assert all(test_tw.buffered_points_gdf['buffer_distance'] == 0.5)

def test_buffer_points_dynamic(sample_data):
    points_gdf, flowlines_gdf = sample_data
    test_tw = TerrainWorksLink(
        points=points_gdf,
        flowlines=flowlines_gdf,
        points_identifier='id',
        flowlines_identifier='flow_id',
        water_name='stream_name',
        flowline_name='flowline_name',
        buffer_m='dynamic_buffer'
    )

    test_tw.buffer_points()
    assert 'buffer_distance' in test_tw.buffered_points_gdf.columns

def test_intersect_points_flowlines(sample_data):
    points_gdf, flowlines_gdf = sample_data
    test_tw = TerrainWorksLink(
        points=points_gdf,
        flowlines=flowlines_gdf,
        points_identifier='id',
        flowlines_identifier='flow_id',
        water_name='stream_name',
        flowline_name='flowline_name',
        buffer_m=0.5
    )

    test_tw.buffer_points()
    test_tw.intersect_points_flowlines()

    assert 'distance_to_point_on_line' in test_tw.hydrolinked_gdf.columns
    assert test_tw.buffered_points_gdf.geometry.geom_type[0] == 'Polygon'


def test_select_closest_flowline(sample_data):
    points_gdf, flowlines_gdf = sample_data
    test_tw = TerrainWorksLink(
        points=points_gdf,
        flowlines=flowlines_gdf,
        points_identifier='id',
        flowlines_identifier='flow_id',
        water_name='stream_name',
        flowline_name='flowline_name',
        buffer_m=0.5
    )

    test_tw.buffer_points()
    test_tw.intersect_points_flowlines()
    test_tw.select_closest_flowline()

    assert 'num_flowlines' in test_tw.hydrolinked_gdf.columns



def test_error_handling(sample_data):
    points_gdf, flowlines_gdf = sample_data
    test_tw = TerrainWorksLink(
        points=points_gdf,
        flowlines=flowlines_gdf,
        points_identifier='id',
        flowlines_identifier='flow_id',
        water_name='stream_name',
        flowline_name='flowline_name'
    )

    test_tw.message = "Test error"
    test_tw.error_handling()
    assert (test_tw.status == 0) and (test_tw.message == "Test error")
