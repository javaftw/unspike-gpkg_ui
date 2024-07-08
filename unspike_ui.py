import os
import warnings
from typing import Tuple, Union

import fiona
import numpy as np
from fiona.errors import FionaDeprecationWarning
from pyproj import CRS
from shapely.geometry import shape, mapping, Polygon, MultiPolygon
from shapely.validation import make_valid
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as MplPolygon
from matplotlib.collections import PatchCollection

# Suppress the Fiona class deprecation warning
warnings.filterwarnings('ignore', category=FionaDeprecationWarning)
# Suppress the PyProj warning about lossy conversion
warnings.filterwarnings('ignore', category=UserWarning)

def calculate_angle(p0: Tuple[float, float], p1: Tuple[float, float], p2: Tuple[float, float]) -> float:
    """
    Calculate the angle between three points in 2D space.

    Args:
        p0 (Tuple[float, float]): Coordinates of the first point.
        p1 (Tuple[float, float]): Coordinates of the second point (vertex).
        p2 (Tuple[float, float]): Coordinates of the third point.

    Returns:
        float: Angle in degrees.
    """
    # Calculate vectors from p1 to p0 and p1 to p2
    v1 = np.array(p0) - np.array(p1)
    v2 = np.array(p2) - np.array(p1)

    # Calculate dot product and norms
    dot_product = np.dot(v1, v2)
    norm_v1 = np.linalg.norm(v1)
    norm_v2 = np.linalg.norm(v2)

    # Handle cases where vectors are extremely small to avoid division by zero
    if norm_v1 < 1e-10 or norm_v2 < 1e-10:
        return 180.0

    # Calculate and return the angle using the dot product formula
    cosine_angle = np.clip(dot_product / (norm_v1 * norm_v2), -1.0, 1.0)
    return np.degrees(np.arccos(cosine_angle))

def filter_polygon(poly: Polygon, min_angle: float) -> Tuple[Polygon, int]:
    """
    Filter a polygon by removing vertices that form angles less than the specified minimum.

    Args:
        poly (Polygon): Input polygon to filter.
        min_angle (float): Minimum angle threshold in degrees.

    Returns:
        Tuple[Polygon, int]: Filtered polygon and number of spikes removed.
    """
    def check_angle(prev, curr, next):
        angle = calculate_angle(prev, curr, next)
        if angle >= min_angle:
            return curr, 0
        return None, 1

    coords = list(poly.exterior.coords)
    new_coords = []
    spikes_removed = 0

    # Iterate through all vertices except the last (which is the same as the first)
    for i in range(len(coords) - 1):
        prev, curr, next = coords[i - 1], coords[i], coords[(i + 1) % (len(coords) - 1)]
        point, spike_removed = check_angle(prev, curr, next)
        if point is not None:
            new_coords.append(point)
        spikes_removed += spike_removed

    # Handle the case where we need to adjust the start/end point
    if len(new_coords) >= 3:
        start_end_point, spike_removed = check_angle(new_coords[-1], new_coords[0], new_coords[1])
        if start_end_point is None:
            # If the start/end point forms a spike, replace it with the midpoint
            midpoint = tuple((np.array(new_coords[-1]) + np.array(new_coords[1])) / 2)
            new_coords[0] = midpoint
            spikes_removed += spike_removed
    else:
        # If we don't have enough points to form a valid polygon, return an empty one
        return Polygon(), spikes_removed

    # Close the polygon by adding the first point at the end
    new_coords.append(new_coords[0])
    filtered_poly = Polygon(new_coords)

    # Ensure the resulting polygon is valid
    if not filtered_poly.is_valid:
        filtered_poly = make_valid(filtered_poly)

    return filtered_poly, spikes_removed

def filter_vertices(geometry: Union[Polygon, MultiPolygon], min_angle: float) -> Tuple[Union[Polygon, MultiPolygon], int]:
    """
    Filter vertices of a geometry (Polygon or MultiPolygon) based on the minimum angle.

    Args:
        geometry (Union[Polygon, MultiPolygon]): Input geometry to filter.
        min_angle (float): Minimum angle threshold in degrees.

    Returns:
        Tuple[Union[Polygon, MultiPolygon], int]: Filtered geometry and total spikes removed.

    Raises:
        ValueError: If the input geometry is not a Polygon or MultiPolygon.
    """
    if isinstance(geometry, Polygon):
        return filter_polygon(geometry, min_angle)
    elif isinstance(geometry, MultiPolygon):
        filtered_polys = []
        total_spikes_removed = 0
        for poly in geometry.geoms:
            filtered_poly, spikes_removed = filter_polygon(poly, min_angle)
            if not filtered_poly.is_empty:
                filtered_polys.extend([filtered_poly] if isinstance(filtered_poly, Polygon) else filtered_poly.geoms)
            total_spikes_removed += spikes_removed

        # Determine the appropriate return type based on the number of filtered polygons
        if len(filtered_polys) > 1:
            return MultiPolygon(filtered_polys), total_spikes_removed
        elif len(filtered_polys) == 1:
            return filtered_polys[0], total_spikes_removed
        else:
            return Polygon(), total_spikes_removed
    else:
        raise ValueError(f"Unsupported geometry type: {type(geometry)}")

def read_gpkg_polygons(file_path, min_angle):
    polygons = []
    filtered_polygons = []
    total_spikes_removed = 0
    try:
        with fiona.open(file_path, 'r') as src:
            for feature in src:
                geom = shape(feature['geometry'])
                if geom.geom_type == 'Polygon' or geom.geom_type == 'MultiPolygon':
                    filtered_geom, spikes_removed = filter_vertices(geom, min_angle)
                    total_spikes_removed += spikes_removed
                    if not filtered_geom.is_empty:
                        polygons.append([list(filtered_geom.exterior.coords)])
                        filtered_polygons.append(filtered_geom)
    except Exception as e:
        print(f"An error occurred while reading the GeoPackage: {e}")
    return polygons, filtered_polygons, total_spikes_removed

def plot_polygons(polygons):
    fig, ax = plt.subplots()
    patches = []
    
    for poly in polygons:
        for ring in poly:
            patches.append(MplPolygon(ring, closed=True))
    
    collection = PatchCollection(patches, facecolors='cyan', edgecolors='blue', linewidths=1.5, alpha=0.5)
    ax.add_collection(collection)
    
    ax.autoscale()
    ax.set_aspect('equal')
    
    # Adding titles and labels
    ax.set_title("Polygon Plot from GeoPackage")
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    
    # Adding gridlines
    ax.grid(True, linestyle='--', linewidth=0.5)
    
    plt.show()

def write_gpkg_polygons(filtered_polygons, input_path, output_path):
    try:
        with fiona.open(input_path, 'r') as src:
            meta = src.meta
            input_crs = CRS(src.crs)

            new_schema = {'geometry': 'MultiPolygon', 'properties': meta['schema']['properties']}

            with fiona.open(output_path, 'w', crs=src.crs, driver='GPKG', schema=new_schema) as dst:
                for feature, geom in zip(src, filtered_polygons):
                    if isinstance(geom, Polygon):
                        geom = MultiPolygon([geom])
                    feature['geometry'] = mapping(geom)
                    dst.write({'geometry': feature['geometry'], 'properties': feature['properties']})
    except Exception as e:
        print(f"An error occurred while writing the GeoPackage: {e}")

def main(gpkg_path, min_angle):
    try:
        if not os.path.exists(gpkg_path):
            raise FileNotFoundError(f"The file '{gpkg_path}' does not exist.")
        if not os.path.isfile(gpkg_path):
            raise ValueError(f"The path '{gpkg_path}' is not a file.")
        
        print(f"Reading GeoPackage file: {gpkg_path}")
        polygons, filtered_polygons, total_spikes_removed = read_gpkg_polygons(gpkg_path, min_angle)
        if not polygons:
            print("No geometries found.")
            return
        plot_polygons(polygons)
        
        output_path = f"{os.path.splitext(gpkg_path)[0]}_unspiked{os.path.splitext(gpkg_path)[1]}"
        write_gpkg_polygons(filtered_polygons, gpkg_path, output_path)
        print(f"Total spikes removed: {total_spikes_removed}")
        print(f"Unspiked GeoPackage written to: {output_path}")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    gpkg_path = "spiky-polygons.gpkg"  # Replace with your actual file path
    min_angle = 10.0  # Set your minimum angle threshold here
    main(gpkg_path, min_angle)
