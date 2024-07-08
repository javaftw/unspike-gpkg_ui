import os
import fiona
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as MplPolygon
from matplotlib.collections import PatchCollection
from shapely.geometry import shape

def read_gpkg_polygons(file_path):
    polygons = []
    try:
        with fiona.open(file_path, 'r') as src:
            for feature in src:
                geom = shape(feature['geometry'])
                if geom.geom_type == 'Polygon':
                    polygons.append([list(geom.exterior.coords)])
                elif geom.geom_type == 'MultiPolygon':
                    for poly in geom.geoms:
                        polygons.append([list(poly.exterior.coords)])
    except Exception as e:
        print(f"An error occurred while reading the GeoPackage: {e}")
    return polygons

def plot_polygons(polygons):
    fig, ax = plt.subplots()
    patches = []
    
    for poly in polygons:
        for ring in poly:
            patches.append(MplPolygon(ring, closed=True))
    
    collection = PatchCollection(patches, facecolors='none', edgecolors='b')
    ax.add_collection(collection)
    
    ax.autoscale()
    ax.set_aspect('equal')
    plt.show()

def main(gpkg_path):
    try:
        if not os.path.exists(gpkg_path):
            raise FileNotFoundError(f"The file '{gpkg_path}' does not exist.")
        if not os.path.isfile(gpkg_path):
            raise ValueError(f"The path '{gpkg_path}' is not a file.")
        
        print(f"Reading GeoPackage file: {gpkg_path}")
        polygons = read_gpkg_polygons(gpkg_path)
        if not polygons:
            print("No geometries found.")
            return
        plot_polygons(polygons)
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    gpkg_path = "spiky-polygons.gpkg"  # Replace with your actual file path
    main(gpkg_path)
