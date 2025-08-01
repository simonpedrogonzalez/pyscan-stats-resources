# Part 1
import geopandas as gpd
import pyscan
import matplotlib.pyplot as plt

def plot_points(ax, pts, c):
    xs = [p[0] for p in pts]
    ys = [p[1] for p in pts]
    ax.scatter(xs, ys, color=c)

def plot_points_traj(ax, pts, c):
    xs = [p[0] for p in pts]
    ys = [p[1] for p in pts]
    ax.plot(xs, ys, color=c)

def plot_approx(ax, traj_pts, core_set_pts):
    plot_points_traj(ax, traj_pts, "g")
    plot_points(ax, core_set_pts, "b")
    ax.set_axis_off()

# Load the shapefile using geopandas
gdf = gpd.read_file("data/cb_2017_us_county_500k.shp")

# Use the first shape only
geom = gdf.geometry.iloc[0]

# Extract and flatten coordinates (handles Polygon and MultiPolygon)
if geom.geom_type == "Polygon":
    parts = [geom.exterior.coords]
elif geom.geom_type == "MultiPolygon":
    parts = [poly.exterior.coords for poly in geom.geoms]
else:
    raise ValueError(f"Unexpected geometry type: {geom.geom_type}")

pts_raw = [(x, y) for part in parts for (x, y) in part]
pts = [pyscan.Point(x, y, 1.0) for x, y in pts_raw]

alpha = 0.02
r_min = 0.05

# Part 2 - Convex Hull
f, ax = plt.subplots(figsize=(16, 12))
core_set_pts = pyscan.hull(pts)
plot_approx(ax, pts, core_set_pts)
plt.show()

# Part 3 - Halfplane Kernel
f, ax = plt.subplots(figsize=(16, 12))
core_set_pts = pyscan.halfplane_kernel(pts, alpha)
plot_approx(ax, pts, core_set_pts)
plt.show()

# Part 4 - DP Compress
f, ax = plt.subplots(figsize=(16, 12))
core_set_pts = pyscan.dp_compress(pts, alpha)
plot_approx(ax, pts, core_set_pts)
plt.show()

# Part 5 - Polygon Grid
core_set_pts = pyscan.polygon_grid(pts, r_min * 2)
f, ax = plt.subplots(figsize=(16, 12))
plot_approx(ax, pts, core_set_pts)
plt.show()

# Part 6 - Polygon Grid Even
core_set_pts = pyscan.polygon_grid_even(pts, r_min * 2, alpha)
f, ax = plt.subplots(figsize=(16, 12))
plot_approx(ax, pts, core_set_pts)
plt.show()

# Part 7 - Polygon Grid Hull
core_set_pts = pyscan.polygon_grid_hull(pts, r_min, alpha)
f, ax = plt.subplots(figsize=(16, 12))
plot_approx(ax, pts, core_set_pts)
plt.show()

