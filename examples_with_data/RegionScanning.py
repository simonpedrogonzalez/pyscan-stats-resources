

# Unzip co-est2019-alldata.zip and cb_2017_us_county_500k.zip
# Data Sources https://www.census.gov/data/datasets/time-series/demo/popest/2010s-counties-total.html#par_textimage_70769902
# https://www.census.gov/geographies/mapping-files/2017/geo/carto-boundary-file.html
import geopandas as gpd
import pandas as pd
import pyscan
import matplotlib.pyplot as plt


def plot_points(ax, pts, c):
    xs = []
    ys = []
    for pt in pts:
        xs.append(pt[0] )
        ys.append(pt[1])
    ax.scatter(xs, ys, color=c, marker='.')

def plot_points_traj(ax, pts, c):
    xs = []
    ys = []
    for pt in pts:
        xs.append(pt[0])
        ys.append(pt[1])
    ax.plot(xs, ys, color=c)

def plot_approx(ax, regions, core_set_pts):
    for reg in regions:
        plot_points_traj(ax, reg, "g")
    plot_points(ax, core_set_pts, "b")
    ax.set_axis_off()
    
    
# Part 2

gdf = gpd.read_file("data/cb_2017_us_county_500k.shp")
# Filter to continental US
gdf = gdf.cx[-124.84:-66.9, 24.396:49.4].copy()

df = pd.read_csv("data/co-est2019-alldata.csv", encoding='latin-1', dtype={'STATE': str, 'COUNTY': str, 'POPESTIMATE2010': int, 'POPESTIMATE2017': int})
# Set index as STATE + COUNTY
df['FIPS'] = df['STATE'] + df['COUNTY']
population2010 = df.set_index('FIPS')['POPESTIMATE2010'].to_dict()
population2017 = df.set_index('FIPS')['POPESTIMATE2017'].to_dict()

regions = []
weights2017 = []
weights2010 = []

for _, row in gdf.iterrows():
    fips = row['GEOID']
    if fips not in population2017 or fips not in population2010:
        continue

    geometry = row.geometry

    # Get all exterior points from Polygon or MultiPolygon
    if geometry.geom_type == "Polygon":
        parts = [geometry.exterior.coords]
    elif geometry.geom_type == "MultiPolygon":
        parts = [part.exterior.coords for part in geometry.geoms]
    else:
        continue

    # Flatten into one list of points
    points = [(x, y) for part in parts for (x, y) in part]

    # Append all at once
    regions.append([pyscan.Point(x, y, 1.0) for x, y in points])
    weights2017.append(population2017[fips])
    weights2010.append(population2010[fips])


alpha = .02
r_min = .05

# Part 4

core_set_pts2010 = pyscan.polygon_sample(regions, weights2010, 10000)
f, ax = plt.subplots(figsize=(16, 12))
plot_approx(ax, regions, core_set_pts2010)
plt.show()

core_set_pts2017 = pyscan.polygon_sample(regions, weights2017, 10000)
_, ax = plt.subplots(figsize=(16, 12))
plot_approx(ax, regions, core_set_pts2017)
plt.show()

# Part 5


net = pyscan.my_sample(core_set_pts2017, 200) + pyscan.my_sample(core_set_pts2010, 200)
disc_f = pyscan.DISC
disk, d_val = pyscan.max_disk(net,
                              [pyscan.WPoint(1.0, p[0], p[1], 1.0) for p in core_set_pts2017],
                              [pyscan.WPoint(1.0, p[0], p[1], 1.0) for p in core_set_pts2010],
                              disc_f)
# Part 6

_, ax = plt.subplots(figsize=(16, 12))
plt.axis('off')
plot_points(ax, core_set_pts2010, "r")
plot_points(ax, core_set_pts2017, "b")
d = plt.Circle(disk.get_origin(), disk.get_radius(), color='g', alpha=.8)
ax.add_artist(d)
plt.show()

# Part 7

net = pyscan.my_sample(core_set_pts2017, 200) + pyscan.my_sample(core_set_pts2010, 200)
disc_f = pyscan.DISC
disk, d_val = pyscan.max_disk_scale(net,
                                  [pyscan.WPoint(1.0, p[0], p[1], 1.0) for p in core_set_pts2017],
                                  [pyscan.WPoint(1.0, p[0], p[1], 1.0) for p in core_set_pts2010],
                                  1,
                                  disc_f)

disk2, d_val = pyscan.max_disk_scale(net,
                                  [pyscan.WPoint(1.0, p[0], p[1], 1.0) for p in core_set_pts2017],
                                  [pyscan.WPoint(1.0, p[0], p[1], 1.0) for p in core_set_pts2010],
                                  .5,
                                  disc_f)
_, ax = plt.subplots(figsize=(16, 12))
plt.axis('off')
plot_points(ax, core_set_pts2010, "r")
plot_points(ax, core_set_pts2017, "b")
d = plt.Circle(disk.get_origin(), disk.get_radius(), color='g', alpha=.8)
ax.add_artist(d)
d = plt.Circle(disk2.get_origin(), disk2.get_radius(), color='k', alpha=.8)
ax.add_artist(d)
plt.show()


