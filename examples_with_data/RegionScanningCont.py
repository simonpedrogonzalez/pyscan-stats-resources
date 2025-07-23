# Part 1: Setup and Data Loading

import geopandas as gpd
import pandas as pd
import matplotlib.pyplot as plt
import pyscan


def plot_points(ax, pts, c):
    xs = [pt[0] for pt in pts]
    ys = [pt[1] for pt in pts]
    ax.scatter(xs, ys, color=c, marker='.')
6
def plot_points_traj(ax, pts, c):
    xs = [pt[0] for pt in pts]
    ys = [pt[1] for pt in pts]
    ax.plot(xs, ys, color=c)

def plot_approx(ax, regions, core_set_pts):
    for reg in regions:
        plot_points_traj(ax, reg, "g")
    plot_points(ax, core_set_pts, "b")
    ax.set_axis_off()


# Load county shapefile
gdf = gpd.read_file("data/cb_2017_us_county_500k.shp")

# Filter to continental US
gdf = gdf.cx[-124.84:-66.9, 24.396:49.4].copy()

# Load population estimates (2010 and 2017)
df = pd.read_csv("data/co-est2019-alldata.csv", encoding='latin-1',
                 dtype={'STATE': str, 'COUNTY': str})
df['FIPS'] = df['STATE'] + df['COUNTY']
df['POPESTIMATE2010'] = pd.to_numeric(df['POPESTIMATE2010'], errors='coerce').fillna(0).astype(int)
df['POPESTIMATE2017'] = pd.to_numeric(df['POPESTIMATE2017'], errors='coerce').fillna(0).astype(int)

population2010 = df.set_index('FIPS')['POPESTIMATE2010'].to_dict()
population2017 = df.set_index('FIPS')['POPESTIMATE2017'].to_dict()

# Part 2: Geometry Extraction and Weight Assignment

regions = []
weights2010 = []
weights2017 = []

for _, row in gdf.iterrows():
    fips = row['GEOID']
    if fips not in population2010 or fips not in population2017:
        continue

    geom = row.geometry

    # Simplify geometry to reduce runtime
    geom = geom.simplify(0.1, preserve_topology=True)

    if geom.geom_type == 'Polygon':
        parts = [geom.exterior.coords]
    elif geom.geom_type == 'MultiPolygon':
        largest = max(geom.geoms, key=lambda g: g.area)
        parts = [largest.exterior.coords]
    else:
        continue

    points = [(x, y) for part in parts for (x, y) in part]

    # Drop counties with any point outside bounds
    if any(not (-124.84 <= x <= -66.9 and 24.396 <= y <= 49.4) for x, y in points):
        continue

    pt_objs = [pyscan.Point(x, y, 1.0) for x, y in points]
    regions.append(pt_objs)
    weights2010.append(population2010[fips])
    weights2017.append(population2017[fips])

# n = 500
# avg_region_size = sum(len(r) for r in regions) / len(regions)
# print(f"Average region size: {avg_region_size:.2f} points")
# regions = regions[:n]
# weights2010 = weights2010[:n]
# weights2017 = weights2017[:n]
# Part 3: Disk Scan (Region Comparison)

disc_f = pyscan.DISC
alpha = 0.2        # spatial error
r_min = 0.4 # 0.4        # min radius (24 miles)
r_max = 4.0        # max radius (240 miles)

print("HERE")
# This takes a while
disk, value = pyscan.max_disk_region(
    regions, regions,
    weights2010,
    regions, weights2017,
    r_min, r_max, alpha, disc_f
)

# Part 4: Plot the Result

_, ax = plt.subplots(figsize=(16, 12))
plt.axis('off')
plot_approx(ax, regions, [])
d = plt.Circle(disk.get_origin(), disk.get_radius(), color='g', alpha=.8)
ax.add_artist(d)
# plt.savefig("disk_region.png", bbox_inches='tight', pad_inches=0.1)
plt.show()
