# Part 1

import pyscan
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import random

# Unzip phila_crime_incidents.zip
# Data source https://www.kaggle.com/datasets/eamonntweedy/philadelphia-crime-incident-data-2006-2022?resource=download

def get_coord(i, lst):
    return [pt[i] for pt in lst]

baseline = []
measured = []

dtype = {
    "text_general_code": str,
    "lon": "float64",
    "lat": "float64",
}

df = pd.read_csv(
    "data/phila_crime_incidents.csv/phila_crime_incidents.csv",
    dtype=dtype
)

df = df[df["lon"].between(-75.3, -74.95) & df["lat"].between(39.8, 40.2)]
df = df[df["text_general_code"].str.contains("thefts", case=False, na=False)]
# Take 50000 random samples
df = df.sample(n=50000, random_state=42)

# Convert the DataFrame to a list of tuples
for index, row in df.iterrows():
    x = row["lon"]
    y = row["lat"]
    point = pyscan.WPoint(1.0, x, y, 1.0)
    baseline.append(point)
    measured.append(point)


# with open("data/phila_crime_incidents.csv/phila_crime_incidents.csv", "r") as f:
#     reader = csv.DictReader(f)
#     numAnomalies = 0
#     for row in reader:
#         try:
#             x = float(row["lon"])
#             y = float(row["lat"])
#             if not (39.8 < y < 40.2 and -76 < x):
#                 continue
#         except:
#             continue

#         point = pyscan.WPoint(1.0, x, y, 1.0)
#         baseline.append(point)
#         if "thefts" in row["text_general_code"].lower():
#             measured.append(point)

# Part 2            

mixed_pts = list(zip(measured, ["r"] * len(measured))) + list(zip(baseline, ["b"] * len(baseline)))
random.shuffle(mixed_pts)
xs = get_coord(0, [p[0] for p in mixed_pts])
ys = get_coord(1, [p[0] for p in mixed_pts])
cs = [p[1] for p in mixed_pts]
f, ax = plt.subplots(figsize=(16, 12))
ax.set_axis_off()
ax.scatter(xs, ys, color=cs, marker='.')
plt.show()

# Part 3

net = pyscan.my_sample(measured, 100) + pyscan.my_sample(baseline, 100)
f, ax = plt.subplots(figsize=(16, 12))
ax.set_axis_off()
ax.scatter(get_coord(0, net), get_coord(1, net), marker='.')
plt.show()


# Part 4

#take the samples
m_sample = pyscan.my_sample(measured, 5000)
b_sample = pyscan.my_sample(baseline, 5000)

#plot the samples
mixed_samples = list(zip(m_sample, ["r"] * len(m_sample))) + list(zip(b_sample, ["b"] * len(b_sample)))
random.shuffle(mixed_samples)
xs = get_coord(0, [p[0] for p in mixed_samples])
ys = get_coord(1, [p[0] for p in mixed_samples])
cs = [p[1] for p in mixed_samples]
_, ax = plt.subplots(figsize=(16, 12))
ax.set_axis_off()
ax.scatter(xs, ys, color=cs, marker='.')
plt.show()

# Part 5



disc_f = pyscan.KULLDORF
halfplane, h_val = pyscan.max_halfplane(net, m_sample, b_sample, disc_f)


# l = (a, b, c) where ax + by + c = 0 => (-ax -c) / b = y
def f(mr, x):
    l = mr.get_coords()
    return -(x * l[0] + l[2]) / l[1]

_, ax = plt.subplots(figsize=(16, 12))
ax.set_axis_off()
ax.scatter(xs, ys, color=cs, marker='.')
ax.plot([-75.3, -74.95], [f(halfplane, -75.2), f(halfplane, -74.95)])
plt.show()

# Part 6


disk, d_val = pyscan.max_disk(net, m_sample, b_sample, disc_f)

# Plot disk
_, ax = plt.subplots(figsize=(16, 12))
ax.set_axis_off()
ax.scatter(xs, ys, color=cs, marker='.')
d = plt.Circle(disk.get_origin(), disk.get_radius(), color='g', alpha=.8)
ax.add_artist(d)
plt.show()


# Part 7

# all of these methods create a grid first and then scan the grid using various methods.
grid = pyscan.Grid(100, m_sample, b_sample)

# a slow exact method
subgrid = pyscan.max_subgrid(grid, disc_f)
rect = grid.toRectangle(subgrid)

_, ax = plt.subplots(figsize=(16, 12))
ax.set_axis_off()
ax.scatter(xs, ys, color=cs, marker='.')
r = plt.Rectangle((rect.lowX(), rect.lowY()), rect.upX() - rect.lowX(), rect.upY() - rect.lowY(),
                  alpha=.8,
                  color='g')
ax.add_artist(r)
plt.show()

# Part 8

# a faster approximate method. Recommend you use the RKULLDORF method for this since the standard Kulldorf
# function can be a bit unstable for small and very large regions.
disc_f = pyscan.RKULLDORF

subgrid = pyscan.max_subgrid_convex(grid, .01, disc_f)
rect = grid.toRectangle(subgrid)

_, ax = plt.subplots(figsize=(16, 12))
ax.set_axis_off()
ax.scatter(xs, ys, color=cs, marker='.')
r = plt.Rectangle((rect.lowX(), rect.lowY()), rect.upX() - rect.lowX(), rect.upY() - rect.lowY(),
                  alpha=.8,
                  color='g')
ax.add_artist(r)
plt.show()


