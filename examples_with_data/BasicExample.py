import pyscan
import matplotlib.pyplot as plt
import numpy as np
import random

def get_coord(i, lst):
    return [pt[i] for pt in lst]


# l = (a, b, c) where ax + by + c = 0 => (-ax -c) / b = y
def f(mr, x):
        l = mr.get_coords()
        return -(x * l[0] + l[2]) / l[1]

def plot_plane(red, blue, max_region, eps):

    n = np.rint(1 / eps)
    s = np.rint(1 / (2 * eps * eps))
    # create a region containing 5% of the points. Inside of this region points are more likely to be red.
    #print(get_coord(1, red))
    _, ax = plt.subplots(figsize=(16, 12))
    ax.scatter(get_coord(0, red), get_coord(1, red), marker=".", c="r")
    ax.scatter(get_coord(0, blue), get_coord(1, blue), marker=".", c="b")

    net = pyscan.my_sample(red, n) + pyscan.my_sample(blue, n)
    red_s = pyscan.my_sample(red, s)
    blue_s = pyscan.my_sample(blue, s)
    approx_region, _ = pyscan.max_halfplane(net, red_s, blue_s, pyscan.KULLDORF)

    ax.scatter(get_coord(0, net), get_coord(1, net), marker="x", c="k")
    ax.plot([0, 1], [f(max_region, 0), f(max_region, 1)], c="k", linewidth=4.0)
    ax.plot([0, 1], [f(approx_region, 0), f(approx_region, 1)], c="g", linewidth=4.0)
    plt.xlim([0, 1])
    plt.ylim([0, 1])
    plt.axis('off')
    plt.show()



# generate some random points
pts = [pyscan.WPoint(1.0, random.random(), random.random(), 1.0) for i in range(10000)]
red, blue, max_region = pyscan.plant_halfplane(pts, .05, .5, .8)
_, ax = plt.subplots(figsize=(18, 12))

ax.scatter(get_coord(0, red), get_coord(1, red), marker=".", c="r")
ax.scatter(get_coord(0, blue), get_coord(1, blue), marker=".", c="b")
ax.plot([0, 1], [f(max_region, 0), f(max_region, 1)], c="k", linewidth=4.0)
ax.set_xlim([0, 1])
ax.set_ylim([0, 1])
plt.axis('off')
plt.show()

plot_plane(red, blue, max_region, .2)
plot_plane(red, blue, max_region, .1)
plot_plane(red, blue, max_region, .05)
plot_plane(red, blue, max_region, .025)
plot_plane(red, blue, max_region, .01)

red, blue, max_region = pyscan.plant_halfplane(pts, .05, .5, .6)

_, ax = plt.subplots(figsize=(16, 12))
plt.scatter(get_coord(0, red), get_coord(1, red), marker=".", c="r")
plt.scatter(get_coord(0, blue), get_coord(1, blue), marker=".", c="b")
plt.plot([0, 1], [f(max_region, 0), f(max_region, 1)], c="k", linewidth=4.0)
plt.xlim([0, 1])
plt.ylim([0, 1])
plt.axis("off")
plt.show()

plot_plane(red, blue, max_region, .025)
plot_plane(red, blue, max_region, .005)
