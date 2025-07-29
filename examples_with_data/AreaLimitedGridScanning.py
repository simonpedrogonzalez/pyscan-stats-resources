import pyscan
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import matplotlib.patches as patches
import time

def xr_to_coords(ds):
    ds_stacked = ds.stack(xy=['x', 'y'])
    ds_non_nan = ds_stacked.dropna('xy')
    xs = ds.x.values
    x_non_nan = ds_non_nan.x.values
    ys = ds.y.values
    y_non_nan = ds_non_nan.y.values
    values_non_nan = ds_non_nan.value.values
    non_nan_array = np.column_stack((x_non_nan, y_non_nan, values_non_nan))
    return non_nan_array, np.unique(xs), np.unique(ys)

def plot_result(b_arr, m_arr, extents, rect=None, rect2=None):
    print('Plotting...')
    non_nan_b = b_arr[~np.isnan(b_arr)]
    non_nan_m = m_arr[~np.isnan(m_arr)]
    vmin = min(non_nan_b.min(), non_nan_m.min(), np.abs(non_nan_b - non_nan_m).min())
    vmax = max(non_nan_b.max(), non_nan_m.max(), np.abs(non_nan_b - non_nan_m).max())
    fig, axs = plt.subplots(1, 2, figsize=(6, 3))

    # Plot Baseline (b)
    im1 = axs[0].imshow(b_arr, cmap='coolwarm', vmin=vmin, vmax=vmax, extent=extents)
    axs[0].set_title('Baseline')

    # Plot Measure (m)
    im2 = axs[1].imshow(m_arr, cmap='coolwarm', vmin=vmin, vmax=vmax, extent=extents)
    axs[1].set_title('Measure')

    if rect is not None:
        for ax in axs:
            rect_patch = patches.Rectangle(
                (rect.lowX(), rect.lowY()), rect.upX() - rect.lowX(), rect.upY() - rect.lowY(),
                linewidth=2, edgecolor='k', facecolor='none'
            )
            ax.add_patch(rect_patch)
        axs[0].text(rect.lowX() -1, rect.lowY() - 1, 'Area Limited \n Subgrid \n (750 cells)', color='k', fontsize=11)

    if rect2 is not None:
        for ax in axs:
            rect_patch = patches.Rectangle(
                (rect2.lowX(), rect2.lowY()), rect2.upX() - rect2.lowX(), rect2.upY() - rect2.lowY(),
                linewidth=2, edgecolor='r', facecolor='none'
            )
            ax.add_patch(rect_patch)
        axs[0].text(rect2.lowX(), rect2.lowY() - 0.5, 'Unlimited Subgrid', color='darkred', fontsize=11)

    plt.tight_layout()
    plt.show()



ds_october = xr.open_dataset(f"data/pm_utah_2016_october.nc")
ds_january = xr.open_dataset(f"data/pm_utah_2016_january.nc")

max_slice = min(len(ds_october.x), len(ds_october.y))
# max_slice = 100
ds_m = ds_october.isel(x=slice(0, max_slice), y=slice(0, max_slice))
ds_b = ds_january.isel(x=slice(0, max_slice), y=slice(0, max_slice))
extents = [ds_october.x.min(), ds_october.x.max(), ds_october.y.min(), ds_october.y.max()]

m_arr, unique_xs, unique_ys = xr_to_coords(ds_m)
b_arr, _, _ = xr_to_coords(ds_b)
    
grid = pyscan.Grid(unique_xs, unique_ys, m_arr, b_arr)

disc_f = pyscan.RKULLDORF

max_area = 750 # given by the number of cells incldued

print('Limited grid size run...')
t0 = time.time()
subgrid = pyscan.max_subgrid_convex(grid, 1e-5, disc_f, max_area)
t1 = time.time()
print('Time in seconds:', (t1-t0))
rect = grid.toRectangle(subgrid)

cell_size_x = ds_october.x[1] - ds_october.x[0]
cell_size_y = ds_october.y[1] - ds_october.y[0]

col0, col1, row0, row1 = subgrid.lowCol(), subgrid.upCol(), subgrid.lowRow(), subgrid.upRow()

width = col1 - col0 + 1
height = row1 - row0 + 1
fvalue = subgrid.fValue()

xB, xU, yB, yU = rect.lowX(), rect.upX(), rect.lowY(), rect.upY()
# Rect size in number of cells
rect_size_area = width * height
print(f'Fvalue: {fvalue}')
print(f'Width: {width} cells, Height: {height} cells')
print(f'Area: {rect_size_area} cells')
print(f'Area limit: {max_area} cells')
print(f'Rectangle: ({col0}, {col1}), ({row0}, {row1})')

# REGULAR
disc_f = pyscan.RKULLDORF
print('Unlimited grid size run...')
t0 = time.time()
subgrid = pyscan.max_subgrid_convex(grid, 1e-5, disc_f)
t1 = time.time()
print('Time in seconds:', (t1-t0))
rect2 = grid.toRectangle(subgrid)

plot_result(ds_january.value.values, ds_october.value.values, extents, rect, rect2)
print('Done.')