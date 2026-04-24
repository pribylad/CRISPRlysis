# %%
# Take png 01_1*, 02_2, 02_6, 03_1_1, 03_1_2, 05_3_DEL_end_summary_norm_violinplot from folders
# fkb12_WT_260423\
# fkb12_ku80_260423\
# fkb12_mre11_260423\
# fkb12_mulam_260423\
# fkb12_pollam_260423\
# fkb12_polmu_260423\
# fkb12_polq_260423\
# fkb12_polqku80_260423\
# and print images from each folder to one line in a single figure, and save the figure as png and pdf
# %%
import matplotlib.pyplot as plt
from pathlib import Path
from PIL import Image

main_dir = Path("/home/lefnui/CRISPRlysis_outputs/Ponce_unpublished/")
# target = 'fkb12'
target = 'rbcs1'
data_sources = [
    f"{target}_WT_260423",
    f"{target}_ku70_260423",
    f"{target}_ku80_260423",
    f"{target}_mre11_260423",
    f"{target}_mulam_260423",
    f"{target}_pollam_260423",
    f"{target}_polmu_260423",
    f"{target}_polq_260423"
]
charts_output_dir = main_dir
output_format = "png"
# output_format = "pdf"
# %%
# Define the list of image filenames to be included in the figure
# Use glob to find the actual image files matching the patterns in each data source folder
image_filenames = []
patterns = [
    "01_1_*.png",
    "02_2_*.png",
    "02_6_*.png",
    "03_1_1_*.png",
    "03_1_2_*.png",
    "05_3_DEL_end_summary_norm_violinplot*.png"
]
for pattern in patterns:
    for data_source in data_sources:
        image_filenames.extend(list((main_dir / data_source).glob(pattern)))
image_filenames
# %%
# Sort the image filenames based on the parent folder name (data source) and then by the filename itself to ensure they are in the correct order

image_filenames.sort(key=lambda x: (x.parent.name, x.name))
image_filenames
# %%
# Create a figure with no of rows equal to the number of data sources and no of columns based on len(image_filenames) for the images
num_rows = len(data_sources)
num_cols = len(patterns)
fig, axes = plt.subplots(num_rows, num_cols, figsize=(num_cols*4, num_rows*4))
# Loop through the sorted image filenames and plot them in the corresponding subplot
for i, image_file in enumerate(image_filenames):
    row = i // num_cols
    col = i % num_cols
    img = Image.open(image_file)
    axes[row, col].imshow(img)
    axes[row, col].axis('off')  # Hide the axes
    # Set the title of the subplot to the parent folder name (data source)
    axes[row, col].set_title(image_file.parent.name, fontsize=10)
# Adjust layout and save the figure
plt.tight_layout()
plt.savefig(charts_output_dir / f'{target}_combined_figures.{output_format}', dpi=300, bbox_inches='tight')
plt.savefig(charts_output_dir / f'{target}_combined_figures.pdf', dpi=300, bbox_inches='tight')
# %%
