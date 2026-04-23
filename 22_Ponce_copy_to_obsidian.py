# pylint: disable=line-too-long
# pylint: disable=invalid-name

# %%
import shutil
from pathlib import Path

# Your paths (unchanged; pathlib handles spaces)
FKB12 = Path('/home/lefnui/CRISPRlysis_outputs/Ponce_unpublished/fkb12_260423/')
RBCS1 = Path('/home/lefnui/CRISPRlysis_outputs/Ponce_unpublished/rbcs1_260423/')
RSP3 = Path('/home/lefnui/CRISPRlysis_outputs/Ponce_unpublished/rsp3_260423/')
FKB12_obsidian = Path('/mnt/g/My Drive/Obsidian/Obsidian/Projects/Edinburgh/02_Cas9_Chlamydomonas/Images/FKB12/')
RBCS1_obsidian = Path('/mnt/g/My Drive/Obsidian/Obsidian/Projects/Edinburgh/02_Cas9_Chlamydomonas/Images/RBCS1/')
RSP3_obsidian = Path('/mnt/g/My Drive/Obsidian/Obsidian/Projects/Edinburgh/02_Cas9_Chlamydomonas/Images/RSP3/')

list_of_prefixes = ['01_1', '01_2', '02_3', '02_5', '02_7', '03_1', '05_1', '05_2', '05_3']

# Create target dirs if missing
for obsidian_dir in [FKB12_obsidian, RBCS1_obsidian, RSP3_obsidian]:
    obsidian_dir.mkdir(parents=True, exist_ok=True)

# Copy using shutil (best: no shell, preserves metadata)
for prefix in list_of_prefixes:
    for file in FKB12.glob(prefix + '*.png'):
        shutil.copy2(file, FKB12_obsidian)
    for file in RBCS1.glob(prefix + '*.png'):
        shutil.copy2(file, RBCS1_obsidian)
    for file in RSP3.glob(prefix + '*.png'):
        shutil.copy2(file, RSP3_obsidian)