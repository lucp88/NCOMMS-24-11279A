"""
Figure 4 — Whole-brain group-level activation map (Regulation > Baseline).

Projects a volumetric SPM contrast (EXP group, Regulation > Baseline,
FWE-corrected at p < .05) onto the fsaverage cortical surface for each
hemisphere and opens interactive 3D views in the browser. The resulting
renderings are used as the basis for the whole-brain panels in Figure 4.

Inputs:
  - reg_bas_FWE_05_exp.nii : second-level contrast from SPM.
    Already hard-thresholded at FWE p < .05 (t >= 7.69).
    Non-significant voxels are NaN; significant voxels carry their t-value.
    Located one directory up from this script (shared across figures).

Outputs:
  - Two interactive browser-based surface renderings (left and right pial
    hemispheres), opened via nilearn's `view_surf`. No files are written
    automatically; use the browser's export/screenshot function to save.

Rendering choices:
  - Mesh:    pial surface (folded cortex, matches the paper figure).
             Swap `fsaverage.pial_*` for `fsaverage.infl_*` to show the
             inflated surface instead (sulci visible as curvature shading).
  - Background: `fsaverage.sulc_*` (continuous sulcal depth). Gives the
             smooth grey gyri/sulci gradient seen in the paper. Using the
             binary curvature sign instead produces a harsher cartoony look.
  - Colormap: 'inferno' (dark purple -> orange -> yellow). Swap to 'hot'
             or 'afmhot' for a warmer red-yellow palette.
  - Threshold: `threshold` is a *display* threshold only — vertices with
             |value| < threshold are rendered transparent so the grey brain
             shows through. It does NOT re-threshold the underlying data.
             Since the input map is already FWE-thresholded at t >= 7.69,
             the display threshold only controls what counts as "background"
             during surface interpolation edge-blending.

Dependencies:
  pip install nilearn nibabel numpy

Data source: reg_bas_FWE_05_exp.nii (one level up, shared across figures)

Author: Lucas Peek
Last updated: April 2026
"""

import os
import numpy as np
import nibabel as nib
from nilearn import plotting, datasets, surface

# ============================================================
#  SETTINGS
# ============================================================

# Path to the NIfTI contrast (one directory above this script)
STAT_FILE = os.path.join('..', 'reg_bas_FWE_05_exp.nii') #reg_bas_FWE_05_exp_zeroed.nii

# Display parameters
THRESHOLD = 0.05        # display threshold (transparent below this value)
VMIN      = 0
VMAX      = 14          # colormap saturation cap (data max is ~14)
CMAP      = 'hot'       # colormap used for activations
BLACK_BG  = True        # black background matches the paper aesthetic

# ============================================================
#  LOAD TEMPLATE & STAT MAP
# ============================================================

# fsaverage surface meshes (pial, inflated, sulcal depth, curvature) for
# both hemispheres, downloaded on first use.
fsaverage = datasets.fetch_surf_fsaverage(mesh='fsaverage')

# Load the volumetric contrast map (MNI152 1mm space).
# Non-significant voxels are NaN; `np.nan_to_num` replaces them with 0 so
# that `surface.vol_to_surf` can safely interpolate near cluster edges.
stat_img = nib.load(STAT_FILE)
data = np.nan_to_num(stat_img.get_fdata(), nan=0)
stat_img = nib.Nifti1Image(data, stat_img.affine, stat_img.header)

# ============================================================
#  PROJECT VOLUME ONTO CORTICAL SURFACE
# ============================================================

# `vol_to_surf` samples the volumetric map along each cortical vertex,
# turning 3D voxel data into a 1D per-vertex texture on the pial surface.
# Linear interpolation is used by default; switch to `interpolation='nearest_most_frequent'`
# for a crisper, less blurred projection at cluster edges.
texture_L = surface.vol_to_surf(stat_img, fsaverage.pial_left, interpolation='linear')
texture_R = surface.vol_to_surf(stat_img, fsaverage.pial_right, interpolation='linear')
# texture_L = surface.vol_to_surf(stat_img, fsaverage.pial_left)
# texture_R = surface.vol_to_surf(stat_img, fsaverage.pial_right)

# ============================================================
#  RENDER: LEFT HEMISPHERE
# ============================================================

view_L = plotting.view_surf(
    surf_mesh=fsaverage.pial_left,     # folded pial surface
    surf_map=texture_L,                # projected stat texture
    bg_map=fsaverage.sulc_left,        # continuous sulcal depth shading
    threshold=THRESHOLD,
    # vmin=VMIN,
    # vmax=VMAX,
    colorbar=True,
    black_bg=BLACK_BG,
    cmap=CMAP,
    symmetric_cmap=False               # one-sided map (no negatives)
)
view_L.open_in_browser()

# ============================================================
#  RENDER: RIGHT HEMISPHERE
# ============================================================

view_R = plotting.view_surf(
    surf_mesh=fsaverage.pial_right,
    surf_map=texture_R,
    bg_map=fsaverage.sulc_right,
    threshold=THRESHOLD,
    # vmin=VMIN,
    # vmax=VMAX,
    colorbar=True,
    black_bg=BLACK_BG,
    cmap=CMAP,
    symmetric_cmap=False
)
view_R.open_in_browser()

# ============================================================
#  RENDER: BOTH HEMISPHERES - statmap matches manuscript exactly - 
# ============================================================
view = plotting.view_img_on_surf(
    stat_img,
    surf_mesh='fsaverage',
    threshold=THRESHOLD,
    cmap=CMAP,
    symmetric_cmap=False,
    black_bg=BLACK_BG,
    # view='left', 
)
view.open_in_browser()
