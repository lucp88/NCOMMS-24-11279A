"""
Figure 7 — Structural Equation Models (SEM) of neurofeedback regulation pathways.

Fits a path model (semopy) to trial-level data for the EXP and CONT groups
separately, testing how self-regulation activity (FFA/OFA) propagates through
task-evoked detection and recognition signals to behavioral performance.

Model specification:
  - Regulation and Training_Day predict self-regulation (regFFA_s, regOFA_s)
  - Self-regulation + Session predict task-evoked detection activity (FFA, OFA)
  - Self-regulation + Session + detection activity predict task-evoked
    recognition activity (FFA, OFA)
  - Self-regulation + task-evoked activity + Session + Training_Day predict
    behavioral performance (detection_frms, recognition_frms)
  - Covariances specified for: within-layer ROI pairs, cross-task within-ROI
    pairs, and detection–recognition performance

Outputs:
  - Path diagram images for EXP and CONT groups (PNG)
  - Model fit statistics printed to console
  - Model reports saved to working directory

Dependencies:
  pip install semopy pandas openpyxl graphviz
  brew install graphviz  (system-level Graphviz for the 'dot' executable)

Data source: sourceData.xlsx, tab "Fig. 7"

Author: Lucas Peek
Last updated: March 2026
"""

import os
import semopy as sem
import pandas as pd
from edit_path_diagrams import edit_path_graph

# ============================================================
#  SETTINGS
# ============================================================

# Path to source data (relative to this script's directory)
SOURCE_FILE = os.path.join('..', 'sourceData.xlsx')
SHEET_NAME = 'Fig. 7'

# Output directory for figures (same as script directory)
# Note: if running interactively (e.g., Spyder), use os.getcwd() instead
OUTPUT_DIR = os.getcwd()

# ============================================================
#  MODEL SPECIFICATION
# ============================================================

model_desc = """
    # Self-regulation activity predicted by session and training day
    regOFA_s, regFFA_s ~ Regulation + Training_Day

    # Task-evoked detection activity predicted by self-regulation and session
    task_detect_FFA_beta_s ~ regOFA_s + regFFA_s + Regulation
    task_detect_OFA_beta_s ~ regOFA_s + regFFA_s + Regulation

    # Task-evoked recognition activity predicted by self-regulation, session,
    # and detection activity (mediation path)
    task_recogn_FFA_beta_s ~ regOFA_s + regFFA_s + Regulation + task_detect_FFA_beta_s
    task_recogn_OFA_beta_s ~ regOFA_s + regFFA_s + Regulation + task_detect_OFA_beta_s

    # Behavioral performance predicted by self-regulation, task activity,
    # session, and training day
    detection_frms ~ regOFA_s + regFFA_s + task_detect_FFA_beta_s + task_detect_OFA_beta_s + Regulation + Training_Day
    recognition_frms ~ regOFA_s + regFFA_s + task_recogn_FFA_beta_s + task_recogn_OFA_beta_s + Regulation + Training_Day

    # Variances
    Regulation ~~ Regulation
    regOFA_s ~~ regOFA_s
    regFFA_s ~~ regFFA_s
    task_detect_FFA_beta_s ~~ task_detect_FFA_beta_s
    task_detect_OFA_beta_s ~~ task_detect_OFA_beta_s
    task_recogn_FFA_beta_s ~~ task_recogn_FFA_beta_s
    task_recogn_OFA_beta_s ~~ task_recogn_OFA_beta_s
    Training_Day ~~ Training_Day
    detection_frms ~~ detection_frms
    recognition_frms ~~ recognition_frms

    # Covariances: within-layer ROI pairs
    regFFA_s ~~ regOFA_s
    task_detect_FFA_beta_s ~~ task_detect_OFA_beta_s
    task_recogn_FFA_beta_s ~~ task_recogn_OFA_beta_s

    # Covariances: cross-task within-ROI (detection ↔ recognition)
    task_detect_FFA_beta_s ~~ task_recogn_FFA_beta_s
    task_detect_OFA_beta_s ~~ task_recogn_OFA_beta_s

    # Covariance: behavioral outcomes
    detection_frms ~~ recognition_frms
"""

# Variables and connections to hide in the path diagram for clarity
vars2remove = ['Training_Day']
connections2remove = [
    'task_detect_FFA_beta_s -> task_recogn_FFA_beta_s',
    'task_detect_OFA_beta_s -> task_recogn_OFA_beta_s',
    'detection_frms -> recognition_frms'
]
changeVarNames = pd.DataFrame([])

# ============================================================
#  LOAD DATA
# ============================================================

data_all = pd.read_excel(SOURCE_FILE, sheet_name=SHEET_NAME)
data_exp = data_all[data_all['Group'] == 'EXP'].copy()
data_cont = data_all[data_all['Group'] == 'CONT'].copy()

# ============================================================
#  FIT & PLOT: EXP GROUP
# ============================================================

print('=' * 60)
print('EXP GROUP')
print('=' * 60)

model_exp = sem.ModelMeans(model_desc, intercepts=False)
res_exp = model_exp.fit(data_exp, groups=['subID'], obj='ML')

# Generate and save path diagram
g_orig_exp, g_new_exp = edit_path_graph(
    model_exp, vars2remove, connections2remove,
    changeVarNames, std_est=True
)
g_new_exp.render(
    os.path.join(OUTPUT_DIR, 'fig7_path_diagram_EXP'),
    format='png', cleanup=True
)

insp = model_exp.inspect(std_est=True)
print(insp[insp['op'] == '~'].to_string())

# ============================================================
#  FIT & PLOT: CONT GROUP
# ============================================================

print('=' * 60)
print('CONT GROUP')
print('=' * 60)

model_cont = sem.ModelMeans(model_desc, intercepts=False)
res_cont = model_cont.fit(data_cont, groups=['subID'], obj='ML')

# Generate and save path diagram
g_orig_cont, g_new_cont = edit_path_graph(
    model_cont, vars2remove, connections2remove,
    changeVarNames, std_est=True
)
g_new_cont.render(
    os.path.join(OUTPUT_DIR, 'fig7_path_diagram_CONT'),
    format='png', cleanup=True
)

print('=' * 60)
print(f'Path diagrams saved to: {OUTPUT_DIR}')
print('  - fig7_path_diagram_EXP.png')
print('  - fig7_path_diagram_CONT.png')
print('=' * 60)
