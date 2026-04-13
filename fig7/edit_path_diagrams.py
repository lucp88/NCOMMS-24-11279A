"""
edit_path_diagrams.py — Helper for customizing semopy path diagrams.

Takes a fitted semopy model and produces a modified Graphviz graph object with:
  - Significant paths (p < 0.05) drawn in black, with line width scaled to
    effect size and style indicating direction (solid = positive, dotted = negative).
    Significance stars replace p-values (* p<.05, ** p<.01, *** p<.001).
  - Non-significant paths (p >= 0.05) drawn as invisible lines (hidden from view).
  - Optional removal of specified variables and/or connections from the diagram.
  - Optional renaming of variable labels for publication-ready figures.

Dependencies:
  - semopy (pip install semopy)
  - graphviz (pip install graphviz)
  - Graphviz system package (brew install graphviz)

Author: Lucas Peek
Created: February 2022
Last updated: March 2026
"""

import re
import semopy as sem


def edit_path_graph(model, vars2exclude, connections2exclude,
                    changeVarNames, std_est=False):
    """
    Generate and customize a semopy path diagram.

    Parameters
    ----------
    model : semopy.Model or semopy.ModelMeans
        A fitted SEM model.
    vars2exclude : list of str
        Variable names to remove from the diagram entirely
        (e.g., ['Training_Day']).
    connections2exclude : list of str
        Path connections to remove, formatted as 'source -> target'
        (e.g., ['detection_frms -> recognition_frms']).
    changeVarNames : pd.DataFrame
        DataFrame with 'old' and 'new' columns for renaming variables.
        Pass an empty DataFrame to skip renaming.
    std_est : bool, optional
        If True, display standardized estimates on paths. Default: False.

    Returns
    -------
    graphObject : graphviz.Digraph
        The original, unmodified path diagram.
    graphObject_new : graphviz.Digraph
        The customized path diagram with styling applied.
    """

    # Generate the base path diagram from semopy
    if std_est:
        graphObject = sem.semplot(model, "tmp.png", std_ests=True)
    else:
        graphObject = sem.semplot(model, "tmp.png", std_ests=False)

    graphObject_new = graphObject.copy()
    body = graphObject_new.body

    # Line width scaling range
    maxLineWidth = 5
    minLineWidth = 0.25

    # --- Step 1: Remove excluded variables ---
    # Any line in the graph body containing an excluded variable is removed
    if vars2exclude:
        blacklist = []
        for var in vars2exclude:
            for it, line in enumerate(body):
                if re.findall(var, line):
                    blacklist.append(it)
        body = [l for i, l in enumerate(body) if i not in blacklist]

    # --- Step 2: Rename variables ---
    if not changeVarNames.empty:
        for idx in range(len(changeVarNames)):
            body = [line.replace(changeVarNames.old[idx],
                                 changeVarNames.new[idx])
                    for line in body]

    # --- Step 3: Remove excluded connections ---
    if connections2exclude:
        blacklist = []
        for con in connections2exclude:
            for it, line in enumerate(body):
                if re.findall(con, line):
                    blacklist.append(it)
        body = [l for i, l in enumerate(body) if i not in blacklist]

    # --- Step 4: Style paths based on significance and effect direction ---

    # Extract all path estimates from the graph body for scaling
    pat = r'(?<=label=")(-?\d+\.\-?\d+)|(-?\d+)(?=\\\\np-val)'
    ests = re.findall(pat, str(body))
    ests = [abs(float(e[0])) for e in ests]

    limitUp = max(ests)
    limitLow = min(ests)

    pattern = r'[*>]'
    for it, line in enumerate(body):
        check1 = re.findall(pattern, line)
        check2 = re.findall(r'\"(.*)\\np', line)

        if check1 and check2:
            effect = float(re.findall(r'\"(.*)\\np', line)[0])
            pVal_str = re.findall(r'\:(.*)"', line)[0]
            pVal = float(pVal_str)

            # Line style: solid for positive, dotted for negative effects
            style = 'solid' if effect > 0 else 'dotted'

            # Scale line width to effect magnitude
            lineWeight = str(
                ((abs(effect) - limitLow) / (limitUp - limitLow)) *
                (maxLineWidth - minLineWidth) + minLineWidth
            )

            if pVal < 0.05:
                # Significant: draw in black with significance stars
                color = 'black'
                if pVal <= 0.001:
                    sigStars = '***'
                elif pVal <= 0.01:
                    sigStars = '**'
                else:
                    sigStars = '*'

                line = line.replace('p-val:' + pVal_str, sigStars)
                endOfLine = line.find(']')
                newLine = (line[:endOfLine] +
                           ', style = ' + style +
                           ', penwidth = ' + lineWeight +
                           ', color = ' + color +
                           line[endOfLine:])
            else:
                # Non-significant: hide the path (invisible)
                color = 'grey'
                style = 'invisible'

                toRemove = re.findall(r'\[(.*)]', line)
                line = line.replace(toRemove[0], '')
                endOfLine = line.find('[')
                newLine = (line[:endOfLine + 1] +
                           'style = ' + style +
                           ', penwidth = ' + lineWeight +
                           ', color = ' + color +
                           line[endOfLine + 1:])

            body[it] = newLine

    graphObject_new.body = body
    return graphObject, graphObject_new
