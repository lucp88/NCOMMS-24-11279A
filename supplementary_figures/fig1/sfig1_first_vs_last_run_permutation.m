%% Supplementary Figure 1 — First vs Last Run Between-Session Permutation
%
% Supports the manuscript statement (NFB Self-regulation Performance):
%   "a focused comparison between the initial and final runs of each
%    training session revealed that FFA activity was significantly more
%    pronounced in the FFA-targeted sessions than in the OFA-targeted ones
%    during the final run, a contrast that was not observed in the initial
%    run."
%
% Replicates the within-group, between-session permutation test used for
% the main Fig. 2 panels E-H, but restricted to a single regulation block
% at a time instead of averaging across runs 1-7. Run 1 (first) and
% Run 7 (last) are tested separately.
%
% Layout: one figure per group (EXP, CONT), each 2 rows x 2 cols.
%     Rows: FFA ROI | OFA ROI
%     Cols: Run 1 (first) | Run 7 (last)
%     Panels labelled a-d per figure.
%
% Contrast direction (same as main Fig. 2):
%   FFA ROI -> FFA session - OFA session
%   OFA ROI -> OFA session - FFA session
%   (i.e., target session minus non-target session for that ROI)
%
% Statistical test: volume-by-volume permutation test on paired t-statistic
% (one-sample ttest on difference scores), with sign-flipping
% (10,000 permutations).
%
% Styling follows Supplementary Figure 7 (STS time-course) for visual
% consistency across supplementary time-course figures.
%
% Dependencies: shadedErrorBar.m (Rob Campbell, 2009) -- located in this dir
%
% Data source: sourceData.xlsx, tab "Fig. 2"
%
% Author: Lucas Peek
% Last updated: April 2026

clear; clc; close all;

%% ============================
%  SETTINGS
%  ============================

% Group sizes
nExpSubs  = 22;
nConSubs  = 20;

% Acquisition parameters
nrVols    = 57;

% Run selections for the two columns
runSelections = {1, 7};
runLabels     = {'Run 1 (first)', 'Run 7 (last)'};
nCols         = numel(runSelections);

% Testing parameters
oneSided_test = 0;
    test_left   = 0;
    test_right  = 1;
twoSided_test = 1;

permutCount = 10000;

% Thresholds
if oneSided_test == 1
    pThresh1 = 0.05;
    pThresh2 = 0.95;
    pThresh3 = 0.10;
    pThresh4 = 0.90;
elseif twoSided_test == 1
    pThresh1 = 0.025;
    pThresh2 = 0.975;
    pThresh3 = 0.05;
    pThresh4 = 0.95;
end

% Epoch timing (volume indices)
basOff1  = 10;
nfOns1   = 11;
nfOff1   = 42;
taskOns1 = 43;
taskOff1 = 57;

% Plot colors (matching Fig. 2 and SFig. 7)
lineColors = [0.28 0.47 1;    % FFA session (blue)
              0.9  0.3  0.3]; % OFA session (red)

% Grey-scale epoch shadings (matching SFig. 7)
prtColors = [
    0.95 0.95 0.95;  % baseline
    0.90 0.90 0.90;  % regulation
    0.85 0.85 0.85   % task
];

roinames   = {'FFA', 'OFA'};
sessnames  = {'FFA', 'OFA'};
groupNames = {'EXP', 'CONT'};

% Per-figure layout: 2 rows (ROIs) x 2 cols (runs), one figure per group
nROIs = numel(roinames);

% Panel letters for each 2x2 figure (reset per group)
panelLetters = {'a','b'; 'c','d'};

% Shared y-axis limits across all panels (comparability)
yl = [-2.5 5.5];

%% ============================
%  LOAD DATA
%  ============================

sourceFile = fullfile('..', '..', 'sourceData.xlsx');
T = readtable(sourceFile, 'Sheet', 'Fig. 2');

volCols = arrayfun(@(v) sprintf('Vol%d', v), 1:nrVols, 'UniformOutput', false);

expIDs = unique(T.SubjectID(strcmp(T.Group, 'EXP')));
conIDs = unique(T.SubjectID(strcmp(T.Group, 'CONT')));
allsubs = {1:numel(expIDs), (numel(expIDs)+1):(numel(expIDs)+numel(conIDs))};

allIDs   = [expIDs; conIDs];
nSubs    = numel(allIDs);
nRunsMax = 7;

RunsColl = cell(1, 2);
for iSess = 1:2
    sessLabel = sessnames{iSess};
    data4D = zeros(nRunsMax, nrVols, 2, nSubs);

    for iROI = 1:2
        roiLabel = roinames{iROI};
        for iSub = 1:nSubs
            subRows = T(strcmp(T.Session, sessLabel) & ...
                        strcmp(T.ROI, roiLabel) & ...
                        T.SubjectID == allIDs(iSub), :);

            subRows = sortrows(subRows, 'Run');
            for iRun = 1:min(nRunsMax, height(subRows))
                data4D(iRun, :, iROI, iSub) = ...
                    table2array(subRows(iRun, volCols));
            end
        end
    end
    RunsColl{iSess} = data4D;
end

%% ============================
%  PERMUTATION TESTS & PLOTTING — one figure per group (2 rows x 2 cols)
%  ============================

% Pre-allocate containers for outputs
h0     = cell(numel(groupNames), nROIs, nCols, 50);
pVal0  = cell(numel(groupNames), nROIs, nCols, 50);
ci0    = cell(numel(groupNames), nROIs, nCols, 50);
stats0 = cell(numel(groupNames), nROIs, nCols, 50);
PPval  = zeros(numel(groupNames), nROIs, nCols, 50);

for iGroup = 1:numel(groupNames)

    figure('Position', [100 100 900 600], 'Color', 'w');

    for iROI = 1:nROIs
        for iCol = 1:nCols

            runs2graph = runSelections{iCol};
            nRunsSel   = numel(runs2graph);

            % --- FFA session (condition A) ---
            tmpA = RunsColl{1}(runs2graph, 1:nrVols, iROI, allsubs{iGroup});
            if nRunsSel == 1
                validIdx = squeeze(tmpA(1,1,1,:)) ~= 0;
                tmpA = tmpA(:,:,:,validIdx);
            else
                tmpA(tmpA == 0) = NaN;
                tmpA = nanmean(tmpA, 1);
            end
            tmpA = squeeze(tmpA);                            % volumes x subjects
            cA      = tmpA - mean(tmpA(1:basOff1, :));       % baseline-corrected
            meanA   = mean(cA, 2)';
            roi_seA = std(cA, 0, 2)' / sqrt(size(cA, 2));

            % --- OFA session (condition B) ---
            tmpB = RunsColl{2}(runs2graph, 1:nrVols, iROI, allsubs{iGroup});
            if nRunsSel == 1
                validIdx = squeeze(tmpB(1,1,1,:)) ~= 0;
                tmpB = tmpB(:,:,:,validIdx);
            else
                tmpB(tmpB == 0) = NaN;
                tmpB = nanmean(tmpB, 1);
            end
            tmpB = squeeze(tmpB);
            cB      = tmpB - mean(tmpB(1:basOff1, :));
            meanB   = mean(cB, 2)';
            roi_seB = std(cB, 0, 2)' / sqrt(size(cB, 2));

            % -----------------------------------------------------------
            %  PLOT
            % -----------------------------------------------------------
            spIdx = (iROI - 1) * nCols + iCol;
            ax = subplot(nROIs, nCols, spIdx); hold on;

            ylim(yl);

            % Grey-scale epoch shading (SFig. 7 style)
            patch([0 basOff1+1 basOff1+1 0], ...
                  [yl(1) yl(1) yl(2) yl(2)], prtColors(1,:), ...
                  'LineStyle', 'none', 'FaceAlpha', 0.5, ...
                  'HandleVisibility', 'off');
            patch([nfOns1 nfOff1+1 nfOff1+1 nfOns1], ...
                  [yl(1) yl(1) yl(2) yl(2)], prtColors(2,:), ...
                  'LineStyle', 'none', 'FaceAlpha', 0.5, ...
                  'HandleVisibility', 'off');
            patch([taskOns1 taskOff1 taskOff1 taskOns1], ...
                  [yl(1) yl(1) yl(2) yl(2)], prtColors(3,:), ...
                  'LineStyle', 'none', 'FaceAlpha', 0.5, ...
                  'HandleVisibility', 'off');

            % Zero baseline (subtle horizontal reference)
            plot([0 nrVols], [0 0], '-', ...
                 'Color', [0.5 0.5 0.5], 'LineWidth', 0.6, ...
                 'HandleVisibility', 'off');

            % SEM band + mean line (shadedErrorBar)
            hA = shadedErrorBar([], meanA, roi_seA, ...
                'lineprops', {'-', 'color', lineColors(1,:), ...
                              'LineWidth', 2});
            hB = shadedErrorBar([], meanB, roi_seB, ...
                'lineprops', {'-', 'color', lineColors(2,:), ...
                              'LineWidth', 2});

            % Epoch boundary lines (dashed black, SFig. 7 style)
            xline(basOff1+1, '--k', 'LineWidth', 1.5, 'HandleVisibility', 'off');
            xline(nfOff1+1,  '--k', 'LineWidth', 1.5, 'HandleVisibility', 'off');

            % In-plot phase labels (SFig. 7 style) — only on top row
            if iROI == 1
                text((basOff1+1)/2, yl(2)*0.92, 'Baseline', ...
                    'FontSize', 11, 'HorizontalAlignment', 'center');
                text((nfOns1 + nfOff1+1)/2, yl(2)*0.92, 'Regulation', ...
                    'FontSize', 11, 'HorizontalAlignment', 'center');
                text((taskOns1 + taskOff1)/2, yl(2)*0.92, 'Task', ...
                    'FontSize', 11, 'HorizontalAlignment', 'center');
            end

            % Axis styling (SFig. 7 style)
            set(ax, 'FontSize', 12, 'FontName', 'Helvetica');
            grid on; box off;
            xlim([0 nrVols]);

            % Axis labels — bottom row x-label, left column y-label
            if iROI == nROIs
                xlabel('Time (TR)', 'FontSize', 12, 'FontName', 'Helvetica');
            end
            if iCol == 1
                ylabel({[roinames{iROI} ' ROI']; 'BOLD Response (a.u.)'}, ...
                       'FontSize', 12, 'FontName', 'Helvetica');
            end

            % Column title only on top row
            if iROI == 1
                title(runLabels{iCol}, 'FontSize', 12, ...
                      'FontName', 'Helvetica', 'FontWeight', 'bold');
            end

            % Panel letter in upper-left corner
            text(0.02, 0.96, panelLetters{iROI, iCol}, ...
                 'Units', 'normalized', 'FontSize', 14, ...
                 'FontName', 'Helvetica', 'FontWeight', 'bold', ...
                 'HorizontalAlignment', 'left', ...
                 'VerticalAlignment', 'top');

            % Legend only in panel a
            if iROI == 1 && iCol == 1
                legend([hA.mainLine, hB.mainLine], ...
                       {'Sess: FFA', 'Sess: OFA'}, ...
                       'Location', 'southeast', 'AutoUpdate', 'off', ...
                       'FontSize', 10);
            end

            % -----------------------------------------------------------
            %  PERMUTATION TEST (paired, sign-flipping)
            % -----------------------------------------------------------
            for vol = 1:50

                if iROI == 1
                    con0 = cA(vol,:) - cB(vol,:);
                elseif iROI == 2
                    con0 = cB(vol,:) - cA(vol,:);
                end

                [h0{iGroup, iROI, iCol, vol}, ...
                 pVal0{iGroup, iROI, iCol, vol}, ...
                 ci0{iGroup, iROI, iCol, vol}, ...
                 stats0{iGroup, iROI, iCol, vol}] = ttest(con0);

                stats_all = zeros(1, permutCount);
                for iPerm = 1:permutCount
                    con_perm = con0 .* sign(randn(1, length(con0)));
                    [~, ~, ~, s] = ttest(con_perm);
                    stats_all(iPerm) = s.tstat;
                end

                PPval(iGroup, iROI, iCol, vol) = ...
                    sum(stats_all <= stats0{iGroup, iROI, iCol, vol}.tstat) ...
                    / permutCount;

                % -------------------------------------------------------
                %  SIGNIFICANCE MARKERS (follow the curves, SFig. 7 style)
                % -------------------------------------------------------
                sems = [roi_seA(vol), roi_seB(vol)];

                if oneSided_test == 1
                    if test_left == 1
                        [v, idx] = min([meanA(vol), meanB(vol)]);
                        if PPval(iGroup, iROI, iCol, vol) < pThresh1
                            text(vol, v - sems(idx) - 0.25, '*', ...
                                 'FontSize', 14, ...
                                 'HorizontalAlignment', 'center');
                        elseif PPval(iGroup, iROI, iCol, vol) <= pThresh3 && ...
                               PPval(iGroup, iROI, iCol, vol) > pThresh1
                            text(vol, v - sems(idx) - 0.5, '.', ...
                                 'FontSize', 12, ...
                                 'HorizontalAlignment', 'center');
                        end
                    elseif test_right == 1
                        [v, idx] = max([meanA(vol), meanB(vol)]);
                        if PPval(iGroup, iROI, iCol, vol) > pThresh2
                            text(vol, v + sems(idx) + 0.25, '*', ...
                                 'FontSize', 14, ...
                                 'HorizontalAlignment', 'center');
                        elseif PPval(iGroup, iROI, iCol, vol) <= pThresh2 && ...
                               PPval(iGroup, iROI, iCol, vol) > pThresh4
                            text(vol, v + sems(idx) + 0.5, '.', ...
                                 'FontSize', 12, ...
                                 'HorizontalAlignment', 'center');
                        end
                    end

                elseif twoSided_test == 1
                    [v, idx] = max([meanA(vol), meanB(vol)]);
                    if PPval(iGroup, iROI, iCol, vol) < pThresh1 || ...
                       PPval(iGroup, iROI, iCol, vol) > pThresh2
                        text(vol, v + sems(idx) + 0.25, '*', ...
                             'FontSize', 14, ...
                             'HorizontalAlignment', 'center');
                    elseif (PPval(iGroup, iROI, iCol, vol) <= pThresh3 && ...
                            PPval(iGroup, iROI, iCol, vol) > pThresh1) || ...
                           (PPval(iGroup, iROI, iCol, vol) <= pThresh2 && ...
                            PPval(iGroup, iROI, iCol, vol) > pThresh4)
                        text(vol, v + sems(idx) + 0.5, '.', ...
                             'FontSize', 12, ...
                             'HorizontalAlignment', 'center');
                    end
                end
            end
        end
    end

    % Per-figure super-title
    sgtitle(sprintf(['Supplementary Figure 1: First vs last regulation ' ...
                     'run — %s group'], groupNames{iGroup}), ...
            'FontSize', 13, 'FontName', 'Helvetica');

    set(gcf, 'PaperPositionMode', 'auto');
end
