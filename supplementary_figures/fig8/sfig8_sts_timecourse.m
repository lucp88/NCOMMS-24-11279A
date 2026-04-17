%% Supplementary Figure 8 — STS Time-Course: Between-Session Comparison
% Compares FFA-session vs OFA-session time-courses in right STS for the
% EXP group. STS is a co-activated but non-targeted face-processing region;
% we test whether its BOLD response differs between training sessions.
%
% Contrast direction:
%   FFA session − OFA session (same convention as Fig. 2 FFA ROI panel)
%
% Statistical test: volume-by-volume permutation test on paired t-statistic
% (one-sample t-test on difference scores), with sign-flipping
% (10,000 permutations). Two-sided.
%
% Dependencies: shadedErrorBar.m (Rob Campbell, 2009)
%
% Data source: sourceData.xlsx, tab "SFig. 8-TC"
% Expected format: each row = one run, columns:
%   SubjectID | Group | Session | ROI | Run | Vol1 | Vol2 | ... | Vol57
%   Group: 'EXP' only
%   Session: 'FFA' or 'OFA'
%   ROI: 'STS'
%
% Paper: ReFOFA (NCOMMS-24-11279)
% Author: Lucas Peek
% Last updated: April 2026

clear; clc; close all;

%% ============================
%  SETTINGS
%  ============================

nExpSubs   = 22;
nrVols     = 57;
runs2graph = 1:7;

% Two-sided permutation test
twoSided_test = 1;
permutCount   = 10000;

pThresh1 = 0.025;   % significant (two-sided α/2)
pThresh2 = 0.975;
pThresh3 = 0.05;    % trending
pThresh4 = 0.95;

% Epoch timing (volume indices)
basOff1  = 10;
nfOns1   = 11;
nfOff1   = 42;
taskOns1 = 43;
taskOff1 = 57;

% Plot colors (matching Fig. 2)
lineColors = [0.28 0.47 1;    % FFA session (blue)
              0.9  0.3  0.3]; % OFA session (red)

prtColors = [
    0.95 0.95 0.95;  % baseline
    0.90 0.90 0.90;  % regulation
    0.85 0.85 0.85   % task
];

%% ============================
%  LOAD DATA
%  ============================

% Read from sourceData.xlsx "SFig. 8-TC" tab
sourceFile = fullfile('..', '..', 'sourceData.xlsx');
T = readtable(sourceFile, 'Sheet', 'SFig. 8-TC');

% Volume column names
volCols = arrayfun(@(v) sprintf('Vol%d', v), 1:nrVols, 'UniformOutput', false);

% Subject IDs
subIDs = unique(T.SubjectID);
nSubs  = numel(subIDs);

sessnames = {'FFA', 'OFA'};
nRuns     = numel(runs2graph);

% Reconstruct data array: RunsColl{session}(runs, volumes, subjects)
RunsColl = cell(1, 2);
for iSess = 1:2
    sessLabel = sessnames{iSess};
    data3D = zeros(nRuns, nrVols, nSubs);

    for iSub = 1:nSubs
        subRows = T(strcmp(T.Session, sessLabel) & ...
                     T.SubjectID == subIDs(iSub), :);
        subRows = sortrows(subRows, 'Run');
        for iRun = 1:min(nRuns, height(subRows))
            data3D(iRun, :, iSub) = table2array(subRows(iRun, volCols));
        end
    end
    RunsColl{iSess} = data3D;
end

%% ============================
%  PERMUTATION TEST & PLOTTING
%  ============================

figure('Position', [100 100 700 400]);
hold on;

% --- FFA session data ---
tmpA = RunsColl{1}(runs2graph, 1:nrVols, :);
tmpA(tmpA == 0) = NaN;
tmpA = squeeze(nanmean(tmpA, 1));   % volumes × subjects

currSMeanA = mean(tmpA, 2)';
currSMeanA = currSMeanA - mean(currSMeanA(1:basOff1));

cA     = tmpA - mean(tmpA(1:basOff1, :));
roi_seA = std(cA, 0, 2)' / sqrt(size(cA, 2));

% --- OFA session data ---
tmpB = RunsColl{2}(runs2graph, 1:nrVols, :);
tmpB(tmpB == 0) = NaN;
tmpB = squeeze(nanmean(tmpB, 1));   % volumes × subjects

currSMeanB = mean(tmpB, 2)';
currSMeanB = currSMeanB - mean(currSMeanB(1:basOff1));

cB     = tmpB - mean(tmpB(1:basOff1, :));
roi_seB = std(cB, 0, 2)' / sqrt(size(cB, 2));

% --- Epoch shading ---
yl = [-1 5.5];
ylim(yl);

patch([0 basOff1+1 basOff1+1 0], ...
      [yl(1) yl(1) yl(2) yl(2)], prtColors(1,:), ...
      'LineStyle', 'none', 'FaceAlpha', 0.5, 'HandleVisibility', 'off');
patch([nfOns1 nfOff1+1 nfOff1+1 nfOns1], ...
      [yl(1) yl(1) yl(2) yl(2)], prtColors(2,:), ...
      'LineStyle', 'none', 'FaceAlpha', 0.5, 'HandleVisibility', 'off');
patch([taskOns1 taskOff1 taskOff1 taskOns1], ...
      [yl(1) yl(1) yl(2) yl(2)], prtColors(3,:), ...
      'LineStyle', 'none', 'FaceAlpha', 0.5, 'HandleVisibility', 'off');

% --- SEM-shaded time-courses ---
shadedErrorBar([], mean(cA, 2), roi_seA, ...
    'lineprops', {'-', 'color', lineColors(1,:), 'LineWidth', 2});
shadedErrorBar([], mean(cB, 2), roi_seB, ...
    'lineprops', {'-', 'color', lineColors(2,:), 'LineWidth', 2});

% --- Epoch boundary lines and labels ---
xline(basOff1+1, '--k', 'LineWidth', 1.5);
xline(nfOff1+1, '--k', 'LineWidth', 1.5);

text((basOff1+1)/2, yl(2)*0.92, 'Baseline', ...
    'FontSize', 12, 'HorizontalAlignment', 'center');
text((nfOns1 + nfOff1+1)/2, yl(2)*0.92, 'Regulation', ...
    'FontSize', 12, 'HorizontalAlignment', 'center');
text((taskOns1 + taskOff1)/2, yl(2)*0.92, 'Task', ...
    'FontSize', 12, 'HorizontalAlignment', 'center');

% --- Formatting ---
set(gca, 'FontSize', 12, 'FontName', 'Helvetica');
grid on; box off;
xlim([0 nrVols]);
xlabel('Time (TR)');
ylabel('BOLD Response (a.u.)');
legend({'Sess: FFA', 'Sess: OFA'}, 'Location', 'southeast', 'AutoUpdate', 'off');

% --- Volume-by-volume permutation test (paired, sign-flipping) ---
% Contrast: FFA session − OFA session
for vol = 1:50
    con0 = cA(vol,:) - cB(vol,:);

    % Observed t-statistic
    [h0{vol}, pVal0{vol}, ci0{vol}, stats0{vol}] = ttest(con0);

    % Sign-flip permutation
    stats_all = zeros(1, permutCount);
    for iPerm = 1:permutCount
        con_perm = con0 .* sign(randn(1, length(con0)));
        [~, ~, ~, s] = ttest(con_perm);
        stats_all(iPerm) = s.tstat;
    end

    % Permutation p-value
    PPval(vol) = sum(stats_all <= stats0{vol}.tstat) / permutCount;

    % Two-sided significance markers
    if PPval(vol) < pThresh1 || PPval(vol) > pThresh2
        tmpSE = [roi_seA(vol) roi_seB(vol)];
        [v, idx] = max([currSMeanA(vol) currSMeanB(vol)]);
        text(vol, v + tmpSE(idx) + 0.25, '*', 'FontSize', 14);
    elseif (PPval(vol) <= pThresh3 && PPval(vol) > pThresh1) || ...
           (PPval(vol) <= pThresh2 && PPval(vol) > pThresh4)
        tmpSE = [roi_seA(vol) roi_seB(vol)];
        [v, idx] = max([currSMeanA(vol) currSMeanB(vol)]);
        text(vol, v + tmpSE(idx) + 0.5, '.', 'FontSize', 12);
    end
end

sgtitle('Supplementary Figure 8: STS Between-Session Comparison');
