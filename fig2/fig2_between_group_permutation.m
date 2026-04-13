%% Figure 2 — Between-Group Permutation Test (Panels A–D)
% Compares EXP vs CONT group time-courses for each Session × ROI.
%   Panel A: FFA session, FFA ROI (target)
%   Panel B: FFA session, OFA ROI (non-target)
%   Panel C: OFA session, FFA ROI (non-target)
%   Panel D: OFA session, OFA ROI (target)
%
% Statistical test: volume-by-volume permutation test on independent-
% samples t-statistic (ttest2), with label-shuffling (5,000 permutations).
%
% Dependencies: shadedErrorBar.m (Rob Campbell, 2009)
%
% Data source: sourceData.xlsx, tab "Fig. 2"
% Expected format: each row = one run, columns:
%   SubjectID | Group | Session | ROI | Run | Vol1 | Vol2 | ... | Vol57
%   Group: 'EXP' or 'CONT'
%   Session: 'FFA' or 'OFA'
%   ROI: 'FFA' or 'OFA'
%
% Author: Lucas Peek
% Last updated: March 2026

clear; clc; close all;

%% ============================
%  SETTINGS (do not change)
%  ============================

% Group indices (after sorting: EXP first, then CONT)
nExpSubs  = 22;
nConSubs  = 20;

% Acquisition parameters
nrVols    = 57;
runs2graph = 1:7;

% Testing parameters
oneSided_test = 0;
    test_left   = 0;  % condition of interest < contrast condition
    test_right  = 1;  % condition of interest > contrast condition
twoSided_test = 1;

permutCount = 10000;

% Thresholds
if oneSided_test == 1
    pThresh1 = 0.05;   % significance threshold
    pThresh2 = 0.95;
    pThresh3 = 0.10;   % trending threshold
    pThresh4 = 0.90;
elseif twoSided_test == 1
    pThresh1 = 0.025;
    pThresh2 = 0.975;
    pThresh3 = 0.05;
    pThresh4 = 0.95;
end

% Epoch timing (volume indices)
basOff1  = 10;   % baseline ends at volume 10
nfOns1   = 11;   % regulation onset
nfOff1   = 42;   % regulation offset
taskOns1 = 43;   % task onset
taskOff1 = 57;   % task offset (= nrVols)

% Plot colors (from checking_dynamic_FBdisp_inclDerivs.m)
lineColors = [0.28 0.47 1;    % EXP group (blue)
              0.9  0.3  0.3]; % CONT group (red)
prtColors  = [0.28 0.47 1;    % baseline phase (blue)
              0.9  0.3  0.3;  % regulation phase (red)
              0.14 0.71 0.37]; % task phase (green)

roinames = {'FFA', 'OFA'};
sessnames = {'FFA', 'OFA'};

%% ============================
%  LOAD DATA
%  ============================

% Read from sourceData.xlsx "Fig. 2" tab
sourceFile = fullfile('..', 'sourceData.xlsx');
T = readtable(sourceFile, 'Sheet', 'Fig. 2');

% Extract volume columns (Vol1 through Vol57)
volCols = arrayfun(@(v) sprintf('Vol%d', v), 1:nrVols, 'UniformOutput', false);

% Get unique subject IDs per group
expIDs = unique(T.SubjectID(strcmp(T.Group, 'EXP')));
conIDs = unique(T.SubjectID(strcmp(T.Group, 'CONT')));

% Reconstruct RunsColl: {session}(runs, volumes, ROI, subjects)
% This matches the original data structure from the .mat file
RunsColl = cell(1, 2);
for iSess = 1:2
    sessLabel = sessnames{iSess};
    allIDs = [expIDs; conIDs];
    nSubs = numel(allIDs);
    nRuns = numel(runs2graph);

    data4D = zeros(nRuns, nrVols, 2, nSubs);

    for iROI = 1:2
        roiLabel = roinames{iROI};
        for iSub = 1:nSubs
            subRows = T(strcmp(T.Session, sessLabel) & ...
                        strcmp(T.ROI, roiLabel) & ...
                        T.SubjectID == allIDs(iSub), :);

            % Sort by run number and select runs2graph
            subRows = sortrows(subRows, 'Run');
            for iRun = 1:min(nRuns, height(subRows))
                data4D(iRun, :, iROI, iSub) = ...
                    table2array(subRows(iRun, volCols));
            end
        end
    end
    RunsColl{iSess} = data4D;
end

% Define subject index ranges (EXP = 1:nExp, CONT = nExp+1:end)
expSubs = 1:numel(expIDs);
conSubs = (numel(expIDs)+1):(numel(expIDs)+numel(conIDs));

%% ============================
%  PERMUTATION TESTS & PLOTTING
%  ============================

graphc = 1;
figure('Position', [100 100 600 900]);

for iSess = 1:2
    for iROI = 1:2

        % --- EXP group data ---
        tmpA = RunsColl{iSess}(runs2graph, 1:nrVols, iROI, expSubs);

        if numel(runs2graph) == 1
            % Single run: remove subjects with missing data (zeros)
            validIdx = squeeze(tmpA(1,1,1,:)) ~= 0;
            tmpA = tmpA(:,:,:,validIdx);
        else
            % Multiple runs: replace zeros with NaN, average across runs
            tmpA(tmpA == 0) = NaN;
            tmpA = nanmean(tmpA, 1);
        end
        tmpA = squeeze(tmpA); % volumes × subjects

        % Baseline-correct (mean of first basOff1 volumes)
        currSMeanA = mean(tmpA, 2)';
        currSMeanA = currSMeanA - mean(currSMeanA(1:basOff1));

        cA     = tmpA - mean(tmpA(1:basOff1, :));
        roi_seA = std(cA, 0, 2)' / sqrt(size(cA, 2));

        % --- CONT group data ---
        tmpB = RunsColl{iSess}(runs2graph, 1:nrVols, iROI, conSubs);

        if numel(runs2graph) == 1
            validIdx = squeeze(tmpB(1,1,1,:)) ~= 0;
            tmpB = tmpB(:,:,:,validIdx);
        else
            tmpB(tmpB == 0) = NaN;
            tmpB = nanmean(tmpB, 1);
        end
        tmpB = squeeze(tmpB);

        currSMeanB = mean(tmpB, 2)';
        currSMeanB = currSMeanB - mean(currSMeanB(1:basOff1));

        cB     = tmpB - mean(tmpB(1:basOff1, :));
        roi_seB = std(cB, 0, 2)' / sqrt(size(cB, 2));

        % --- Combined mean for plotting ---
        currSMean = [currSMeanA; currSMeanB];

        % --- Plot ---
        subplot(4, 1, graphc); hold on;
        h = plot(currSMean', 'LineWidth', 2);
        set(h, {'color'}, {lineColors(1,:); lineColors(2,:)});

        title(sprintf('%s session — %s ROI', sessnames{iSess}, roinames{iROI}));
        xticks([5 28 50]);
        xticklabels({'Baseline', 'Regulation', 'Task'});
        legend('EXP', 'CONT', 'AutoUpdate', 'off');

        ylimits = ylim;

        % Phase shading
        patch([0 basOff1+1 basOff1+1 0], ...
              [ylimits(1) ylimits(1) ylimits(2) ylimits(2)], ...
              prtColors(1,:), 'LineStyle', 'none');
        patch([nfOns1 nfOff1+1 nfOff1+1 nfOns1], ...
              [ylimits(1) ylimits(1) ylimits(2) ylimits(2)], ...
              prtColors(2,:), 'LineStyle', 'none');
        patch([taskOns1 taskOff1 taskOff1 taskOns1], ...
              [ylimits(1) ylimits(1) ylimits(2) ylimits(2)], ...
              prtColors(3,:), 'LineStyle', 'none');
        alpha(0.2);

        % SEM shading
        shadedErrorBar([], mean(cA, 2), roi_seA, ...
            'lineprops', {'-', 'color', lineColors(1,:), 'LineWidth', 2});
        shadedErrorBar([], mean(cB, 2), roi_seB, ...
            'lineprops', {'-', 'color', lineColors(2,:), 'LineWidth', 2});

        grid on;
        ylim([ylimits(1) 4]);
        xlim([0 taskOff1]);

        % --- Volume-by-volume permutation test ---
        for vol = 1:50

            % Observed t-statistic (independent samples)
            [h0{graphc, vol}, pVal0{graphc, vol}, ...
             ci0{graphc, vol}, stats0{graphc, vol}] = ttest2(cA(vol,:), cB(vol,:));

            % Build pooled data and shuffle group labels
            group1 = cA(vol,:);  % EXP
            group2 = cB(vol,:);  % CONT
            pool   = [group1, group2];
            labels = [ones(1, length(group1)), 2*ones(1, length(group2))];

            stats_all = zeros(1, permutCount);
            for iPerm = 1:permutCount
                rand_labels = labels(randperm(numel(labels)));
                grp1_perm = pool(rand_labels == 1);
                grp2_perm = pool(rand_labels == 2);
                [~, ~, ~, s] = ttest2(grp1_perm, grp2_perm);
                stats_all(iPerm) = s.tstat;
            end
            totStats{vol} = stats_all;

            % Permutation p-value: proportion of null distribution
            % at or below the observed t-statistic
            PPval(graphc, vol) = sum(stats_all <= stats0{graphc, vol}.tstat) / permutCount;

            % Mark significance on plot
            if oneSided_test == 1
                if test_left == 1
                    if PPval(graphc, vol) < pThresh1
                        tmpSE = [roi_seA(vol), roi_seB(vol)];
                        [v, idx] = min(currSMean(:, vol));
                        text(vol, v - tmpSE(idx) - 0.25, '*', 'FontSize', 14);
                    elseif PPval(graphc, vol) <= pThresh3 && PPval(graphc, vol) > pThresh1
                        tmpSE = [roi_seA(vol), roi_seB(vol)];
                        [v, idx] = min(currSMean(:, vol));
                        text(vol, v - tmpSE(idx) - 0.5, '.', 'FontSize', 12);
                    end

                elseif test_right == 1
                    if PPval(graphc, vol) > pThresh2
                        tmpSE = [roi_seA(vol), roi_seB(vol)];
                        [v, idx] = max(currSMean(:, vol));
                        text(vol, v + tmpSE(idx) + 0.25, '*', 'FontSize', 14);
                        alpha(0.2);
                    elseif PPval(graphc, vol) <= pThresh2 && PPval(graphc, vol) > pThresh4
                        tmpSE = [roi_seA(vol), roi_seB(vol)];
                        [v, idx] = max(currSMean(:, vol));
                        text(vol, v + tmpSE(idx) + 0.5, '.', 'FontSize', 14, 'FontWeight', 'normal');
                        alpha(0.1);
                    end
                end

            elseif twoSided_test == 1
                if PPval(graphc, vol) < pThresh1 || PPval(graphc, vol) > pThresh2
                    tmpSE = [roi_seA(vol), roi_seB(vol)];
                    [v, idx] = max(currSMean(:, vol));
                    text(vol, v + tmpSE(idx) + 0.25, '*', 'FontSize', 14);
                    alpha(0.1);
                elseif (PPval(graphc, vol) <= pThresh3 && PPval(graphc, vol) > pThresh1) || ...
                       (PPval(graphc, vol) <= pThresh2 && PPval(graphc, vol) > pThresh4)
                    tmpSE = [roi_seA(vol), roi_seB(vol)];
                    [v, idx] = max(currSMean(:, vol));
                    text(vol, v + tmpSE(idx) + 0.5, '.', 'FontSize', 12);
                    alpha(0.1);
                end
            end
        end

        graphc = graphc + 1;
    end
end

sgtitle('Figure 2: Between-Group Comparison (EXP vs CONT)');
