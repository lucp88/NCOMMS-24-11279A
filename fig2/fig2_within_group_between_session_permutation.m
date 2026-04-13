%% Figure 2 — Within-Group, Between-Session Permutation Test (Panels E–H)
% Compares FFA-session vs OFA-session time-courses within each Group × ROI.
%   Panel E: EXP group, OFA ROI  (OFA session vs FFA session)
%   Panel F: EXP group, FFA ROI  (no differential effect expected)
%   Panel G: CONT group, OFA ROI (no modulation expected)
%   Panel H: CONT group, FFA ROI (no modulation expected)
%
% Contrast direction:
%   If ROI = FFA: FFA session − OFA session
%   If ROI = OFA: OFA session − FFA session
%   (i.e., target session − non-target session for that ROI)
%
% Statistical test: volume-by-volume permutation test on paired t-statistic
% (one-sample ttest on difference scores), with sign-flipping
% (5,000 permutations).
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

% Group sizes
nExpSubs  = 22;
nConSubs  = 20;

% Acquisition parameters
nrVols    = 57;
runs2graph = 1:7;

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

% Plot colors (from checking_dynamic_FBdisp_inclDerivs.m)
lineColors = [0.28 0.47 1;    % Session: FFA (blue)
              0.9  0.3  0.3]; % Session: OFA (red)
prtColors  = [0.28 0.47 1;    % baseline phase (blue)
              0.9  0.3  0.3;  % regulation phase (red)
              0.14 0.71 0.37]; % task phase (green)

roinames   = {'FFA', 'OFA'};
sessnames  = {'FFA', 'OFA'};
groupNames = {'EXP', 'CONT'};

%% ============================
%  LOAD DATA
%  ============================

% Read from sourceData.xlsx "Fig. 2" tab
sourceFile = fullfile('..', 'sourceData.xlsx');
T = readtable(sourceFile, 'Sheet', 'Fig. 2');

% Extract volume columns
volCols = arrayfun(@(v) sprintf('Vol%d', v), 1:nrVols, 'UniformOutput', false);

% Get unique subject IDs per group
expIDs = unique(T.SubjectID(strcmp(T.Group, 'EXP')));
conIDs = unique(T.SubjectID(strcmp(T.Group, 'CONT')));
allsubs = {1:numel(expIDs), (numel(expIDs)+1):(numel(expIDs)+numel(conIDs))};

% Reconstruct RunsColl: {session}(runs, volumes, ROI, subjects)
allIDs = [expIDs; conIDs];
nSubs  = numel(allIDs);
nRuns  = numel(runs2graph);

RunsColl = cell(1, 2);
for iSess = 1:2
    sessLabel = sessnames{iSess};
    data4D = zeros(nRuns, nrVols, 2, nSubs);

    for iROI = 1:2
        roiLabel = roinames{iROI};
        for iSub = 1:nSubs
            subRows = T(strcmp(T.Session, sessLabel) & ...
                        strcmp(T.ROI, roiLabel) & ...
                        T.SubjectID == allIDs(iSub), :);

            subRows = sortrows(subRows, 'Run');
            for iRun = 1:min(nRuns, height(subRows))
                data4D(iRun, :, iROI, iSub) = ...
                    table2array(subRows(iRun, volCols));
            end
        end
    end
    RunsColl{iSess} = data4D;
end

%% ============================
%  PERMUTATION TESTS & PLOTTING
%  ============================

graphc = 1;
figure('Position', [100 100 600 900]);

for iGroup = 1:numel(groupNames)
    for iROI = 1:2

        % --- FFA session data (condition A) ---
        tmpA = RunsColl{1}(runs2graph, 1:nrVols, iROI, allsubs{iGroup});

        if numel(runs2graph) == 1
            validIdx = squeeze(tmpA(1,1,1,:)) ~= 0;
            tmpA = tmpA(:,:,:,validIdx);
        else
            tmpA(tmpA == 0) = NaN;
            tmpA = nanmean(tmpA, 1);
        end
        tmpA = squeeze(tmpA); % volumes × subjects

        currSMeanA = mean(tmpA, 2)';
        currSMeanA = currSMeanA - mean(currSMeanA(1:basOff1));

        cA     = tmpA - mean(tmpA(1:basOff1, :));
        roi_seA = std(cA, 0, 2)' / sqrt(size(cA, 2));

        % --- OFA session data (condition B) ---
        tmpB = RunsColl{2}(runs2graph, 1:nrVols, iROI, allsubs{iGroup});

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

        title(sprintf('%s group — %s ROI', groupNames{iGroup}, roinames{iROI}));
        xticks([5 28 50]);
        xticklabels({'Baseline', 'Regulation', 'Task'});
        legend('Sess: FFA', 'Sess: OFA', 'AutoUpdate', 'off');

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
        xlim([0 nrVols]);

        % --- Volume-by-volume permutation test (paired, sign-flipping) ---
        for vol = 1:50

            % Compute difference score with directional contrast:
            %   FFA ROI → FFA session − OFA session
            %   OFA ROI → OFA session − FFA session
            % (target session minus non-target session)
            if iROI == 1
                con0 = cA(vol,:) - cB(vol,:);  % FFA sess − OFA sess
            elseif iROI == 2
                con0 = cB(vol,:) - cA(vol,:);  % OFA sess − FFA sess
            end

            % Observed t-statistic (one-sample t-test on difference)
            [h0{graphc, vol}, pVal0{graphc, vol}, ...
             ci0{graphc, vol}, stats0{graphc, vol}] = ttest(con0);

            % Sign-flip permutation
            stats_all = zeros(1, permutCount);
            for iPerm = 1:permutCount
                con_perm = con0 .* sign(randn(1, length(con0)));
                [~, ~, ~, s] = ttest(con_perm);
                stats_all(iPerm) = s.tstat;
            end
            totStats{vol} = stats_all;

            % Permutation p-value
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
                        text(vol, v - tmpSE(idx) - 0.5, '.', 'FontSize', 14);
                    end

                elseif test_right == 1
                    if PPval(graphc, vol) > pThresh2
                        tmpSE = [roi_seA(vol), roi_seB(vol)];
                        [v, idx] = max(currSMean(:, vol));
                        text(vol, v + tmpSE(idx) + 0.25, '*', 'FontSize', 14);
                    elseif PPval(graphc, vol) <= pThresh2 && PPval(graphc, vol) > pThresh4
                        tmpSE = [roi_seA(vol), roi_seB(vol)];
                        [v, idx] = max(currSMean(:, vol));
                        text(vol, v + tmpSE(idx) + 0.5, '.', 'FontSize', 14);
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

sgtitle('Figure 2: Within-Group Between-Session Comparison');
