% =========================================================================
% PoseAnalysisMaster.m  (CSV edition -- combined_summary.csv format)
%
% Scatter/fit analysis comparing experimental pose-landing frequencies
% against the six CSA/CRSA algorithm variants, reading the flat
% combined_summary.csv (or a single *_sim_summary.csv) format written by
% chute_drop_batch.py's write_condition_summary() -- the same file format
% consumed by plot_pose_probabilities.m.
%
%   CSA_Dyn   = CSA,  transitioning poses zeroed, remaining renormalized  ("before")  <- geom_CSA_B_pct
%   CSA_Stat  = CSA,  transitioning poses redistribute to destination pose ("after")  <- geom_CSA_A_pct
%   CSA_Nan   = CSA,  no filtering at all (raw CSA weight normalized)     ("none")   <- geom_CSA_N_pct
%   CRSA_Dyn  = CRSA, transitioning poses zeroed, remaining renormalized  ("before")  <- geom_CRSA_B_pct
%   CRSA_Stat = CRSA, transitioning poses redistribute to destination pose ("after")  <- geom_CRSA_A_pct
%   CRSA_Nan  = CRSA, no filtering at all (raw CRSA weight normalized)     ("none")   <- geom_CRSA_N_pct
%
% UNMATCHED POSES (pose_id == -1) are always excluded from this analysis:
% there is no algorithmic value for them to scatter against, so they
% cannot contribute to a per-pose algo-vs-experimental comparison.
%
% FLAGGED-UNSTABLE POSES (flagged_unstable_by_geometry == true, i.e. the
% geometric analysis noted TRANSITIONS / FLOOR-UNSTABLE for that pose) are
% handled specially:
%   - EXCLUDED from every statistic (MSE, Pearson r, R^2, the best-fit
%     line itself, and all per-part breakdowns) -- a pose the geometry
%     itself flagged as not a legitimate final resting pose shouldn't be
%     used to judge how well the algorithm predicts legitimate resting
%     poses.
%   - STILL PLOTTED on the scatter (so you can see where they land), but
%     rendered in light grey instead of their part color, to visually
%     separate them from the points that drove the fit/stats.
%
% OUTPUT FILES (all saved to OutputDir):
%   PoseAnalysisMaster_Report.txt                    -- scatter/fit stats table (overall + per-part)
%   PoseAnalysisMaster_ScatterFit_<AlgoName>.png (x6) -- one standalone figure per
%                                                        algorithm variant (algo vs exp probability)
%
% POINT COLOR CONVENTION (scatter/fit figures):
%   Each part gets its own color, applied consistently across all six
%   figures, at low opacity (0.25) with small markers -- EXCEPT poses
%   flagged unstable by geometry, which are always light grey regardless
%   of part, and excluded from all statistics (see above).
%   The gray dashed line is the y=x calibration reference; the red solid
%   line is the pooled least-squares best fit (fit to stable poses only).
%
% LEGENDS (per figure):
%   - Top-left: data legend (calibration line, scatter note, best-fit
%     equation, and the flagged-unstable marker note if any are present)
%   - Bottom-right: part-color legend, drawn individually on every figure
%     (no shared/top legend across figures)
%
% --------------------------------------------------------------------
% USAGE
% --------------------------------------------------------------------
%   PoseAnalysisMaster('combined_summary.csv')
%   PoseAnalysisMaster('combined_summary.csv', 'OutputDir', 'analysis_out')
%   PoseAnalysisMaster('combined_summary.csv', 'ApplyModifier', true)
%   PoseAnalysisMaster('combined_summary.csv', 'Parts', {'Df4a','Qf4i'})
%
% --------------------------------------------------------------------
% NAME-VALUE OPTIONS
% --------------------------------------------------------------------
%   'OutputDir'     Folder to save the report + figures into. Created if
%                   it doesn't exist. Default: a 'PoseAnalysisMaster_out'
%                   subfolder next to summaryPath.
%   'Parts'         Cell array of part names, in the order you want them
%                   colored/legended. Default: all unique values in the
%                   'part' column, sorted alphabetically.
%   'PartColors'    Nx3 RGB matrix, one row per entry in Parts (in the
%                   same order). Default: a built-in palette, cycled if
%                   there are more parts than palette rows.
%   'ApplyModifier' true/false. If true, applies a per-pose multiplicative
%                   modifier to EACH of the six geom_*_pct columns,
%                   independently, within each (part, alpha_deg, beta_deg)
%                   condition, and renormalizes -- exactly the same
%                   post-hoc renormalization used in
%                   plot_pose_probabilities.m, since geom_*_pct is already
%                   normalized per condition and the normalization
%                   constant cancels algebraically:
%                     P_i * mod_i / sum_j(P_j * mod_j)
%                       = raw_i * mod_i / sum_j(raw_j * mod_j)
%                   Default modifier: mod_i = 1/mean(ratio_wall_i, ratio_floor_i).
%                   Default: false.
%   'ModifierFcn'   Optional function handle @(ratioWall, ratioFloor) -> mod_i.
%                   Ignored if ApplyModifier is false.
%                   Default: @(rw, rf) 1 ./ mean([rw, rf], 2, 'omitnan')
% =========================================================================

function PoseAnalysisMaster(summaryPath, varargin)

close all; clc;

% -------------------------------------------------------------------------
% OPTIONS
% -------------------------------------------------------------------------
p = inputParser;
addRequired(p, 'summaryPath', @(x) ischar(x) || isstring(x));
addParameter(p, 'OutputDir', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'Parts', {}, @iscell);
addParameter(p, 'PartColors', [], @(x) isnumeric(x) && (isempty(x) || size(x,2)==3));
addParameter(p, 'ApplyModifier', false, @islogical);
addParameter(p, 'ModifierFcn', @(rw, rf) 1 ./ mean([rw, rf], 2, 'omitnan'), ...
    @(x) isa(x, 'function_handle'));
parse(p, summaryPath, varargin{:});
opt = p.Results;

UNSTABLE_COLOR = [0.75 0.75 0.75];   % light grey for flagged-unstable poses

if ~isfile(opt.summaryPath)
    error('PoseAnalysisMaster:fileNotFound', 'Could not find file: %s', opt.summaryPath);
end

if isempty(opt.OutputDir)
    [srcDir, ~, ~] = fileparts(opt.summaryPath);
    if isempty(srcDir), srcDir = pwd; end
    resultsRoot = fullfile(srcDir, 'PoseAnalysisMaster_out');
else
    resultsRoot = char(opt.OutputDir);
end
if ~exist(resultsRoot, 'dir'), mkdir(resultsRoot); end

% Six CSA/CRSA variants. Order drives figure/save order and matches the
% original xlsx-based convention (B->Dyn "before", A->Stat "after",
% N->Nan "none").
algoNames  = {'CSA_Dyn','CSA_Stat','CSA_Nan','CRSA_Dyn','CRSA_Stat','CRSA_Nan'};
geomCols   = {'geom_CSA_B_pct','geom_CSA_A_pct','geom_CSA_N_pct', ...
              'geom_CRSA_B_pct','geom_CRSA_A_pct','geom_CRSA_N_pct'};
nAlgos     = 6;

% -------------------------------------------------------------------------
% LOAD DATA
% -------------------------------------------------------------------------
fprintf('Loading %s...\n', opt.summaryPath);
T = readtable(opt.summaryPath, 'TextType', 'string');

requiredCols = {'part','alpha_deg','beta_deg','pose_id','sim_freq_pct'};
for i = 1:numel(requiredCols)
    if ~ismember(requiredCols{i}, T.Properties.VariableNames)
        error('PoseAnalysisMaster:missingColumn', ...
            'Expected column "%s" not found in %s. Is this a combined_summary.csv / *_sim_summary.csv written by chute_drop_batch.py?', ...
            requiredCols{i}, opt.summaryPath);
    end
end
for i = 1:nAlgos
    if ~ismember(geomCols{i}, T.Properties.VariableNames)
        error('PoseAnalysisMaster:missingGeomColumn', ...
            'Expected column "%s" not found in %s.', geomCols{i}, opt.summaryPath);
    end
end
if opt.ApplyModifier
    for c = {'ratio_wall','ratio_floor'}
        if ~ismember(c{1}, T.Properties.VariableNames)
            error('PoseAnalysisMaster:missingModifierColumn', ...
                'ApplyModifier=true requires column "%s", not found in %s.', c{1}, opt.summaryPath);
        end
    end
end

% Drop unmatched poses (pose_id == -1): no algorithmic value to scatter
% them against.
nBefore = height(T);
T = T(T.pose_id ~= -1, :);
fprintf('Dropped %d unmatched (pose_id == -1) row(s); %d row(s) remain.\n', ...
    nBefore - height(T), height(T));

if isempty(T)
    disp('No matched rows to analyze. Exiting.');
    return;
end

% -------------------------------------------------------------------------
% FLAGGED-UNSTABLE MASK (kept in the plot, excluded from stats)
% -------------------------------------------------------------------------
unstableAll = false(height(T), 1);
if ismember('flagged_unstable_by_geometry', T.Properties.VariableNames)
    fu = T.flagged_unstable_by_geometry;
    if iscell(fu) || isstring(fu)
        unstableAll = strcmpi(string(fu), 'true') | strcmpi(string(fu), '1');
    else
        unstableAll = logical(fu);
    end
else
    warning('PoseAnalysisMaster:noUnstableColumn', ...
        '"flagged_unstable_by_geometry" column not found -- treating all poses as stable.');
end
fprintf('Flagged-unstable rows: %d / %d (excluded from stats, still plotted in grey).\n\n', ...
    sum(unstableAll), height(T));

% -------------------------------------------------------------------------
% PARTS / COLORS
% -------------------------------------------------------------------------
if isempty(opt.Parts)
    partNames = cellstr(unique(T.part, 'sorted'));
else
    partNames = cellfun(@char, opt.Parts, 'UniformOutput', false);
end

defaultPalette = [ ...
    0.85 0.10 0.10;   % red
    0.95 0.55 0.10;   % orange
    0.90 0.80 0.10;   % yellow
    0.20 0.65 0.25;   % green
    0.15 0.40 0.85;   % blue
    0.55 0.30 0.75;   % purple
    0.10 0.65 0.65;   % teal
    0.75 0.30 0.45 ]; % magenta/rose

if isempty(opt.PartColors)
    nP = numel(partNames);
    partColors = defaultPalette(mod((0:nP-1)', size(defaultPalette,1)) + 1, :);
else
    if size(opt.PartColors, 1) ~= numel(partNames)
        error('PoseAnalysisMaster:partColorMismatch', ...
            'PartColors must have one row per entry in Parts (%d rows expected, got %d).', ...
            numel(partNames), size(opt.PartColors, 1));
    end
    partColors = opt.PartColors;
end

% Restrict T (and unstableAll in lock-step) to requested parts (in case
% 'Parts' was given as a subset)
partSel = ismember(T.part, string(partNames));
T           = T(partSel, :);
unstableAll = unstableAll(partSel);
if isempty(T)
    error('PoseAnalysisMaster:noRowsAfterPartFilter', ...
        'No rows remain after filtering to Parts = {%s}.', strjoin(partNames, ', '));
end

fprintf('Parts: %s\n\n', strjoin(partNames, ', '));

% -------------------------------------------------------------------------
% BUILD FLAT ROW ARRAYS (already pose-aligned by the CSV -- no matching
% or parsing needed, unlike the old xlsx + MASTER_summary.txt pipeline)
% -------------------------------------------------------------------------
expAll  = double(T.sim_freq_pct);
algoMat = zeros(height(T), nAlgos);
for a = 1:nAlgos
    v = double(T.(geomCols{a}));
    v(isnan(v)) = 0;
    algoMat(:, a) = v;
end

% ------------------------------------------------------------------------
% POST-HOC MODIFIER (optional) -- applied per (part, alpha_deg, beta_deg)
% condition group, independently to each of the six algorithm columns,
% then renormalized within that group. Since geom_*_pct is already
% normalized per condition, this exactly reproduces what re-running
% ChutePoseAnalysis.m with the modifier baked into the raw weights would
% have produced -- see math note in the header comment.
%
% Flagged-unstable poses ARE included in this renormalization (their raw
% weight still physically exists and affects the condition's total mass),
% they are just excluded downstream from the stats/fit -- the same way
% they're excluded from stats whether or not a modifier is applied.
% ------------------------------------------------------------------------
if opt.ApplyModifier
    fprintf('Applying post-hoc modifier (mod_i = 1/avg(ratio_wall, ratio_floor)) per condition...\n');
    rw = double(T.ratio_wall);
    rf = double(T.ratio_floor);
    modVals = opt.ModifierFcn(rw, rf);
    modVals(isnan(modVals) | isinf(modVals)) = 0;

    [~, ~, condIdxAll] = unique(T(:, {'part','alpha_deg','beta_deg'}), 'rows');
    nCond = max(condIdxAll);
    for a = 1:nAlgos
        colVals = algoMat(:, a);
        for gc = 1:nCond
            sel = (condIdxAll == gc);
            weighted = colVals(sel) .* modVals(sel);
            denom = sum(weighted);
            if denom > 0
                colVals(sel) = 100 * weighted / denom;
            else
                colVals(sel) = 0;
            end
        end
        algoMat(:, a) = colVals;
    end
    fprintf('  done.\n\n');
end

% Per-row part index (into partNames / partColors), for point coloring
partIdxAll = zeros(height(T), 1);
for i = 1:height(T)
    partIdxAll(i) = find(strcmpi(partNames, T.part(i)), 1);
end

fprintf('Total rows collected: %d\n\n', height(T));

% =========================================================================
% SHARED EXCLUSION MASK
%
%   A row is dropped from EVERYTHING (plot and stats alike) when
%   exp == 0 AND at least half (>= 3 of 6) of the algorithms also
%   predict 0 for that row. This is independent of, and applied before,
%   the flagged-unstable exclusion below. Because the mask is computed
%   once from all six algorithm columns and applied identically to every
%   figure, all six scatter plots are built from the same retained set
%   of rows.
% =========================================================================
zeroCnt  = sum(algoMat == 0, 2);
keepMask = ~(expAll == 0 & zeroCnt >= (nAlgos / 2));

fprintf('Rows retained after shared exclusion filter: %d / %d\n', ...
    sum(keepMask), height(T));

% statsMask: rows used for every statistic/fit (plotted AND stable).
% plotMask (== keepMask): rows drawn on the scatter at all, stable or not.
statsMask = keepMask & ~unstableAll;
fprintf('Of those, %d are flagged-unstable (plotted grey, excluded from stats); %d remain for stats.\n\n', ...
    sum(keepMask & unstableAll), sum(statsMask));

% =========================================================================
% SCATTER/FIT: ALGORITHM PROBABILITY vs EXPERIMENTAL PROBABILITY
%
%   One standalone figure per algorithm (six total), saved individually.
%   Pooled across every part, condition, and pose, using the shared
%   keepMask above for what gets PLOTTED, and statsMask for what feeds
%   every statistic (MSE, Pearson r, R^2, and the best-fit line itself).
%   For each algorithm:
%     - small, low-opacity (alpha=0.25) scatter points: colored by part
%       for stable poses, light grey for flagged-unstable poses
%     - a gray dashed y=x line (calibration reference)
%     - a red solid least-squares best-fit line, computed via linear
%       regression (polyfit) pooled across all parts, STABLE POSES ONLY
%     - MSE and R^2 are computed relative to the y=x line (STABLE POSES
%       ONLY), since the question of interest is calibration/accuracy
%       (how close is the algorithm to truth), not merely linear
%       association. Pearson r and the best-fit-line R^2 are also
%       computed and reported (also stable-only), but as secondary
%       numbers describing linear association/bias only (R^2_fit ==
%       pearsonR^2 for a simple OLS fit).
%     - a top-left legend showing the calibration line, scatter note,
%       best-fit equation, and (if any are plotted) the flagged-unstable
%       marker note
%     - a bottom-right legend mapping part -> color, drawn individually
%       on this figure
% =========================================================================
fprintf('Building scatter/fit figures...\n');

yPlot        = expAll(keepMask);
partIdxPlot  = partIdxAll(keepMask);
unstablePlot = unstableAll(keepMask);

colorsPlot = partColors(partIdxPlot, :);
colorsPlot(unstablePlot, :) = repmat(UNSTABLE_COLOR, sum(unstablePlot), 1);

scatterXPlot = cell(nAlgos, 1);   % all plotted points (stable + unstable)
for a = 1:nAlgos
    xAllA           = algoMat(:, a);
    scatterXPlot{a} = xAllA(keepMask);
end

% Stats-only arrays (stable poses only)
yStats       = expAll(statsMask);
partIdxStats = partIdxAll(statsMask);
scatterXStats = cell(nAlgos, 1);
for a = 1:nAlgos
    xAllA            = algoMat(:, a);
    scatterXStats{a} = xAllA(statsMask);
end

fitStats = struct('mse', nan(1,nAlgos), 'pearsonR', nan(1,nAlgos), ...
                   'r2_yx', nan(1,nAlgos), 'r2_fit', nan(1,nAlgos), ...
                   'slope', nan(1,nAlgos), 'intercept', nan(1,nAlgos), ...
                   'n', zeros(1,nAlgos));

% Per-part, per-algorithm stats (for the report only; stable poses only)
partStats = struct('mse', nan(nAlgos, numel(partNames)), ...
                    'r2_yx', nan(nAlgos, numel(partNames)), ...
                    'pearsonR', nan(nAlgos, numel(partNames)), ...
                    'n', zeros(nAlgos, numel(partNames)));

outFilesFit = cell(nAlgos, 1);

for a = 1:nAlgos
    xPlot = scatterXPlot{a};
    xFitData = scatterXStats{a};
    yFitData = yStats;

    hFig = figure('Visible', 'off', 'Color', 'w', 'Units', 'inches', ...
        'Position', [0 0 6.5 6]);
    ax = axes('Position', [0.13, 0.11, 0.83, 0.80]);
    hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');

    axMax = max([xPlot; yPlot], [], 'omitnan') * 1.05;
    if isempty(axMax) || axMax <= 0, axMax = 1; end

    % --- y = x reference line (gray), drawn first ---
    hYX = plot(ax, [0 axMax], [0 axMax], '--', 'Color', [0.6 0.6 0.6], ...
        'LineWidth', 1.4, 'DisplayName', 'Calibration line (y = x)');

    % --- scatter of ALL plotted pose probabilities: part-colored for
    %     stable poses, light grey for flagged-unstable poses ---
    hSc = scatter(ax, xPlot, yPlot, 20, colorsPlot, 'filled', ...
        'MarkerFaceAlpha', 0.25, 'MarkerEdgeColor', 'none', ...
        'DisplayName', 'Pose probabilities (colored by part)');

    hUn = gobjects(0);
    if any(unstablePlot)
        hUn = scatter(ax, xPlot(unstablePlot), yPlot(unstablePlot), 20, UNSTABLE_COLOR, ...
            'filled', 'MarkerFaceAlpha', 0.45, 'MarkerEdgeColor', [0.45 0.45 0.45], ...
            'LineWidth', 0.5, 'DisplayName', 'Flagged unstable (excluded from stats)');
    end

    % --- best-fit line (red, solid), from linear regression on STABLE
    %     (statsMask) data only ---
    if numel(xFitData) >= 2 && std(xFitData) > 0
        pf   = polyfit(xFitData, yFitData, 1);
        xFit = [0 axMax];
        yFit = polyval(pf, xFit);
        hFit = plot(ax, xFit, yFit, '-', 'Color', [0.80 0.10 0.10], ...
            'LineWidth', 2.0, 'DisplayName', ...
            sprintf('Best fit (y = %.3fx + %.3f), stable poses only', pf(1), pf(2)));
        fitStats.slope(a)     = pf(1);
        fitStats.intercept(a) = pf(2);
    else
        pf   = [NaN NaN];
        hFit = gobjects(0);
    end

    % --- overall statistics (relative to y=x and to best-fit line),
    %     STABLE POSES ONLY ---
    resid_yx = yFitData - xFitData;
    mse_val  = mean(resid_yx.^2, 'omitnan');
    ssRes_yx = sum(resid_yx.^2, 'omitnan');
    ssTot    = sum((yFitData - mean(yFitData, 'omitnan')).^2, 'omitnan');
    if ssTot > 0
        r2_yx = 1 - ssRes_yx / ssTot;
    else
        r2_yx = nan;
    end

    if numel(xFitData) >= 2 && std(xFitData) > 0 && std(yFitData) > 0
        rMat     = corrcoef(xFitData, yFitData);
        pearsonR = rMat(1,2);
    else
        pearsonR = nan;
    end

    if all(isfinite(pf))
        yHatFit   = polyval(pf, xFitData);
        ssRes_fit = sum((yFitData - yHatFit).^2, 'omitnan');
        if ssTot > 0
            r2_fit = 1 - ssRes_fit / ssTot;
        else
            r2_fit = nan;
        end
    else
        r2_fit = nan;
    end

    fitStats.mse(a)      = mse_val;
    fitStats.pearsonR(a) = pearsonR;
    fitStats.r2_yx(a)    = r2_yx;
    fitStats.r2_fit(a)   = r2_fit;
    fitStats.n(a)        = numel(xFitData);

    % --- per-part statistics (relative to y=x only), STABLE POSES ONLY,
    %     for the report ---
    for pIdx = 1:numel(partNames)
        selP = (partIdxStats == pIdx);
        if sum(selP) < 2, continue; end
        xP = xFitData(selP); yP = yFitData(selP);
        residP = yP - xP;
        ssResP = sum(residP.^2, 'omitnan');
        ssTotP = sum((yP - mean(yP, 'omitnan')).^2, 'omitnan');
        partStats.mse(a, pIdx)   = mean(residP.^2, 'omitnan');
        partStats.n(a, pIdx)     = numel(xP);
        if ssTotP > 0
            partStats.r2_yx(a, pIdx) = 1 - ssResP / ssTotP;
        end
        if std(xP) > 0 && std(yP) > 0
            rMatP = corrcoef(xP, yP);
            partStats.pearsonR(a, pIdx) = rMatP(1,2);
        end
    end

    xlim(ax, [0, axMax]);
    ylim(ax, [0, axMax]);
    axis(ax, 'square');
    xlabel(ax, sprintf('%s probability (%%)', strrep(algoNames{a}, '_', '\_')), 'FontSize', 10);
    ylabel(ax, 'Experimental probability (%)', 'FontSize', 10);
    title(ax, sprintf('%s:  MSE=%.3f | r=%.3f | R^2_{y=x}=%.3f  (stable poses only, n=%d)', ...
        strrep(algoNames{a}, '_', '\_'), mse_val, pearsonR, r2_yx, numel(xFitData)), ...
        'FontSize', 10.5, 'FontWeight', 'bold');

    % Data legend, top-left: calibration line, scatter note, unstable
    % note (if present), best-fit equation
    legHandles = [hYX, hSc];
    if ~isempty(hUn), legHandles = [legHandles, hUn]; end %#ok<AGROW>
    if ~isempty(hFit), legHandles = [legHandles, hFit]; end %#ok<AGROW>
    legend(ax, legHandles, 'Location', 'northwest', 'FontSize', 7, 'Box', 'on');

    % --- part-color legend, bottom-right corner of the plot itself,
    %     drawn individually on this figure (no shared/top legend across
    %     figures).
    axPos = ax.Position;
    legCornerX = axPos(1) + axPos(3)-0.05;
    legCornerY = axPos(2)+0.15;
    legAxPart = axes('Position', [legCornerX, legCornerY, 0.001, 0.001], 'Visible', 'off');
    hold(legAxPart, 'on');
    dummyH = gobjects(numel(partNames), 1);
    for pIdx = 1:numel(partNames)
        dummyH(pIdx) = scatter(legAxPart, NaN, NaN, 40, partColors(pIdx, :), ...
            'filled', 'DisplayName', partNames{pIdx});
    end
    legend(legAxPart, dummyH, partNames, 'Location', 'southeast', ...
        'Box', 'on', 'FontSize', 7);

    modTag = '';
    if opt.ApplyModifier, modTag = '_mod'; end
    outFileFit = fullfile(resultsRoot, ...
        sprintf('PoseAnalysisMaster_ScatterFit_%s%s.png', algoNames{a}, modTag));
    exportgraphics(hFig, outFileFit, 'Resolution', 150);
    delete(hFig);
    outFilesFit{a} = outFileFit;
    fprintf('  Saved: %s\n', outFileFit);
end

% =========================================================================
% WRITE OUTPUT TEXT FILE -- scatter/fit stats, overall + per-part
% =========================================================================
reportFile = fullfile(resultsRoot, 'PoseAnalysisMaster_Report.txt');
fid = fopen(reportFile, 'w');
if fid < 0
    error('Cannot open report file: %s', reportFile);
end

writeLine(fid, repmat('=', 1, 92));
writeLine(fid, '  SCATTER/FIT: ALGORITHM vs EXPERIMENTAL PROBABILITY');
writeLine(fid, '  Six CSA/CRSA variants (Dyn/Stat/Nan), pooled across all parts,');
writeLine(fid, '  conditions, and poses. Each variant is saved as its own figure.');
fprintf(fid, '  Source: %s\n', opt.summaryPath);
if opt.ApplyModifier
    writeLine(fid, '  Post-hoc modifier APPLIED: mod_i = 1/avg(ratio_wall, ratio_floor),');
    writeLine(fid, '  renormalized within each (part, alpha_deg, beta_deg) condition.');
end
writeLine(fid, '');
writeLine(fid, '  GLOSSARY');
writeLine(fid, '  CSA_Dyn / CRSA_Dyn   : transitioning poses zeroed, remaining renormalized ("before")');
writeLine(fid, '  CSA_Stat / CRSA_Stat : transitioning poses redistribute to destination pose ("after")');
writeLine(fid, '  CSA_Nan / CRSA_Nan   : no filtering at all, raw weight normalized ("none")');
writeLine(fid, '');
writeLine(fid, '  Unmatched poses (pose_id == -1) were excluded before any of the');
writeLine(fid, '  statistics below -- there is no algorithmic value to scatter them');
writeLine(fid, '  against.');
writeLine(fid, '');
writeLine(fid, '  Poses flagged unstable by geometry (flagged_unstable_by_geometry ==');
writeLine(fid, '  true) are EXCLUDED from every statistic and from the best-fit line');
writeLine(fid, '  below, but ARE still shown on each scatter figure, in light grey,');
writeLine(fid, '  so their landing location remains visible.');
fprintf(fid, '  Flagged-unstable rows (post shared-exclusion-filter): %d\n', sum(keepMask & unstableAll));
writeLine(fid, '');
writeLine(fid, '  Shared exclusion filter (applied identically to all six');
writeLine(fid, '  algorithms/figures, independent of the flagged-unstable exclusion');
writeLine(fid, '  above): a row is dropped entirely (plot and stats) when experimental');
writeLine(fid, '  probability == 0 AND at least half of the six algorithms also');
writeLine(fid, '  predict 0 for that row.');
fprintf(fid, '  Rows retained after shared exclusion filter: %d / %d\n', sum(keepMask), height(T));
fprintf(fid, '  Rows used for statistics (stable + retained): %d / %d\n', sum(statsMask), height(T));
writeLine(fid, '');
writeLine(fid, '  MSE and R^2_{y=x} are computed relative to the y=x (perfect');
writeLine(fid, '  estimation / calibration) line: they measure calibration/accuracy');
writeLine(fid, '  against truth, not just linear association.');
writeLine(fid, '  Pearson r and R^2_{fit} are computed relative to the least-squares');
writeLine(fid, '  best-fit line (linear regression via polyfit) and measure linear');
writeLine(fid, '  association/bias only (R^2_fit == PearsonR^2 for a simple OLS fit).');
writeLine(fid, repmat('=', 1, 92));
writeLine(fid, '');
writeLine(fid, '  OVERALL (pooled across all parts, stable poses only)');
fprintf(fid, '  %-10s  %6s  %10s  %10s  %12s  %12s  %8s  %10s\n', ...
    'Algo', 'N', 'MSE', 'PearsonR', 'R2 (y=x)', 'R2 (fit)', 'Slope', 'Intercept');
fprintf(fid, '  %s\n', repmat('-', 1, 84));
for a = 1:nAlgos
    fprintf(fid, '  %-10s  %6d  %10.4f  %10.4f  %12.4f  %12.4f  %8.4f  %10.4f\n', ...
        algoNames{a}, fitStats.n(a), fitStats.mse(a), fitStats.pearsonR(a), ...
        fitStats.r2_yx(a), fitStats.r2_fit(a), fitStats.slope(a), fitStats.intercept(a));
end
writeLine(fid, '');
writeLine(fid, repmat('-', 1, 92));
writeLine(fid, '');
writeLine(fid, '  PER-PART BREAKDOWN (relative to y=x only, stable poses only; parts');
writeLine(fid, '  with < 2 retained rows for a given algorithm are omitted and shown');
writeLine(fid, '  as N/A)');
writeLine(fid, '');
for a = 1:nAlgos
    fprintf(fid, '  %s\n', algoNames{a});
    fprintf(fid, '    %-8s  %6s  %10s  %10s  %12s\n', 'Part', 'N', 'MSE', 'PearsonR', 'R2 (y=x)');
    fprintf(fid, '    %s\n', repmat('-', 1, 52));
    for pIdx = 1:numel(partNames)
        nP = partStats.n(a, pIdx);
        if nP < 2
            fprintf(fid, '    %-8s  %6d  %10s  %10s  %12s\n', partNames{pIdx}, nP, 'N/A', 'N/A', 'N/A');
        else
            fprintf(fid, '    %-8s  %6d  %10.4f  %10.4f  %12.4f\n', partNames{pIdx}, nP, ...
                partStats.mse(a, pIdx), partStats.pearsonR(a, pIdx), partStats.r2_yx(a, pIdx));
        end
    end
    writeLine(fid, '');
end
writeLine(fid, repmat('=', 1, 92));
fclose(fid);
fprintf('Saved report: %s\n\n', reportFile);

disp('Done.');
end % PoseAnalysisMaster


% =========================================================================
%  writeLine -- write a string + newline to fid
% =========================================================================
function writeLine(fid, str)
fprintf(fid, '%s\n', str);
end