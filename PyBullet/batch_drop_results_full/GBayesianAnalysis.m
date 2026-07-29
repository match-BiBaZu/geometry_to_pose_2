% =========================================================================
% GBayesianAnalysis.m  (CSV edition -- combined_summary.csv format)
%
% Statistical goodness-of-fit and model-comparison analysis for the six
% CSA/CRSA pose-prediction algorithm variants, reading the flat
% combined_summary.csv (or a single *_sim_summary.csv) format written by
% chute_drop_batch.py's write_condition_summary() -- the same file format
% consumed by plot_pose_probabilities.m and PoseAnalysisMaster.m.
%
%   CSA_Dyn   = CSA,  transitioning poses zeroed, remaining renormalized  ("before")  <- geom_CSA_B_pct
%   CSA_Stat  = CSA,  transitioning poses redistribute to destination pose ("after")  <- geom_CSA_A_pct
%   CSA_Nan   = CSA,  no filtering at all (raw CSA weight normalized)     ("none")   <- geom_CSA_N_pct
%   CRSA_Dyn  = CRSA, transitioning poses zeroed, remaining renormalized  ("before")  <- geom_CRSA_B_pct
%   CRSA_Stat = CRSA, transitioning poses redistribute to destination pose ("after")  <- geom_CRSA_A_pct
%   CRSA_Nan  = CRSA, no filtering at all (raw CRSA weight normalized)     ("none")   <- geom_CRSA_N_pct
%
% DATA INGESTION -- MUCH SIMPLER THAN THE OLD XLSX/TXT PIPELINE:
% the CSV already stores one row per (part, alpha_deg, beta_deg, pose_id)
% with the experimental count (sim_count) and all six algorithmic
% predictions (geom_*_pct) on the SAME row. There is no separate
% experimental-workbook parsing, no *_MASTER_summary.txt text parsing,
% and no pose-number matching step -- chute_drop_batch.py already did
% that matching when it wrote the file. Each condition is just
% "group rows by (part, alpha_deg, beta_deg)".
%
% alpha_deg / beta_deg are treated as roll / pitch respectively, matching
% the convention already used in plot_pose_probabilities.m and
% PoseAnalysisMaster.m (title strings there read
% "roll=%g, pitch=%g" from alpha_deg, beta_deg).
%
% UNMATCHED POSES (pose_id == -1, i.e. drop trials whose settled
% orientation matched no catalog pose within quat_match_tol): by default
% these ARE included as their own explicit category per condition, with
% every algorithm's predicted probability floored to the same epsilon
% used for any zero-probability cell (no algorithm predicts "unmatched").
% This keeps each condition's multinomial self-consistent -- counts sum
% to exactly N, the total trials actually run in that condition -- and
% treats "the algorithm didn't anticipate a trial landing nowhere in the
% catalog" as a genuine (if usually small) miss for every algorithm
% alike, rather than silently discarding those trials. Set
% 'IncludeUnmatched' to false to drop those trials entirely instead
% (N shrinks to just the matched trials for that condition).
%
% FLAGGED-UNSTABLE POSES (geom_notes == "TRANSITIONS" /
% flagged_unstable_by_geometry == TRUE in the CSV): these are catalog
% poses that ChutePoseAnalysis.m has flagged, ONCE, at cataloging time,
% as geometrically transient equilibria -- a FIXED, per-pose-id property
% that does not vary by condition (roll/pitch) or by algorithm. This is
% distinct from, and upstream of, the per-condition zeroing that CSA_B/
% CRSA_B ("before") and CSA_A/CRSA_A ("after") already apply -- that
% zeroing is a per-condition, per-algorithm decision about a specific
% trial's transition state, whereas flagged_unstable_by_geometry is a
% fixed catalog-level label. Because the flag is identical across all
% six algorithms, it is safe to use it to redefine the category universe
% for ALL of them at once without giving any algorithm an unfair,
% self-selected pass -- unlike post-hoc removal of zero-cell events
% (which would be circular, since it's the model's own residual driving
% the removal). Set 'ExcludeFlaggedUnstable' to true to drop these rows
% from the category set for every algorithm, renormalize each of the six
% geom_*_pct columns over the remaining poses per condition, and shrink N
% for that condition to the trial count over the remaining (non-flagged)
% poses only -- the same "shrink N to the surviving universe" pattern
% already used by 'IncludeUnmatched' = false. Default: false.
%
% using:
%
%   (1) G-TEST (log-likelihood-ratio test) -- per algorithm, ABSOLUTE
%       goodness of fit: "is this algorithm's predicted distribution
%       statistically consistent with what was actually observed?"
%       G is computed per (part,alpha_deg,beta_deg) condition from raw
%       trial counts, then summed across all conditions (G is additive
%       for independent multinomial draws).
%
%       Three quantities are reported for this test:
%         - G/df: the G statistic normalized by total degrees of
%           freedom. This is reported ABOVE the p-value in the summary
%           graphic as a quick-glance severity index (G/df near 1
%           roughly corresponds to what pure sampling noise would
%           produce; much larger than 1 signals systematic misfit).
%         - p_asymptotic: the standard chi^2(df) approximation. This is
%           UNRELIABLE here because many pose categories have very low
%           expected counts (sparse multinomial cells inflate df without
%           contributing proportional signal to G), which is why it can
%           saturate near 1 even for a well-fitting model.
%         - p_bootstrap: a parametric bootstrap p-value. For each
%           algorithm, we simulate many synthetic datasets assuming that
%           algorithm's predicted distribution IS the true generating
%           distribution, compute G on each synthetic dataset, and see
%           where the real G falls in that simulated null distribution.
%           This does not rely on the chi^2 asymptotic approximation and
%           is the trustworthy number to report/present.
%
%   (2) BAYESIAN POSTERIOR -- RELATIVE model comparison: "given all the
%       data, how should belief be split across the six CSA/CRSA
%       variants?" Reuses the same per-condition multinomial
%       log-likelihoods computed for the G-test (logL_a = logL_sat -
%       G_a/2), pooled with a flat prior and converted to posterior
%       weights via softmax.
%       NOTE: this is a RELATIVE statement among these six candidates
%       only. It says nothing about whether the "winning" algorithm is
%       absolutely accurate -- see (3) for that.
%
%   (3) CALIBRATION (reliability diagram) -- ABSOLUTE accuracy check per
%       algorithm, independent of the other five: for every (condition,
%       pose), bin by the algorithm's predicted probability, and check
%       whether the observed frequency in that bin actually matches. A
%       well-calibrated algorithm's points fall on the y=x diagonal.
%       This is the tool that answers "does this model predict real
%       probabilities with X accuracy" rather than "which model is best
%       relative to its competitors."
%
%       Bins are EQUAL-WEIGHT (quantile) bins, not fixed-width. With a
%       modest number of conditions/trials backing this analysis,
%       fixed-width bins across the predicted-probability range leave
%       the high-probability tail almost empty while the near-zero end
%       is oversaturated; equal-weight bins instead cut the sorted,
%       trial-weighted data so every bin carries a comparable amount of
%       evidence. Bin count is chosen adaptively (see pickNBins) from
%       the actual amount of data available, separately for each of the
%       two variants below, and edges are shared across all six
%       algorithms within a variant so the panels stay comparable.
%
%       TWO variants are computed:
%         - ALL POSES: every (condition, pose) pair, regardless of
%           whether that pose ever occurred. This is dominated, by raw
%           trial-weight, by the large number of poses that correctly
%           never occur (predicted ~0%, observed ~0%), which inflates
%           the apparent calibration quality without testing whether the
%           algorithm is any good at sizing the probability of poses
%           that actually happen.
%         - OBSERVED-ONLY: restricted to (condition, pose) pairs where
%           that pose was experimentally observed at least once in that
%           specific condition (sim_count > 0). This is the harder, more
%           meaningful test: "when a pose is real, does this algorithm's
%           predicted probability for it actually track how often it
%           happened?" No tunable threshold is involved -- the filter is
%           simply "did this pose occur here" -- so it carries none of
%           the researcher-degrees-of-freedom risk of an arbitrary
%           probability cutoff. It also has far less underlying data
%           than the all-poses variant, hence fewer bins.
%
% All quantities are computed on RAW TRIAL COUNTS (sim_count, not
% normalized percentages), so conditions with more physical trials
% correctly carry more statistical weight than sparse ones.
%
% COLOR CONVENTION (summary graphic): the three CSA variants are shown
% in shades of BLUE (light -> dark = Dyn -> Stat -> Nan) and the three
% CRSA variants are shown in shades of RED (light -> dark = Dyn -> Stat
% -> Nan), so the two algorithm families are visually distinguishable
% at a glance and the within-family ordering is consistent everywhere.
%
% OUTPUT FILES (saved to OutputDir):
%   GBayesianAnalysis_Report.txt              -- G/df, p-value/posterior/calibration tables
%   GBayesianAnalysis_Calibration_<Algo>.png  -- one calibration image per algorithm (6 files)
%   GBayesianAnalysis_AlgorithmSupport.png    -- algorithm support across conditions
%   GBayesianAnalysis_BayesianPosterior.png   -- headline relative-comparison chart
%   GBayesianAnalysis_GTest.png               -- absolute goodness-of-fit chart
%
% --------------------------------------------------------------------
% USAGE
% --------------------------------------------------------------------
%   GBayesianAnalysis('combined_summary.csv')
%   GBayesianAnalysis('combined_summary.csv', 'OutputDir', 'analysis_out')
%   GBayesianAnalysis('combined_summary.csv', 'IncludeUnmatched', false)
%   GBayesianAnalysis('combined_summary.csv', 'ApplyModifier', true)
%   GBayesianAnalysis('combined_summary.csv', 'ExcludeFlaggedUnstable', true)
%   GBayesianAnalysis('combined_summary.csv', 'Parts', {'Df4a','Qf4i'})
%
% --------------------------------------------------------------------
% NAME-VALUE OPTIONS
% --------------------------------------------------------------------
%   'OutputDir'       Folder to save the report + figures into. Created
%                     if it doesn't exist. Default: a
%                     'GBayesianAnalysis_out' subfolder next to
%                     summaryPath.
%   'Parts'           Cell array of part names to include. Default: all
%                     unique values in the 'part' column.
%   'IncludeUnmatched' true/false. If true (default), pose_id == -1 rows
%                     are kept as their own category per condition, with
%                     every algorithm's predicted probability floored to
%                     epsilon (see header note above). If false, those
%                     trials are dropped entirely and N shrinks
%                     accordingly for that condition.
%   'ExcludeFlaggedUnstable' true/false. If true, catalog poses with
%                     flagged_unstable_by_geometry == TRUE (equivalently
%                     geom_notes == "TRANSITIONS") are dropped from the
%                     category set for ALL SIX algorithms, the six
%                     geom_*_pct columns are renormalized over the
%                     remaining poses within each (part, alpha_deg,
%                     beta_deg) condition, and N for that condition
%                     shrinks to the trial count over the remaining
%                     (non-flagged) poses only. See header note above for
%                     why this is a legitimate, non-circular restriction
%                     (fixed catalog-level label, identical across
%                     algorithms) as opposed to post-hoc zero-cell
%                     removal. Requires a 'flagged_unstable_by_geometry'
%                     column in the CSV. Default: false.
%   'ApplyModifier'   true/false. If true, applies a per-pose
%                     multiplicative modifier to EACH of the six
%                     geom_*_pct columns, independently, within each
%                     (part, alpha_deg, beta_deg) condition, and
%                     renormalizes -- the same post-hoc renormalization
%                     used in plot_pose_probabilities.m /
%                     PoseAnalysisMaster.m, since geom_*_pct is already
%                     normalized per condition and the normalization
%                     constant cancels algebraically. This changes what
%                     ALL THREE analyses (G-test, Bayesian posterior,
%                     calibration) see as each algorithm's prediction,
%                     since all three are computed downstream of the
%                     (possibly modified) predicted-probability vector.
%                     Default modifier: mod_i = 1/mean(ratio_wall_i, ratio_floor_i).
%                     Default: false.
%   'ModifierFcn'     Optional function handle @(ratioWall, ratioFloor) -> mod_i.
%                     Ignored if ApplyModifier is false.
%                     Default: @(rw, rf) 1 ./ mean([rw, rf], 2, 'omitnan')
% =========================================================================

function GBayesianAnalysis(summaryPath, varargin)

close all; clc;

% -------------------------------------------------------------------------
% OPTIONS
% -------------------------------------------------------------------------
p = inputParser;
addRequired(p, 'summaryPath', @(x) ischar(x) || isstring(x));
addParameter(p, 'OutputDir', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'Parts', {}, @iscell);
addParameter(p, 'IncludeUnmatched', true, @islogical);
addParameter(p, 'ExcludeFlaggedUnstable', false, @islogical);
addParameter(p, 'ApplyModifier', false, @islogical);
addParameter(p, 'ModifierFcn', @(rw, rf) 1 ./ mean([rw, rf], 2, 'omitnan'), ...
    @(x) isa(x, 'function_handle'));
addParameter(p, 'ModifierMaxBoost', 20, @(x) isnumeric(x) && isscalar(x) && x > 0);
parse(p, summaryPath, varargin{:});
opt = p.Results;

if ~isfile(opt.summaryPath)
    error('GBayesianAnalysis:fileNotFound', 'Could not find file: %s', opt.summaryPath);
end

if isempty(opt.OutputDir)
    [srcDir, ~, ~] = fileparts(opt.summaryPath);
    if isempty(srcDir), srcDir = pwd; end
    resultsRoot = fullfile(srcDir, 'GBayesianAnalysis_out');
else
    resultsRoot = char(opt.OutputDir);
end
if ~exist(resultsRoot, 'dir'), mkdir(resultsRoot); end

% Six CSA/CRSA variants. Order matters: it drives both panel placement
% (CSA row, then CRSA row) and the color assignment below.
algoNames = {'CSA_Dyn','CSA_Stat','CSA_Nan','CRSA_Dyn','CRSA_Stat','CRSA_Nan'};
geomCols  = {'geom_CSA_B_pct','geom_CSA_A_pct','geom_CSA_N_pct', ...
             'geom_CRSA_B_pct','geom_CRSA_A_pct','geom_CRSA_N_pct'};
nAlgos    = 6;

% Flat prior over the six algorithms (no reason to favor any a priori)
priorVec = ones(1, nAlgos) / nAlgos;

% Floor probability substituted when an algorithm predicts EXACTLY 0
% for a pose that was actually observed (n_i > 0), OR for the unmatched
% category when IncludeUnmatched is true. Prevents G = Inf and flags
% the event for reporting.
epsilonFrac = 1e-6;

% Bootstrap settings
nBootstrap = 2000;
randSeed   = 42;

% -------------------------------------------------------------------------
% LOAD DATA
% -------------------------------------------------------------------------
fprintf('Loading %s...\n', opt.summaryPath);
T = readtable(opt.summaryPath, 'TextType', 'string');

requiredCols = {'part','alpha_deg','beta_deg','pose_id','sim_count','sim_freq_pct'};
for i = 1:numel(requiredCols)
    if ~ismember(requiredCols{i}, T.Properties.VariableNames)
        error('GBayesianAnalysis:missingColumn', ...
            'Expected column "%s" not found in %s. Is this a combined_summary.csv / *_sim_summary.csv written by chute_drop_batch.py?', ...
            requiredCols{i}, opt.summaryPath);
    end
end
for i = 1:nAlgos
    if ~ismember(geomCols{i}, T.Properties.VariableNames)
        error('GBayesianAnalysis:missingGeomColumn', ...
            'Expected column "%s" not found in %s.', geomCols{i}, opt.summaryPath);
    end
end
if opt.ApplyModifier
    for c = {'ratio_wall','ratio_floor'}
        if ~ismember(c{1}, T.Properties.VariableNames)
            error('GBayesianAnalysis:missingModifierColumn', ...
                'ApplyModifier=true requires column "%s", not found in %s.', c{1}, opt.summaryPath);
        end
    end
end
if opt.ExcludeFlaggedUnstable
    if ~ismember('flagged_unstable_by_geometry', T.Properties.VariableNames)
        error('GBayesianAnalysis:missingUnstableColumn', ...
            'ExcludeFlaggedUnstable=true requires column "flagged_unstable_by_geometry", not found in %s.', ...
            opt.summaryPath);
    end
end

if isempty(opt.Parts)
    partNames = cellstr(unique(T.part, 'sorted'));
else
    partNames = cellfun(@char, opt.Parts, 'UniformOutput', false);
end
T = T(ismember(T.part, string(partNames)), :);
if isempty(T)
    error('GBayesianAnalysis:noRowsAfterPartFilter', ...
        'No rows remain after filtering to Parts = {%s}.', strjoin(partNames, ', '));
end
fprintf('Parts: %s\n', strjoin(partNames, ', '));

if ~opt.IncludeUnmatched
    nBeforeUnm = height(T);
    T = T(T.pose_id ~= -1, :);
    fprintf('IncludeUnmatched=false: dropped %d unmatched row(s); %d row(s) remain.\n', ...
        nBeforeUnm - height(T), height(T));
end

% ------------------------------------------------------------------------
% EXCLUDE FLAGGED-UNSTABLE POSES (optional) -- a FIXED, per-pose-id
% catalog label (identical across all conditions and all six algorithms;
% see header note), NOT a per-condition/per-algorithm transition
% classification. Because the label is the same for every algorithm,
% dropping these rows for all six at once and renormalizing does not
% give any algorithm a self-selected pass -- every algorithm is scored
% against the same reduced category universe. This intentionally mirrors
% the IncludeUnmatched=false pattern above: rows are removed from T
% before any per-condition grouping happens below, so N naturally
% shrinks to the surviving (non-flagged) trial count per condition.
%
% NOTE on interpretation: unlike unmatched trials (which by construction
% have sim_count folded into a synthetic pose_id==-1 row and zero
% predicted probability from every algorithm), a flagged-unstable pose
% can still have sim_count > 0 -- i.e. a trial's settled orientation DID
% match a catalog pose that geometry flags as transient. Excluding those
% rows also discards those observed hits from N. If a flagged pose has
% a non-trivial hit count, it's worth checking those specific drop
% trajectories (loosely-set quat_match_tol catching a mid-transition
% pause, vs. a genuine rare metastable state) before leaning on this
% option -- this flag changes the SCORE, not the underlying question of
% whether the flag itself is set correctly upstream.
% ------------------------------------------------------------------------
if opt.ExcludeFlaggedUnstable
    isFlagged = false(height(T), 1);
    flagCol = T.flagged_unstable_by_geometry;
    if isnumeric(flagCol) || islogical(flagCol)
        isFlagged = logical(flagCol);
    else
        isFlagged = strcmpi(strtrim(string(flagCol)), 'TRUE') | (string(flagCol) == "1");
    end

    nFlaggedRows   = sum(isFlagged);
    nFlaggedHits   = sum(T.sim_count(isFlagged));
    if nFlaggedHits > 0
        fprintf(['ExcludeFlaggedUnstable=true: %d flagged-unstable row(s) removed, ', ...
            'of which %d row(s) had sim_count > 0 (total %d observed trial(s) discarded from N -- ', ...
            'worth auditing those specific hits before trusting this fully).\n'], ...
            nFlaggedRows, sum(isFlagged & T.sim_count > 0), nFlaggedHits);
    else
        fprintf('ExcludeFlaggedUnstable=true: %d flagged-unstable row(s) removed (all sim_count == 0).\n', ...
            nFlaggedRows);
    end

    T = T(~isFlagged, :);
    fprintf('  %d row(s) remain.\n', height(T));

    % Renormalize each of the six geom_*_pct columns to sum to 100 over
    % the surviving (non-flagged) poses within each (part, alpha_deg,
    % beta_deg) condition. Conditions where an algorithm's remaining
    % mass sums to 0 are left at 0 (all remaining predicted probability
    % was concentrated on now-excluded poses); such a case would itself
    % be worth flagging separately, but is left as-is here rather than
    % silently redistributed.
    [~, ~, condIdxRenorm] = unique(T(:, {'part','alpha_deg','beta_deg'}), 'rows');
    nCondRenorm = max(condIdxRenorm);
    for a = 1:nAlgos
        colVals = double(T.(geomCols{a}));
        for gc = 1:nCondRenorm
            sel = (condIdxRenorm == gc);
            s = sum(colVals(sel), 'omitnan');
            if s > 0
                colVals(sel) = colVals(sel) / s * 100;
            end
        end
        T.(geomCols{a}) = colVals;
    end
    fprintf('  Renormalized all six geom_*_pct columns over remaining poses per condition.\n');
end
fprintf('\n');

% Coerce geom_*_pct to numeric fractions (0-1), treating NaN (the
% "" that chute_drop_batch.py writes for the synthetic unmatched row) as
% 0 -- no algorithm predicts "unmatched".
algoFracAll = zeros(height(T), nAlgos);
for a = 1:nAlgos
    v = double(T.(geomCols{a}));
    v(isnan(v)) = 0;
    algoFracAll(:, a) = v / 100;
end

% ------------------------------------------------------------------------
% POST-HOC MODIFIER (optional) -- applied per (part, alpha_deg, beta_deg)
% condition group, independently to each of the six algorithm columns,
% then renormalized within that group. Since geom_*_pct is already
% normalized per condition, this exactly reproduces what re-running
% ChutePoseAnalysis.m with the modifier baked into the raw weights would
% have produced -- see math note in the header comment. Applied here,
% before any grouping below, so the G-test / Bayesian / calibration
% analyses all see the modified predictions consistently.
%
% BUG FIX vs. earlier revision: mod_i = 1/avg(ratio_wall, ratio_floor) is
% mathematically undefined (Inf) whenever a pose has NO wall contact at
% all (very common -- plenty of legitimate resting poses never touch a
% wall) or a negligible floor moment arm, since computeTransitionRatios
% leaves those ratios at exactly 0. Mapping that Inf straight to
% modVal = 0 (as an earlier version of this script did) SILENTLY KILLS
% the probability of exactly the most stable, most physically common
% poses -- the opposite of what "boost poses with a small ratio" should
% mean. Undefined ratios (both NaN, e.g. truly no data) now leave the
% pose's weight UNCHANGED (mod = 1); a divide-by-zero Inf is capped at
% ModifierMaxBoost instead of collapsing to 0; and any other unusually
% large finite modVal is capped the same way, so one outlier ratio can't
% dominate a condition's renormalization.
%
% origAlgoFracAll (the pre-modifier predictions) and modValsUsed are
% retained so the zero-cell diagnostic below can report WHICH mechanism
% produced each "predicted 0, observed >0" event: the synthetic
% unmatched category, a zero already baked into the CSV by
% ChutePoseAnalysis.m's own transition/floor-instability zeroing, or this
% modifier step.
% ------------------------------------------------------------------------
origAlgoFracAll = algoFracAll;   % snapshot before any modifier, for diagnostics
[~, ~, condIdxPre] = unique(T(:, {'part','alpha_deg','beta_deg'}), 'rows');

modValsUsed = ones(height(T), 1);   % 1 = "no adjustment"; overwritten below if ApplyModifier
nClampedHigh = 0;

if opt.ApplyModifier
    fprintf('Applying post-hoc modifier (mod_i = 1/avg(ratio_wall, ratio_floor)) per condition...\n');
    rw = double(T.ratio_wall);
    rf = double(T.ratio_floor);
    modVals = opt.ModifierFcn(rw, rf);

    modVals(isnan(modVals)) = 1;                 % no ratio data at all -> leave unchanged
    isHigh = isinf(modVals) | (modVals > opt.ModifierMaxBoost);
    nClampedHigh = sum(isHigh);
    modVals(isHigh) = opt.ModifierMaxBoost;       % cap instead of exploding OR collapsing to 0
    modVals(modVals < 0) = 0;                     % guard against a pathological custom ModifierFcn

    if nClampedHigh > 0
        fprintf('  Note: %d row(s) had an undefined/huge ratio and were capped at ModifierMaxBoost=%g ', ...
            nClampedHigh, opt.ModifierMaxBoost);
        fprintf('(previously these were incorrectly zeroed -- see header note).\n');
    end

    nCondPre = max(condIdxPre);
    for a = 1:nAlgos
        colVals = algoFracAll(:, a);
        for gc = 1:nCondPre
            sel = (condIdxPre == gc);
            weighted = colVals(sel) .* modVals(sel);
            denom = sum(weighted);
            if denom > 0
                colVals(sel) = weighted / denom;
            else
                colVals(sel) = 0;
            end
        end
        algoFracAll(:, a) = colVals;
    end
    modValsUsed = modVals;
    fprintf('  done.\n\n');
end

% -------------------------------------------------------------------------
% PITCH/ROLL GRID (auto-detected from data; alpha_deg = roll, beta_deg =
% pitch, matching the convention already used in plot_pose_probabilities.m
% and PoseAnalysisMaster.m)
% -------------------------------------------------------------------------
rollVals  = unique(T.alpha_deg)';   % roll levels
pitchVals = unique(T.beta_deg)';    % pitch levels
nR = numel(rollVals);
nP = numel(pitchVals);

% =========================================================================
% PASS 1 -- walk every (part, alpha_deg, beta_deg) condition present in
% the table, compute per-condition G and log-likelihood for each
% algorithm directly from sim_count. Also retain the raw counts and
% probability matrices themselves (needed for the bootstrap and
% calibration steps below). Unlike the old xlsx/txt pipeline, no
% cross-source pose matching is needed here: each row already carries
% both the experimental count and the six algorithmic predictions for
% the same pose.
% =========================================================================
fprintf('Computing per-condition statistics...\n');

[condKeys, ~, condIdxAll] = unique(T(:, {'part','alpha_deg','beta_deg'}), 'rows');
nCondsTotal = height(condKeys);

condPart     = cell(nCondsTotal, 1);
condRoll     = zeros(nCondsTotal, 1);
condPitch    = zeros(nCondsTotal, 1);
condN        = zeros(nCondsTotal, 1);
condG        = nan(nCondsTotal, nAlgos);
condDf       = nan(nCondsTotal, 1);
condLogL     = nan(nCondsTotal, nAlgos);
condAlgoFrac = cell(nCondsTotal, 1);
condCounts   = cell(nCondsTotal, 1);
condUsed     = false(nCondsTotal, 1);

logLPitchRoll = zeros(nP, nR, nAlgos);
havePitchRoll = false(nP, nR);

zeroCellCount = zeros(1, nAlgos);  % diagnostic: algorithm predicted 0,
                                    % experiment observed > 0

for c = 1:nCondsTotal
    sel = (condIdxAll == c);
    rows = T(sel, :);

    totalTrials = sum(rows.sim_count);
    if totalTrials <= 0
        continue;
    end

    nPoses = height(rows);
    if nPoses < 2
        continue;   % need at least 2 categories for a meaningful df
    end

    nCounts  = double(rows.sim_count);
    algoFrac = algoFracAll(sel, :);   % nPoses x nAlgos, row-aligned to rows/nCounts

    % ---- multinomial log-likelihoods ----
    logL_sat = 0;
    logL_a   = zeros(1, nAlgos);
    for i = 1:nPoses
        ni = nCounts(i);
        if ni <= 0, continue; end   % contributes 0 to every logL
        logL_sat = logL_sat + ni * log(ni / totalTrials);
        for a = 1:nAlgos
            pia = algoFrac(i, a);
            if pia <= 0
                pia = epsilonFrac;
                zeroCellCount(a) = zeroCellCount(a) + 1;
            end
            logL_a(a) = logL_a(a) + ni * log(pia);
        end
    end

    G_a  = 2 * (logL_sat - logL_a);
    df_c = nPoses - 1;

    partName = char(condKeys.part(c));
    rollV    = condKeys.alpha_deg(c);
    pitchV   = condKeys.beta_deg(c);

    condPart{c}     = partName;
    condRoll(c)     = rollV;
    condPitch(c)    = pitchV;
    condN(c)        = totalTrials;
    condG(c, :)     = G_a;
    condDf(c)       = df_c;
    condLogL(c, :)  = logL_a;
    condAlgoFrac{c} = algoFrac;
    condCounts{c}   = nCounts;
    condUsed(c)     = true;

    ip = find(pitchVals == pitchV, 1);
    ir = find(rollVals  == rollV,  1);
    logLPitchRoll(ip, ir, :) = squeeze(logLPitchRoll(ip, ir, :))' + logL_a;
    havePitchRoll(ip, ir) = true;

    fprintf('  [OK] %s  roll=%d  pitch=%d  N=%d  nPoses=%d\n', ...
        partName, round(rollV), round(pitchV), totalTrials, nPoses);
end

% Compact down to used conditions only
condPart     = condPart(condUsed);
condRoll     = condRoll(condUsed);
condPitch    = condPitch(condUsed); %#ok<NASGU> % kept for potential future use/debugging
condN        = condN(condUsed);
condG        = condG(condUsed, :);
condDf       = condDf(condUsed);
condLogL     = condLogL(condUsed, :);
condAlgoFrac = condAlgoFrac(condUsed);
condCounts   = condCounts(condUsed);

nConds = numel(condN);
if nConds == 0
    disp('No usable conditions found. Exiting.');
    return;
end
fprintf('\nTotal conditions used: %d\n\n', nConds);

% =========================================================================
% AGGREGATE -- G-TEST (absolute goodness of fit, per algorithm)
% =========================================================================
G_total   = sum(condG, 1);
df_total  = sum(condDf);
Gdf_total = G_total / df_total;   % quick-glance severity index (G/df)
pValueAsymptotic = zeros(1, nAlgos);
for a = 1:nAlgos
    % Upper-tail chi^2 p-value via regularized incomplete gamma function
    % (no Statistics Toolbox dependency). UNRELIABLE for this dataset --
    % see header notes on sparse-category df inflation. Kept only for
    % transparency/comparison against the bootstrap p-value below.
    pValueAsymptotic(a) = gammainc(G_total(a) / 2, df_total / 2, 'upper');
end

% =========================================================================
% AGGREGATE -- BAYESIAN POSTERIOR (relative comparison, pooled across all
% conditions). Reuses the same per-condition logL already computed above.
% =========================================================================
logL_total = sum(condLogL, 1);
w = logL_total + log(priorVec);
posterior = exp(w - max(w));
posterior = posterior / sum(posterior);

% Per-(pitch,roll) pooled posterior, for the severity-trend panel
postPitchRoll = nan(nP, nR, nAlgos);
for ip = 1:nP
    for ir = 1:nR
        if ~havePitchRoll(ip, ir), continue; end
        wpr = squeeze(logLPitchRoll(ip, ir, :))' + log(priorVec);
        ppr = exp(wpr - max(wpr));
        ppr = ppr / sum(ppr);
        postPitchRoll(ip, ir, :) = ppr;
    end
end

% =========================================================================
% PARAMETRIC BOOTSTRAP -- proper (non-asymptotic) p-value per algorithm.
%
% For algorithm a: assume its predicted distribution IS the true
% generating distribution. Simulate nBootstrap synthetic datasets from
% that assumption (same N per condition as the real data), compute G for
% each synthetic dataset against that same predicted distribution, and
% build an empirical null distribution of G_total. The bootstrap p-value
% is the fraction of simulated G_total values that are >= the REAL,
% observed G_total for that algorithm -- i.e. "if this algorithm's
% predictions were exactly correct, how often would we see a G this
% large or larger just from sampling noise?"
% =========================================================================
fprintf('Running parametric bootstrap (%d replicates per algorithm)...\n', nBootstrap);
rng(randSeed);

GbootTotal = zeros(nAlgos, nBootstrap);
for c = 1:nConds
    Nc = condN(c);
    for a = 1:nAlgos
        Gvals = simulateNullG(Nc, condAlgoFrac{c}(:, a), nBootstrap, epsilonFrac);
        GbootTotal(a, :) = GbootTotal(a, :) + Gvals;
    end
    if mod(c, 10) == 0
        fprintf('  ...condition %d / %d\n', c, nConds);
    end
end

pValueBoot = zeros(1, nAlgos);
for a = 1:nAlgos
    pValueBoot(a) = mean(GbootTotal(a, :) >= G_total(a));
end
fprintf('Bootstrap complete.\n\n');

% =========================================================================
% CALIBRATION (RELIABILITY DIAGRAM) -- absolute accuracy per algorithm.
%
% Pool every (condition, pose) predicted probability and its observed
% frequency (raw count / condition trial count), bin by predicted
% probability, and compute the trial-count-weighted mean predicted
% probability and observed frequency in each bin. A well-calibrated
% algorithm's points fall on the y = x diagonal.
%
% Bins are EQUAL-WEIGHT (quantile), not fixed-width -- see pickNBins and
% weightedQuantileEdges below.
%
% TWO variants are computed here:
%   calibX/calibY/calibW/calibErr           -- ALL (condition, pose) pairs.
%   calibX_obs/calibY_obs/calibW_obs/calibErr_obs
%                                            -- restricted to (condition,
%     pose) pairs where that pose was EXPERIMENTALLY OBSERVED at least
%     once in that specific condition (sim_count > 0).
% =========================================================================
fprintf('Computing calibration curves...\n');

% ---- equal-weight bin count + edges ----
targetCondsPerBin = 3;
minRowsPerBin     = 4;
minBins           = 4;
maxBins           = 15;

meanCondN = mean(condN);

poolP_all = []; poolN_all = [];
poolP_obs = []; poolN_obs = [];
for a = 1:nAlgos
    for c = 1:nConds
        pVec = condAlgoFrac{c}(:, a);
        nVec = condCounts{c};
        Nc   = condN(c);
        poolP_all = [poolP_all; pVec]; %#ok<AGROW>
        poolN_all = [poolN_all; repmat(Nc, numel(pVec), 1)]; %#ok<AGROW>
        obsSel = nVec > 0;
        poolP_obs = [poolP_obs; pVec(obsSel)]; %#ok<AGROW>
        poolN_obs = [poolN_obs; repmat(Nc, sum(obsSel), 1)]; %#ok<AGROW>
    end
end

nBinsAll = pickNBins(poolP_all, poolN_all, meanCondN, targetCondsPerBin, ...
    minRowsPerBin, minBins, maxBins);
nBinsObs = pickNBins(poolP_obs, poolN_obs, meanCondN, targetCondsPerBin, ...
    minRowsPerBin, minBins, maxBins);

binEdgesAll = weightedQuantileEdges(poolP_all, poolN_all, nBinsAll);
binEdgesObs = weightedQuantileEdges(poolP_obs, poolN_obs, nBinsObs);

nBinsAll = numel(binEdgesAll) - 1;   % may shrink if tied values collapsed edges
nBinsObs = numel(binEdgesObs) - 1;

fprintf('  All-poses calibration:      %d equal-weight bins (target >= %d conditions'' worth of trials/bin)\n', ...
    nBinsAll, targetCondsPerBin);
fprintf('  Observed-only calibration:  %d equal-weight bins (target >= %d conditions'' worth of trials/bin)\n', ...
    nBinsObs, targetCondsPerBin);

calibX = nan(nAlgos, nBinsAll);   % weighted mean predicted probability
calibY = nan(nAlgos, nBinsAll);   % weighted mean observed frequency
calibW = zeros(nAlgos, nBinsAll); % total weight (trials) in bin
calibErr = zeros(1, nAlgos);      % weighted mean |predicted - observed|

calibX_obs = nan(nAlgos, nBinsObs);
calibY_obs = nan(nAlgos, nBinsObs);
calibW_obs = zeros(nAlgos, nBinsObs);
calibErr_obs = zeros(1, nAlgos);
nObsPoints = zeros(1, nAlgos);    % how many (condition,pose) pairs qualify

for a = 1:nAlgos
    allP = [];
    allN = [];
    allCount = [];
    for c = 1:nConds
        pVec = condAlgoFrac{c}(:, a);
        nVec = condCounts{c};
        Nc   = condN(c);
        allP = [allP; pVec]; %#ok<AGROW>
        allN = [allN; repmat(Nc, numel(pVec), 1)]; %#ok<AGROW>
        allCount = [allCount; nVec]; %#ok<AGROW>
    end

    % ---- variant 1: all poses, equal-weight bins ----
    binIdx = discretize(allP, binEdgesAll);
    for b = 1:nBinsAll
        sel = (binIdx == b);
        if ~any(sel), continue; end
        sumN = sum(allN(sel));
        if sumN <= 0, continue; end
        calibX(a, b) = sum(allP(sel) .* allN(sel)) / sumN;
        calibY(a, b) = sum(allCount(sel)) / sumN;
        calibW(a, b) = sumN;
    end
    validBins = ~isnan(calibX(a, :)) & calibW(a, :) > 0;
    if any(validBins)
        calibErr(a) = sum(calibW(a, validBins) .* abs(calibX(a, validBins) - calibY(a, validBins))) ...
                      / sum(calibW(a, validBins));
    else
        calibErr(a) = NaN;
    end

    % ---- variant 2: observed-only (pose actually occurred in that
    % condition, i.e. sim_count > 0), equal-weight bins. No tunable
    % threshold is involved. ----
    obsSel = allCount > 0;
    nObsPoints(a) = sum(obsSel);
    Pobs = allP(obsSel); Nobs = allN(obsSel); Cobs = allCount(obsSel);
    binIdxObs = discretize(Pobs, binEdgesObs);
    for b = 1:nBinsObs
        sel = (binIdxObs == b);
        if ~any(sel), continue; end
        sumN = sum(Nobs(sel));
        if sumN <= 0, continue; end
        calibX_obs(a, b) = sum(Pobs(sel) .* Nobs(sel)) / sumN;
        calibY_obs(a, b) = sum(Cobs(sel)) / sumN;
        calibW_obs(a, b) = sumN;
    end
    validBinsObs = ~isnan(calibX_obs(a, :)) & calibW_obs(a, :) > 0;
    if any(validBinsObs)
        calibErr_obs(a) = sum(calibW_obs(a, validBinsObs) .* abs(calibX_obs(a, validBinsObs) - calibY_obs(a, validBinsObs))) ...
                           / sum(calibW_obs(a, validBinsObs));
    else
        calibErr_obs(a) = NaN;
    end
end
fprintf('Calibration complete.\n\n');

% =========================================================================
% GRAPHICS -- one image per graph (9 files total)
%   6 files : calibration reliability diagram, one per algorithm
%             (CSA family in shades of BLUE, CRSA family in shades of RED)
%   1 file  : algorithm support across conditions (pooled posterior per
%             pitch/roll), plain-language title, alpha/beta axis labels
%   1 file  : Bayesian posterior (headline relative comparison)
%   1 file  : G-test (absolute fit, G/df above bootstrap p-value)
% =========================================================================
fprintf('Building summary figures...\n');

% Color convention: CSA family = shades of BLUE (light -> dark for
% Dyn -> Stat -> Nan); CRSA family = shades of RED (light -> dark for
% Dyn -> Stat -> Nan). Order must match algoNames above.
algoColors = [0.68 0.85 0.95;   % CSA_Dyn   - light blue
              0.25 0.55 0.85;   % CSA_Stat  - medium blue
              0.02 0.20 0.45;   % CSA_Nan   - dark blue
              0.97 0.75 0.72;   % CRSA_Dyn  - light red
              0.85 0.30 0.25;   % CRSA_Stat - medium red
              0.50 0.03 0.03];  % CRSA_Nan  - dark red

modTag = '';
if opt.ApplyModifier, modTag = [modTag '_mod']; end
if opt.ExcludeFlaggedUnstable, modTag = [modTag '_stableOnly']; end

% ---- shared calibration axis scale, computed once across all six
% algorithms and both variants, so every calibration image below uses
% the identical scale for direct visual comparison ----
validAll    = ~isnan(calibX) & calibW > 0;
validObsAll = ~isnan(calibX_obs) & calibW_obs > 0;
axMaxCalib = max([calibX(validAll); calibY(validAll); ...
                   calibX_obs(validObsAll); calibY_obs(validObsAll)], [], 'omitnan') * 100 * 1.15;
if isempty(axMaxCalib) || ~isfinite(axMaxCalib) || axMaxCalib <= 0
    axMaxCalib = 100;
end

% ---- Files 1-6: one calibration image per algorithm (CSA in blue,
% CRSA in red; order follows algoNames) ----
for a = 1:nAlgos
    hFigC = figure('Visible', 'off', 'Color', 'w', 'Units', 'inches', ...
        'Position', [0 0 5 5]);
    axC = axes(hFigC);
    hold(axC, 'on'); grid(axC, 'on'); box(axC, 'on');

    plot(axC, [0 axMaxCalib], [0 axMaxCalib], '--', 'Color', [0.6 0.6 0.6], ...
        'LineWidth', 1.2, 'DisplayName', 'Perfect calibration');

    validBins = ~isnan(calibX(a, :)) & calibW(a, :) > 0;
    plot(axC, calibX(a, validBins)*100, calibY(a, validBins)*100, '-o', ...
        'Color', algoColors(a,:), 'MarkerFaceColor', algoColors(a,:), ...
        'LineWidth', 1.6, 'MarkerSize', 5, 'DisplayName', 'All poses');

    validBinsObs = ~isnan(calibX_obs(a, :)) & calibW_obs(a, :) > 0;
    plot(axC, calibX_obs(a, validBinsObs)*100, calibY_obs(a, validBinsObs)*100, '--s', ...
        'Color', algoColors(a,:), 'MarkerFaceColor', 'w', 'MarkerEdgeColor', algoColors(a,:), ...
        'LineWidth', 1.4, 'MarkerSize', 6, 'DisplayName', 'Observed-only');

    xlim(axC, [0 axMaxCalib]); ylim(axC, [0 axMaxCalib]);
    axis(axC, 'square');
    xlabel(axC, 'Predicted probability (%)');
    ylabel(axC, 'Observed frequency (%)');
    title(axC, sprintf('%s Calibration  (err: all=%.2f pp, obs=%.2f pp)', ...
        strrep(algoNames{a}, '_', '\_'), calibErr(a), calibErr_obs(a)), ...
        'FontWeight', 'bold', 'Color', algoColors(a,:), 'FontSize', 10);
    legend(axC, 'Location', 'northoutside', 'Orientation', 'horizontal', 'Box', 'off', 'FontSize', 7);

    outFileCalib = fullfile(resultsRoot, sprintf('GBayesianAnalysis_Calibration_%s%s.png', algoNames{a}, modTag));
    exportgraphics(hFigC, outFileCalib, 'Resolution', 200);
    delete(hFigC);
    fprintf('  Saved: %s\n', outFileCalib);
end

% ---- File 7: algorithm support across conditions ----
% x-axis uses alpha (roll) / beta (pitch) notation, e.g. "a15b25" for
% roll=15 deg, pitch=25 deg, slanted 45 degrees for readability.
hFigSupport = figure('Visible', 'off', 'Color', 'w', 'Units', 'inches', ...
    'Position', [0 0 13 6]);
ax3 = axes(hFigSupport);
hold(ax3, 'on'); grid(ax3, 'on'); box(ax3, 'on');
xLabels = {};
xVals   = [];
xi = 0;
plotData = nan(nP*nR, nAlgos);
for ip = 1:nP
    for ir = 1:nR
        if ~havePitchRoll(ip, ir), continue; end
        xi = xi + 1;
        xVals(end+1) = xi; %#ok<AGROW>
        xLabels{end+1} = sprintf('\\alpha%d\\beta%d', round(rollVals(ir)), round(pitchVals(ip))); %#ok<AGROW>
        plotData(xi, :) = squeeze(postPitchRoll(ip, ir, :))' * 100;
    end
end
plotData = plotData(1:xi, :);
for a = 1:nAlgos
    plot(ax3, xVals, plotData(:, a), '-o', 'Color', algoColors(a,:), ...
        'MarkerFaceColor', algoColors(a,:), 'LineWidth', 1.6, ...
        'MarkerSize', 4, 'DisplayName', strrep(algoNames{a}, '_', '\_'));
end
set(ax3, 'XTick', xVals, 'XTickLabel', xLabels, 'XTickLabelRotation', 45);
xlabel(ax3, '\alpha = roll angle, \beta = pitch angle (degrees)');
ylabel(ax3, 'Posterior probability (%)');
title(ax3, {'Algorithm Support Across Test Conditions', ...
    '(pooled Bayesian posterior per roll/pitch condition; blue = CSA family, red = CRSA family)'}, ...
    'FontWeight', 'bold');
ylim(ax3, [0, 100]);
legend(ax3, 'Location', 'northoutside', 'Orientation', 'horizontal', 'Box', 'off', 'FontSize', 8);

outFileSupport = fullfile(resultsRoot, sprintf('GBayesianAnalysis_AlgorithmSupport%s.png', modTag));
exportgraphics(hFigSupport, outFileSupport, 'Resolution', 200);
delete(hFigSupport);
fprintf('  Saved: %s\n', outFileSupport);

% ---- File 8: Bayesian posterior (headline relative-comparison result) ----
hFigPost = figure('Visible', 'off', 'Color', 'w', 'Units', 'inches', ...
    'Position', [0 0 7 5.5]);
ax1 = axes(hFigPost);
b1 = bar(ax1, posterior * 100, 'FaceColor', 'flat');
b1.CData = algoColors;
set(ax1, 'XTickLabel', strrep(algoNames, '_', '\_'), 'XTickLabelRotation', 30);
xlabel(ax1, 'Algorithm');
ylabel(ax1, 'Posterior probability (%)');
title(ax1, {'Bayesian Model Comparison', ...
    '(relative support among the 6 variants, flat prior)', ' '}, ...
    'FontWeight', 'bold');
ylim(ax1, [0, 100]);
grid(ax1, 'on');
for a = 1:nAlgos
    text(ax1, a, posterior(a)*100 + 3, sprintf('%.1f%%', posterior(a)*100), ...
        'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 9);
end

outFilePost = fullfile(resultsRoot, sprintf('GBayesianAnalysis_BayesianPosterior%s.png', modTag));
exportgraphics(hFigPost, outFilePost, 'Resolution', 200);
delete(hFigPost);
fprintf('  Saved: %s\n', outFilePost);

% ---- File 9: G-test with G/df + bootstrap p-value (absolute goodness
% of fit) ----
hFigG = figure('Visible', 'off', 'Color', 'w', 'Units', 'inches', ...
    'Position', [0 0 10 5.5]);
ax2 = axes(hFigG);
b2 = bar(ax2, G_total, 'FaceColor', 'flat');
b2.CData = algoColors;
set(ax2, 'XTickLabel', strrep(algoNames, '_', '\_'));
xlabel(ax2, 'Algorithm');
ylabel(ax2, sprintf('G statistic (df = %d total)', df_total));
title(ax2, {'G-Test Goodness of Fit', ...
    '(G/df severity index above bootstrap p-value; vs. observed pose counts)'}, ...
    'FontWeight', 'bold');
grid(ax2, 'on');
yMaxG = max(G_total) * 1.4;
if ~isfinite(yMaxG) || yMaxG <= 0, yMaxG = 1; end
ylim(ax2, [0, yMaxG]);
for a = 1:nAlgos
    if pValueBoot(a) < (1 / nBootstrap)
        pLabel = sprintf('p_{boot} < %.4f', 1/nBootstrap);
    else
        pLabel = sprintf('p_{boot} = %.4f', pValueBoot(a));
    end
    labelStr = sprintf('G/df = %.2f\n%s', Gdf_total(a), pLabel);
    text(ax2, a, G_total(a) + 0.03*yMaxG, labelStr, ...
        'HorizontalAlignment', 'center', 'FontSize', 8, 'VerticalAlignment', 'bottom');
end

outFileG = fullfile(resultsRoot, sprintf('GBayesianAnalysis_GTest%s.png', modTag));
exportgraphics(hFigG, outFileG, 'Resolution', 200);
delete(hFigG);
fprintf('  Saved: %s\n', outFileG);

% =========================================================================
% REPORT
% =========================================================================
reportFile = fullfile(resultsRoot, 'GBayesianAnalysis_Report.txt');
fid = fopen(reportFile, 'w');
if fid < 0
    error('Cannot open report file: %s', reportFile);
end

writeLine(fid, repmat('=', 1, 96));
writeLine(fid, '  G-TEST, BAYESIAN MODEL COMPARISON, AND CALIBRATION: CSA/CRSA (Dyn/Stat/Nan)');
writeLine(fid, '  Computed from RAW TRIAL COUNTS (sim_count, not normalized percentages),');
writeLine(fid, '  pooled across all parts and roll/pitch conditions.');
fprintf(fid, '  Source: %s\n', opt.summaryPath);
fprintf(fid, '  Conditions used: %d   |   Total df: %d   |   Bootstrap reps: %d\n', ...
    nConds, df_total, nBootstrap);
if opt.IncludeUnmatched
    writeLine(fid, '  Unmatched (pose_id == -1) trials INCLUDED as their own category per');
    writeLine(fid, '  condition, with every algorithm floored to epsilon probability for it.');
else
    writeLine(fid, '  Unmatched (pose_id == -1) trials EXCLUDED entirely (IncludeUnmatched=false).');
end
if opt.ExcludeFlaggedUnstable
    writeLine(fid, '  Flagged-unstable poses (flagged_unstable_by_geometry == TRUE, a FIXED');
    writeLine(fid, '  catalog-level label identical across all six algorithms) EXCLUDED from');
    writeLine(fid, '  the category set for every algorithm; geom_*_pct renormalized over the');
    writeLine(fid, '  remaining poses per condition; N shrunk to the surviving trial count.');
end
if opt.ApplyModifier
    writeLine(fid, '  Post-hoc modifier APPLIED: mod_i = 1/avg(ratio_wall, ratio_floor),');
    writeLine(fid, '  renormalized within each (part, alpha_deg, beta_deg) condition, BEFORE');
    writeLine(fid, '  the G-test / Bayesian / calibration analyses below.');
end
writeLine(fid, repmat('=', 1, 96));
writeLine(fid, '');
writeLine(fid, '  GLOSSARY');
writeLine(fid, '  CSA_Dyn / CRSA_Dyn   : transitioning poses zeroed, remaining renormalized ("before")');
writeLine(fid, '  CSA_Stat / CRSA_Stat : transitioning poses redistribute to destination pose ("after")');
writeLine(fid, '  CSA_Nan / CRSA_Nan   : no filtering at all, raw weight normalized ("none")');
writeLine(fid, '  G/df                 : G-statistic normalized by degrees of freedom (severity index)');
writeLine(fid, '');
writeLine(fid, repmat('-', 1, 96));
writeLine(fid, '');
writeLine(fid, '  G-TEST (absolute goodness of fit vs. observed data)');
writeLine(fid, '  G/df is a quick-glance severity index (values near 1 are broadly consistent with');
writeLine(fid, '  sampling noise; much larger values signal systematic misfit).');
writeLine(fid, '  p_asymptotic uses the standard chi^2(df) approximation and is UNRELIABLE');
writeLine(fid, '  here (sparse pose categories inflate df without proportional signal).');
writeLine(fid, '  p_bootstrap is a parametric-bootstrap p-value and is the trustworthy one:');
writeLine(fid, '  it simulates data under "this algorithm is exactly correct" and checks');
writeLine(fid, '  where the real G falls in that simulated distribution.');
writeLine(fid, '');
fprintf(fid, '  %-10s  %12s  %8s  %10s  %14s  %14s  %14s\n', ...
    'Algo', 'G_total', 'df', 'G/df', 'p_asymptotic', 'p_bootstrap', 'logL_total');
fprintf(fid, '  %s\n', repmat('-', 1, 86));
for a = 1:nAlgos
    fprintf(fid, '  %-10s  %12.3f  %8d  %10.3f  %14.4g  %14.4g  %14.2f\n', ...
        algoNames{a}, G_total(a), df_total, Gdf_total(a), pValueAsymptotic(a), pValueBoot(a), logL_total(a));
end
writeLine(fid, '');
writeLine(fid, repmat('-', 1, 96));
writeLine(fid, '');
writeLine(fid, '  BAYESIAN POSTERIOR (relative support among the 6 CSA/CRSA variants, flat prior)');
writeLine(fid, '  Reuses the per-condition log-likelihoods above; NOT an absolute');
writeLine(fid, '  correctness claim -- only meaningful as a comparison among these');
writeLine(fid, '  six candidates.');
writeLine(fid, '');
fprintf(fid, '  %-10s  %14s\n', 'Algo', 'Posterior (%)');
fprintf(fid, '  %s\n', repmat('-', 1, 28));
for a = 1:nAlgos
    fprintf(fid, '  %-10s  %14.2f\n', algoNames{a}, posterior(a)*100);
end
writeLine(fid, '');
writeLine(fid, repmat('-', 1, 96));
writeLine(fid, '');
writeLine(fid, '  CALIBRATION BINNING');
writeLine(fid, '  Bins are EQUAL-WEIGHT (quantile) bins, not fixed-width: edges are chosen');
writeLine(fid, '  so each bin carries roughly the same trial-count weight, rather than');
writeLine(fid, '  spacing edges evenly across the predicted-probability range.');
fprintf(fid, '  %d conditions and ~%.0f trials/condition (%.0f trials total) back this analysis.\n', ...
    nConds, meanCondN, sum(condN));
fprintf(fid, '  All-poses: %d bins. Observed-only: %d bins (fewer, since it has far less\n', ...
    nBinsAll, nBinsObs);
writeLine(fid, '  underlying data). Edges are shared across all six algorithms within each');
writeLine(fid, '  variant so the panels stay directly comparable.');
writeLine(fid, '');
writeLine(fid, '  CALIBRATION -- ALL POSES (absolute accuracy, independent of the other');
writeLine(fid, '  algorithms). For each algorithm, all (condition, pose) predicted');
writeLine(fid, '  probabilities are binned and compared to the trial-count-weighted');
writeLine(fid, '  observed frequency in that bin. calib_error is the weighted mean');
writeLine(fid, '  absolute gap between predicted and observed, in percentage points.');
writeLine(fid, '  CAUTION: this variant is dominated by the large number of poses that');
writeLine(fid, '  correctly never occur (predicted ~0%, observed ~0%), which inflates');
writeLine(fid, '  apparent calibration quality. See the OBSERVED-ONLY variant below for');
writeLine(fid, '  the harder, more meaningful test.');
writeLine(fid, '');
fprintf(fid, '  %-10s  %20s\n', 'Algo', 'calib_error (pp)');
fprintf(fid, '  %s\n', repmat('-', 1, 34));
for a = 1:nAlgos
    fprintf(fid, '  %-10s  %20.3f\n', algoNames{a}, calibErr(a));
end
writeLine(fid, '');
writeLine(fid, repmat('-', 1, 96));
writeLine(fid, '');
writeLine(fid, '  CALIBRATION -- OBSERVED-ONLY (the poses that actually happened)');
writeLine(fid, '  Restricted to (condition, pose) pairs where that pose was');
writeLine(fid, '  experimentally observed at least once in that specific condition');
writeLine(fid, '  (sim_count > 0). No tunable probability threshold is used -- the');
writeLine(fid, '  filter is simply "did this pose occur here" -- so this introduces no');
writeLine(fid, '  researcher-degrees-of-freedom risk. This answers: "when a pose is');
writeLine(fid, '  real, does this algorithm size its probability correctly?" -- a much');
writeLine(fid, '  harder and more relevant test than the all-poses variant above.');
writeLine(fid, '');
fprintf(fid, '  %-10s  %20s  %16s\n', 'Algo', 'calib_error (pp)', 'n data points');
fprintf(fid, '  %s\n', repmat('-', 1, 52));
for a = 1:nAlgos
    fprintf(fid, '  %-10s  %20.3f  %16d\n', algoNames{a}, calibErr_obs(a), nObsPoints(a));
end
writeLine(fid, '');
writeLine(fid, repmat('-', 1, 96));
writeLine(fid, '');
writeLine(fid, '  ZERO-CELL DIAGNOSTIC');
writeLine(fid, '  Count of (condition, pose) events where an algorithm predicted');
writeLine(fid, '  EXACTLY 0 probability for a pose that was experimentally observed');
writeLine(fid, '  (n_i > 0), including the synthetic unmatched category when');
writeLine(fid, '  IncludeUnmatched is true (no algorithm predicts "unmatched"). A floor');
writeLine(fid, '  of eps was substituted to keep G finite; a high count for a given');
writeLine(fid, '  algorithm is a red flag worth auditing directly.');
writeLine(fid, '');
fprintf(fid, '  %-10s  %14s\n', 'Algo', 'Zero-cell events');
fprintf(fid, '  %s\n', repmat('-', 1, 28));
for a = 1:nAlgos
    fprintf(fid, '  %-10s  %14d\n', algoNames{a}, zeroCellCount(a));
end
writeLine(fid, '');
writeLine(fid, repmat('=', 1, 96));
fclose(fid);
fprintf('Saved report: %s\n\n', reportFile);

disp('Done.');
end % GBayesianAnalysis


% =========================================================================
%  simulateNullG -- vectorized parametric bootstrap draws + G statistic
%
%  Given N trials and a true probability vector p (K x 1, fractions),
%  draws nBoot independent synthetic multinomial samples and computes G
%  (comparing each synthetic sample back against the SAME p that
%  generated it) for every replicate at once. Returns a 1 x nBoot vector.
% =========================================================================
function Gvals = simulateNullG(N, p, nBoot, epsilonFrac)
K = numel(p);
pFloor = p(:);
pFloor(pFloor <= 0) = epsilonFrac;
pFloor = pFloor / sum(pFloor);   % renormalize after flooring zeros

% Build sampling edges from pFloor (not raw p) so that every category
% has a strictly positive width -- discretize requires strictly
% increasing edges, and raw p can have exact-zero categories (algorithms
% that hard-exclude poses) which would otherwise produce duplicate,
% non-increasing edge values.
edges = [0, cumsum(pFloor')];
edges(end) = 1;   % guard against floating-point drift

r = rand(N, nBoot);
binIdx = discretize(r, edges);   % N x nBoot, values in 1..K

counts = zeros(K, nBoot);
for k = 1:K
    counts(k, :) = sum(binIdx == k, 1);
end

nz = counts > 0;

logTermSat = zeros(K, nBoot);
logTermSat(nz) = counts(nz) .* log(counts(nz) / N);

pMat = repmat(pFloor, 1, nBoot);
logTermAlgo = zeros(K, nBoot);
logTermAlgo(nz) = counts(nz) .* log(pMat(nz));

Gvals = 2 * (sum(logTermSat, 1) - sum(logTermAlgo, 1));
end


% =========================================================================
%  pickNBins -- adaptively choose an equal-weight bin count
%
%  Picks the number of calibration bins from whichever of two
%  constraints is more restrictive:
%    (a) trial-weight per bin: total weight / (targetCondsPerBin *
%        meanCondN) -- each bin should represent at least
%        targetCondsPerBin conditions' worth of trials, so a bin's
%        calibration point isn't driven by one condition's noise.
%    (b) row count per bin: numel(p) / minRowsPerBin -- each bin should
%        contain at least minRowsPerBin distinct (condition,pose) rows.
%  The result is clamped to [minBins, maxBins].
% =========================================================================
function nBins = pickNBins(p, w, meanCondN, targetCondsPerBin, minRowsPerBin, minBins, maxBins)
if isempty(p)
    nBins = minBins;
    return;
end
totalW    = sum(w);
nRows     = numel(p);
nByWeight = floor(totalW / (targetCondsPerBin * meanCondN));
nByRows   = floor(nRows / minRowsPerBin);
nBins = min([nByWeight, nByRows]);
nBins = max(nBins, minBins);
nBins = min(nBins, maxBins);
end


% =========================================================================
%  weightedQuantileEdges -- equal-weight (quantile) bin edges
%
%  Sorts (x,w) by x, walks the cumulative weight, and cuts wherever it
%  crosses each 1/nBins fraction of total weight, so every resulting bin
%  carries roughly equal trial-count weight rather than equal probability
%  width. Falls back gracefully if repeated x-values or a small sample
%  force duplicate cut points (fewer, wider bins than requested).
% =========================================================================
function edges = weightedQuantileEdges(x, w, nBins)
x = x(:); w = w(:);
[xs, ord] = sort(x);
ws = w(ord);
cw = cumsum(ws);
totalW = cw(end);

edges = zeros(1, nBins + 1);
edges(1) = 0;
targets = totalW * (1:nBins-1) / nBins;
for i = 1:numel(targets)
    j = find(cw >= targets(i), 1, 'first');
    edges(i+1) = xs(j);
end
edges(end) = max(x) * 1.0001 + eps;   % guard so the max value falls inside the last bin

edges = unique(edges, 'stable');      % collapse any duplicate cut points
edges = sort(edges);
if numel(edges) < 3
    % Degenerate case (too little spread in the data): fall back to a
    % single bin spanning the full range.
    edges = [0, max(x) * 1.0001 + eps];
end
end


% =========================================================================
%  writeLine
% =========================================================================
function writeLine(fid, str)
fprintf(fid, '%s\n', str);
end