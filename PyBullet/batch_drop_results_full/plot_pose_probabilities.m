function plot_pose_probabilities(summaryPath, varargin)
%PLOT_POSE_PROBABILITIES  Grouped bar chart: experimental (simulated) vs
%   one or more algorithmic (geometric) pose-frequency distributions.
%
%   Reads a "*_sim_summary.csv" (single condition) or "combined_summary.csv"
%   (many conditions) file written by chute_drop_batch.py's
%   write_condition_summary(), and plots each pose's experimental frequency
%   (sim_freq_pct, from the Monte-Carlo drop trials) next to one or more
%   matching algorithmic/geometric frequency columns (from
%   ChutePoseAnalysis.m, e.g. geom_CRSA_N_pct) as a grouped bar chart --
%   one figure per (part, alpha, beta) condition found in the file.
%
% --------------------------------------------------------------------
% USAGE
% --------------------------------------------------------------------
%   % single-condition file:
%   plot_pose_probabilities('Ql4i_25R_20P_sim_summary.csv')
%
%   % combined file covering multiple conditions -- makes one figure per
%   % (part, alpha, beta) combination found:
%   plot_pose_probabilities('combined_summary.csv')
%
%   % choose which algorithmic column to compare against (default: 'CRSA_N'):
%   plot_pose_probabilities('combined_summary.csv', 'GeomColumn', 'CSA_N')
%
%   % plot MULTIPLE algorithmic columns side by side in the same figure,
%   % all next to the same experimental bars:
%   plot_pose_probabilities('combined_summary.csv', ...
%       'GeomColumn', {'CSA_N','CRSA_N','CRSA_A'})
%
%   % apply a post-hoc per-pose modifier to the algorithmic percentages
%   % WITHOUT recalculating from raw weights. Because geom_*_pct is
%   % already normalized (sums to 100 per condition), renormalizing with
%   % a multiplicative modifier is mathematically identical to applying
%   % the modifier before normalization -- the original normalization
%   % constant cancels out:
%   %     P_i' = (P_i * mod_i) / sum_j(P_j * mod_j)
%   %          = (raw_i * mod_i) / sum_j(raw_j * mod_j)
%   % Default modifier: mod_i = 1 / mean(ratio_wall_i, ratio_floor_i).
%   % Applied independently to EACH selected GeomColumn (each is its own
%   % per-condition normalization, so each gets renormalized on its own).
%   plot_pose_probabilities('combined_summary.csv', 'ApplyModifier', true)
%
%   % save each figure as a PNG into a folder instead of just displaying:
%   plot_pose_probabilities('combined_summary.csv', 'SaveDir', 'figs')
%
%   % restrict to specific parts/conditions:
%   plot_pose_probabilities('combined_summary.csv', 'Parts', {'Ql4i'})
%
% --------------------------------------------------------------------
% NAME-VALUE OPTIONS
% --------------------------------------------------------------------
%   'GeomColumn'    Which geometric column(s) to plot against sim_freq_pct.
%                   A single string, OR a cell array of strings, each one
%                   of: 'CSA_B','CSA_A','CSA_N','CRSA_B','CRSA_A','CRSA_N'
%                   (matches the geom_<GeomColumn>_pct column names written
%                   by chute_drop_batch.py). Default: 'CRSA_N'.
%                       B = "before" merging by CSA weight
%                       A = "after"  merging
%                       N = normalized (sums to 100% across poses) -- this
%                           is almost always the right one to compare
%                           directly against sim_freq_pct, which is also a
%                           normalized (of all trials) percentage.
%                   When multiple columns are given, each becomes its own
%                   bar group in EVERY figure, all sharing the same
%                   experimental bars for that condition.
%   'IncludeUnmatched' true/false. If true, includes the pose_id == -1 row
%                   (drop trials whose settled orientation matched NO
%                   catalog pose within quat_match_tol) as its own bar
%                   category, with algorithmic value(s) forced to 0 since
%                   the geometric analysis has no equivalent "unmatched"
%                   bucket. Default: true (it's useful to SEE if unmatched
%                   trials are a nontrivial fraction -- if so, revisit
%                   --quat-match-tol or the pose catalog rather than
%                   silently dropping them from the comparison).
%   'ApplyModifier' true/false. If true, applies a per-pose multiplicative
%                   modifier to each selected algorithmic column and
%                   renormalizes within each condition, computed purely
%                   from columns already present in the summary CSV -- no
%                   recalculation against ChutePoseAnalysis.m raw weights
%                   is needed (see math note above). Default modifier is
%                   mod_i = 1/mean(ratio_wall_i, ratio_floor_i), using this
%                   CSV's 'ratio_wall'/'ratio_floor' columns. The unmatched
%                   (pose_id == -1) row has no ratios and is left at 0
%                   regardless. Default: false.
%   'ModifierFcn'   Optional function handle @(ratioWall, ratioFloor) -> mod_i,
%                   applied elementwise per pose, in case you want a
%                   different modifier than the 1/avg(Rf,Rw) default.
%                   Ignored if ApplyModifier is false.
%                   Default: @(rw, rf) 1 ./ mean([rw, rf], 2, 'omitnan')
%   'Parts'         Cell array of part-name strings to include (default:
%                   all parts found in the file).
%   'SortBy'        'pose_id' (default) or 'sim_freq_pct' (descending,
%                   useful for putting the biggest discrepancies/most
%                   common poses first). When multiple GeomColumns are
%                   given, 'sim_freq_pct' sort still uses the experimental
%                   column only (so all bar groups share one pose order).
%   'SaveDir'       If given, each figure is saved as PNG into this folder
%                   (created if it doesn't exist) instead of just being
%                   left on screen. Default: '' (no saving).
%   'ShowStats'     true/false. If true, annotates each figure with the
%                   Pearson correlation and RMSE between the experimental
%                   values and EACH algorithmic column's values, one line
%                   per column (computed over MATCHED poses only, i.e.
%                   excluding the unmatched bucket even if
%                   IncludeUnmatched is true). Default: true.
%
% --------------------------------------------------------------------
% EXPECTED COLUMNS (as written by chute_drop_batch.py)
% --------------------------------------------------------------------
%   part, alpha_deg, beta_deg, pose_id, sim_count, sim_freq_pct,
%   geom_CSA_B_pct, geom_CSA_A_pct, geom_CSA_N_pct,
%   geom_CRSA_B_pct, geom_CRSA_A_pct, geom_CRSA_N_pct,
%   ratio_wall, thresh_wall, ratio_floor, thresh_floor,
%   geom_notes, flagged_unstable_by_geometry
%
% Rows flagged flagged_unstable_by_geometry == true (i.e. the geometric
% analysis noted TRANSITIONS / FLOOR-UNSTABLE for that pose) are marked
% with a red outline on EVERY algorithmic bar for that pose, since a
% simulated hit there is landing on a pose the geometry considered
% invalid/transitional, regardless of which algorithmic column is shown.

    p = inputParser;
    addRequired(p, 'summaryPath', @(x) ischar(x) || isstring(x));
    addParameter(p, 'GeomColumn', 'CRSA_N', @(x) ischar(x) || isstring(x) || iscell(x));
    addParameter(p, 'IncludeUnmatched', true, @islogical);
    addParameter(p, 'ApplyModifier', false, @islogical);
    addParameter(p, 'ModifierFcn', @(rw, rf) 1 ./ mean([rw, rf], 2, 'omitnan'), ...
        @(x) isa(x, 'function_handle'));
    addParameter(p, 'Parts', {}, @iscell);
    addParameter(p, 'SortBy', 'pose_id', @(x) any(strcmpi(x, {'pose_id','sim_freq_pct'})));
    addParameter(p, 'SaveDir', '', @(x) ischar(x) || isstring(x));
    addParameter(p, 'ShowStats', true, @islogical);
    parse(p, summaryPath, varargin{:});
    opt = p.Results;

    if ~isfile(opt.summaryPath)
        error('plot_pose_probabilities:fileNotFound', ...
            'Could not find file: %s', opt.summaryPath);
    end

    % ------------------------------------------------------------------
    % Normalize GeomColumn into a cellstr, one or more entries
    % ------------------------------------------------------------------
    validGeomNames = {'CSA_B','CSA_A','CSA_N','CRSA_B','CRSA_A','CRSA_N'};
    if ischar(opt.GeomColumn) || isstring(opt.GeomColumn)
        geomNames = cellstr(opt.GeomColumn);
    else
        geomNames = cellfun(@char, opt.GeomColumn, 'UniformOutput', false);
    end
    geomNames = geomNames(:)';
    for i = 1:numel(geomNames)
        if ~any(strcmpi(geomNames{i}, validGeomNames))
            error('plot_pose_probabilities:invalidGeomColumn', ...
                'GeomColumn entry "%s" is invalid. Must be one of: %s', ...
                geomNames{i}, strjoin(validGeomNames, ', '));
        end
        geomNames{i} = upper(geomNames{i});
    end
    if numel(unique(geomNames)) ~= numel(geomNames)
        error('plot_pose_probabilities:duplicateGeomColumn', ...
            'GeomColumn contains duplicate entries: %s', strjoin(geomNames, ', '));
    end
    nGeom = numel(geomNames);

    T = readtable(opt.summaryPath, 'TextType', 'string');

    requiredCols = {'part','alpha_deg','beta_deg','pose_id','sim_freq_pct'};
    for i = 1:numel(requiredCols)
        if ~ismember(requiredCols{i}, T.Properties.VariableNames)
            error('plot_pose_probabilities:missingColumn', ...
                'Expected column "%s" not found in %s. Is this a *_sim_summary.csv or combined_summary.csv written by chute_drop_batch.py?', requiredCols{i}, opt.summaryPath);
        end
    end

    geomCols = strings(1, nGeom);
    for i = 1:nGeom
        gc = "geom_" + geomNames{i} + "_pct";
        if ~ismember(gc, T.Properties.VariableNames)
            error('plot_pose_probabilities:missingGeomColumn', ...
                'Column "%s" not found. Available geom_*_pct columns: %s', ...
                gc, strjoin(T.Properties.VariableNames( ...
                    startsWith(T.Properties.VariableNames, 'geom_') & ...
                    endsWith(T.Properties.VariableNames, '_pct')), ', '));
        end
        geomCols(i) = gc;
    end

    if opt.ApplyModifier
        modReqCols = {'ratio_wall','ratio_floor'};
        for i = 1:numel(modReqCols)
            if ~ismember(modReqCols{i}, T.Properties.VariableNames)
                error('plot_pose_probabilities:missingModifierColumn', ...
                    'ApplyModifier=true requires column "%s", not found in %s.', ...
                    modReqCols{i}, opt.summaryPath);
            end
        end
    end

    if ~isempty(opt.Parts)
        T = T(ismember(T.part, string(opt.Parts)), :);
        if isempty(T)
            error('plot_pose_probabilities:noRowsAfterPartFilter', ...
                'No rows remain after filtering to Parts = {%s}. Check part names match the "part" column exactly.', ...
                strjoin(string(opt.Parts), ', '));
        end
    end

    if ~opt.IncludeUnmatched
        T = T(T.pose_id ~= -1, :);
    end

    if ~isempty(opt.SaveDir) && ~isfolder(opt.SaveDir)
        mkdir(opt.SaveDir);
    end

    % Color palette for algorithmic bar groups (cycled if nGeom > palette size)
    algPalette = [ ...
        0.85 0.55 0.10; ...   % orange
        0.30 0.65 0.30; ...   % green
        0.55 0.30 0.75; ...   % purple
        0.10 0.65 0.65; ...   % teal
        0.75 0.30 0.45; ...   % magenta/rose
        0.45 0.45 0.10 ];     % olive
    expColor = [0.20 0.45 0.85];

    % --- one figure per (part, alpha_deg, beta_deg) condition ----------
    [conditions, ~, condIdx] = unique(T(:, {'part','alpha_deg','beta_deg'}), 'rows');

    for c = 1:height(conditions)
        rows = T(condIdx == c, :);

        switch lower(opt.SortBy)
            case 'pose_id'
                rows = sortrows(rows, 'pose_id');
            case 'sim_freq_pct'
                rows = sortrows(rows, 'sim_freq_pct', 'descend');
        end

        poseLabels = string(rows.pose_id);
        poseLabels(rows.pose_id == -1) = "unmatched";

        expVals = rows.sim_freq_pct;
        matchedMask = (rows.pose_id ~= -1);

        % algValsAll: nPoses x nGeom
        algValsAll = zeros(height(rows), nGeom);
        for gi = 1:nGeom
            v = double(rows.(geomCols(gi)));
            v(isnan(v)) = 0;
            algValsAll(:, gi) = v;
        end

        % ------------------------------------------------------------
        % POST-HOC MODIFIER (no recalculation against raw weights needed)
        %
        % geom_*_pct is already normalized within this condition:
        %   P_i = raw_i / sum(raw)
        % Renormalizing the *_pct values directly with a multiplicative
        % modifier is algebraically identical to applying the modifier
        % to the raw weights before normalizing, since the normalization
        % constant sum(raw) cancels:
        %   P_i * mod_i / sum_j(P_j * mod_j) = raw_i*mod_i / sum_j(raw_j*mod_j)
        % Applied independently per GeomColumn (each has its own raw
        % weight / normalization).
        % ------------------------------------------------------------
        if opt.ApplyModifier
            rw = double(rows.ratio_wall);
            rf = double(rows.ratio_floor);
            modVals = opt.ModifierFcn(rw, rf);
            modVals(~matchedMask) = 1;   % irrelevant: those rows are 0 anyway
            modVals(isnan(modVals) | isinf(modVals)) = 0;  % no ratio data -> can't weight it in

            for gi = 1:nGeom
                weighted = algValsAll(:, gi) .* modVals;
                denom = sum(weighted(matchedMask));
                if denom > 0
                    algValsAll(:, gi) = 100 * weighted / denom;
                else
                    algValsAll(:, gi) = zeros(height(rows), 1);
                end
            end
        end

        unstableFlag = false(height(rows), 1);
        if ismember('flagged_unstable_by_geometry', rows.Properties.VariableNames)
            fu = rows.flagged_unstable_by_geometry;
            if iscell(fu) || isstring(fu)
                unstableFlag = strcmpi(string(fu), 'true') | strcmpi(string(fu), '1');
            else
                unstableFlag = logical(fu);
            end
        end

        fig = figure('Name', sprintf('%s  alpha=%g  beta=%g', ...
            rows.part(1), rows.alpha_deg(1), rows.beta_deg(1)), ...
            'Color', 'w');
        ax = axes(fig); %#ok<LAXES>

        barData = [expVals, algValsAll];
        b = bar(ax, barData, 'grouped');
        b(1).FaceColor = expColor;
        b(1).DisplayName = 'Experimental (simulated)';
        for gi = 1:nGeom
            col = algPalette(mod(gi-1, size(algPalette,1)) + 1, :);
            b(gi+1).FaceColor = col;
            dispName = strrep(geomNames{gi}, '_', '\_');
            if opt.ApplyModifier
                b(gi+1).DisplayName = sprintf('Algorithmic (%s, modified)', dispName);
            else
                b(gi+1).DisplayName = sprintf('Algorithmic (%s)', dispName);
            end
        end

        % Outline every algorithmic bar (across all geom columns) for any
        % pose the geometric analysis flagged as unstable/transitional.
        if any(unstableFlag)
            hold(ax, 'on');
            firstOutline = true;
            for gi = 1:nGeom
                xPos = b(gi+1).XEndPoints;
                yVal = algValsAll(:, gi);
                hOut = plot(ax, xPos(unstableFlag), yVal(unstableFlag), ...
                    'rs', 'MarkerSize', 14, 'LineWidth', 1.5);
                if firstOutline
                    hOut.DisplayName = 'Flagged unstable/transitional (geometry)';
                    firstOutline = false;
                else
                    hOut.Annotation.LegendInformation.IconDisplayStyle = 'off';
                end
            end
            hold(ax, 'off');
        end

        ax.XTick = 1:height(rows);
        ax.XTickLabel = poseLabels;
        ax.XTickLabelRotation = 0;
        xlabel(ax, 'Pose ID');
        ylabel(ax, 'Frequency (%)');
        titleStr = sprintf('%s -- roll=%g\\circ, pitch=%g\\circ: experimental vs algorithmic pose probability', ...
            rows.part(1), rows.alpha_deg(1), rows.beta_deg(1));
        if opt.ApplyModifier
            titleStr = [titleStr, sprintf('  [modified: mod_i = 1/avg(R_f,R_w)]')]; %#ok<AGROW>
        end
        title(ax, titleStr);
        legend(ax, 'Location', 'best');
        grid(ax, 'on');
        box(ax, 'on');

        if opt.ShowStats
            e = expVals(matchedMask);
            statLines = cell(1, nGeom);
            for gi = 1:nGeom
                a = algValsAll(matchedMask, gi);
                if numel(e) >= 2 && std(e) > 0 && std(a) > 0
                    R = corrcoef(e, a);
                    rVal = R(1, 2);
                else
                    rVal = NaN;
                end
                rmse = sqrt(mean((e - a).^2));
                statLines{gi} = sprintf('%s: n=%d, r=%.3f, RMSE=%.2f pp', ...
                    geomNames{gi}, numel(e), rVal, rmse);
            end
            statStr = strjoin(statLines, newline);
            text(ax, 0.98, 0.95, statStr, 'Units', 'normalized', ...
                'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
                'BackgroundColor', [1 1 1 0.7], 'EdgeColor', [0.6 0.6 0.6]);
        end

        if ~isempty(opt.SaveDir)
            modTag = '';
            if opt.ApplyModifier, modTag = '_mod'; end
            geomTag = strjoin(geomNames, '-');
            fname = sprintf('%s_%gR_%gP_%s%s.png', rows.part(1), ...
                rows.alpha_deg(1), rows.beta_deg(1), geomTag, modTag);
            saveas(fig, fullfile(opt.SaveDir, fname));
            fprintf('Saved %s\n', fullfile(opt.SaveDir, fname));
        end
    end
end