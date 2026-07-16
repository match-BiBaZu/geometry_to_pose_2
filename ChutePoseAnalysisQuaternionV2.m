% =========================================================================
% ChutePoseAnalysis.m  —  Quaternion Edition (v4 — CSA/CRSA raw-weight parity)
% =========================================================================
%
% Top-level pipeline:
%   1. User picks a folder of STL parts and a roll/pitch angle grid.
%   2. For each part:
%        - load geometry, compute centroid, convex hull, candidate
%          floor-resting planes.
%   3. For each (roll, pitch) chute orientation:
%        - build the chute rotation (Rchute) and gravity/slide/wall
%          directions in the chute frame.
%        - enumerate candidate resting poses (theta sweep per floor plane).
%        - run primary + secondary stability checks.
%        - floor instability check: zero Qs where l_A/l_f >= THRESH_FLOOR.
%        - transition check (wall only): transitions(si) = ratioWall >= THRESH_WALL.
%          NOTE (fix): this check is PURELY GEOMETRIC (based on ratioWall
%          alone) and is no longer gated on Qs>0. The floor-instability
%          check and the wall-transition check are independent geometric
%          tests; each of CSA/CRSA applies the shared "transitions" mask
%          to its own raw weight afterward. Previously "transitions" was
%          forced false for any pose whose CSA weight (Qs) had already
%          been zeroed by the floor check, which silently suppressed the
%          wall-transition filtering for CRSA_B/CRSA_A as well (since they
%          share this same mask) whenever a pose's CSA weight was zero.
%        - merge duplicate poses, compute CSA and CRSA scores.
%        - evaluate six probability algorithms (CSA_B, CSA_A, CSA_N, CRSA_B, CRSA_A, CRSA_N).
%        - export per-pose PDF figures and text summaries:
%            _GEOMETRIC_summary.txt  — all poses after mergeByTheta
%            _MASTER_summary.txt     — stable poses with full analysis
%        - export critical value scatter figures to "critical value metrics" folder:
%            transitionAnalysis  — ratioWall vs pose, red bg = transitioning
%            stabilityAnalysis   — ratioFloor vs pose, red bg = unstable
%
% Critical ratio thresholds (geometry-derived, no friction):
%   THRESH_WALL  = 2.731  — l_A/l_w wall ratio above which pose transitions
%   THRESH_FLOOR = 2.296  — l_A/l_f floor ratio above which pose is unstable
%
% Six sub-algorithms (v4 — CSA/CRSA RAW-WEIGHT PARITY):
%   CSA_B  = CSA,  transitioning poses zeroed, remaining renormalized  ("before")
%   CSA_A  = CSA,  transitioning poses redistribute to destination     ("after")
%   CSA_N  = CSA,  no filtering at all (raw CSA weight normalized)     ("none")
%   CRSA_B = CRSA, transitioning poses zeroed, remaining renormalized  ("before")
%   CRSA_A = CRSA, transitioning poses redistribute to destination     ("after")
%   CRSA_N = CRSA, no filtering at all (raw CRSA weight normalized)    ("none")
%
%   IMPORTANT (fix in this version): CSA_N/CSA_B/CSA_A are now computed
%   from QsRaw — a snapshot of the CSA weight taken immediately after
%   computeCSA(), BEFORE the secondary-stability check and the floor-
%   instability check zero anything out. This mirrors exactly how
%   CRSA_N/CRSA_B/CRSA_A are computed from crsaRaw, which is never touched
%   by those two checks. Previously, CSA_N/CSA_B/CSA_A were all computed
%   from the same (already-zeroed) "Qs", and because the wall-transition
%   condition (ratioWall >= THRESH_WALL) is geometrically almost redundant
%   with the secondary-stability / floor-instability checks (both are
%   frictionless-tipping criteria, just formulated differently), nearly
%   every pose with transitions(si)==true already had Qs(si)==0 by the
%   time CSA_B/CSA_A were computed — making the "zero transitioning poses"
%   and "redistribute transitioning poses" steps no-ops, so CSA_B == CSA_A
%   == CSA_N in practice. QsRaw fixes this by giving CSA_B/CSA_A something
%   non-zero to actually zero/redistribute, exactly like CRSA_B/CRSA_A do
%   with crsaRaw.
%
%   Qs itself (the secondary/floor-zeroed version) is UNCHANGED in this
%   version and continues to drive pose merging (mergeByCSAValues),
%   sumQ, the PDF figures, the critical-value scatter plots, and the
%   MASTER summary "Q=0 / FLOOR-UNSTABLE" notes — so none of that
%   diagnostic machinery changes behavior. Only which weight vector feeds
%   probCSA_N/probCSA_B/probCSA_A has changed.
% =========================================================================

function ChutePoseAnalysis()

close all; clc;

% -------------------------------------------------------------------------
% CRITICAL RATIO THRESHOLDS
% -------------------------------------------------------------------------
THRESH_WALL  = 2.731;   % l_A/l_w  — wall transition threshold
THRESH_FLOOR = 2.296;   % l_A/l_f  — floor instability threshold

% -------------------------------------------------------------------------
% TOLERANCES
% -------------------------------------------------------------------------
wallTol       = 0.1;
planeTol      = 0.08;
thetaTol_deg  = 2;
planeMergeTol = 0.05;
omegaTol      = 0.005;
hTol          = 0.5;
quatMatchTol  = 0.05;

% =========================================================================
% STEP 1 — SELECT PARTS FOLDER
% =========================================================================
partsFolder = uigetdir(pwd, 'Select folder containing STL files');
if isequal(partsFolder, 0)
    disp('No folder selected. Exiting.'); return;
end

stlFiles = dir(fullfile(partsFolder, '*.stl'));
if isempty(stlFiles)
    errordlg('No STL files found in selected folder.', 'Error'); return;
end

fprintf('Found %d STL file(s).\n', numel(stlFiles));

% =========================================================================
% STEP 2 — ANGLE INPUT DIALOGS
% =========================================================================
rollChoice = questdlg( ...
    'Roll angle (about X-axis): static value or range?', ...
    'Roll Input', 'Static', 'Range', 'Static');
if isempty(rollChoice), disp('Cancelled.'); return; end

if strcmp(rollChoice, 'Static')
    ans_r = inputdlg({'Roll angle (deg):'}, 'Roll - Static', 1, {'0'});
    if isempty(ans_r), disp('Cancelled.'); return; end
    rollAngles = str2double(ans_r{1});
else
    ans_r = inputdlg( ...
        {'Start (deg):', 'End (deg):', 'Increment (deg):'}, ...
        'Roll - Range', 1, {'0','30','10'});
    if isempty(ans_r), disp('Cancelled.'); return; end
    rollAngles = str2double(ans_r{1}) : str2double(ans_r{3}) : str2double(ans_r{2});
end

pitchChoice = questdlg( ...
    'Pitch angle (about Y-axis): static value or range?', ...
    'Pitch Input', 'Static', 'Range', 'Static');
if isempty(pitchChoice), disp('Cancelled.'); return; end

if strcmp(pitchChoice, 'Static')
    ans_p = inputdlg({'Pitch angle (deg):'}, 'Pitch - Static', 1, {'0'});
    if isempty(ans_p), disp('Cancelled.'); return; end
    pitchAngles = str2double(ans_p{1});
else
    ans_p = inputdlg( ...
        {'Start (deg):', 'End (deg):', 'Increment (deg):'}, ...
        'Pitch - Range', 1, {'0','30','10'});
    if isempty(ans_p), disp('Cancelled.'); return; end
    pitchAngles = str2double(ans_p{1}) : str2double(ans_p{3}) : str2double(ans_p{2});
end

fprintf('Roll angles:  %s deg\n', mat2str(round(rollAngles)));
fprintf('Pitch angles: %s deg\n', mat2str(round(pitchAngles)));

% =========================================================================
% STEP 3 — MAIN LOOP: STL x ROLL x PITCH
% =========================================================================
for fi = 1:numel(stlFiles)

    filePath = fullfile(partsFolder, stlFiles(fi).name);
    [~, partName, ~] = fileparts(filePath);

    fprintf('\n=================================================================\n');
    fprintf('Part: %s\n', partName);
    fprintf('=================================================================\n');

    % ------------------------------------------------------------------
    % LOAD GEOMETRY
    % ------------------------------------------------------------------
    partGeometry = fegeometry(filePath);
    partMesh     = generateMesh(partGeometry);

    facesAndVertices = triangulation(partMesh);
    partFaces        = facesAndVertices.ConnectivityList;
    partVertices     = facesAndVertices.Points;

    partSpan  = max(partVertices) - min(partVertices);
    partScale = max(partSpan);
    hTol_part = partScale * 0.03;

    convexHullFaces = convhull(partVertices, 'Simplify', true);
    chullVertexIdx  = unique(convexHullFaces(:));

    centroidCoordinates = centroidOfPolyhedron(partVertices, partFaces);
    fprintf('Centroid: [%s]\n', num2str(centroidCoordinates));

    % ------------------------------------------------------------------
    % FIND FLOOR PLANES
    % ------------------------------------------------------------------
    [restingPlaneVertices, restingPlaneEquations] = ...
        findFloorPlanes(partVertices, convexHullFaces, planeTol);

    fprintf('Found %d candidate floor resting planes.\n', numel(restingPlaneVertices));

    % ------------------------------------------------------------------
    % ANGLE GRID LOOP
    % ------------------------------------------------------------------
    for ri = 1:numel(rollAngles)
        for pi_idx = 1:numel(pitchAngles)

            clc;
            fprintf('Part: %s\n', partName);
            fprintf('Roll angles:  %s deg\n', mat2str(round(rollAngles)));
            fprintf('Pitch angles: %s deg\n', mat2str(round(pitchAngles)));

            chuteRoll_deg  = rollAngles(ri);
            chutePitch_deg = pitchAngles(pi_idx);

            folderName  = sprintf('%s_%dR_%dP', partName, ...
                round(chuteRoll_deg), round(chutePitch_deg));
            summaryRoot = fullfile(partsFolder, 'results', partName, folderName);
            if ~exist(summaryRoot, 'dir'), mkdir(summaryRoot); end

            % Critical value metrics output folder
            metricsRoot = fullfile(partsFolder, 'results', partName, ...
                folderName, 'critical value metrics');
            if ~exist(metricsRoot, 'dir'), mkdir(metricsRoot); end

            runTag = folderName;
            fprintf('\n--- %s ---\n', runTag);

            % ----------------------------------------------------------
            % CHUTE FRAME
            % ----------------------------------------------------------
            roll  = -deg2rad(chuteRoll_deg);
            pitch =  deg2rad(chutePitch_deg);

            q_roll  = q_fromAxisAngle([1,0,0], roll);
            q_pitch = q_fromAxisAngle([0,1,0], pitch);
            q_chute = q_compose(q_roll, q_pitch);

            Rchute     = q_toRotm(q_chute);
            gravityDir = (Rchute' * [0; 0; -1])';

            slideDir    = normaliseVec((Rchute * [1; 0; 0])');
            floorNorm_c = normaliseVec((Rchute * [0; 0; 1])');
            wallNorm_c  = normaliseVec((Rchute * [0; 1; 0])');

            % ----------------------------------------------------------
            % POSE ENUMERATION
            % ----------------------------------------------------------
            [rawCandidates, planeQuats] = thetaFromHullEdges( ...
                partVertices, convexHullFaces, chullVertexIdx, ...
                restingPlaneVertices, restingPlaneEquations, ...
                centroidCoordinates, wallTol, planeTol);

            fprintf('Raw candidates: %d\n', numel(rawCandidates));

            allRestingPositions = mergeByTheta( ...
                rawCandidates, thetaTol_deg, planeMergeTol, partVertices);

            numCandidates = numel(allRestingPositions);
            fprintf('After theta merging: %d\n', numCandidates);

            clear rawCandidates;

            % ----------------------------------------------------------
            % GEOMETRIC SUMMARY — all poses after mergeByTheta
            % ----------------------------------------------------------
            writeGeometricSummary( ...
                allRestingPositions, planeQuats, ...
                summaryRoot, folderName, partName, ...
                chuteRoll_deg, chutePitch_deg);

            % ----------------------------------------------------------
            % PRIMARY STABILITY CHECK
            % ----------------------------------------------------------
            gd = gravityDir(:) / norm(gravityDir);
            if abs(gd(1)) < 0.9, tmp = [1;0;0]; else, tmp = [0;1;0]; end
            e1g = cross(gd, tmp);  e1g = e1g / norm(e1g);
            e2g = cross(gd, e1g);  e2g = e2g / norm(e2g);

            stablePositions = false(numCandidates, 1);
            for ci = 1:numCandidates
                pos    = allRestingPositions(ci);
                vertsW = pos.verticesWorld;
                centW  = pos.centroidWorld(:)';
                chVerts = vertsW(chullVertexIdx, :);
                pts2D   = [chVerts * e1g, chVerts * e2g];
                cent2D  = [centW  * e1g,  centW  * e2g];
                try
                    hIdx_    = convhull(pts2D(:,1), pts2D(:,2));
                    hullPoly = pts2D(hIdx_, :);
                    stablePositions(ci) = inpolygon( ...
                        cent2D(1), cent2D(2), hullPoly(:,1), hullPoly(:,2));
                catch
                    stablePositions(ci) = false;
                end
            end

            numStable = sum(stablePositions);
            fprintf('Stable positions: %d\n', numStable);

            if numStable == 0
                fprintf('No stable positions for %s. Skipping.\n', runTag);
                continue;
            end

            stableIdx = find(stablePositions);

            % ----------------------------------------------------------
            % CSA CALCULATION
            % ----------------------------------------------------------
            gWorld = [0;0;-1];
            [Qs, omegas, heights, ~, hullProjs] = computeCSA( ...
                allRestingPositions, stableIdx, chullVertexIdx, Rchute, gWorld);

            % Raw CSA weight, captured BEFORE any validity zeroing below,
            % so CSA_N/CSA_B/CSA_A can be computed from an untouched weight
            % source exactly the way CRSA_N/CRSA_B/CRSA_A are computed from
            % untouched crsaRaw (see v4 header note). "Qs" continues to be
            % zeroed by the secondary-stability and floor-instability
            % checks below and remains the weight used for pose merging,
            % display/plots, sumQ, and MASTER summary notes.
            QsRaw = Qs;

            % ----------------------------------------------------------
            % SECONDARY STABILITY CHECK
            % ----------------------------------------------------------
            secondaryUnstable = false(numStable, 1);
            for si = 1:numStable
                pi_   = stableIdx(si);
                pos   = allRestingPositions(pi_);
                C0    = pos.centroidWorld(:)';
                vertsW_si = (Rchute * (pos.verticesWorld' - C0'))' + C0;
                centW_si  = C0;

                supportIdx_si = union(pos.floorContactVertIdx, pos.wallContactVertIdx);
                supportVerts  = vertsW_si(supportIdx_si, :);

                if size(supportVerts, 1) < 3
                    secondaryUnstable(si) = true;
                else
                    gravCoords = supportVerts * gWorld;
                    minGrav    = min(gravCoords);
                    planePoint = minGrav * gWorld';
                    supportProj = zeros(size(supportVerts));
                    for kk = 1:size(supportVerts,1)
                        p_  = supportVerts(kk,:);
                        d_  = dot((p_ - planePoint), gWorld');
                        supportProj(kk,:) = p_ - d_ * gWorld';
                    end
                    supportProj = uniquetol(supportProj, 1e-5, 'ByRows', true);

                    if size(supportProj, 1) < 3
                        secondaryUnstable(si) = true;
                    else
                        if abs(gWorld(1)) < 0.9, tmp2 = [1;0;0]; else, tmp2 = [0;1;0]; end
                        e1b = cross(gWorld, tmp2); e1b = e1b/norm(e1b);
                        e2b = cross(gWorld, e1b);  e2b = e2b/norm(e2b);
                        support2D = [supportProj*e1b, supportProj*e2b];
                        cent2D_si = [dot(centW_si, e1b), dot(centW_si, e2b)];

                        try
                            hIdx2D = convhull(support2D(:,1), support2D(:,2));
                            hPoly2D = support2D(hIdx2D, :);
                            inHull = inpolygon(cent2D_si(1), cent2D_si(2), ...
                                hPoly2D(:,1), hPoly2D(:,2));
                            if ~inHull
                                secondaryUnstable(si) = true;
                            end
                        catch
                            secondaryUnstable(si) = true;
                        end
                    end
                end

                if secondaryUnstable(si)
                    Qs(si) = 0;
                end
            end

            numSecondaryFailed = sum(secondaryUnstable);
            if numSecondaryFailed > 0
                fprintf('Secondary stability zeroed: %d pose(s)\n', numSecondaryFailed);
            end

            % ----------------------------------------------------------
            % TRANSITION MOMENT ARM RATIOS
            % Must be computed BEFORE floor instability and transition
            % checks so ratioFloor and ratioWall are available.
            % ----------------------------------------------------------
            [ratioWall, ratioFloor, ~, momentArmGeo] = computeTransitionRatios( ...
                allRestingPositions, stableIdx, Rchute, ...
                slideDir, floorNorm_c, wallNorm_c, planeTol);

            % ----------------------------------------------------------
            % FLOOR INSTABILITY CHECK (geometry-based)
            % If l_A/l_f >= THRESH_FLOOR the pose is geometrically
            % unstable — zero its Qs weight.
            % ----------------------------------------------------------
            floorUnstable = false(numStable, 1);
            for si = 1:numStable
                if ratioFloor(si) >= THRESH_FLOOR
                    floorUnstable(si) = true;
                    Qs(si) = 0;
                    fprintf('  Floor instability (rF=%.3f >= %.3f): zeroed Pose %d\n', ...
                        ratioFloor(si), THRESH_FLOOR, si);
                end
            end

            numFloorUnstable = sum(floorUnstable);
            if numFloorUnstable > 0
                fprintf('Floor instability check zeroed: %d pose(s)\n', numFloorUnstable);
            end

            % ----------------------------------------------------------
            % WALL TRANSITION CHECK (geometry-based, wall only)
            % FIX: this is a PURE geometric test on ratioWall only.
            % It must NOT be gated on Qs(si), since Qs is a CSA-specific
            % weight already zeroed by the (independent) floor-instability
            % check above. Gating on Qs==0 previously suppressed the wall
            % transition flag for exactly the poses that most needed it,
            % which in turn suppressed CRSA_B/CRSA_A filtering for those
            % poses since CRSA shares this same "transitions" mask.
            % ----------------------------------------------------------
            transitions = false(numStable, 1);
            for si = 1:numStable
                transitions(si) = (ratioWall(si) >= THRESH_WALL);
            end

            % ----------------------------------------------------------
            % MERGE DUPLICATE POSES
            % (QsRaw is merged in lock-step with Qs — same groupings,
            %  same summing — so probCSA_N/B/A stay index-aligned with
            %  stableIdx/ratioWall/ratioFloor/crsaRaw after this point.)
            % ----------------------------------------------------------
            [stableIdx, Qs, omegas, heights, mergedGroups, QsRaw] = mergeByCSAValues( ...
                stableIdx, Qs, omegas, heights, omegaTol, hTol_part, ...
                allRestingPositions, Rchute, QsRaw);

            % Re-index ratioWall, ratioFloor, transitions to merged set
            % (mergeByCSAValues keeps the best-Q representative per group,
            %  whose index position in the pre-merge arrays is bestGlobal)
            numStable = numel(stableIdx);
            sumQ      = sum(Qs);

            fprintf('After CSA value merging: %d\n', numStable);

            % Recompute ratios on the merged stable set
            [ratioWall, ratioFloor, ~, momentArmGeo] = computeTransitionRatios( ...
                allRestingPositions, stableIdx, Rchute, ...
                slideDir, floorNorm_c, wallNorm_c, planeTol);

            % Re-apply floor-instability check after merge (in case merge
            % changed indices). This only affects Qs (the display/merge
            % weight), NOT QsRaw, and is independent of the wall-transition
            % mask below. floorUnstableMerged records WHICH poses failed
            % this check, for use in the _B/_A instability logic below.
            floorUnstableMerged = false(numStable, 1);
            for si = 1:numStable
                if ratioFloor(si) >= THRESH_FLOOR
                    Qs(si) = 0;
                    floorUnstableMerged(si) = true;
                end
            end
            sumQ = sum(Qs);

            % Re-apply wall-transition check after merge — again purely
            % geometric (ratioWall only), NOT gated on Qs.
            transitions = false(numStable, 1);
            for si = 1:numStable
                transitions(si) = (ratioWall(si) >= THRESH_WALL);
            end

            % ----------------------------------------------------------
            % COMBINED INSTABILITY FLAG (for _B) + DOMINANT-FAILURE-MODE
            % CLASSIFICATION (for _A)
            %
            % A pose is treated as "not a legitimate final resting pose"
            % if it fails EITHER the wall-transition check OR the
            % floor-instability check. This combined flag drives _B,
            % which does not attempt to trace a destination — it just
            % needs to know a pose is invalid, then diffuses its mass
            % proportionally across whatever remains.
            %
            % _A is more particular: it can only trace a destination for
            % wall-transitioning poses (computeTransitionChains only
            % follows wall-tip rotations — there is currently no floor-tip
            % chain follower). For a pose that fails BOTH checks, we ask
            % which failure is more severe — using normalized excess
            % (ratio/threshold - 1) so wall and floor margins are
            % comparable on the same scale — and let that dominant mode
            % decide how _A treats the pose:
            %   wall dominant  -> trace the wall-tip chain to a specific
            %                     destination pose (transitionDest)
            %   floor dominant -> no known destination; zero the pose's
            %                     weight with no redistribution target
            %                     (behaves like _B for that single pose)
            % ----------------------------------------------------------
            instabilityFlag = transitions | floorUnstableMerged;

            wallMargin  = (ratioWall  ./ THRESH_WALL)  - 1;
            floorMargin = (ratioFloor ./ THRESH_FLOOR) - 1;

            wallDominant  = false(numStable, 1);
            floorDominant = false(numStable, 1);

            for si = 1:numStable
                if transitions(si) && floorUnstableMerged(si)
                    % Both checks failed — the larger normalized excess
                    % identifies the dominant (more severe) failure mode.
                    if wallMargin(si) >= floorMargin(si)
                        wallDominant(si)  = true;
                    else
                        floorDominant(si) = true;
                    end
                elseif transitions(si)
                    wallDominant(si)  = true;
                elseif floorUnstableMerged(si)
                    floorDominant(si) = true;
                end
            end

            % ----------------------------------------------------------
            % CRSA RAW SCORES
            % ----------------------------------------------------------
            crsaRaw = computeCRSA( ...
                allRestingPositions, stableIdx, Rchute, gWorld, wallNorm_c);

            % ----------------------------------------------------------
            % SIX SUB-ALGORITHMS: CSA_B, CSA_A, CSA_N, CRSA_B, CRSA_A, CRSA_N
            %
            % CSA_N/CSA_B/CSA_A are computed from QsRaw (unzeroed), and
            % CRSA_N/CRSA_B/CRSA_A are computed from crsaRaw (unzeroed) —
            % same helper functions, same masks, applied to each
            % algorithm's own untouched raw weight.
            %
            % _B uses the COMBINED instabilityFlag (wall OR floor) since
            % it never attempts to trace a destination — it just needs to
            % know a pose is invalid.
            %
            % _A uses the dominant-failure-mode classification above:
            % wall-dominant poses get the traced-chain treatment via
            % transitionDest; floor-dominant poses get zeroed with no
            % redistribution target, since no floor-tip chain exists yet.
            % ----------------------------------------------------------

            % CSA_N — CSA, no filtering at all
            probCSA_N = normalizeWeights(QsRaw);

            % CSA_B — CSA, invalid poses (wall OR floor) zeroed, renormalized
            probCSA_B = zeroAndNormalize(QsRaw, instabilityFlag);

            % CRSA_N — CRSA, no filtering at all
            probCRSA_N = normalizeWeights(crsaRaw);

            % CRSA_B — CRSA, invalid poses (wall OR floor) zeroed, renormalized
            probCRSA_B = zeroAndNormalize(crsaRaw, instabilityFlag);

            % Transition chains (pure geometry — wall-tip tracing only,
            % computed for every wall-transitioning pose regardless of
            % whether it's wall- or floor-dominant; the dominant-mode
            % masks below decide whether _A actually uses this trace)
            [transitionDest, transitionRotatedVerts, transitionRotatedCent, ...
                refQuats, refQuatsAll, chainLog] = ...
                computeTransitionChains( ...
                allRestingPositions, stableIdx, transitions, ...
                ratioWall, THRESH_WALL, ...
                Rchute, q_chute, quatMatchTol, centroidCoordinates, ...
                planeQuats, mergedGroups);

            % For _A: floor-dominant poses have no traceable destination
            % (no floor-tip chain follower exists), so their raw weight
            % is zeroed pre-emptively and their transitionDest is forced
            % to 0 so redistributeWeights doesn't (mis)redirect them via
            % a wall-trace that isn't actually their dominant failure mode.
            QsRaw_forA          = QsRaw;
            crsaRaw_forA        = crsaRaw;
            QsRaw_forA(floorDominant)   = 0;
            crsaRaw_forA(floorDominant) = 0;

            transitionDest_forA = transitionDest;
            transitionDest_forA(floorDominant) = 0;

            % CSA_A — wall-dominant poses redistribute to traced destination;
            % floor-dominant poses zeroed with no redistribution target
            probCSA_A = redistributeWeights(QsRaw_forA, transitionDest_forA);

            % CRSA_A — same logic, CRSA raw weight
            probCRSA_A = redistributeWeights(crsaRaw_forA, transitionDest_forA);

            % ----------------------------------------------------------
            % CONSOLE OUTPUT
            % ----------------------------------------------------------
            fprintf('\n%-4s %-12s %-12s  %-9s %-9s %-9s %-9s %-9s %-9s %-10s %-10s\n', ...
                'Pos','Floor Plane','Theta (deg)', ...
                'CSA_B(%)','CSA_A(%)','CSA_N(%)','CRSA_B(%)','CRSA_A(%)','CRSA_N(%)', ...
                'rWall','rFloor');
            fprintf('%s\n', repmat('-',1,100));
            for si = 1:numStable
                ci_ = stableIdx(si);
                pos = allRestingPositions(ci_);
                fprintf('%-4d %-12d %-12.0f  %-9.2f %-9.2f %-9.2f %-9.2f %-9.2f %-9.2f %-10.3f %-10.3f\n', ...
                    si, pos.floorPlaneIdx, round(rad2deg(pos.theta)), ...
                    probCSA_B(si), probCSA_A(si), probCSA_N(si), ...
                    probCRSA_B(si), probCRSA_A(si), probCRSA_N(si), ...
                    ratioWall(si), ratioFloor(si));
            end

            % ----------------------------------------------------------
            % PDF + MASTER SUMMARY
            % ----------------------------------------------------------
            pdfPath = fullfile(summaryRoot, [folderName, '.pdf']);

            saveFiguresAndPDF( ...
                allRestingPositions, stableIdx, ...
                partVertices, partFaces, convexHullFaces, chullVertexIdx, ...
                Rchute, Qs, omegas, heights, sumQ, ...
                probCSA_B, probCSA_A, probCSA_N, probCRSA_B, probCRSA_A, probCRSA_N, ...
                transitionDest, transitions, ...
                chuteRoll_deg, chutePitch_deg, ...
                partName, pdfPath, ...
                slideDir, floorNorm_c, wallNorm_c, ...
                ratioWall, ratioFloor, THRESH_WALL, THRESH_FLOOR, planeTol, ...
                transitionRotatedVerts, transitionRotatedCent, ...
                refQuats, refQuatsAll, quatMatchTol, chainLog, momentArmGeo);

            writeMasterSummary( ...
                allRestingPositions, stableIdx, Qs, omegas, heights, sumQ, ...
                probCSA_B, probCSA_A, probCSA_N, probCRSA_B, probCRSA_A, probCRSA_N, ...
                transitionDest, transitions, ...
                ratioWall, ratioFloor, THRESH_WALL, THRESH_FLOOR, ...
                refQuats, refQuatsAll, ...
                chuteRoll_deg, chutePitch_deg, partName, summaryRoot, folderName);

            % ----------------------------------------------------------
            % CRITICAL VALUE METRICS FIGURES
            % ----------------------------------------------------------
            saveCriticalValueFigures( ...
                allRestingPositions, stableIdx, Qs, ...
                ratioWall, ratioFloor, ...
                THRESH_WALL, THRESH_FLOOR, ...
                partName, chuteRoll_deg, chutePitch_deg, metricsRoot, folderName);

            % ----------------------------------------------------------
            % CLEAR RUN-LEVEL VARIABLES
            % ----------------------------------------------------------
            clear allRestingPositions planeQuats stableIdx stablePositions;
            clear Qs QsRaw omegas heights crsaRaw hullProjs;
            clear ratioWall ratioFloor transitions momentArmGeo;
            clear probCSA_B probCSA_A probCSA_N probCRSA_B probCRSA_A probCRSA_N;
            clear transitionDest transitionRotatedVerts transitionRotatedCent;
            clear refQuats refQuatsAll chainLog mergedGroups;
            clear secondaryUnstable supportVerts supportProj support2D;
            clear vertsW_si centW_si;
            clear floorUnstable floorUnstableMerged instabilityFlag;
            clear wallMargin floorMargin wallDominant floorDominant;
            clear QsRaw_forA crsaRaw_forA transitionDest_forA;

        end % pitch
    end % roll

    clear partGeometry partMesh facesAndVertices;
    clear partFaces partVertices convexHullFaces chullVertexIdx;
    clear restingPlaneVertices restingPlaneEquations;
    clear centroidCoordinates;

end % STL files

disp('=================================================================');
disp('Batch run complete.');
disp('=================================================================');

end % ChutePoseAnalysis


% =========================================================================
%  CRITICAL VALUE METRICS FIGURES
% =========================================================================
function saveCriticalValueFigures( ...
    allRestingPositions, stableIdx, Qs, ...
    ratioWall, ratioFloor, ...
    THRESH_WALL, THRESH_FLOOR, ...
    partName, chuteRoll_deg, chutePitch_deg, metricsRoot, folderName)

numStable = numel(stableIdx);
poseNums  = (1:numStable)';

COL_NONZERO = [0.10, 0.65, 0.20];   % green circle  — Qs > 0
COL_ZEROED  = [0.85, 0.10, 0.10];   % red triangle  — Qs == 0

BG_GREEN = [0.72, 0.93, 0.72];      % safe region background
BG_RED   = [0.98, 0.75, 0.75];      % critical region background
BG_ALPHA = 0.55;

isNonZero = Qs > 0;

    function buildFig(ratioVec, thresh, figTitle, yLabel, pdfName)

        yMax = max([ratioVec(:); thresh * 1.35]) * 1.08 + 0.05;
        yMin = 0;
        xMin = 0.5;
        xMax = numStable + 0.5;

        fig = figure('Color','w','Visible','off', ...
            'Units','inches','Position',[0 0 10 6], ...
            'PaperUnits','inches','PaperSize',[10 6], ...
            'PaperPosition',[0 0 10 6]);

        ax = axes(fig);
        hold(ax, 'on');
        grid(ax, 'on');
        box(ax,  'on');
        ax.FontSize  = 11;
        ax.GridAlpha = 0.20;
        ax.Layer     = 'top';
        ax.XTick     = 1:numStable;

        xlim(ax, [xMin, xMax]);
        ylim(ax, [yMin, yMax]);

        patch(ax, [xMin xMax xMax xMin], ...
            [yMin yMin min(thresh, yMax) min(thresh, yMax)], ...
            BG_GREEN, 'FaceAlpha', BG_ALPHA, 'EdgeColor', 'none', ...
            'HandleVisibility', 'off');

        if thresh < yMax
            patch(ax, [xMin xMax xMax xMin], ...
                [thresh thresh yMax yMax], ...
                BG_RED, 'FaceAlpha', BG_ALPHA, 'EdgeColor', 'none', ...
                'HandleVisibility', 'off');
        end

        xline(ax, NaN);
        hl = yline(ax, thresh, '-', ...
            'Color', [0.30 0.30 0.30], 'LineWidth', 2.0, ...
            'HandleVisibility', 'off');

        axPos   = ax.Position;
        yNorm   = (thresh - yMin) / (yMax - yMin);
        textY   = yNorm + 0.03;
        textY   = min(textY, 0.97);

        annotation(fig, 'textbox', ...
            'Units', 'normalized', ...
            'Position', [axPos(1) + axPos(3)*0.70, ...
                         axPos(2) + axPos(4)*textY, ...
                         0.20, 0.06], ...
            'String', sprintf('critical = %.3f', thresh), ...
            'FontSize', 10, ...
            'FontName', 'Courier', ...
            'FontWeight', 'bold', ...
            'BackgroundColor', [1 1 1], ...
            'EdgeColor', [0.30 0.30 0.30], ...
            'LineWidth', 1.0, ...
            'Margin', 3, ...
            'FitBoxToText', 'on', ...
            'HorizontalAlignment', 'center', ...
            'Interpreter', 'none');

        hNZ = []; hZ = [];
        xNZ = poseNums(isNonZero);  yNZ = ratioVec(isNonZero);
        xZ  = poseNums(~isNonZero); yZ  = ratioVec(~isNonZero);

        if ~isempty(xNZ)
            hNZ = scatter(ax, xNZ, yNZ, 80, COL_NONZERO, 'o', 'filled', ...
                'MarkerEdgeColor', COL_NONZERO * 0.5, 'LineWidth', 1.0);
        end
        if ~isempty(xZ)
            hZ = scatter(ax, xZ, yZ, 100, COL_ZEROED, '^', 'filled', ...
                'MarkerEdgeColor', COL_ZEROED * 0.5, 'LineWidth', 1.0);
        end

        xlabel(ax, 'Pose index', 'FontSize', 12);
        ylabel(ax, yLabel,       'FontSize', 12, 'Interpreter', 'tex');
        title(ax, sprintf('%s  |  %s  |  Roll=%d° Pitch=%d°', ...
            figTitle, partName, round(chuteRoll_deg), round(chutePitch_deg)), ...
            'FontSize', 13, 'FontWeight', 'bold');

        lgdH = []; lgdL = {};
        if ~isempty(hNZ), lgdH(end+1) = hNZ; lgdL{end+1} = 'Non-zero (Qs > 0)'; end %#ok<AGROW>
        if ~isempty(hZ),  lgdH(end+1) = hZ;  lgdL{end+1} = 'Zeroed (Qs = 0)';   end %#ok<AGROW>
        if ~isempty(lgdH)
            legend(ax, lgdH, lgdL, 'Location', 'northwest', ...
                'FontSize', 10, 'Box', 'on');
        end

        pdfPath = fullfile(metricsRoot, [folderName, '_', pdfName, '.pdf']);
        exportgraphics(fig, pdfPath, 'Resolution', 200);
        fprintf('  Saved: %s\n', pdfPath);
        close(fig);
    end

buildFig(ratioWall, THRESH_WALL, ...
    'transitionAnalysis', 'l_A / l_w  (wall ratio)', ...
    'transitionAnalysis');

buildFig(ratioFloor, THRESH_FLOOR, ...
    'stabilityAnalysis', 'l_A / l_f  (floor ratio)', ...
    'stabilityAnalysis');

end % saveCriticalValueFigures


% =========================================================================
%  QUATERNION PRIMITIVES
% =========================================================================

function q = q_fromAxisAngle(axis, phi)
axis = axis(:)' / norm(axis(:));
s    = sin(phi/2);
q    = [cos(phi/2), s*axis(1), s*axis(2), s*axis(3)];
end

function q = q_compose(q1, q2)
w1=q1(1); x1=q1(2); y1=q1(3); z1=q1(4);
w2=q2(1); x2=q2(2); y2=q2(3); z2=q2(4);
q = [w2*w1 - x2*x1 - y2*y1 - z2*z1, ...
    w2*x1 + x2*w1 + y2*z1 - z2*y1, ...
    w2*y1 - x2*z1 + y2*w1 + z2*x1, ...
    w2*z1 + x2*y1 - y2*x1 + z2*w1];
end

function [Vout, cout] = q_rotateCloud(q, V, pivot, cent)
Rmat = q_toRotm(q);
p    = pivot(:)';
Vout = (Rmat * (V    - p)')' + p;
cout = (Rmat * (cent(:)' - p)')' + p;
end

function R = q_toRotm(q)
w=q(1); x=q(2); y=q(3); z=q(4);
R = [1-2*(y^2+z^2),   2*(x*y-z*w),   2*(x*z+y*w);
    2*(x*y+z*w), 1-2*(x^2+z^2),   2*(y*z-x*w);
    2*(x*z-y*w),   2*(y*z+x*w), 1-2*(x^2+y^2)];
end

function q = q_fromRotm(R)
tr = R(1,1)+R(2,2)+R(3,3);
if tr > 0
    s = 0.5/sqrt(tr+1);
    q = [0.25/s, (R(3,2)-R(2,3))*s, (R(1,3)-R(3,1))*s, (R(2,1)-R(1,2))*s];
elseif (R(1,1)>R(2,2)) && (R(1,1)>R(3,3))
    s = 2*sqrt(1+R(1,1)-R(2,2)-R(3,3));
    q = [(R(3,2)-R(2,3))/s, 0.25*s, (R(1,2)+R(2,1))/s, (R(1,3)+R(3,1))/s];
elseif R(2,2)>R(3,3)
    s = 2*sqrt(1+R(2,2)-R(1,1)-R(3,3));
    q = [(R(1,3)-R(3,1))/s, (R(1,2)+R(2,1))/s, 0.25*s, (R(2,3)+R(3,2))/s];
else
    s = 2*sqrt(1+R(3,3)-R(1,1)-R(2,2));
    q = [(R(2,1)-R(1,2))/s, (R(1,3)+R(3,1))/s, (R(2,3)+R(3,2))/s, 0.25*s];
end
q = q / norm(q);
end

function d = q_geodesic(q1, q2)
dp = abs(dot(q1(:), q2(:)));
dp = min(dp, 1.0);
d  = 2*acos(dp);
end

function q = q_poseQuat(vertsC, centC)
hIdx = p4_hullIdx(vertsC);
pts  = vertsC(hIdx, :) - centC;
if size(pts,1) < 3
    q = [1,0,0,0]; return;
end
[~,~,V] = svd(pts, 'econ');
R = V;
if det(R) < 0, R(:,3) = -R(:,3); end
q = q_fromRotm(R);
[~, mi] = max(abs(q));
if q(mi) < 0, q = -q; end
end


% =========================================================================
%  POSE ENUMERATION — HULL EDGE THETA
% =========================================================================
function [rawCandidates, planeQuats] = thetaFromHullEdges( ...
    partVertices, convexHullFaces, chullVertexIdx, ...
    restingPlaneVertices, restingPlaneEquations, ...
    centroidCoordinates, wallTol, planeTol)

rawCandidates = struct( ...
    'floorPlaneIdx',            {}, ...
    'wallSide',                 {}, ...
    'wallContactVertIdx',       {}, ...
    'floorContactVertIdx',      {}, ...
    'theta',                    {}, ...
    'verticesWorld',            {}, ...
    'centroidWorld',            {}, ...
    'floorZ',                   {}, ...
    'centroidHeightAboveFloor', {} );

numFloorPlanes = numel(restingPlaneVertices);
planeQuats = repmat({[1,0,0,0]}, 1, numFloorPlanes);

edgeTol         = wallTol;
neighbourDeg    = 5;
neighbourN      = 20;
fallbackSamples = 720;

for planeIdx = 1:numFloorPlanes

    planeEq       = restingPlaneEquations(planeIdx,:);
    planeNormal   = planeEq(1:3);
    planeNormUnit = planeNormal / norm(planeNormal);
    supportIdx    = restingPlaneVertices{planeIdx};
    if numel(supportIdx) < 2, continue; end

    faceNorm = planeNormUnit(:);
    targetN  = [0;0;-1];
    dotFT    = dot(faceNorm, targetN);

    if abs(dotFT - 1) < 1e-8
        q_align = [1,0,0,0];
    elseif abs(dotFT + 1) < 1e-8
        perp    = null(faceNorm');
        ax      = perp(:,1) / norm(perp(:,1));
        q_align = q_fromAxisAngle(ax', pi);
    else
        ax  = cross(faceNorm, targetN);
        ax  = ax / norm(ax);
        ang = acos(max(-1, min(1, dotFT)));
        q_align = q_fromAxisAngle(ax', ang);
    end

    planeQuats{planeIdx} = q_align;

    R_floor  = q_toRotm(q_align);
    C0       = centroidCoordinates(:)';
    vertsRot = (R_floor * (partVertices' - C0'))' + C0;
    centRot  = C0;

    suppRot  = vertsRot(supportIdx, :);
    floorZ   = mean(suppRot(:, 3));
    centroidHeightAboveFloor = centRot(3) - floorZ;
    if centroidHeightAboveFloor < 1e-6, continue; end

    hullVerts2D = vertsRot(chullVertexIdx, 1:2);

    try
        ord2D = convhull(hullVerts2D(:,1), hullVerts2D(:,2), 'Simplify', false);
    catch
        ord2D = [];
    end

    exactThetas = [];

    if ~isempty(ord2D)
        edgePts = hullVerts2D(ord2D(1:end-1), :);
        M = size(edgePts, 1);

        for ei = 1:M
            p0 = edgePts(ei, :);
            p1 = edgePts(mod(ei, M) + 1, :);
            dXY = p1 - p0;
            edgeLen = norm(dXY);
            if edgeLen < 1e-10, continue; end
            dXY_n = dXY / edgeLen;

            dp     = hullVerts2D - p0;
            d_perp = abs(dp(:,1) * dXY_n(2) - dp(:,2) * dXY_n(1));
            onEdge = find(d_perp <= edgeTol);
            if numel(onEdge) < 2, continue; end

            ex = dXY_n(1);
            ey = dXY_n(2);
            theta_exact = atan2(ex, ey);

            for tDelta = [0, pi]
                th = mod(theta_exact + tDelta, 2*pi);
                exactThetas(end+1) = th; %#ok<AGROW>
            end
        end
    end

    thetaCandidates = exactThetas(:)';

    dNeigh = deg2rad(neighbourDeg);
    for et = exactThetas
        neigh = et + linspace(-dNeigh, dNeigh, 2*neighbourN+1);
        thetaCandidates = [thetaCandidates, neigh]; %#ok<AGROW>
    end

    thetaFallback   = linspace(0, 2*pi, fallbackSamples+1);
    thetaFallback   = thetaFallback(1:end-1);
    thetaCandidates = [thetaCandidates, thetaFallback];

    thetaCandidates = mod(thetaCandidates, 2*pi);
    thetaCandidates = sort(unique(round(thetaCandidates / deg2rad(0.01)) * deg2rad(0.01)));

    thetaWallKeys = {};

    for ti = 1:numel(thetaCandidates)
        theta = thetaCandidates(ti);

        allVertsWorld = rotatePtsAroundZ(vertsRot, centRot, theta);

        floorZcurrent = min(allVertsWorld(:,3));
        floorResid    = abs(allVertsWorld(:,3) - floorZcurrent);
        onFloor       = find(floorResid < planeTol);
        if numel(onFloor) < 2, continue; end

        floorPts2D = allVertsWorld(onFloor, 1:2);
        centXY     = [centRot(1), centRot(2)];

        if numel(onFloor) >= 3
            try
                hIdx_    = convhull(floorPts2D(:,1), floorPts2D(:,2));
                hPoly    = floorPts2D(hIdx_, :);
                floorOK  = inpolygon(centXY(1), centXY(2), hPoly(:,1), hPoly(:,2));
            catch
                floorOK = false;
            end
        else
            floorOK = (centXY(1) >= min(floorPts2D(:,1)) - wallTol) && ...
                (centXY(1) <= max(floorPts2D(:,1)) + wallTol);
        end
        if ~floorOK, continue; end

        projPts      = allVertsWorld(chullVertexIdx, 1:2);
        yWall        = max(projPts(:, 2));
        onWall_local = find(abs(projPts(:, 2) - yWall) < wallTol);

        if numel(onWall_local) >= 2
            onWall      = chullVertexIdx(onWall_local);
            wallPtsXY   = allVertsWorld(onWall, 1:2);
            wallPtsXY_u = uniquetol(wallPtsXY, wallTol, 'ByRows', true);

            if size(wallPtsXY_u, 1) >= 2
                wallPtsX = wallPtsXY_u(:, 1);
                centX    = centRot(1);

                if centX >= min(wallPtsX) - wallTol && ...
                        centX <= max(wallPtsX) + wallTol

                    wc      = sort(onWall(:)');
                    fc      = sort(onFloor(:)');
                    topoKey = sprintf('W_%s__F_%s', mat2str(wc), mat2str(fc));

                    if ~any(strcmp(thetaWallKeys, topoKey))
                        thetaWallKeys{end+1} = topoKey; %#ok<AGROW>

                        pos.floorPlaneIdx            = planeIdx;
                        pos.wallSide                 = 1;
                        pos.wallContactVertIdx       = wc;
                        pos.floorContactVertIdx      = fc;
                        pos.theta                    = theta;
                        pos.verticesWorld            = allVertsWorld;
                        pos.centroidWorld            = centRot;
                        pos.floorZ                   = floorZ;
                        pos.centroidHeightAboveFloor = centroidHeightAboveFloor;

                        rawCandidates(end+1) = pos; %#ok<AGROW>
                    end
                end
            end
        end

    end % theta candidates

    clear vertsRot thetaCandidates thetaWallKeys exactThetas;

end % planeIdx
end % thetaFromHullEdges


% =========================================================================
%  PROBABILITY-ALGORITHM WEIGHT HELPERS
%  (shared normalization / filtering logic used identically by both the
%   CSA weight source (QsRaw) and the CRSA weight source (crsaRaw), so that
%   CSA_B/CSA_A/CSA_N and CRSA_B/CRSA_A/CRSA_N are computed with exactly
%   the same math, just fed a different raw-weight vector. Both QsRaw and
%   crsaRaw are UNTOUCHED by the secondary-stability / floor-instability
%   validity checks — that is what lets the "_B"/"_A" filtering actually
%   change something relative to "_N" for both algorithms alike.)
% =========================================================================
function probOut = normalizeWeights(weights)
% "_N" (none) — no filtering at all, just normalize the raw weight.
n = numel(weights);
s = sum(weights);
if s > 0
    probOut = 100 * weights / s;
else
    probOut = zeros(n, 1);
end
end

function probOut = zeroAndNormalize(weights, zeroMask)
% "_B" (before) — zero out transitioning poses, then renormalize the
% remaining weight (removed mass is redistributed proportionally across
% all surviving poses, not sent to any specific destination pose).
w = weights;
w(zeroMask) = 0;
probOut = normalizeWeights(w);
end

function probOut = redistributeWeights(weights, transitionDest)
% "_A" (after) — zero out transitioning poses, but hand their weight to
% the specific destination pose reached by following the physical
% tip-over chain (transitionDest), then renormalize.
n = numel(weights);
wOut = weights;
for si = 1:n
    d = transitionDest(si);
    if d > 0 && d ~= si
        wOut(d)  = wOut(d) + weights(si);
        wOut(si) = 0;
    end
end
probOut = normalizeWeights(wOut);
end


% =========================================================================
%  TRANSITION CHAINS  (pure geometry — weight-agnostic)
%
%  Determines, for every transitioning pose, which destination pose it
%  physically tips over into (transitionDest), by following the wall-tip
%  rotation chain and matching quaternions against every stable pose's
%  reference orientation. This geometric result does not depend on which
%  raw weight vector (QsRaw or crsaRaw) will later be redistributed using it
%  — redistributeWeights() applies it afterward to whichever weight source
%  is being evaluated.
% =========================================================================
function [transitionDest, transitionRotatedVerts, transitionRotatedCent, ...
    refQuats, refQuatsAll, chainLog] = ...
    computeTransitionChains( ...
    allRestingPositions, stableIdx, transitions, ...
    ratioWall, threshWall, ...
    Rchute, q_chute, quatMatchTol, centroidCoordinates, ...
    planeQuats, mergedGroups) %#ok<INUSL,INUSD>

numStable = numel(stableIdx);

transitionDest         = zeros(numStable, 1);
transitionRotatedVerts = cell(numStable, 1);
transitionRotatedCent  = cell(numStable, 1);
chainLog               = cell(numStable, 1);

% ------------------------------------------------------------------
% BUILD refQuats AND refQuatsAll
% ------------------------------------------------------------------
refQuats    = zeros(numStable, 4);
refQuatsAll = cell(numStable, 1);

for si = 1:numStable

    pos      = allRestingPositions(stableIdx(si));
    q_align  = planeQuats{pos.floorPlaneIdx};
    q_theta  = q_fromAxisAngle([0,0,1], pos.theta);
    q_ref    = q_compose(q_align, q_theta);
    [~, mi]  = max(abs(q_ref));
    if q_ref(mi) < 0, q_ref = -q_ref; end
    refQuats(si,:) = q_ref;

    if nargin >= 11 && ~isempty(mergedGroups) && numel(mergedGroups) >= si
        grpIdx = mergedGroups{si};
    else
        grpIdx = stableIdx(si);
    end
    grpIdx = grpIdx(:)';

    qList = zeros(numel(grpIdx), 4);
    for gk = 1:numel(grpIdx)
        pos_g   = allRestingPositions(grpIdx(gk));
        qa      = planeQuats{pos_g.floorPlaneIdx};
        qt      = q_fromAxisAngle([0,0,1], pos_g.theta);
        qg      = q_compose(qa, qt);
        [~, mi] = max(abs(qg));
        if qg(mi) < 0, qg = -qg; end
        qList(gk,:) = qg;
    end

    if isempty(qList)
        qList = q_ref;
    end
    refQuatsAll{si} = qList;
end

maxIter = 5;

for si = 1:numStable
    if ~transitions(si), continue; end

    log = {};

    [vertsC, centC] = p4_initPose(si, allRestingPositions, stableIdx, Rchute);

    pos0      = allRestingPositions(stableIdx(si));
    floorIdx0 = pos0.floorContactVertIdx(:);
    wallIdx0  = pos0.wallContactVertIdx(:);

    % Wall-only transitions
    transType = 'wall';

    firstRotDone = false;
    vertsC_cur   = vertsC;
    centC_cur    = centC;
    wallIdx_cur  = wallIdx0;
    floorIdx_cur = floorIdx0;
    destSi       = 0;

    q_accumulated = refQuats(si,:);

    for iter = 1:maxIter

        logEntry = struct();
        logEntry.iter      = iter;
        logEntry.transType = transType;

        [vertsC_rot, centC_rot, phi_q, pivotUsed, axisUsed, q_rot] = ...
            p4_wallRotation_q(vertsC_cur, centC_cur, wallIdx_cur);

        logEntry.phi_deg   = rad2deg(phi_q);
        logEntry.pivot     = pivotUsed;
        logEntry.axis      = axisUsed;

        if ~isempty(q_rot)
            q_accumulated = q_compose(q_accumulated, q_rot);
            [~, mi] = max(abs(q_accumulated));
            if q_accumulated(mi) < 0, q_accumulated = -q_accumulated; end
        end

        logEntry.q_accumulated = q_accumulated;

        if isempty(vertsC_rot)
            [vertsC_zs, centC_zs, ~, ~] = p4_reseat(vertsC_cur, centC_cur);
            logEntry.vertsC_rot = vertsC_zs;
            logEntry.centC_rot  = centC_zs;

            if ~firstRotDone
                transitionRotatedVerts{si} = vertsC_zs;
                transitionRotatedCent{si}  = centC_zs;
                firstRotDone = true;
            end

            [matchSi_zero, matchDist_zero] = p4_matchQuatComposed(q_accumulated, refQuatsAll, quatMatchTol);
            logEntry.matchSi   = matchSi_zero;
            logEntry.matchDist = matchDist_zero;

            if matchSi_zero > 0
                destSi = matchSi_zero;
            end
            log{end+1} = logEntry; %#ok<AGROW>
            break;
        end

        if ~firstRotDone
            transitionRotatedVerts{si} = vertsC_rot;
            transitionRotatedCent{si}  = centC_rot;
            firstRotDone = true;
        end

        logEntry.vertsC_rot = vertsC_rot;
        logEntry.centC_rot  = centC_rot;

        [vertsC_seat, centC_seat, newFloorIdx, newWallIdx] = ...
            p4_reseat(vertsC_rot, centC_rot);

        [matchSi, matchDist] = p4_matchQuatComposed(q_accumulated, refQuatsAll, quatMatchTol);
        logEntry.matchSi   = matchSi;
        logEntry.matchDist = matchDist;

        if matchSi > 0
            transitionRotatedVerts{si} = vertsC_seat;
            transitionRotatedCent{si}  = centC_seat;
            destSi = matchSi;
            log{end+1} = logEntry; %#ok<AGROW>
            break;
        end

        isStable = false;
        if ~isempty(newWallIdx)
            wallXvals = vertsC_seat(newWallIdx, 1);
            centX     = centC_seat(1);
            isStable  = (centX >= min(wallXvals) - 1e-6) && ...
                (centX <= max(wallXvals) + 1e-6);
        else
            isStable = true;
        end

        logEntry.isStable = isStable;
        log{end+1} = logEntry; %#ok<AGROW>

        if isStable
            transitionRotatedVerts{si} = vertsC_seat;
            transitionRotatedCent{si}  = centC_seat;
            break;
        end

        vertsC_cur   = vertsC_seat;
        centC_cur    = centC_seat;
        wallIdx_cur  = newWallIdx;
        floorIdx_cur = newFloorIdx;

        clear vertsC_rot centC_rot vertsC_seat centC_seat;

    end % iter

    transitionDest(si) = destSi;
    chainLog{si}       = log;

    clear vertsC centC vertsC_cur centC_cur log;

end % si

end % computeTransitionChains


% =========================================================================
%  p4_matchQuatComposed
% =========================================================================
function [matchSi, bestDist] = p4_matchQuatComposed(q_query, refQuatsAll, quatMatchTol)

matchSi  = 0;
bestDist = inf;
bestSj   = 0;

for sj = 1:numel(refQuatsAll)
    qList = refQuatsAll{sj};
    if isempty(qList), continue; end
    for kk = 1:size(qList, 1)
        d = q_geodesic(q_query, qList(kk, :));
        if d < bestDist
            bestDist = d;
            bestSj   = sj;
        end
    end
end

if bestDist <= quatMatchTol
    matchSi = bestSj;
end

end


% =========================================================================
%  p4_wallRotation_q
% =========================================================================
function [vertsR, centR, phi, pivotOut, axisOut, q_rot] = ...
    p4_wallRotation_q(vertsC, centC, wallIdx)

vertsR   = [];
centR    = [];
phi      = 0;
pivotOut = [0,0,0];
axisOut  = [0,0,1];
q_rot    = [];

if isempty(wallIdx), return; end

wallVertsC = vertsC(wallIdx, :);
[~, pivLoc] = max(wallVertsC(:, 1));
pivot = wallVertsC(pivLoc, :);
pivotOut = pivot;

hullIdx3D = p4_hullIdx(vertsC);
if numel(hullIdx3D) < 3, return; end

mergeTol = max(range(vertsC(:,1)), range(vertsC(:,2))) * 1e-3;
hullXY_raw = vertsC(hullIdx3D, 1:2);
hullXY = uniquetol(hullXY_raw, mergeTol, 'ByRows', true, 'DataScale', 1);

if size(hullXY, 1) < 3, return; end

try
    ord2D = convhull(hullXY(:,1), hullXY(:,2));
catch
    return;
end

ord2D = ord2D(1:end-1);
pts2D = hullXY(ord2D, :);
M     = size(pts2D, 1);
if M < 2, return; end

pivXY = pivot(1:2);
dists = vecnorm(pts2D - pivXY, 2, 2);
[~, pivLoc2D] = min(dists);

prevLoc = mod(pivLoc2D - 2, M) + 1;
nextLoc = mod(pivLoc2D,     M) + 1;
prevXY  = pts2D(prevLoc, :);
nextXY  = pts2D(nextLoc, :);

if nextXY(1) >= prevXY(1)
    neighXY = nextXY;
else
    neighXY = prevXY;
end

dX  = neighXY(1) - pivXY(1);
dY  = neighXY(2) - pivXY(2);
phi = atan2(-dY, dX);

if abs(phi) < 1e-8, return; end

q_rot  = q_fromAxisAngle([0,0,1], phi);
axisOut = [0,0,1];

[vertsR, centR] = q_rotateCloud(q_rot, vertsC, pivot, centC);

end


% =========================================================================
%  p4_floorRotation_q
% =========================================================================
function [vertsR, centR, phi, pivotOut, axisOut, q_rot] = ...
    p4_floorRotation_q(vertsC, centC, floorIdx, ...
    allRestingPositions, stableIdx, srcSi, Rchute)

vertsR   = [];
centR    = [];
phi      = 0;
pivotOut = [0,0,0];
axisOut  = [0,-1,0];
q_rot    = [];

numStable = numel(stableIdx);

floorVertsC = vertsC(floorIdx, :);
[~, pivLoc] = max(floorVertsC(:, 1));
pivot = floorVertsC(pivLoc, :);
pivotOut = pivot;

bestPhi = Inf;

for sj = 1:numStable
    if sj == srcSi, continue; end

    pos_sj = allRestingPositions(stableIdx(sj));
    if isempty(pos_sj.floorContactVertIdx), continue; end

    [vertsJ, ~] = p4_initPose(sj, allRestingPositions, stableIdx, Rchute);
    floorIdxJ   = pos_sj.floorContactVertIdx(:);
    [nJ, ~]     = p4_floorNormal(vertsJ, floorIdxJ);

    phi_j = atan2(-nJ(1), nJ(3));

    if phi_j > 1e-6 && phi_j < bestPhi
        bestPhi = phi_j;
    end

    clear vertsJ;
end

if isinf(bestPhi), return; end

phi = bestPhi;

q_rot   = q_fromAxisAngle([0,1,0], -phi);
axisOut = [0,-1,0];

[vertsR, centR] = q_rotateCloud(q_rot, vertsC, pivot, centC);

end


% =========================================================================
%  p4_initPose
% =========================================================================
function [vertsW, centW] = p4_initPose(si, allRestingPositions, stableIdx, Rchute) %#ok<INUSD>
pos    = allRestingPositions(stableIdx(si));
vertsW = pos.verticesWorld;
centW  = pos.centroidWorld(:)';

floorIdx = pos.floorContactVertIdx(:);
wallIdx  = pos.wallContactVertIdx(:);

floorZ        = min(vertsW(floorIdx, 3));
vertsW(:,3)   = vertsW(:,3) - floorZ;
centW(3)      = centW(3)    - floorZ;

wallY         = max(vertsW(wallIdx, 2));
vertsW(:,2)   = vertsW(:,2) - wallY;
centW(2)      = centW(2)    - wallY;
end


% =========================================================================
%  p4_reseat
% =========================================================================
function [vertsOut, centOut, floorIdx, wallIdx] = p4_reseat(vertsC, centC)

span       = max(max(vertsC,[],1) - min(vertsC,[],1));
contactTol = max(0.05, span * 0.01);

hullIdx = p4_hullIdx(vertsC);

zVals       = vertsC(hullIdx, 3);
zMin        = min(zVals);
vertsC(:,3) = vertsC(:,3) - zMin;
centC(3)    = centC(3)    - zMin;

yVals = vertsC(hullIdx, 2);
yMax  = max(yVals);
if ~isnan(yMax)
    vertsC(:,2) = vertsC(:,2) - yMax;
    centC(2)    = centC(2)    - yMax;
end

zVals2   = vertsC(hullIdx, 3);
yVals2   = vertsC(hullIdx, 2);
floorIdx = hullIdx(abs(zVals2 - 0) <= contactTol);
wallIdx  = hullIdx(abs(yVals2 - 0) <= contactTol);

vertsOut = vertsC;
centOut  = centC;

end


% =========================================================================
%  p4_floorNormal, p4_hullIdx
% =========================================================================
function [n, meanPt] = p4_floorNormal(vertsC, floorIdx)
fv     = vertsC(floorIdx(:), :);
meanPt = mean(fv, 1);
if size(fv, 1) >= 3
    [~, ~, V] = svd(fv - meanPt);
    n = V(:,3)';
else
    n = [0, 0, 1];
end
if n(3) < 0, n = -n; end
n = n / norm(n);
end

function idx = p4_hullIdx(vertsC)
try
    hf  = convhull(vertsC, 'Simplify', true);
    idx = unique(hf(:));
catch
    idx = (1:size(vertsC,1))';
end
end


% =========================================================================
%  CRSA RAW SCORES  (v3 — dual floor/wall hull, blocked-edge filtering)
%
%  Implements: P_i ∝ sum_k sum_j (1/m_i) * (Q_i - Q'_ij) / h_ti   [Eq. 6]
%
%  Q_i    = solid angle from centroid over the COMBINED floor+wall
%           support hull (unchanged from before).
%  Q'_ij  = solid angle from a CRITICAL point over that same hull.
%           The critical point keeps the centroid's original in-plane
%           (X,Y) or (X,Z) coordinates and replaces only the coordinate
%           normal to whichever surface (floor or wall) the transition
%           edge belongs to, with R = perpendicular distance from the
%           centroid to that edge's axis LINE.
%  m_i    = count of FREE (non-blocked) transition edges across both
%           the floor hull and the wall hull.
%
%  Local frame convention (matches p4_initPose / p4_reseat):
%     floor at Z=0, wall at Y=0, body occupies Z>=0, Y<=0.
%     n_f = [0,0,1]   (floor -> into body)
%     n_w = [0,-1,0]  (wall  -> into body)
%
%  Edge classification:
%     For each edge of a hull, the OUTWARD in-plane normal is the one
%     pointing away from that hull's own 2D centroid (standard convex
%     polygon convention == "rotate away from the middle of the part")
%     An edge is BLOCKED (excluded from the sum and from m_i) if its
%     outward normal points toward the OTHER support surface, since
%     tipping that way immediately re-engages that surface rather than
%     freeing the part (would require infinite energy to actually tip).
%
%  Degenerate 2-point hulls: project both points onto the floor/wall
%  intersection line (the local X-axis) and add those as two more
%  points before hulling, so a 2-point "line" becomes a proper quad
%  bounded by the floor-wall seam.
% =========================================================================
function crsaRaw = computeCRSA( ...
    allRestingPositions, stableIdx, Rchute, gWorld, wallNorm_c) %#ok<INUSL>

numStable = numel(stableIdx);
crsaRaw   = zeros(numStable, 1);
gWorld    = gWorld(:) / norm(gWorld(:));

NEG_CLAMP_TOL = -1e-6;   % numerical-noise guard only, not a hard max(0,...)
BLOCK_TOL     = 1e-6;    % dot-product threshold for "points toward other surface"

for si = 1:numStable
    pi_  = stableIdx(si);
    pos  = allRestingPositions(pi_);

    % ------------------------------------------------------------------
    % Build the LOCAL frame: floor at Z=0, wall at Y=0, body Z>=0,Y<=0
    % (same convention as p4_initPose)
    % ------------------------------------------------------------------
    vertsC = pos.verticesWorld;
    centC  = pos.centroidWorld(:)';

    floorContactIdx = pos.floorContactVertIdx(:);
    wallContactIdx  = pos.wallContactVertIdx(:);

    if isempty(floorContactIdx) || size(vertsC,1) < 3
        continue;
    end

    floorZ       = min(vertsC(floorContactIdx, 3));
    vertsC(:,3)  = vertsC(:,3) - floorZ;
    centC(3)     = centC(3)    - floorZ;

    if ~isempty(wallContactIdx)
        wallY        = max(vertsC(wallContactIdx, 2));
        vertsC(:,2)  = vertsC(:,2) - wallY;
        centC(2)     = centC(2)    - wallY;
    end

    n_f = [0, 0, 1];
    n_w = [0, -1, 0];

    % ------------------------------------------------------------------
    % Combined support hull (for Q_i) — SAME as before, world/Rchute frame
    % ------------------------------------------------------------------
    C0     = pos.centroidWorld(:)';
    vertsW = (Rchute * (pos.verticesWorld' - C0'))' + C0;
    centW  = C0;

    supportIdx   = union(floorContactIdx, wallContactIdx);
    supportVerts = vertsW(supportIdx, :);
    if size(supportVerts,1) < 3, continue; end

    % ------------------------------------------------------------------
    % HEIGHT (h): dual floor/wall ray-cast along gravity — IDENTICAL to
    % computeCSA's height logic. h = min(t_floor, t_wall): whichever
    % support plane the centroid's gravity ray reaches first. This is
    % copied verbatim from computeCSA (same SVD-normal extraction, same
    % denom/t_floor/t_wall guards, no normal re-orientation flip — CSA
    % doesn't flip its floor normal either) so CRSA's denominator can
    % never diverge from CSA's reported heights(si) for the same pose.
    % ------------------------------------------------------------------
    rayDir = [0,0,-1];

    t_floor = Inf;
    if numel(floorContactIdx) >= 3
        floorVertsW = vertsW(floorContactIdx, :);
        floorMeanW  = mean(floorVertsW, 1);
        [~,~,Vf]    = svd(floorVertsW - floorMeanW);
        floorNormalW = Vf(:,3)';
        floorNormalW = floorNormalW / norm(floorNormalW);

        denomF = dot(rayDir, floorNormalW);
        if abs(denomF) > 1e-10
            tF = dot(floorMeanW - centW, floorNormalW) / denomF;
            if tF >= 0
                t_floor = tF;
            end
        end
    end

    t_wall = Inf;
    if numel(wallContactIdx) >= 3
        wallVertsW = vertsW(wallContactIdx, :);
        wallMeanW  = mean(wallVertsW, 1);
        [~,~,Vw]   = svd(wallVertsW - wallMeanW);
        wallNormalW = Vw(:,3)';
        wallNormalW = wallNormalW / norm(wallNormalW);

        denomW = dot(rayDir, wallNormalW);
        if abs(denomW) > 1e-10
            tW = dot(wallMeanW - centW, wallNormalW) / denomW;
            if tW >= 0
                t_wall = tW;
            end
        end
    end

    if isinf(t_floor) && isinf(t_wall), continue; end
    h = min(t_floor, t_wall);
    if h < 1e-10, continue; end

    % ------------------------------------------------------------------
    % (unchanged) support-hull solid angle Q_i
    % ------------------------------------------------------------------
    gravCoords  = supportVerts * gWorld;
    minGrav     = min(gravCoords);
    planePoint  = minGrav * gWorld';
    supportProj = zeros(size(supportVerts));
    for kk = 1:size(supportVerts,1)
        p_  = supportVerts(kk,:);
        d_  = dot((p_ - planePoint), gWorld');
        supportProj(kk,:) = p_ - d_ * gWorld';
    end
    supportProj = uniquetol(supportProj, 1e-5, 'ByRows', true);
    if size(supportProj,1) < 3, continue; end

    if abs(gWorld(1)) < 0.9, tmp2 = [1;0;0]; else, tmp2 = [0;1;0]; end
    e1b = cross(gWorld, tmp2); e1b = e1b/norm(e1b);
    e2b = cross(gWorld, e1b);  e2b = e2b/norm(e2b);
    support2D = [supportProj*e1b, supportProj*e2b];

    try
        hullIdx_c = convhull(support2D(:,1), support2D(:,2));
    catch
        continue
    end

    hull3D = supportProj(hullIdx_c, :);   % closed loop, last==first
    Nsup   = size(hull3D,1) - 1;
    if Nsup < 3, continue; end

    Q_i = solidAnglePolygon(hull3D, centW);
    if Q_i <= 0, continue; end

    spanScale    = max([range(vertsC(:,1)), range(vertsC(:,2)), range(vertsC(:,3))]);
    edgePointTol = max(spanScale * 1e-4, 1e-6);
    
    [floorEdgesFree, floorCritLocal] = buildSubHullEdges( ...
        vertsC, floorContactIdx, [1 2], 3, n_f, centC, BLOCK_TOL, true,  edgePointTol);
    
    [wallEdgesFree, wallCritLocal] = buildSubHullEdges( ...
        vertsC, wallContactIdx, [1 3], 2, n_w, centC, BLOCK_TOL, false, edgePointTol);

    % ------------------------------------------------------------------
    % Accumulate deltaSum over all FREE edges (floor + wall)
    % ------------------------------------------------------------------
    deltaSum = 0;
    mCount   = 0;

    allCritLocal = [floorCritLocal; wallCritLocal];
    allFree      = [floorEdgesFree; wallEdgesFree];

    for kk = 1:size(allCritLocal,1)
        if ~allFree(kk), continue; end

        critLocal = allCritLocal(kk,:);

        critWorld = localToWorldCrit(critLocal, pos, floorContactIdx, ...
            wallContactIdx, Rchute);

        Qcrit = solidAnglePolygon(hull3D, critWorld);
        delta = Q_i - Qcrit;
        if delta < NEG_CLAMP_TOL
            delta = 0;
        else
            delta = max(0, delta);
        end

        deltaSum = deltaSum + delta;
        mCount   = mCount + 1;
    end

    if mCount == 0, continue; end

    crsaRaw(si) = deltaSum / (mCount * h);
end
end % computeCRSA


% =========================================================================
%  buildSubHullEdges
%
%  Builds the 2D convex hull for one support surface (floor or wall) in
%  the LOCAL frame, classifies each edge as free/blocked, and returns
%  the critical point (in LOCAL 3D coords) for each edge in hull order.
% =========================================================================
function [isFree, critLocal] = buildSubHullEdges( ...
    vertsC, contactIdx, planeIdx2D, normalIdx3, n_self, centC, blockTol, isFloor, pointTol)

isFree    = false(0,1);
critLocal = zeros(0,3);

if nargin < 9 || isempty(pointTol)
    pointTol = 1e-6;
end

if numel(contactIdx) < 2
    return;
end

pts3 = vertsC(contactIdx, :);

if numel(contactIdx) == 2
    seamPts = pts3;
    seamPts(:,2) = 0;
    seamPts(:,3) = 0;
    pts3 = [pts3; seamPts];
end

pts2D = pts3(:, planeIdx2D);
pts2D = uniquetol(pts2D, pointTol, 'ByRows', true, 'DataScale', 1);
if size(pts2D,1) < 3
    return;
end

try
    ord = convhull(pts2D(:,1), pts2D(:,2));
catch
    return;
end
ord    = ord(1:end-1);
poly2D = pts2D(ord, :);

% NEW: collapse duplicate/near-collinear hull edges into one before
% classifying/critical-pointing them.
poly2D = mergeCollinearEdges(poly2D, 1e-3);

M = size(poly2D,1);
if M < 3
    return;
end

hullCentroid2D = mean(poly2D, 1);

isFree    = false(M,1);
critLocal = zeros(M,3);

for ei = 1:M
    % ... rest of the function body is UNCHANGED from here down ...
    p0_2D = poly2D(ei, :);
    p1_2D = poly2D(mod(ei,M)+1, :);
    tangent2D = p1_2D - p0_2D;
    tl = norm(tangent2D);
    if tl < 1e-10
        continue;
    end
    tangent2D = tangent2D / tl;

    % outward in-plane normal: perpendicular to tangent, pointing
    % away from this hull's own 2D centroid ("rotate away from the
    % middle of the part")
    perp2D = [-tangent2D(2), tangent2D(1)];
    edgeMid2D = 0.5*(p0_2D + p1_2D);
    if dot(perp2D, edgeMid2D - hullCentroid2D) < 0
        perp2D = -perp2D;
    end

    % embed outward normal + edge endpoints back into local 3D
    p0_3 = zeros(1,3); p0_3(planeIdx2D) = p0_2D; p0_3(normalIdx3) = 0;
    p1_3 = zeros(1,3); p1_3(planeIdx2D) = p1_2D; p1_3(normalIdx3) = 0;
    outward3 = zeros(1,3); outward3(planeIdx2D) = perp2D;

    % ------------------------------------------------------------
    % Blocked test: does the outward normal point toward the OTHER
    % surface?
    %   floor edge blocked if outward3 . [0,1,0]  > blockTol  (toward wall)
    %   wall  edge blocked if outward3 . [0,0,-1] > blockTol  (toward floor)
    % ------------------------------------------------------------
    if isFloor
        blocked = dot(outward3, [0,1,0]) > blockTol;
    else
        blocked = dot(outward3, [0,0,-1]) > blockTol;
    end

    if blocked
        isFree(ei) = false;
        continue;
    end
    isFree(ei) = true;

    % ------------------------------------------------------------
    % R = perpendicular distance from centroid to the edge AXIS LINE
    % ------------------------------------------------------------
    edgeTangent3 = p1_3 - p0_3;
    etl = norm(edgeTangent3);
    if etl < 1e-10
        isFree(ei) = false;
        continue;
    end
    edgeTangent3 = edgeTangent3 / etl;

    t_foot = dot(centC - p0_3, edgeTangent3);
    footPt = p0_3 + t_foot * edgeTangent3;
    R = norm(centC - footPt);

    % critPt: same in-plane (X,Y) or (X,Z) as centroid, only the
    % surface-normal coordinate replaced with R (signed along n_self)
    critPt = footPt;
    critPt(normalIdx3) = footPt(normalIdx3) + R * n_self(normalIdx3);
    if n_self(normalIdx3) == 0
        critPt(normalIdx3) = R;
    end

    critLocal(ei, :) = critPt;
end
end % buildSubHullEdges

function polyOut = mergeCollinearEdges(polyIn, angTol)
% Drop hull vertices whose two adjacent edges are (numerically) parallel,
% i.e. merge two "edges" that are really one straight edge split by an
% extra near-collinear vertex (from the 2-point seam augmentation, or
% near-duplicate contact points that survived point-level dedup).
M = size(polyIn, 1);
if M <= 3
    polyOut = polyIn;
    return;
end

poly    = polyIn;
changed = true;
guard   = 0;

while changed && guard < 10
    changed = false;
    guard   = guard + 1;
    Mc = size(poly,1);
    if Mc <= 3, break; end
    keep = true(Mc,1);

    for i = 1:Mc
        pPrev = poly(mod(i-2,Mc)+1, :);
        pCur  = poly(i, :);
        pNext = poly(mod(i,Mc)+1, :);

        v1 = pCur  - pPrev;
        v2 = pNext - pCur;
        l1 = norm(v1); l2 = norm(v2);

        if l1 < 1e-12 || l2 < 1e-12
            keep(i) = false;   % coincident point
            continue;
        end

        v1n = v1/l1; v2n = v2/l2;
        crossMag = abs(v1n(1)*v2n(2) - v1n(2)*v2n(1));

        if crossMag < angTol
            keep(i) = false;   % pCur lies on the line through its neighbors
        end
    end

    if sum(keep) < 3, break; end
    if any(~keep)
        poly = poly(keep, :);
        changed = true;
    end
end

polyOut = poly;
end


% =========================================================================
%  localToWorldCrit
%
%  Maps a critical point expressed in the LOCAL frame (floor at Z=0,
%  wall at Y=0, offsets subtracted from pos.verticesWorld/centroidWorld)
%  back into the WORLD/Rchute-rotated frame used for hull3D/centW, by
%  re-applying the same floor/wall offset subtraction used to build the
%  local frame, then the same Rchute rotation used to build vertsW/centW.
% =========================================================================
function ptWorld = localToWorldCrit(ptLocal, pos, floorContactIdx, wallContactIdx, Rchute)

vertsC0 = pos.verticesWorld;

floorZ  = min(vertsC0(floorContactIdx, 3));
offsetZ = floorZ;

if ~isempty(wallContactIdx)
    wallY   = max(vertsC0(wallContactIdx, 2));
    offsetY = wallY;
else
    offsetY = 0;
end

ptPreOffset = ptLocal;
ptPreOffset(3) = ptPreOffset(3) + offsetZ;
ptPreOffset(2) = ptPreOffset(2) + offsetY;

C0 = pos.centroidWorld(:)';
ptWorld = (Rchute * (ptPreOffset(:) - C0(:)))' + C0;

end % localToWorldCrit


% =========================================================================
%  TRANSITION MOMENT ARM RATIOS
%  Returns ratioWall (l_A/l_w) and ratioFloor (l_A/l_f) for every stable
%  pose. transitions output is placeholder only — actual transition logic
%  lives in the main loop using the hardcoded thresholds.
% =========================================================================
function [ratioWall, ratioFloor, transitions, momentArmGeo] = computeTransitionRatios( ...
    allRestingPositions, stableIdx, Rchute, ...
    slideDir, floorNorm_c, wallNorm_c, planeTol) %#ok<INUSL>

numStable   = numel(stableIdx);
ratioWall   = zeros(numStable,1);
ratioFloor  = zeros(numStable,1);
transitions = false(numStable,1);   % not used directly; kept for signature compat

emptyGeo = struct( ...
    'pivotWall_c',   [],  ...
    'pivotFloor_c',  [],  ...
    'centC',         [],  ...
    'l_A_wall',      0,   ...
    'l_w',           0,   ...
    'l_A_floor',     0,   ...
    'l_f',           0    );
momentArmGeo = repmat(emptyGeo, numStable, 1);

for si = 1:numStable
    pi_  = stableIdx(si);
    pos  = allRestingPositions(pi_);

    vertsC          = pos.verticesWorld;
    centC           = pos.centroidWorld(:)';
    floorContactIdx = pos.floorContactVertIdx;
    wallContactIdx  = pos.wallContactVertIdx;

    floorZ       = min(vertsC(floorContactIdx, 3));
    vertsC(:,3)  = vertsC(:,3) - floorZ;
    centC(3)     = centC(3)    - floorZ;

    wallY        = max(vertsC(wallContactIdx, 2));
    vertsC(:,2)  = vertsC(:,2) - wallY;
    centC(2)     = centC(2)    - wallY;

    momentArmGeo(si).centC = centC;

    floorVertsC = vertsC(floorContactIdx, :);
    wallVertsC  = vertsC(wallContactIdx,  :);

    if ~isempty(wallVertsC)
        [~, idxW]    = max(wallVertsC(:, 1));
        pivotWall_c  = wallVertsC(idxW, :);

        l_A_wall = abs(centC(2));
        l_w      = abs(centC(1) - pivotWall_c(1));

        momentArmGeo(si).pivotWall_c = pivotWall_c;
        momentArmGeo(si).l_A_wall    = l_A_wall;
        momentArmGeo(si).l_w         = l_w;

        if l_w > planeTol
            ratioWall(si) = l_A_wall / l_w;
        end
    end

    if ~isempty(floorVertsC)
        [~, idxF]    = max(floorVertsC(:, 1));
        pivotFloor_c = floorVertsC(idxF, :);

        l_A_floor = abs(centC(3));
        l_f       = abs(centC(1) - pivotFloor_c(1));

        momentArmGeo(si).pivotFloor_c = pivotFloor_c;
        momentArmGeo(si).l_A_floor    = l_A_floor;
        momentArmGeo(si).l_f          = l_f;

        if l_f > planeTol
            ratioFloor(si) = l_A_floor / l_f;
        end
    end
end
end


% =========================================================================
%  SAVE FIGURES AND PDF
% =========================================================================
function saveFiguresAndPDF( ...
    allRestingPositions, stableIdx, ...
    partVertices, partFaces, convexHullFaces, chullVertexIdx, ...
    Rchute, Qs, omegas, heights, sumQ, ...
    probCSA_B, probCSA_A, probCSA_N, probCRSA_B, probCRSA_A, probCRSA_N, ...
    transitionDest, transitions, ...
    chuteRoll_deg, chutePitch_deg, ...
    partName, pdfPath, ...
    slideDir, floorNorm_c, wallNorm_c, ...
    ratioWall, ratioFloor, threshWall, threshFloor, planeTol, ...
    transitionRotatedVerts, transitionRotatedCent, ...
    refQuats, refQuatsAll, quatMatchTol, chainLog, momentArmGeo) %#ok<INUSD>

numStable  = numel(stableIdx);
partSpan   = max(partVertices) - min(partVertices);
chuteLen   = max(partSpan) * 5;
chuteDepth = max(partSpan) * 3;
wallHeight = max(partSpan) * 3;
Xh         = chuteLen / 2;
axisExt    = max(partSpan) * 0.55;

p4ReceivedFrom = zeros(numStable,1);
for si = 1:numStable
    d = transitionDest(si);
    if d > 0, p4ReceivedFrom(d) = si; end
end

wallTransAxisW  = normaliseVec((Rchute * [0;0;1])');
floorTransAxisW = normaliseVec((Rchute * [0;1;0])');

for si = 1:numStable

    pi_  = stableIdx(si);
    pos  = allRestingPositions(pi_);

    vertsC          = pos.verticesWorld;
    centC           = pos.centroidWorld(:)';
    floorContactIdx = pos.floorContactVertIdx;
    wallContactIdx  = pos.wallContactVertIdx;

    floorZ         = min(vertsC(floorContactIdx,3));
    vertsC(:,3)    = vertsC(:,3) - floorZ;
    centC(3)       = centC(3)    - floorZ;

    wallY          = max(vertsC(wallContactIdx,2));
    vertsC(:,2)    = vertsC(:,2) - wallY;
    centC(2)       = centC(2)    - wallY;

    vertsWorld = (Rchute * vertsC')';
    centWorld  = (Rchute * centC(:))';

    floorCorners = [-Xh,-chuteDepth,0; Xh,-chuteDepth,0; Xh,0,0; -Xh,0,0];
    wallCorners  = [-Xh,0,0; Xh,0,0; Xh,0,wallHeight; -Xh,0,wallHeight];
    floorW = (Rchute * floorCorners')';
    wallW  = (Rchute * wallCorners')';

    supportIdx    = union(floorContactIdx, wallContactIdx);
    supportVerts  = vertsWorld(supportIdx,:);
    allVertIdx    = (1:size(vertsWorld,1))';
    floorSet      = floorContactIdx(:);
    wallSet       = wallContactIdx(:);
    nonContactSet = setdiff(allVertIdx, supportIdx);
    floorOnlySet  = setdiff(floorSet, wallSet);
    wallOnlySet   = setdiff(wallSet, floorSet);
    bothSet       = intersect(floorSet, wallSet);

    floorVertsW = vertsWorld(floorContactIdx, :);
    floorMeanW  = mean(floorVertsW, 1);
    if size(floorVertsW,1) >= 3
        [~,~,Vsvd] = svd(floorVertsW - floorMeanW);
        floorNW = Vsvd(:,3)';
    else
        floorNW = [0,0,1];
    end
    if dot(floorNW, [0,0,-1]) > 0, floorNW = -floorNW; end
    floorNW = floorNW / norm(floorNW);

    rayDir = [0,0,-1];
    denom  = dot(rayDir, floorNW);
    if abs(denom) > 1e-10
        t_viz         = dot(floorMeanW - centWorld, floorNW) / denom;
        centProjFloor = centWorld + t_viz * rayDir;
    else
        centProjFloor    = centWorld;
        centProjFloor(3) = min(floorVertsW(:,3));
    end
    planeZ = centProjFloor(3);

    supportProj      = supportVerts;
    supportProj(:,3) = planeZ;
    supportProj      = uniquetol(supportProj, 1e-6, 'ByRows', true);

    hullProj3D = [];
    try
        hullIdx2   = convhull(supportProj(:,1), supportProj(:,2));
        hullProj3D = supportProj(hullIdx2(1:end-1),:);
    catch
    end

    isInsideHull = false;
    if ~isempty(hullProj3D) && size(hullProj3D,1) >= 3
        try
            isInsideHull = inpolygon( ...
                centProjFloor(1), centProjFloor(2), ...
                hullProj3D(:,1),  hullProj3D(:,2));
        catch
            isInsideHull = false;
        end
    end
    if isInsideHull
        stabilityStr = 'STABLE';
        stabilityCol = [0.05 0.55 0.05];
    else
        stabilityStr = 'UNSTABLE';
        stabilityCol = [0.85 0.10 0.10];
    end

    ov_vertsWorld = [];
    ov_centWorld  = [];
    ov_faceCol    = [0 0 0];
    ov_edgeCol    = [0 0 0];
    if transitions(si)
        vertsC_ov = transitionRotatedVerts{si};
        centC_ov  = transitionRotatedCent{si};
        if ~isempty(vertsC_ov)
            ov_vertsWorld = (Rchute * vertsC_ov')';
            ov_centWorld  = (Rchute * centC_ov(:))';
            if transitionDest(si) > 0
                ov_faceCol = [0.20, 0.50, 1.00];
                ov_edgeCol = [0.05, 0.25, 0.70];
            else
                ov_faceCol = [1.00, 0.25, 0.20];
                ov_edgeCol = [0.70, 0.08, 0.05];
            end
        end
    end

    bNote    = '';
    aSrcNote = '';
    aDstNote = '';
    if transitions(si),        bNote    = ' [B:removed]'; end
    if transitionDest(si) > 0, aSrcNote = sprintf(' [A:->Pose%d]', transitionDest(si)); end
    if p4ReceivedFrom(si) > 0, aDstNote = sprintf(' [A:+Pose%d]', p4ReceivedFrom(si)); end
    tNote = [bNote, aSrcNote, aDstNote];

    geo = momentArmGeo(si);

    hasWallGeo = ~isempty(geo.pivotWall_c);
    if hasWallGeo
        pivW_w   = (Rchute * geo.pivotWall_c(:))';
        centC_   = geo.centC;
        centWall_c    = centC_;
        centWall_c(3) = geo.pivotWall_c(3);
        armEndW_w     = (Rchute * centWall_c(:))';
        axWall_A = pivW_w - axisExt * wallTransAxisW;
        axWall_B = pivW_w + axisExt * wallTransAxisW;
    end

    hasFloorGeo = ~isempty(geo.pivotFloor_c);
    if hasFloorGeo
        pivF_w   = (Rchute * geo.pivotFloor_c(:))';
        centC_   = geo.centC;
        centFloor_c    = centC_;
        centFloor_c(2) = geo.pivotFloor_c(2);
        armEndF_w      = (Rchute * centFloor_c(:))';
        axFloor_A = pivF_w - axisExt * floorTransAxisW;
        axFloor_B = pivF_w + axisExt * floorTransAxisW;
    end

    quatStr = buildQuatAnnotation(si, transitions, chainLog, refQuats, refQuatsAll, quatMatchTol);

    hFig = figure('Visible','off', ...
        'Color','w','Units','inches','Position',[0 0 18 8], ...
        'PaperUnits','inches','PaperSize',[18 8],'PaperPosition',[0 0 18 8]);

    tl = tiledlayout(hFig, 1, 2, 'TileSpacing','compact','Padding','compact');

    title(tl, sprintf( ...
        '%s  |  Roll=%d\xb0  Pitch=%d\xb0  |  Pose %d (Plane %d, \\theta=%d\xb0)  |  CSA_B=%.1f%% CSA_A=%.1f%% CSA_N=%.1f%%  |  CRSA_B=%.1f%% CRSA_A=%.1f%% CRSA_N=%.1f%%%s', ...
        partName, round(chuteRoll_deg), round(chutePitch_deg), si, ...
        pos.floorPlaneIdx, round(rad2deg(pos.theta)), ...
        probCSA_B(si), probCSA_A(si), probCSA_N(si), ...
        probCRSA_B(si), probCRSA_A(si), probCRSA_N(si), tNote), ...
        'FontSize', 10, 'FontWeight', 'bold', 'Interpreter', 'tex');

    % ── DISPLAY-ONLY mirror + rotate (both ax1 "part in chute" and ax2 ────
    % "CSA geometry" panels). Mirrors ChuteTransitionSolver.m's
    % renderPoseAxes convention: flip the rendered geometry about the
    % yz-plane (x -> -x). Purely a viewer transform — only these *_m
    % copies feed the plots; vertsWorld/centWorld/supportProj/hullProj3D/
    % centProjFloor as originally computed (and all upstream CSA/CRSA/
    % ratio math) are untouched. Mirroring reverses triangle winding, so
    % partFaces columns are flipped to keep lighting correct, and the
    % camera azimuth is negated on both panels to compensate for the
    % handedness flip.
    vertsWorld_m = vertsWorld; vertsWorld_m(:,1) = -vertsWorld_m(:,1);
    centWorld_m  = centWorld;  centWorld_m(1)    = -centWorld_m(1);
    floorW_m = floorW; floorW_m(:,1) = -floorW_m(:,1);
    wallW_m  = wallW;  wallW_m(:,1)  = -wallW_m(:,1);
    if ~isempty(partFaces), partFaces_m = partFaces(:, [1 3 2]); else, partFaces_m = partFaces; end
    ov_vertsWorld_m = ov_vertsWorld; ov_centWorld_m = ov_centWorld;
    if ~isempty(ov_vertsWorld)
        ov_vertsWorld_m(:,1) = -ov_vertsWorld_m(:,1);
        ov_centWorld_m(1)    = -ov_centWorld_m(1);
    end
    if hasWallGeo
        pivW_w_m    = pivW_w;    pivW_w_m(1)    = -pivW_w_m(1);
        armEndW_w_m = armEndW_w; armEndW_w_m(1) = -armEndW_w_m(1);
        axWall_A_m  = axWall_A;  axWall_A_m(1)  = -axWall_A_m(1);
        axWall_B_m  = axWall_B;  axWall_B_m(1)  = -axWall_B_m(1);
        pivOnWall_c    = geo.pivotWall_c; pivOnWall_c(2) = 0;
        pivOnWall_w_m  = (Rchute * pivOnWall_c(:))'; pivOnWall_w_m(1) = -pivOnWall_w_m(1);
    end
    if hasFloorGeo
        pivF_w_m    = pivF_w;    pivF_w_m(1)    = -pivF_w_m(1);
        armEndF_w_m = armEndF_w; armEndF_w_m(1) = -armEndF_w_m(1);
        axFloor_A_m = axFloor_A; axFloor_A_m(1) = -axFloor_A_m(1);
        axFloor_B_m = axFloor_B; axFloor_B_m(1) = -axFloor_B_m(1);
        pivOnFloor_c    = geo.pivotFloor_c; pivOnFloor_c(3) = 0;
        pivOnFloor_w_m  = (Rchute * pivOnFloor_c(:))'; pivOnFloor_w_m(1) = -pivOnFloor_w_m(1);
    end
    supportProj_m   = supportProj;   supportProj_m(:,1)   = -supportProj_m(:,1);
    centProjFloor_m = centProjFloor; centProjFloor_m(1)   = -centProjFloor_m(1);
    if ~isempty(hullProj3D)
        hullProj3D_m = hullProj3D; hullProj3D_m(:,1) = -hullProj3D_m(:,1);
    else
        hullProj3D_m = hullProj3D;
    end

    ax1 = nexttile(tl, 1);
    hold(ax1,'on'); grid(ax1,'on');
    view(ax1,-40,24);
    lighting(ax1,'gouraud'); camlight(ax1,'headlight');
    xlabel(ax1,'X (slide)'); ylabel(ax1,'Y (wall\perp)'); zlabel(ax1,'Z (floor\perp)');

    panelTitle1 = sprintf('Pose %d | %s | P_{CSA_N}=%.2f%%', ...
        si, stabilityStr, probCSA_N(si));
    title(ax1, panelTitle1, 'FontSize', 8, 'FontWeight', 'bold', ...
        'Interpreter', 'tex', 'Color', stabilityCol);

    patch(ax1, floorW_m(:,1), floorW_m(:,2), floorW_m(:,3), ...
        [0.72 0.90 0.72], 'FaceAlpha', 0.25, 'EdgeColor', [0.20 0.55 0.20], 'LineWidth', 0.8);
    patch(ax1, wallW_m(:,1), wallW_m(:,2), wallW_m(:,3), ...
        [0.72 0.90 0.72], 'FaceAlpha', 0.25, 'EdgeColor', [0.20 0.55 0.20], 'LineWidth', 0.8);

    trisurf(partFaces_m, vertsWorld_m(:,1), vertsWorld_m(:,2), vertsWorld_m(:,3), ...
        'Parent', ax1, 'FaceColor', [1.00 0.92 0.23], 'FaceAlpha', 0.50, ...
        'EdgeColor', [0.40 0.35 0.00], 'EdgeAlpha', 0.20);

    scatter3(ax1, vertsWorld_m(floorContactIdx,1), vertsWorld_m(floorContactIdx,2), ...
        vertsWorld_m(floorContactIdx,3), 30, [0.10 0.60 0.10], 'filled', 'MarkerFaceAlpha', 0.9);
    scatter3(ax1, vertsWorld_m(wallContactIdx,1), vertsWorld_m(wallContactIdx,2), ...
        vertsWorld_m(wallContactIdx,3), 30, [0.10 0.60 0.10], 'filled', 'MarkerFaceAlpha', 0.9);

    scatter3(ax1, centWorld_m(1), centWorld_m(2), centWorld_m(3), 70, 'k', 'filled');

    if ~isempty(ov_vertsWorld_m)
        trisurf(partFaces_m, ov_vertsWorld_m(:,1), ov_vertsWorld_m(:,2), ov_vertsWorld_m(:,3), ...
            'Parent', ax1, 'FaceColor', ov_faceCol, 'FaceAlpha', 0.45, ...
            'EdgeColor', ov_edgeCol, 'EdgeAlpha', 0.20);
        scatter3(ax1, ov_centWorld_m(1), ov_centWorld_m(2), ov_centWorld_m(3), ...
            80, ov_faceCol, 'filled', 'MarkerEdgeColor', ov_edgeCol, 'LineWidth', 1.5);
    end

    if hasWallGeo
        wallTrip  = (isfinite(threshWall) && ratioWall(si) >= threshWall);
        wallLineW = 1.8 + 1.2 * wallTrip;
        wallCol   = [1.0 0.45 0.0];
        scatter3(ax1, pivW_w_m(1), pivW_w_m(2), pivW_w_m(3), 100, wallCol, 'd', 'filled');
        plot3(ax1, [pivW_w_m(1), armEndW_w_m(1)], [pivW_w_m(2), armEndW_w_m(2)], ...
            [pivW_w_m(3), armEndW_w_m(3)], '-', 'Color', wallCol, 'LineWidth', wallLineW);
        plot3(ax1, [pivW_w_m(1), pivOnWall_w_m(1)], [pivW_w_m(2), pivOnWall_w_m(2)], ...
            [pivW_w_m(3), pivOnWall_w_m(3)], '-.', 'Color', wallCol, 'LineWidth', wallLineW);
        plot3(ax1, [axWall_A_m(1), axWall_B_m(1)], [axWall_A_m(2), axWall_B_m(2)], ...
            [axWall_A_m(3), axWall_B_m(3)], '--', 'Color', [0.85 0.20 0.85], 'LineWidth', 1.8);
        midW = 0.5*(pivW_w_m + armEndW_w_m);
        text(ax1, midW(1), midW(2), midW(3), ...
            sprintf('  l_A=%.3f\n  l_w=%.3f\n  r=%.3f', geo.l_A_wall, geo.l_w, ratioWall(si)), ...
            'Color', wallCol, 'FontSize', 7, 'FontWeight', 'bold', 'Interpreter', 'none');
    end

    if hasFloorGeo
        floorTrip  = (isfinite(threshFloor) && ratioFloor(si) >= threshFloor);
        floorLineW = 1.8 + 1.2 * floorTrip;
        floorCol   = [0.0 0.65 0.65];
        scatter3(ax1, pivF_w_m(1), pivF_w_m(2), pivF_w_m(3), 100, floorCol, 'd', 'filled');
        plot3(ax1, [pivF_w_m(1), armEndF_w_m(1)], [pivF_w_m(2), armEndF_w_m(2)], ...
            [pivF_w_m(3), armEndF_w_m(3)], '-', 'Color', floorCol, 'LineWidth', floorLineW);
        plot3(ax1, [pivF_w_m(1), pivOnFloor_w_m(1)], [pivF_w_m(2), pivOnFloor_w_m(2)], ...
            [pivF_w_m(3), pivOnFloor_w_m(3)], '-.', 'Color', floorCol, 'LineWidth', floorLineW);
        plot3(ax1, [axFloor_A_m(1), axFloor_B_m(1)], [axFloor_A_m(2), axFloor_B_m(2)], ...
            [axFloor_A_m(3), axFloor_B_m(3)], '--', 'Color', [0.20 0.70 0.20], 'LineWidth', 1.8);
        midF = 0.5*(pivF_w_m + armEndF_w_m);
        text(ax1, midF(1), midF(2), midF(3), ...
            sprintf('  l_A=%.3f\n  l_f=%.3f\n  r=%.3f', geo.l_A_floor, geo.l_f, ratioFloor(si)), ...
            'Color', floorCol, 'FontSize', 7, 'FontWeight', 'bold', 'Interpreter', 'none');
    end

    threshStr = buildThresholdAnnotation( ...
        ratioWall(si), ratioFloor(si), threshWall, threshFloor, ...
        transitions(si), stabilityStr);
    annotation(hFig, 'textbox', 'Units', 'normalized', 'Position', [0.01, 0.68, 0.24, 0.24], ...
        'String', threshStr, 'FontSize', 7.5, 'FontName', 'Courier', ...
        'BackgroundColor', [0.97 0.97 0.97], 'EdgeColor', [0.4 0.4 0.4], ...
        'LineWidth', 0.8, 'Margin', 4, 'FitBoxToText', 'on', 'Interpreter', 'none');

    annotation(hFig, 'textbox', 'Units', 'normalized', 'Position', [0.27, 0.68, 0.24, 0.24], ...
        'String', quatStr, 'FontSize', 7.5, 'FontName', 'Courier', ...
        'BackgroundColor', [0.94 0.97 1.00], 'EdgeColor', [0.30 0.45 0.70], ...
        'LineWidth', 0.8, 'Margin', 4, 'FitBoxToText', 'on', 'Interpreter', 'none');

    zoomPts = vertsWorld_m;
    if ~isempty(ov_vertsWorld_m), zoomPts = [zoomPts; ov_vertsWorld_m]; end
    spanPts = range(zoomPts);
    padAmt  = max(0.20 * max(spanPts) * ones(1,3), 1e-3 * ones(1,3));
    pMin    = min(zoomPts) - padAmt;
    pMax    = max(zoomPts) + padAmt;
    axis(ax1, 'equal');
    xlim(ax1, [pMin(1), pMax(1)]);
    ylim(ax1, [pMin(2), pMax(2)]);
    zlim(ax1, [pMin(3), pMax(3)]);

    ax2 = nexttile(tl, 2);
    hold(ax2,'on'); grid(ax2,'on'); axis(ax2,'equal');
    view(ax2,-40,24);
    xlabel(ax2,'X'); ylabel(ax2,'Y'); zlabel(ax2,'Z');

    panelTitle2 = sprintf('CSA Geometry | \\omega=%.4f sr | h=%.4f | rW=%.3f/%.3f rF=%.3f/%.3f', ...
        omegas(si), heights(si), ratioWall(si), threshWall, ratioFloor(si), threshFloor);
    title(ax2, panelTitle2, 'FontSize', 8, 'FontWeight', 'normal', 'Interpreter', 'tex');

    trisurf(partFaces_m, vertsWorld_m(:,1), vertsWorld_m(:,2), vertsWorld_m(:,3), ...
        'Parent', ax2, 'FaceColor', [0.8 0.8 0.8], 'FaceAlpha', 0.12, ...
        'EdgeColor', [0.7 0.7 0.7], 'EdgeAlpha', 0.08);
    if ~isempty(nonContactSet)
        scatter3(ax2, vertsWorld_m(nonContactSet,1), vertsWorld_m(nonContactSet,2), ...
            vertsWorld_m(nonContactSet,3), 8, [0.65 0.65 0.65], 'filled', 'MarkerFaceAlpha', 0.4);
    end
    if ~isempty(floorOnlySet)
        scatter3(ax2, vertsWorld_m(floorOnlySet,1), vertsWorld_m(floorOnlySet,2), ...
            vertsWorld_m(floorOnlySet,3), 16, [0.15 0.45 0.90], 'filled');
    end
    if ~isempty(wallOnlySet)
        scatter3(ax2, vertsWorld_m(wallOnlySet,1), vertsWorld_m(wallOnlySet,2), ...
            vertsWorld_m(wallOnlySet,3), 16, [0.85 0.15 0.15], 'filled');
    end
    if ~isempty(bothSet)
        scatter3(ax2, vertsWorld_m(bothSet,1), vertsWorld_m(bothSet,2), ...
            vertsWorld_m(bothSet,3), 16, [0.80 0.10 0.80], 'filled');
    end
    for kk = 1:size(supportProj_m,1)
        p_ = supportProj_m(kk,:);
        plot3(ax2, [p_(1),p_(1)], [p_(2),p_(2)], [p_(3),planeZ], ...
            ':', 'Color', [0.6 0.6 0.6 0.5], 'LineWidth', 0.8);
    end
    if ~isempty(hullProj3D_m)
        patch(ax2, hullProj3D_m(:,1), hullProj3D_m(:,2), hullProj3D_m(:,3), ...
            [0.2 0.6 1.0], 'FaceAlpha', 0.30, ...
            'EdgeColor', [0.1 0.35 0.80], 'LineWidth', 1.5);
        for kk = 1:size(hullProj3D_m,1)
            plot3(ax2, [centWorld_m(1), hullProj3D_m(kk,1)], ...
                [centWorld_m(2), hullProj3D_m(kk,2)], ...
                [centWorld_m(3), hullProj3D_m(kk,3)], ...
                'Color', [0.3 0.3 0.3], 'LineWidth', 1.2);
        end
    end
    scatter3(ax2, supportProj_m(:,1), supportProj_m(:,2), supportProj_m(:,3), ...
        10, [0.1 0.35 0.80], 'filled', 'MarkerFaceAlpha', 0.7);
    scatter3(ax2, centWorld_m(1), centWorld_m(2), centWorld_m(3), 55, 'k', 'filled');
    scatter3(ax2, centProjFloor_m(1), centProjFloor_m(2), centProjFloor_m(3), ...
        55, [0.05 0.60 0.05], 'filled');
    plot3(ax2, [centWorld_m(1), centProjFloor_m(1)], ...
        [centWorld_m(2), centProjFloor_m(2)], ...
        [centWorld_m(3), centProjFloor_m(3)], ...
        '--', 'Color', [0.05 0.60 0.05], 'LineWidth', 1.8);

    legendStr = sprintf('\\color[rgb]{0.15,0.45,0.90}Floor  \\color[rgb]{0.85,0.15,0.15}Wall  \\color[rgb]{0.8,0.1,0.8}Both  \\color{black}Centroid  \\color[rgb]{0.05,0.60,0.05}CentProj');
    text(ax2, 0.02, 0.02, legendStr, 'Units','normalized', ...
        'FontSize', 7, 'Interpreter', 'tex', 'VerticalAlignment', 'bottom');
    drawnow;
    if si == 1
        exportgraphics(hFig, pdfPath, 'Resolution', 150, 'Append', false);
    else
        exportgraphics(hFig, pdfPath, 'Resolution', 150, 'Append', true);
    end
    close(hFig);
    drawnow;

    clear vertsC vertsWorld centWorld floorCorners wallCorners floorW wallW;
    clear supportVerts supportProj hullProj3D zoomPts;
    clear floorVertsW ov_vertsWorld ov_centWorld;

end % si
end % saveFiguresAndPDF


% =========================================================================
%  buildQuatAnnotation
% =========================================================================
function lines = buildQuatAnnotation(si, transitions, chainLog, refQuats, refQuatsAll, quatMatchTol)

q_ref = refQuats(si, :);
refStr = sprintf('refQuat: [%6.3f %6.3f %6.3f %6.3f]', ...
    q_ref(1), q_ref(2), q_ref(3), q_ref(4));

mergeNote = '';
if nargin >= 5 && ~isempty(refQuatsAll) && numel(refQuatsAll) >= si
    nAlt = size(refQuatsAll{si}, 1);
    if nAlt > 1
        mergeNote = sprintf('  (+%d symmetric orientation(s) merged in)', nAlt - 1);
    end
end

if ~transitions(si) || isempty(chainLog{si})
    lines = {[refStr, mergeNote]; 'No transition attempted'};
    return;
end

log     = chainLog{si};
lastEntry = log{end};

if isfield(lastEntry, 'q_accumulated')
    q_acc = lastEntry.q_accumulated;
    accStr = sprintf('accumQuat:[%6.3f %6.3f %6.3f %6.3f]', ...
        q_acc(1), q_acc(2), q_acc(3), q_acc(4));
else
    accStr = 'accumQuat: unavailable';
    q_acc  = [];
end

if ~isempty(q_acc)
    dp   = abs(dot(q_ref(:), q_acc(:)));
    dp   = min(dp, 1.0);
    dist = 2 * acos(dp);
    selfDistStr = sprintf('selfDist : %.4f rad (vs own refQuat)', dist);
else
    selfDistStr = 'selfDist : unavailable';
end

if isfield(lastEntry, 'matchDist')
    matchDistStr = sprintf('matchDist: %.4f rad (tol=%.4f)', lastEntry.matchDist, quatMatchTol);
else
    matchDistStr = 'matchDist: unavailable';
end

if isfield(lastEntry, 'matchSi') && lastEntry.matchSi > 0
    matchStr = sprintf('match    : Pose %d  [MATCHED]', lastEntry.matchSi);
elseif isfield(lastEntry, 'matchDist')
    matchStr = sprintf('match    : none (best=%.4f rad)', lastEntry.matchDist);
else
    matchStr = 'match    : none';
end

numIters = numel(log);
iterStr  = sprintf('chain len: %d step(s)', numIters);

lines = {[refStr, mergeNote]; accStr; selfDistStr; matchDistStr; matchStr; iterStr};
end


% =========================================================================
%  buildThresholdAnnotation
% =========================================================================
function lines = buildThresholdAnnotation(rW, rF, tW, tF, isTransition, stabilityStr)

wallStr  = sprintf('Wall  rW=%.3f / tW=%.3f [%s]', rW, tW, sel(rW >= tW, 'TRIP', 'ok'));
floorStr = sprintf('Floor rF=%.3f / tF=%.3f [%s]', rF, tF, sel(rF >= tF, 'UNSTABLE', 'ok'));

if isTransition
    transStr = 'Wall transition: YES';
else
    transStr = 'Wall transition: no';
end

lines = {wallStr; floorStr; transStr; ['Stability: ', stabilityStr]};
end

function s = sel(cond, a, b)
if cond, s = a; else, s = b; end
end


% =========================================================================
%  WRITE GEOMETRIC SUMMARY
% =========================================================================
function writeGeometricSummary( ...
    allRestingPositions, planeQuats, ...
    summaryRoot, folderName, partName, ...
    chuteRoll_deg, chutePitch_deg)

numPoses = numel(allRestingPositions);
fpath    = fullfile(summaryRoot, [folderName, '_GEOMETRIC_summary.txt']);
fid      = fopen(fpath, 'w');

fprintf(fid, '=========================================================================\n');
fprintf(fid, ' GEOMETRIC SUMMARY  (all poses after theta merge — pre-stability)\n');
fprintf(fid, ' Part  : %s\n', partName);
fprintf(fid, ' Chute : Roll = %d deg   Pitch = %d deg\n', ...
    round(chuteRoll_deg), round(chutePitch_deg));
fprintf(fid, ' Total poses: %d\n', numPoses);
fprintf(fid, '=========================================================================\n\n');

fprintf(fid, '%-4s  %-6s  %-8s  %-38s  %-s\n', ...
    'Pose', 'Plane', 'Theta', 'refQuat [w x y z]', 'FloorVerts | WallVerts');
fprintf(fid, '%s\n', repmat('-', 1, 100));

for ci = 1:numPoses
    pos = allRestingPositions(ci);

    q_align = planeQuats{pos.floorPlaneIdx};
    q_theta = q_fromAxisAngle([0,0,1], pos.theta);
    q_ref   = q_compose(q_align, q_theta);
    [~, mi] = max(abs(q_ref));
    if q_ref(mi) < 0, q_ref = -q_ref; end

    quatStr   = sprintf('[%7.4f %7.4f %7.4f %7.4f]', ...
        q_ref(1), q_ref(2), q_ref(3), q_ref(4));
    floorStr  = mat2str(pos.floorContactVertIdx);
    wallStr   = mat2str(pos.wallContactVertIdx);
    contactStr = sprintf('%s | %s', floorStr, wallStr);

    fprintf(fid, '%-4d  %-6d  %-8.2f  %-38s  %-s\n', ...
        ci, pos.floorPlaneIdx, rad2deg(pos.theta), quatStr, contactStr);
end

fprintf(fid, '\n');
fprintf(fid, 'Column key:\n');
fprintf(fid, '  Pose       = pose index in the merged candidate list\n');
fprintf(fid, '  Plane      = floor plane index (from findFloorPlanes)\n');
fprintf(fid, '  Theta      = rotation about floor-normal axis (deg)\n');
fprintf(fid, '  refQuat    = composed quaternion q_align * q_theta [w x y z]\n');
fprintf(fid, '  FloorVerts = vertex indices in floor contact set\n');
fprintf(fid, '  WallVerts  = vertex indices in wall contact set\n');

fclose(fid);
fprintf('  Saved geometric summary: %s\n', fpath);
end


% =========================================================================
%  WRITE MASTER SUMMARY
% =========================================================================
function writeMasterSummary( ...
    allRestingPositions, stableIdx, Qs, omegas, heights, sumQ, ...
    probCSA_B, probCSA_A, probCSA_N, probCRSA_B, probCRSA_A, probCRSA_N, ...
    transitionDest, transitions, ...
    ratioWall, ratioFloor, threshWall, threshFloor, ...
    refQuats, refQuatsAll, ...
    chuteRoll_deg, chutePitch_deg, partName, summaryRoot, folderName)

numStable = numel(stableIdx);

p4ReceivedFrom = zeros(numStable,1);
for si = 1:numStable
    d = transitionDest(si);
    if d > 0, p4ReceivedFrom(d) = si; end
end

fpath = fullfile(summaryRoot, [folderName, '_MASTER_summary.txt']);
fid   = fopen(fpath, 'w');

fprintf(fid, '=========================================================================\n');
fprintf(fid, ' MASTER SUMMARY  (stable poses — post CSA analysis)\n');
fprintf(fid, ' Part  : %s\n', partName);
fprintf(fid, ' Chute : Roll = %d deg   Pitch = %d deg\n', ...
    round(chuteRoll_deg), round(chutePitch_deg));
fprintf(fid, '=========================================================================\n\n');

fprintf(fid, 'Stable poses                  : %d\n', numStable);
fprintf(fid, 'Sum of raw CSA weights (sumQ) : %.6f\n', sumQ);
fprintf(fid, 'Wall  transition threshold    : %.4f  (l_A/l_w >= this -> transitions)\n', threshWall);
fprintf(fid, 'Floor instability threshold   : %.4f  (l_A/l_f >= this -> Qs zeroed)\n\n', threshFloor);

fprintf(fid, '%-4s %-6s %-7s %-38s %-10s %-9s %-8s %-8s %-9s %-9s %-9s %-9s %-9s %-9s %-7s %-7s %-6s %-s\n', ...
    'Pos', 'Plane', 'Theta', 'refQuat [w x y z]', 'Omega(sr)', 'Height', ...
    'rWall', 'rFloor', 'CSA_B(%)', 'CSA_A(%)', 'CSA_N(%)', ...
    'CRSA_B(%)', 'CRSA_A(%)', 'CRSA_N(%)', '->Dest', '+From', 'nMerge', 'Notes');
fprintf(fid, '%s\n', repmat('-', 1, 175));

for si = 1:numStable
    ci_ = stableIdx(si);
    pos = allRestingPositions(ci_);

    q_ref   = refQuats(si, :);
    quatStr = sprintf('[%7.4f %7.4f %7.4f %7.4f]', ...
        q_ref(1), q_ref(2), q_ref(3), q_ref(4));

    notes = '';
    if transitions(si),          notes = [notes, 'TRANSITIONS ']; end
    if Qs(si) == 0 && ~transitions(si)
        if ratioFloor(si) >= threshFloor
            notes = [notes, 'FLOOR-UNSTABLE '];
        else
            notes = [notes, 'Q=0 '];
        end
    end
    notes = strtrim(notes);
    if isempty(notes), notes = '-'; end

    destStr = '-';
    if transitionDest(si) > 0, destStr = sprintf('%d', transitionDest(si)); end
    fromStr = '-';
    if p4ReceivedFrom(si) > 0, fromStr = sprintf('%d', p4ReceivedFrom(si)); end

    nMerge = 1;
    if ~isempty(refQuatsAll) && numel(refQuatsAll) >= si
        nMerge = size(refQuatsAll{si}, 1);
    end

    fprintf(fid, '%-4d %-6d %-7.0f %-38s %-10.4f %-9.4f %-8.3f %-8.3f %-9.2f %-9.2f %-9.2f %-9.2f %-9.2f %-9.2f %-7s %-7s %-6d %-s\n', ...
        si, pos.floorPlaneIdx, round(rad2deg(pos.theta)), quatStr, ...
        omegas(si), heights(si), ratioWall(si), ratioFloor(si), ...
        probCSA_B(si), probCSA_A(si), probCSA_N(si), ...
        probCRSA_B(si), probCRSA_A(si), probCRSA_N(si), ...
        destStr, fromStr, nMerge, notes);
end

fprintf(fid, '%s\n', repmat('-', 1, 175));
fprintf(fid, '%-4s %-6s %-7s %-38s %-10s %-9s %-8s %-8s %-9.2f %-9.2f %-9.2f %-9.2f %-9.2f %-9.2f\n', ...
    'TOT', '', '', '', '', '', '', '', ...
    sum(probCSA_B), sum(probCSA_A), sum(probCSA_N), ...
    sum(probCRSA_B), sum(probCRSA_A), sum(probCRSA_N));

fprintf(fid, '\nAlgorithm key:\n');
fprintf(fid, '  CSA_B  = CSA,  transitioning poses zeroed, remaining renormalized  ("before")\n');
fprintf(fid, '  CSA_A  = CSA,  transitioning poses redistribute to destination pose ("after")\n');
fprintf(fid, '  CSA_N  = CSA,  no filtering at all (raw CSA weight normalized)     ("none")\n');
fprintf(fid, '  CRSA_B = CRSA, transitioning poses zeroed, remaining renormalized  ("before")\n');
fprintf(fid, '  CRSA_A = CRSA, transitioning poses redistribute to destination pose ("after")\n');
fprintf(fid, '  CRSA_N = CRSA, no filtering at all (raw CRSA weight normalized)     ("none")\n');
fprintf(fid, '\n  NOTE: CSA_B/CSA_A/CSA_N are computed from an UNZEROED CSA weight\n');
fprintf(fid, '  (QsRaw), just like CRSA_B/CRSA_A/CRSA_N are computed from an\n');
fprintf(fid, '  unzeroed crsaRaw. The "Qs"/"Q=0"/"FLOOR-UNSTABLE" notes above\n');
fprintf(fid, '  reflect a SEPARATE validity-zeroed weight used only for pose\n');
fprintf(fid, '  merging, sumQ, and the diagnostic plots/figures.\n');
fprintf(fid, '\n  Transition  : wall only, l_A/l_w >= %.4f\n', threshWall);
fprintf(fid, '  Instability : floor only, l_A/l_f >= %.4f  (Qs zeroed)\n', threshFloor);
fprintf(fid, '  rWall  = l_A/l_w moment-arm ratio (wall-tip transition driver)\n');
fprintf(fid, '  rFloor = l_A/l_f moment-arm ratio (floor instability indicator)\n');

fclose(fid);
fprintf('  Saved master summary: %s\n', fpath);
end


% =========================================================================
%  SOLID ANGLE UTILITIES
% =========================================================================
function omega = solidAnglePolygon(hullPts, pt)
verts = hullPts(1:end-1, :);
n     = size(verts,1);
if n < 3; omega = 0; return; end
v0    = verts(1,:);
omega = 0;
for i = 2:n-1
    omega = omega + solidAngleTriangle(v0, verts(i,:), verts(i+1,:), pt);
end
end

function sa = solidAngleTriangle(a, b, c, p)
a = a-p; b = b-p; c = c-p;
la = norm(a); lb = norm(b); lc = norm(c);
num = abs(det([a;b;c]));
den = la*lb*lc + dot(a,b)*lc + dot(b,c)*la + dot(c,a)*lb;
sa  = 2*atan2(num, den);
end


% =========================================================================
%  GEOMETRY HELPERS
% =========================================================================
function v = normaliseVec(v)
n = norm(v);
if n > 0; v = v / n; end
end


% =========================================================================
%  SHARED GEOMETRY FUNCTIONS
% =========================================================================

function [myCentroid] = centroidOfPolyhedron(vertex, faces)
vector1 = vertex(faces(:,2),:) - vertex(faces(:,1),:);
vector2 = vertex(faces(:,3),:) - vertex(faces(:,1),:);
triangAreasTmp = 0.5 * cross(vector1, vector2);
triangAreas(:,1) = sqrt( triangAreasTmp(:,1).^2 + ...
    triangAreasTmp(:,2).^2 + ...
    triangAreasTmp(:,3).^2 );
totArea = sum(triangAreas);
point1 = vertex(faces(:,1),:);
point2 = vertex(faces(:,2),:);
point3 = vertex(faces(:,3),:);
centroidTriangles = (1/3) .* (point1 + point2 + point3);
mg(:,1) = triangAreas(:,1) .* centroidTriangles(:,1);
mg(:,2) = triangAreas(:,1) .* centroidTriangles(:,2);
mg(:,3) = triangAreas(:,1) .* centroidTriangles(:,3);
myCentroid = sum(mg) ./ totArea;
end

function [planeVerts, planeEqs] = findFloorPlanes( ...
    partVertices, convexHullFaces, planeTol)
planeVerts = {};
planeEqs   = [];
for faceNumber = 1:size(convexHullFaces,1)
    v1 = partVertices(convexHullFaces(faceNumber,1),:);
    v2 = partVertices(convexHullFaces(faceNumber,2),:);
    v3 = partVertices(convexHullFaces(faceNumber,3),:);
    normalVector = cross(v2-v1, v3-v1);
    if norm(normalVector) < 1e-10, continue; end
    A = normalVector(1); B = normalVector(2); C = normalVector(3);
    D = dot(normalVector, v1);
    residuals = partVertices * normalVector' - D;
    onPlane   = find(abs(residuals) < planeTol);
    if isempty(onPlane), continue; end
    sortedOnPlane = sort(onPlane(:))';
    isDuplicate = false;
    for k = 1:numel(planeVerts)
        if isequal(sort(planeVerts{k}), sortedOnPlane)
            isDuplicate = true; break;
        end
    end
    if ~isDuplicate
        planeVerts{end+1} = sortedOnPlane;   %#ok<AGROW>
        planeEqs(end+1,:) = [A, B, C, D];   %#ok<AGROW>
    end
end
end

function pts_out = rotatePtsAroundZ(pts_in, pivot, theta)
s = sin(theta); c = cos(theta);
R = [c,-s,0; s,c,0; 0,0,1];
pts_c   = pts_in - pivot;
pts_out = (R * pts_c')' + pivot;
end

function merged = mergeByTheta(rawCandidates, thetaTol_deg, ...
    planeMergeTol, partVertices) %#ok<INUSD>

merged = struct( ...
    'floorPlaneIdx',            {}, ...
    'wallSide',                 {}, ...
    'wallContactVertIdx',       {}, ...
    'floorContactVertIdx',      {}, ...
    'theta',                    {}, ...
    'verticesWorld',            {}, ...
    'centroidWorld',            {}, ...
    'floorZ',                   {}, ...
    'centroidHeightAboveFloor', {} );

if isempty(rawCandidates), return; end

floorIds  = [rawCandidates.floorPlaneIdx];
wallSides = [rawCandidates.wallSide];
thetas    = [rawCandidates.theta];
groups    = unique([floorIds(:), wallSides(:)], 'rows');

for gi = 1:size(groups,1)
    fp  = groups(gi,1);
    ws  = groups(gi,2);
    sel = find(floorIds == fp & wallSides == ws);
    if isempty(sel), continue; end

    groupThetas = thetas(sel(:));
    groupThetas = groupThetas(:);
    thetaTol    = deg2rad(thetaTol_deg);
    used        = false(size(sel));

    for ii = 1:numel(sel)
        if used(ii), continue; end
        theta0  = rawCandidates(sel(ii)).theta;
        dtheta  = angle(exp(1i*(groupThetas - theta0)));
        nearby  = abs(dtheta) < thetaTol;
        clIdx   = sel(nearby);
        used(nearby) = true;

        nWall    = arrayfun(@(i) numel(rawCandidates(i).wallContactVertIdx), clIdx);
        [~,best] = max(nWall);
        rep      = rawCandidates(clIdx(best));

        allWall  = [];
        allFloor = [];
        for kk = 1:numel(clIdx)
            allWall  = union(allWall,  rawCandidates(clIdx(kk)).wallContactVertIdx);
            allFloor = union(allFloor, rawCandidates(clIdx(kk)).floorContactVertIdx);
        end
        rep.wallContactVertIdx  = allWall;
        rep.floorContactVertIdx = allFloor;
        rep.theta = meanCircularAngle([rawCandidates(clIdx).theta]);

        merged(end+1) = rep; %#ok<AGROW>
    end
end
end

function mu = meanCircularAngle(angles)
mu = atan2(mean(sin(angles)), mean(cos(angles)));
if mu < 0, mu = mu + 2*pi; end
end

function [stableIdxOut, QsOut, omegasOut, heightsOut, mergedGroupsOut, QsRawOut] = mergeByCSAValues( ...
    stableIdx, Qs, omegas, heights, omegaTol, hTol, ...
    allRestingPositions, Rchute, QsRaw) %#ok<INUSD>

cloudTol = hTol;
n        = numel(stableIdx);
clouds   = cell(n,1);

for si = 1:n
    pi_    = stableIdx(si);
    pos    = allRestingPositions(pi_);
    C0     = pos.centroidWorld(:)';
    vertsW = (Rchute * (pos.verticesWorld' - C0'))' + C0;
    clouds{si} = vertsW - mean(vertsW,1);
end

used         = false(n,1);
stableIdxOut = [];
QsOut        = [];
omegasOut    = [];
heightsOut   = [];
mergedGroupsOut = {};

haveQsRaw = (nargin >= 9) && ~isempty(QsRaw);
QsRawOut  = [];

for ii = 1:n
    if used(ii), continue; end
    groupMask     = false(n,1);
    groupMask(ii) = true;
    for jj = ii+1:n
        if used(jj), continue; end
        if cloudsMatch(clouds{ii}, clouds{jj}, cloudTol)
            groupMask(jj) = true;
        end
    end
    groupIdx        = find(groupMask);
    used(groupMask) = true;
    [~, bestLocal]  = max(Qs(groupIdx));
    bestGlobal      = groupIdx(bestLocal);
    stableIdxOut(end+1) = stableIdx(bestGlobal); %#ok<AGROW>
    QsOut(end+1)        = sum(Qs(groupIdx));     %#ok<AGROW>
    omegasOut(end+1)    = omegas(bestGlobal);    %#ok<AGROW>
    heightsOut(end+1)   = heights(bestGlobal);   %#ok<AGROW>

    if haveQsRaw
        QsRawOut(end+1) = sum(QsRaw(groupIdx)); %#ok<AGROW>
    end

    mergedGroupsOut{end+1} = stableIdx(groupIdx(:)'); %#ok<AGROW>
end

stableIdxOut = stableIdxOut(:);
QsOut        = QsOut(:);
omegasOut    = omegasOut(:);
heightsOut   = heightsOut(:);
if haveQsRaw
    QsRawOut = QsRawOut(:);
else
    QsRawOut = QsOut;  % fallback: behave like the old single-weight call
end

if numel(stableIdx) ~= numel(stableIdxOut)
    fprintf('  CSA cloud merge: %d -> %d position(s)\n', ...
        numel(stableIdx), numel(stableIdxOut));
end
end

function tf = cloudsMatch(A, B, tol)
if isempty(A) || isempty(B)
    tf = false; return;
end
tf = allPointsMatched(A, B, tol) && allPointsMatched(B, A, tol);
end

function tf = allPointsMatched(src, tgt, tol)
tf = true;
for k = 1:size(src,1)
    if min(vecnorm(tgt - src(k,:), 2, 2)) > tol
        tf = false; return;
    end
end
end

function [Qs, omegas, heights, projCentroids, hullProjs] = computeCSA( ...
    allRestingPositions, stableIdx, chullVertexIdx, Rchute, gWorld) %#ok<INUSD>

numStable     = numel(stableIdx);
Qs            = zeros(numStable,1);
omegas        = zeros(numStable,1);
heights       = zeros(numStable,1);
projCentroids = NaN(numStable,3);
hullProjs     = cell(numStable, 1);

gWorld = gWorld(:) / norm(gWorld(:));
rayDir = [0,0,-1];

for si = 1:numStable
    pi_  = stableIdx(si);
    pos  = allRestingPositions(pi_);

    C0     = pos.centroidWorld(:)';
    vertsW = (Rchute * (pos.verticesWorld' - C0'))' + C0;
    centW  = C0;

    supportIdx   = union(pos.floorContactVertIdx, pos.wallContactVertIdx);
    supportVerts = vertsW(supportIdx,:);
    if size(supportVerts,1) < 3, continue; end

    % ------------------------------------------------------------
    % HEIGHT: ray-cast centroid along gravity onto BOTH the floor
    % plane and the wall plane; the pose's height is whichever
    % plane the ray hits FIRST (smaller t), not the floor alone.
    % ------------------------------------------------------------
    floorContactIdx = pos.floorContactVertIdx;
    wallContactIdx  = pos.wallContactVertIdx;

    t_floor = Inf;
    if numel(floorContactIdx) >= 3
        floorVerts = vertsW(floorContactIdx, :);
        floorMean  = mean(floorVerts, 1);
        [~,~,Vf]   = svd(floorVerts - floorMean);
        floorNormal = Vf(:,3)';
        floorNormal = floorNormal / norm(floorNormal);

        denomF = dot(rayDir, floorNormal);
        if abs(denomF) > 1e-10
            tF = dot(floorMean - centW, floorNormal) / denomF;
            if tF >= 0
                t_floor = tF;
            end
        end
    end

    t_wall = Inf;
    if numel(wallContactIdx) >= 3
        wallVerts = vertsW(wallContactIdx, :);
        wallMean  = mean(wallVerts, 1);
        [~,~,Vw]  = svd(wallVerts - wallMean);
        wallNormal = Vw(:,3)';
        wallNormal = wallNormal / norm(wallNormal);

        denomW = dot(rayDir, wallNormal);
        if abs(denomW) > 1e-10
            tW = dot(wallMean - centW, wallNormal) / denomW;
            if tW >= 0
                t_wall = tW;
            end
        end
    end

    if isinf(t_floor) && isinf(t_wall), continue; end
    h = min(t_floor, t_wall);
    if h < 1e-10, continue; end

    centProjFloor       = centW + h * rayDir;
    heights(si)         = h;
    projCentroids(si,:) = centProjFloor;

    % ------------------------------------------------------------
    % (unchanged) support-hull solid angle calculation
    % ------------------------------------------------------------
    gravCoords  = supportVerts * gWorld;
    minGrav     = min(gravCoords);
    planePoint  = minGrav * gWorld';
    supportProj = zeros(size(supportVerts));
    for kk = 1:size(supportVerts,1)
        p_ = supportVerts(kk,:);
        d_ = dot((p_ - planePoint), gWorld');
        supportProj(kk,:) = p_ - d_ * gWorld';
    end
    supportProj = uniquetol(supportProj, 1e-5, 'ByRows', true);
    if size(supportProj,1) < 3, continue; end

    if abs(gWorld(1)) < 0.9, tmp2 = [1;0;0]; else, tmp2 = [0;1;0]; end
    e1b = cross(gWorld, tmp2); e1b = e1b/norm(e1b);
    e2b = cross(gWorld, e1b);  e2b = e2b/norm(e2b);
    support2D = [supportProj*e1b, supportProj*e2b];

    try
        hullIdx_c = convhull(support2D(:,1), support2D(:,2));
    catch
        continue
    end

    hull3D = supportProj(hullIdx_c(1:end-1),:);
    numH   = size(hull3D,1);
    if numH < 3, continue; end

    Vs = hull3D - centW;
    norms_Vs = vecnorm(Vs,2,2);
    if any(norms_Vs < 1e-10), continue; end
    Vs = Vs ./ norms_Vs;

    omega = 0;
    for k = 1:numH-2
        a_ = Vs(1,:); b_ = Vs(k+1,:); c_ = Vs(k+2,:);
        num_val = abs(dot(a_, cross(b_,c_)));
        den_val = 1 + dot(a_,b_) + dot(b_,c_) + dot(a_,c_);
        if abs(den_val) < 1e-14, continue; end
        omega = omega + 2*atan2(num_val, den_val);
    end
    omegas(si) = omega;
    Qs(si)     = omega / h;
    hullProjs{si} = hull3D;
end
end
