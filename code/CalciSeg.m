function [granules_labeled, summary_stats] = CalciSeg(stack, varargin)
% [granules_labeled, summary_stats] = CalciSeg(stack, Name, Value)
% segments a rows-by-columns-by-time stack into spatial regions (granules).
% Initial segmentation uses extrema-seeded Voronoi cells, correlation-based
% region growing, or reconstruction ICA. An existing label map can be resumed.
% Boundary reassignment and connectivity checks can use parallel workers.
%
% TOOLBOXES:
%   - Image Processing Toolbox        (imregionalmax/imregionalmin, imdilate, imerode,
%                                      bwulterode, bwdist, watershed, bwlabel,
%                                      regionprops, graythresh, imgradient,
%                                      medfilt2, multithresh, imfill, imgaussfilt,
%                                      imboxfilt, strel)
%   - Statistics and Machine Learning Toolbox (pca, rica, knnsearch)
%   - Parallel Computing Toolbox (optional for parallel execution)
%
% With 'keep_border' false, Voronoi initialization assigns all cells seeded
% on the image edge to background (ID 0). Refinement can reassign background.
%
% PRECONDITIONS (caller's responsibility):
%   - The stack should be MOTION-CORRECTED; residual motion inflates variance
%     and correlations and corrupts the segmentation.
%   - Correlation-based steps ('corr' projection/segmentation, refinement
%     'corr', pool_fragments) are sensitive to drift and photobleaching.
%     'linear' removes a linear trend. 'dff' rescales around a fixed baseline.
%
% METHOD SELECTION (rough guide):
%   - Sparse, blob-like sources ....... init_seg_method 'corr' (region growing)
%   - Dense fields / tiling ........... 'voronoi' (optionally 'spatial_smooth')
%   - Overlapping temporal sources .... 'rICA'
% The 'voronoi' default over-segments by design; increase 'spatial_smooth' or
% enforce 'limitPixCount' to obtain fewer, larger granules.
%
% EXAMPLES:
%   % Default Voronoi segmentation with summary statistics
%   [L, stats] = CalciSeg(stack);
%
%   % Correlation-based growing with baseline normalization and no progress display
%   L = CalciSeg(stack, 'init_seg_method','corr', 'detrend','dff', ...
%                'corr_thresh',0.8, 'verbose',false);
%
%   % Fewer, larger granules via smoothing and a size floor, 5 refinements
%   L = CalciSeg(stack, 'spatial_smooth',2, 'limitPixCount',[25 inf], 'n_rep',5);
%
%
% OPTIONAL INPUT PARAMETERS:
% Parameter name        Values and description
% =========================================================================
%
% 'projection_method' : Method for calculating the projection
% (string)              across time.
%                       Default: 'std'
%                       - 'std'    : Standard deviation projection
%                       - 'mean'   : Mean intensity projection
%                       - 'median' : Median intensity projection
%                       - 'max'    : Maximum intensity projection
%                       - 'min'    : Minimum intensity projection
%                       - 'pca'    : Principal component projection. For
%                                    this, "stack" is projected into PC
%                                    space.
%                       - 'corr'   : Pearson correlation of each pixel trace
%                                    with the mean trace of its neighbors.
%                       - 'none'   : Applicable when "stack" has no third
%                                    dimension or when to take the first
%                                    slice of the 3-D stack
%
% 'init_seg_method'   : Method for the initial segmentation.
% (string)              Default: 'voronoi'
%                       - 'voronoi' : Assign pixels to the nearest extrema seeds.
%                       - 'corr'    : Local growing based on correlation.
%                                     This will set 'projection_method' to
%                                     be 'corr'
%                       - 'rICA'    : Reconstruction independent component
%                                     analysis. Use parameter 'rICA_thresh'
%                                     to limit segmentation to regions with
%                                     high loadings.
%
% 'regmax_method'     : Method for determining how to identify local
% (string)              extrema in the intensity projection.
%                       Default: 'raw'
%                       - 'raw'      : Simply apply imregionalmax/-min on
%                                      the projection.
%                       - 'filtered' : Detect maxima after dilation and
%                                      minima after erosion in Voronoi mode. In correlation
%                                      and rICA modes, spatially average the seed map or
%                                      each loading map over a square neighborhood
%                                      of width 2*filter_size+1 (replicated edges).
%                       - 'both'     : Combine both above-mentioned
%                                      methods. Note, this is only
%                                      applicable when init_seg_method is
%                                      set to be 'voronoi'.
%
% 'n_rep'             : Number of iterations for refining the regions.
% (integer)             Default: 0
%
% 'refinement_method' : Measure to fine-tune granule assignment during
% (string)              refinement iterations.
%                       Default: 'rmse'
%                       - 'rmse'   : Root MEDIAN square error, i.e.
%                                    sqrt(median((a-b).^2)). Note this is
%                                    median- not mean-based despite the name.
%                       - 'rmedse' : Alias for 'rmse' (clearer name).
%                       - 'corr'   : Pearson correlation coefficient
%
% 'rICA_thresh'       : Number of Otsu thresholds used to exclude weakly
%                       associated pixels following reconstruction
%                       independent component analysis as initial
%                       segmentation.
%                       - integer 1-20 : apply threshold
%                       - 0            : skip
%
% 'limitPixCount'     : Limits the pixel area per granule
% (matrix or string)    Default:  [1 inf]
%                       - [a, b] : The minimum (a) and maximum (b) number
%                                  of pixels that can be assigned to a
%                                  granule. Note that, the minimum size
%                                  affects the filter size for the
%                                  regmax_method input.
%                       - 'auto' : An automatic assessment based on the
%                                  distribution of granule sizes before
%                                  refinement, using the 'auto_quantiles'
%                                  (default 5th/95th). The filter size for
%                                  regmax_method is set to 1.
%                       Regions below a are merged or dropped if isolated.
%                       Splitting is attempted above b. The minimum is enforced
%                       again after the final connectivity check.
%                       Max-size enforcement is best-effort: shapes that no
%                       strategy can split are left oversized and listed in
%                       summary_stats.oversized_ids.
%
% 'corr_thresh'       : Default seed and growth threshold for 'corr'
% (number)              initialization, and strict threshold for fragment pooling.
%                       Default: 0.85
%
% 'n_rICA'            : Number of features to extract during reconstruction
% (integer)             independent component analysis.
%                       Zero (default) selects a count using cumulative PCA
%                       variance and 'n_PCA'. Positive counts specify the
%                       number of learned features and may exceed the PCA rank.
%
% 'n_PCA'             : Percentage of the total variance of the data that
% (number)              PCA projection should retain, also used to choose the
%                       feature count for rICA when 'n_rICA' is zero.
%                       Range: (0, 100]. Default: 100.
%
% 'fillmissing'       : Fill NaNs and Infs along rows, columns, then time.
% (logical)             Each axis uses linear filling followed by nearest filling.
%                       Default: false
%
% 'resume_CalciSeg'    : Map with granule IDs from a previous segmentation
% (2-D matrix)          session. Supplying a map skips initial segmentation.
%                       Preprocessing and size enforcement still run.
%                       Default: []
%
% 'pool_fragments'    : Merge neighboring granules with correlated mean traces.
% (logical)             The algorithm calculates the Pearson correlation coefficient
%                       between neighboring granules and pools them if the
%                       correlation exceeds the threshold set using input
%                       'corr_thresh' (default: 0.85). The neighbor
%                       comparison is repeated until no granule pair
%                       exceeds the threshold. Note that this setting
%                       ignores the upper limit for pixel area per granule
%                       since it is applied after the refinement step.
%                       Default: false
%
% 'temporal_bin'      : Temporal downsampling factor for faster
% (integer)             segmentation. The stack is binned by averaging
%                       each complete group of N frames. If fewer than two
%                       complete bins remain, use all frames without binning.
%                       Otherwise, segmentation skips any final partial bin.
%                       Full-resolution time courses are used for the
%                       final summary statistics.
%                       Default: 1 (no binning)
%
% 'check_update'      : When true, check GitHub once per MATLAB session for a
% (logical)             newer CalciSeg version and warn if one is available.
%                       Requires an internet connection; fails silently when
%                       offline. No network call is made when false.
%                       Default: false
%
% 'verbose'           : When true, show progress bars and print status
% (logical)             messages (projection info, size estimates, solver
%                       progress). Set to false to suppress progress output.
%                       Warnings remain enabled. Progress bars are suppressed
%                       automatically when no graphical desktop is available.
%                       Default: true
%
% 'seed'              : Random-number generator seed used for the reproducible
% (integer)             rICA decomposition (init_seg_method='rICA'). The RNG
%                       state is saved and restored around the call, so the
%                       caller's global stream is left unchanged.
%                       Default: 1
%
% 'detrend'           : Temporal detrending applied to the data used for
% (string)              SEGMENTATION only (summary statistics stay on the raw
%                       traces). Choose linear drift removal or baseline normalization.
%                       - 'none'   : no detrending (default)
%                       - 'linear' : remove a per-pixel linear trend
%                       - 'dff'    : (F-F0)/F0, using each pixel's 10th percentile.
%                                    Nonpositive baselines are set to eps('single').
%
% 'spatial_smooth'    : Gaussian sigma (in pixels) applied to the projection
% (number)              before extrema detection. Larger values suppress
%                       spurious local maxima and reduce over-segmentation.
%                       Default: 0 (no smoothing)
%
% 'seed_on'           : Which extrema seed the Voronoi tessellation.
% (string)              'both' (default) | 'max' | 'min'.
%
% 'keep_border'       : When true, use only extrema seeds in Voronoi mode
% (logical)             and retain their cells, including cells touching the
%                       edge. When false, add edge seeds and discard all
%                       edge-seeded cells during initialization. Default: false.
%
% 'pixel_size'        : [dy dx] pixel size used to scale coordinates for
% (2-vector)            anisotropic Voronoi tessellation. Default: [1 1].
%
% 'seed_thresh'       : Correlation-growing seed/grow thresholds. When empty
% 'grow_thresh'         they default to 'corr_thresh'; set them separately to
% (number or [])        decouple where a region may start vs. how it grows.
%                       Default: [] (use corr_thresh)
%
% 'rica_sd_threshold' : With [] (default), orient components by spatial skewness.
% (number or [])        A nonnegative value compares mean excursions beyond
%                       the positive and negative thresholds in standard deviations.
%                       A tie keeps the fitted sign.
%
% 'auto_quantiles'    : [lo hi] quantiles used by limitPixCount 'auto'.
% (2-vector)            Default: [0.05 0.95]
%
% 'use_parallel'      : Allow parfor for the split check and refinement. When
% (logical)             false (or the workload is small) they run serially with
%                       identical results. Default: true
%
%
% OUTPUT
% Parameter name      Values and description
% =========================================================================
%
% granules_labeled    : Map with granule IDs. The size of granules_labeled
%                       corresponds to the first two dimensions of stack.
%
% summary_stats       : Summary statistics for each granule and additional
%                       information on the segmentation process saved as
%                       structures.
%                      - projection           : Preprocessed projection image.
%                                               rICA uses its dimensions only.
%                      - avgTCs               : Each granule's mean time
%                                               course
%                      - granule_Corr         : Average within-granule
%                                               correlation
%                      - granule_Avg_img      : Within-granule average
%                                               activity
%                      - granule_Std_img      : Each granule's standard
%                                               deviation over time
%                      - granule_Max_img      : Each granule's maximum
%                                               value over time
%                      - granule_Min_img      : Each granule's minimum
%                                               value over time
%                      - granule_Median_img   : Each granule's median
%                                               value over time
%                      - granule_Corr_img     : Within-granule correlation
%                      - active_region.map    : Binary map of active
%                                               regions based on Otsu's
%                                               method
%                      - active_region.method : Metric used for the activity mask.
%                                               PCA uses 'std'.
%                                               For 'none', the mask is empty.
%                      - quality.SNR          : Per-granule signal-to-noise
%                                               ratio (peak minus 10th percentile,
%                                               divided by difference-based noise).
%                      - quality.compactness  : Per-granule spatial
%                                               compactness (4*pi*area /
%                                               perimeter^2; 1 = circle)
%                      - area                 : Per-granule pixel count.
%                      - centroid             : Per-granule [x y] centroid.
%                      - bbox                 : Per-granule bounding box
%                                               [x y width height].
%                      - pixel_idx            : Cell array of each granule's
%                                               linear pixel indices.
%                      - oversized_ids        : IDs still exceeding the maximum
%                                               size (max enforcement is
%                                               best-effort; empty if none).
%                      - parameters           : Validated options and derived size
%                                               settings. refine_iterations records
%                                               completed passes when n_rep > 0. For rICA,
%                                               n_rICA_requested stores the input and
%                                               n_rICA the fitted feature count.
%
%
% Reproducibility: the rICA path is seeded via the 'seed' input and restores
% the caller's RNG state afterwards. All other steps are deterministic for a
% given input, though results may change across MATLAB/toolbox versions
% (e.g. knnsearch tie-breaking, watershed).
%
%
% Version 1.4 | 23-Nov-26 (R2024a)
%
% =========================================================================

% CalciSeg version tag (compared against GitHub only when 'check_update' is
% requested). Keep this in sync with the "% Version" line in the header.
CALCISEG_VERSION = '1.3';

% For higher-dimensional input, collapse singleton dimensions before validation.
if ndims(stack) > 3
    stack = squeeze(stack);
end
opt = validateInputs(stack, varargin);

% Fail fast with an actionable message if a required toolbox is unavailable.
checkDependencies();

% Optional, opt-in update check. Off by default so the function never makes a
% network call unless asked; runs at most once per MATLAB session and fails
% silently when offline.
if opt.check_update
    checkVersion(CALCISEG_VERSION, ...
        'https://raw.githubusercontent.com/YannickGuenzel/CalciSeg/main/code/CalciSeg.m');
end

% Reduce precision of stack to be single
if ~isa(stack, 'single')
    stack = single(stack);
end%if not single precision

% Check for missing data (single pass over the stack)
if ~all(isfinite(stack), 'all')
    if opt.fillmissing
        stack(isinf(stack)) = NaN;
        % Fill along rows, then columns, then time, with a nearest-value fallback.
        stack = fillmissing(stack, "linear", 1); stack = fillmissing(stack, "nearest", 1);
        stack = fillmissing(stack, "linear", 2); stack = fillmissing(stack, "nearest", 2);
        if ndims(stack) == 3
            stack = fillmissing(stack, "linear", 3);
            stack = fillmissing(stack, "nearest", 3);
        end
        % Reject values still missing after every available axis has been tried.
        if any(~isfinite(stack), 'all')
            error('CalciSeg:fillFailed', ...
                'fillmissing could not remove all NaN/Inf values (e.g. a completely missing pixel or frame). Please clean the data manually.');
        end
    else
        error('CalciSeg:missingData', 'The stack contains NaN or inf values. Please check your data or set the input "fillmissing" to true.')
    end%if fill
end%if missing data

% Temporal downsampling for faster segmentation
reshaped_stack_fullres = [];
if opt.temporal_bin > 1
    [x_s, y_s, t_s] = size(stack);
    t_binned = floor(t_s / opt.temporal_bin);
    % Retain at least two time samples for segmentation.
    if t_binned < 2
        warning('CalciSeg:temporalBin', ...
            'Temporal binning (bin=%d, frames=%d) would leave %d frame(s). Disabling binning.', ...
            opt.temporal_bin, t_s, t_binned);
        opt.temporal_bin = 1;
    else
        reshaped_stack_fullres = reshape(stack, [x_s*y_s, t_s]);
        stack = reshape(stack(:,:,1:t_binned*opt.temporal_bin), [x_s, y_s, opt.temporal_bin, t_binned]);
        % Preserve both spatial axes, including singleton rows or columns.
        stack = reshape(mean(stack, 3), [x_s, y_s, t_binned]);
        if opt.verbose
            disp(['Temporal binning: ', num2str(t_s), ' -> ', num2str(t_binned), ' frames (bin factor: ', num2str(opt.temporal_bin), ')'])
        end
    end
end%if temporal binning

% Define size for local neighborhood
if isnumeric(opt.limitPixCount)
    opt.filter_size = ceil(sqrt(opt.limitPixCount(1)/pi)) + double(opt.limitPixCount(1)==0);
    % Make sure it is an odd number
    opt.filter_size = opt.filter_size + double((mod(opt.filter_size,2)==0));
elseif strcmp(opt.limitPixCount, 'auto')
    opt.filter_size = 1;
end

% Initial preprocessing
[reshaped_stack, projection, reshaped_stack_raw] = initialPreprocessing(stack, opt);

% Perform initial segmentation (or use a previous segmentation result)
if ~isempty(opt.resume_CalciSeg)
    % Compact labels in their original numeric type before any conversion.
    [ids, ~, group] = unique(opt.resume_CalciSeg(:));
    granules_labeled = reshape(double(group) - double(ids(1) == 0), size(opt.resume_CalciSeg));
else
    [granules_labeled, opt] = initialSegmentation(reshaped_stack, projection, opt);
end

% Automatically estimate min/max granule size from the initial segmentation
if (isstring(opt.limitPixCount) || ischar(opt.limitPixCount)) && strcmp(opt.limitPixCount, 'auto')
    % Get pixel count for each region (exclude background ID 0)
    [granuleList, ~, C] = unique(granules_labeled);
    C = histcounts(C, 1:length(granuleList)+1);
    C = [C(:), granuleList(:)];
    C = C(C(:,2) ~= 0, :);
    if isempty(C)
        % No granules were found - quantile() would return NaN and silently
        % disable both size constraints. Fall back to an unconstrained range.
        warning('CalciSeg:autoLimit', ...
            'No granules found; cannot auto-estimate size limits. Using [1 Inf].');
        opt.limitPixCount = [1 inf];
    else
        opt.limitPixCount = [floor(quantile(C(:,1), opt.auto_quantiles(1))), floor(quantile(C(:,1), opt.auto_quantiles(2)))];
        if opt.verbose
            disp(['automatic assessment of the min./max. granule size: [', num2str(opt.limitPixCount),']'])
        end
    end
end%if auto min size

% Enforce granule size constraints (split oversized, merge undersized).
% Note: removeHugeRegions runs first because splitting may create tiny
% fragments that removeTinyRegions then merges.
granules_labeled = removeHugeRegions(granules_labeled, opt.limitPixCount(2));
granules_labeled = removeTinyRegions(reshaped_stack, granules_labeled, opt);
if opt.n_rep>0
    [granules_labeled, opt.refine_iterations] = refineSegmentation(reshaped_stack, granules_labeled, opt);
end%if refine

% Optionally pool fragmented granules based on neighbor correlation
if opt.pool_fragments
    granules_labeled = poolFragmentedGranules(reshaped_stack, granules_labeled, opt);
end%if pool fragments

% Normalize connectivity: guarantee every output granule is a single
% 4-connected component regardless of which path produced it (e.g. watershed
% splitting or a resumed map). No-op when the invariant already holds.
granules_labeled = checkSplits(granules_labeled, opt.use_parallel);

% Connectivity normalization can expose components below the size minimum.
if opt.limitPixCount(1) > 1
    granules_labeled = removeTinyRegions(reshaped_stack, granules_labeled, opt);
end

% Relabel granules with contiguous IDs ordered by descending size
granules_labeled = finalRefinement(granules_labeled);

% Sanity hint: flag likely over-segmentation (only when verbose).
if opt.verbose
    nG_final = numel(unique(granules_labeled(granules_labeled>0)));
    if nG_final > numel(granules_labeled)/4
        warning('CalciSeg:manyGranules', ...
            ['Segmentation produced %d granules (> 1 per 4 pixels), which often ', ...
            'indicates over-segmentation. Consider increasing "spatial_smooth" or ', ...
            'the minimum "limitPixCount".'], nG_final);
    end
end

% Use unbinned fluorescence for statistics, including any imputed samples.
% Detrending and PCA affect the segmentation data only.
if nargout>1
    if opt.temporal_bin > 1
        summary_stats = regionStats(reshaped_stack_fullres, granules_labeled, opt);
    else
        summary_stats = regionStats(reshaped_stack_raw, granules_labeled, opt);
    end
    % Keep the projection at the temporal resolution used for segmentation.
    summary_stats.projection = projection;
    % Include resolved size limits and, for rICA runs, the fitted feature count.
    summary_stats.parameters = opt;
end%if summary stats

end%FCN:CalciSeg

% -------------------------------------------------------------------------

function checkVersion(localVersion, url)
% Compare the local version tag against the "% Version" line of the file
% published on GitHub. Warn only when the remote version is newer. The
% persistent flag limits the check to one attempt while the function is
% loaded, including when the network request fails.
persistent alreadyChecked
if ~isempty(alreadyChecked)
    return
end
alreadyChecked = true;
try
    remoteText = webread(url, weboptions('Timeout', 2));
    tok = regexp(remoteText, 'Version\s+([0-9]+\.[0-9]+)', 'tokens', 'once');
    % Only warn when the remote version is strictly NEWER than the local one
    % (a numeric compare avoids "1.10 < 1.9" string-ordering pitfalls and
    % never nags when the local copy is ahead of the published one).
    if ~isempty(tok) && versionLessThan(localVersion, tok{1})
        warning('CalciSeg:version', ...
            ['A newer CalciSeg version (%s) is available; you are running %s.\n', ...
            'Latest version can be found at: https://github.com/YannickGuenzel/CalciSeg'], ...
            tok{1}, localVersion);
    end
catch
    % No internet connection or GitHub unreachable - skip silently.
end
end%FCN:checkVersion

function tf = versionLessThan(a, b)
% True if version string a (e.g. '1.2') is strictly less than b ('1.3').
va = sscanf(a, '%d.%d')';
vb = sscanf(b, '%d.%d')';
va(end+1:2) = 0; vb(end+1:2) = 0;   % pad to [major minor]
tf = (va(1) < vb(1)) || (va(1) == vb(1) && va(2) < vb(2));
end%FCN:versionLessThan

function checkDependencies()
% Verify the required toolbox functions are available and error early with an
% actionable message if a toolbox is missing (rather than failing deep inside
% a long run). The Parallel Computing Toolbox is optional (parfor degrades to
% serial), so it is not required here.
reqIP = {'imregionalmax','imregionalmin','bwlabel','bwulterode','bwdist', ...
    'watershed','regionprops','imgradient','imgaussfilt','strel'};
reqStats = {'pca','knnsearch','prctile','quantile'};
missIP = reqIP(cellfun(@(f) isempty(which(f)), reqIP));
missStats = reqStats(cellfun(@(f) isempty(which(f)), reqStats));
if ~isempty(missIP)
    error('CalciSeg:missingToolbox', ...
        'CalciSeg requires the Image Processing Toolbox (missing: %s).', strjoin(missIP, ', '));
end
if ~isempty(missStats)
    error('CalciSeg:missingToolbox', ...
        'CalciSeg requires the Statistics and Machine Learning Toolbox (missing: %s).', strjoin(missStats, ', '));
end
end%FCN:checkDependencies

function closeIfValid(h)
% Close a graphics handle if it still exists (used by onCleanup guards so a
% progress bar is removed on both normal and error exits).
if ~isempty(h) && ishandle(h)
    close(h)
end
end%FCN:closeIfValid

function opt = validateInputs(stack, varargin)
% Parse and validate all optional Name-Value input arguments and return
% them in the options struct 'opt'. Also validates the stack dimensions
% and enforces parameter dependencies.

% Set all parameters to their default value
opt.projection_method = 'std';
opt.init_seg_method = 'voronoi';
opt.regmax_method = 'raw';
opt.n_rep = 0;
opt.refinement_method = 'rmse';
opt.rICA_thresh = 0;
opt.limitPixCount = [1 inf];
opt.corr_thresh = 0.85;
opt.n_rICA = 0;
opt.n_PCA = 100;
opt.fillmissing = false;
opt.resume_CalciSeg = [];
opt.pool_fragments = false;
opt.temporal_bin = 1;
opt.check_update = false;
opt.verbose = true;
opt.seed = 1;
opt.detrend = 'none';          % 'none' | 'linear' | 'dff' (segmentation data)
opt.seed_on = 'both';          % Voronoi seeds: 'both' | 'max' | 'min'
opt.spatial_smooth = 0;        % Gaussian sigma applied to the projection (0=off)
opt.keep_border = false;       % keep edge-seeded granules in Voronoi mode
opt.seed_thresh = [];          % correlation-growing seed threshold ([]=corr_thresh)
opt.grow_thresh = [];          % correlation-growing growth threshold ([]=corr_thresh)
opt.pixel_size = [1 1];        % [dy dx] pixel size for anisotropic Voronoi
opt.rica_sd_threshold = [];    % skewness by default, tail comparison when supplied
opt.auto_quantiles = [0.05 0.95]; % quantiles for limitPixCount 'auto'
opt.use_parallel = true;       % allow parfor (falls back to serial otherwise)

% Check optional input arguments (Name-Value pairs). Track which names the
% caller supplied so inert parameters can be flagged afterwards.
supplied = {};
varargin = varargin{1};
if ~isempty(varargin)
    if mod(length(varargin),2)>0
        error('CalciSeg:badNargin', 'Wrong number of arguments (optional inputs must be Name-Value pairs).')
    else
        % Iterate over all input arguments
        for iArg = 1:2:length(varargin)
            if ischar(varargin{iArg}) || isstring(varargin{iArg})
                supplied{end+1} = char(varargin{iArg}); %#ok<AGROW>
            end
            switch varargin{iArg}
                % Projection method -------------------------------------------
                case 'projection_method'
                    valid_projection_methods = {'std', 'mean', 'median', 'max', 'min', 'pca', 'corr', 'none'};
                    if ~ismember(varargin{iArg+1}, valid_projection_methods)
                        error('CalciSeg:invalidInput', 'Invalid "projection_method". Valid options are: %s.', strjoin(valid_projection_methods, ', '));
                    else
                        opt.projection_method = varargin{iArg+1};
                    end
                    % Initial segmentation method ---------------------------------
                case 'init_seg_method'
                    valid_init_seg_methods = {'voronoi', 'corr', 'rICA'};
                    if ~ismember(varargin{iArg+1}, valid_init_seg_methods)
                        error('CalciSeg:invalidInput', 'Invalid "init_seg_method". Valid options are: %s.', strjoin(valid_init_seg_methods, ', '));
                    else
                        opt.init_seg_method = varargin{iArg+1};
                    end
                    % Thresholding after rICA -----------------------------
                case 'rICA_thresh'
                    if ~isnumeric(varargin{iArg+1}) || ~isreal(varargin{iArg+1}) || ~isscalar(varargin{iArg+1}) || ~isfinite(varargin{iArg+1}) || varargin{iArg+1} < 0 || varargin{iArg+1} > 20 || mod(varargin{iArg+1},1) ~= 0
                        error('CalciSeg:invalidInput', 'Input "rICA_thresh" must be an integer between 0 and 20.');
                    else
                        opt.rICA_thresh = double(varargin{iArg+1});
                    end
                    % Method for smoothing regional maxima ----------------
                case 'regmax_method'
                    valid_regmax_methods = {'raw', 'filtered', 'both'};
                    if ~ismember(varargin{iArg+1}, valid_regmax_methods)
                        error('CalciSeg:invalidInput', 'Invalid "regmax_method". Valid options are: %s.', strjoin(valid_regmax_methods, ', '));
                    else
                        opt.regmax_method = varargin{iArg+1};
                    end
                    % Number of refinement iterations ----------------------------
                case 'n_rep'
                    if ~isnumeric(varargin{iArg+1}) || ~isreal(varargin{iArg+1}) || ~isscalar(varargin{iArg+1}) || varargin{iArg+1} < 0 || mod(varargin{iArg+1},1) ~= 0
                        error('CalciSeg:invalidInput', 'Input "n_rep" must be a non-negative scalar integer.');
                    else
                        opt.n_rep = double(varargin{iArg+1});
                    end
                    % Distance metric during refinement -------------------
                case 'refinement_method'
                    % 'rmse' here is root-MEDIAN-square error; 'rmedse' is
                    % accepted as a clearer alias and normalized to 'rmse'.
                    valid_refinement_methods = {'corr', 'rmse', 'rmedse'};
                    if ~ismember(varargin{iArg+1}, valid_refinement_methods)
                        error('CalciSeg:invalidInput', 'Invalid "refinement_method". Valid options are: %s.', strjoin(valid_refinement_methods, ', '));
                    elseif strcmp(varargin{iArg+1}, 'rmedse')
                        opt.refinement_method = 'rmse';
                    else
                        opt.refinement_method = varargin{iArg+1};
                    end
                    % Minimum and maximum number of pixels per region -------------
                case 'limitPixCount'
                    val = varargin{iArg+1};
                    if ~isnumeric(val) && ~((ischar(val) && isrow(val) || isstring(val) && isscalar(val)) && strcmp(val, 'auto'))
                        error('CalciSeg:invalidInput', 'Input "limitPixCount" must be a non-negative pair [min max] or "auto".');
                    elseif isnumeric(val)
                        if ~isreal(val) || ~isvector(val) || numel(val) ~= 2 || isnan(val(2)) || any(val < 0) || val(1) > val(2) || ~isfinite(val(1)) || any(mod(val(isfinite(val)),1) ~= 0)
                            error('CalciSeg:invalidInput', 'Input "limitPixCount" must be [min max] with non-negative integers and finite min (max may be inf).');
                        end
                        val = double(val(:)');
                    end
                    opt.limitPixCount = val;
                    % Correlation threshold ---------------------------------------
                case 'corr_thresh'
                    if ~isnumeric(varargin{iArg+1}) || ~isreal(varargin{iArg+1}) || ~isscalar(varargin{iArg+1}) || ~isfinite(varargin{iArg+1}) || varargin{iArg+1} < 0 || varargin{iArg+1} > 1
                        error('CalciSeg:invalidInput', 'Input "corr_thresh" must be a scalar number between 0 and 1.');
                    else
                        opt.corr_thresh = double(varargin{iArg+1});
                    end
                    % Number of independent components ----------------------------
                case 'n_rICA'
                    if ~isnumeric(varargin{iArg+1}) || ~isreal(varargin{iArg+1}) || ~isscalar(varargin{iArg+1}) || varargin{iArg+1} < 0 || mod(varargin{iArg+1},1) ~= 0
                        error('CalciSeg:invalidInput', 'Input "n_rICA" must be a non-negative scalar integer.');
                    else
                        opt.n_rICA = double(varargin{iArg+1});
                    end
                case 'n_PCA'
                    if ~isnumeric(varargin{iArg+1}) || ~isreal(varargin{iArg+1}) || ~isscalar(varargin{iArg+1}) || ~isfinite(varargin{iArg+1}) || varargin{iArg+1} <= 0 || varargin{iArg+1} > 100
                        error('CalciSeg:invalidInput', 'Input "n_PCA" must be a scalar number larger than 0 and not exceeding 100.');
                    else
                        opt.n_PCA = double(varargin{iArg+1});
                    end
                case 'fillmissing'
                    if ~islogical(varargin{iArg+1}) || ~isscalar(varargin{iArg+1})
                        error('CalciSeg:invalidInput', 'Input "fillmissing" must be a scalar logical (true or false).');
                    else
                        opt.fillmissing = varargin{iArg+1};
                    end
                case 'resume_CalciSeg'
                    val = varargin{iArg+1};
                    if isnumeric(val) && isequal(size(val), [0 0])
                        opt.resume_CalciSeg = [];
                    elseif ~isnumeric(val) || ~isreal(val) || ~ismatrix(val) || size(val,1) ~= size(stack,1) || size(val,2) ~= size(stack,2)
                        error('CalciSeg:invalidInput', 'Input "resume_CalciSeg" must be a 2-D matrix with the same spatial size as "stack".');
                    elseif any(val(:) < 0) || any(~isfinite(val(:))) || any(mod(val(:),1) ~= 0)
                        error('CalciSeg:invalidInput', 'Input "resume_CalciSeg" must contain only finite, non-negative integers.');
                    else
                        opt.resume_CalciSeg = val;
                    end
                case 'pool_fragments'
                    if ~islogical(varargin{iArg+1}) || ~isscalar(varargin{iArg+1})
                        error('CalciSeg:invalidInput', 'Input "pool_fragments" must be a scalar logical (true or false).');
                    else
                        opt.pool_fragments = varargin{iArg+1};
                    end
                case 'temporal_bin'
                    if ~isnumeric(varargin{iArg+1}) || ~isreal(varargin{iArg+1}) || ~isscalar(varargin{iArg+1}) || varargin{iArg+1} < 1 || mod(varargin{iArg+1},1) ~= 0
                        error('CalciSeg:invalidInput', 'Input "temporal_bin" must be a positive integer.');
                    else
                        opt.temporal_bin = double(varargin{iArg+1});
                    end
                case 'check_update'
                    if ~islogical(varargin{iArg+1}) || ~isscalar(varargin{iArg+1})
                        error('CalciSeg:invalidInput', 'Input "check_update" must be a scalar logical (true or false).');
                    else
                        opt.check_update = varargin{iArg+1};
                    end
                case 'verbose'
                    if ~islogical(varargin{iArg+1}) || ~isscalar(varargin{iArg+1})
                        error('CalciSeg:invalidInput', 'Input "verbose" must be a scalar logical (true or false).');
                    else
                        opt.verbose = varargin{iArg+1};
                    end
                case 'seed'
                    if ~isnumeric(varargin{iArg+1}) || ~isreal(varargin{iArg+1}) || ~isscalar(varargin{iArg+1}) || ~isfinite(varargin{iArg+1}) || varargin{iArg+1} < 0 || mod(varargin{iArg+1},1) ~= 0 || varargin{iArg+1} >= 2^32
                        error('CalciSeg:invalidInput', 'Input "seed" must be an integer in [0, 2^32-1].');
                    else
                        opt.seed = double(varargin{iArg+1});
                    end
                case 'detrend'
                    valid_detrend = {'none','linear','dff'};
                    if ~ismember(varargin{iArg+1}, valid_detrend)
                        error('CalciSeg:invalidInput', 'Input "detrend" must be one of: %s.', strjoin(valid_detrend, ', '));
                    else
                        opt.detrend = varargin{iArg+1};
                    end
                case 'seed_on'
                    valid_seed_on = {'both','max','min'};
                    if ~ismember(varargin{iArg+1}, valid_seed_on)
                        error('CalciSeg:invalidInput', 'Input "seed_on" must be one of: %s.', strjoin(valid_seed_on, ', '));
                    else
                        opt.seed_on = varargin{iArg+1};
                    end
                case 'spatial_smooth'
                    if ~isnumeric(varargin{iArg+1}) || ~isreal(varargin{iArg+1}) || ~isscalar(varargin{iArg+1}) || varargin{iArg+1} < 0 || ~isfinite(varargin{iArg+1})
                        error('CalciSeg:invalidInput', 'Input "spatial_smooth" must be a non-negative finite scalar (Gaussian sigma).');
                    else
                        opt.spatial_smooth = double(varargin{iArg+1});
                    end
                case 'keep_border'
                    if ~islogical(varargin{iArg+1}) || ~isscalar(varargin{iArg+1})
                        error('CalciSeg:invalidInput', 'Input "keep_border" must be a scalar logical (true or false).');
                    else
                        opt.keep_border = varargin{iArg+1};
                    end
                case 'seed_thresh'
                    v = varargin{iArg+1};
                    if ~isnumeric(v) || ~isreal(v) || (~isempty(v) && (~isscalar(v) || ~isfinite(v) || v < 0 || v > 1))
                        error('CalciSeg:invalidInput', 'Input "seed_thresh" must be [] or a scalar between 0 and 1.');
                    else
                        opt.seed_thresh = v;
                    end
                case 'grow_thresh'
                    v = varargin{iArg+1};
                    if ~isnumeric(v) || ~isreal(v) || (~isempty(v) && (~isscalar(v) || ~isfinite(v) || v < 0 || v > 1))
                        error('CalciSeg:invalidInput', 'Input "grow_thresh" must be [] or a scalar between 0 and 1.');
                    else
                        opt.grow_thresh = v;
                    end
                case 'pixel_size'
                    v = varargin{iArg+1};
                    if ~isnumeric(v) || ~isreal(v) || ~isvector(v) || numel(v) ~= 2 || any(v <= 0) || any(~isfinite(v))
                        error('CalciSeg:invalidInput', 'Input "pixel_size" must be a 2-element positive vector [dy dx].');
                    else
                        opt.pixel_size = double(v(:)');
                    end
                case 'rica_sd_threshold'
                    v = varargin{iArg+1};
                    if ~isnumeric(v) || ~isreal(v) || (~isempty(v) && (~isscalar(v) || ~isfinite(v) || v < 0))
                        error('CalciSeg:invalidInput', 'Input "rica_sd_threshold" must be [] or a non-negative finite scalar.');
                    end
                    opt.rica_sd_threshold = double(v);
                case 'auto_quantiles'
                    v = varargin{iArg+1};
                    if ~isnumeric(v) || ~isreal(v) || ~isvector(v) || numel(v) ~= 2 || any(~isfinite(v)) || any(v < 0) || any(v > 1) || v(1) >= v(2)
                        error('CalciSeg:invalidInput', 'Input "auto_quantiles" must be [lo hi] with 0<=lo<hi<=1.');
                    else
                        opt.auto_quantiles = double(v(:)');
                    end
                case 'use_parallel'
                    if ~islogical(varargin{iArg+1}) || ~isscalar(varargin{iArg+1})
                        error('CalciSeg:invalidInput', 'Input "use_parallel" must be a scalar logical (true or false).');
                    else
                        opt.use_parallel = varargin{iArg+1};
                    end
                otherwise
                    if ischar(varargin{iArg}) || (isstring(varargin{iArg}) && isscalar(varargin{iArg}))
                        error('CalciSeg:unknownArgument', 'Unknown argument "%s".', varargin{iArg});
                    else
                        error('CalciSeg:invalidArgumentName', ...
                            'Expected a parameter name (char/string) at argument position %d, but received a %s.', ...
                            iArg, class(varargin{iArg}));
                    end
            end%switch
        end%iArg
    end%if mod
end%if custom arguments

% Minimum size cannot exceed the number of pixels
if isnumeric(opt.limitPixCount) && opt.limitPixCount(1) > size(stack,1)*size(stack,2)
    warning('CalciSeg:sizeLimit', 'limitPixCount(1)=%d exceeds the pixel count (%d); no granule can satisfy it.', opt.limitPixCount(1), size(stack,1)*size(stack,2));
end

% Validate 'stack' dimensions and enforce parameter dependencies
if ndims(stack) > 3
    error('CalciSeg:invalidStack', 'Number of dimensions of the input "stack" cannot exceed 3.')
elseif ismatrix(stack)
    warning('CalciSeg:noTemporalComponent', 'Input "stack" has no temporal component. Inputs are adjusted accordingly (projection_method="none"; init_seg_method="voronoi"; refinement_method="rmse"; temporal_bin=1; detrend="none").');
    opt.projection_method = 'none';
    opt.init_seg_method = 'voronoi';
    opt.refinement_method = 'rmse';
    opt.temporal_bin = 1;
    opt.detrend = 'none';
end%if invalid 'stack' input
% For the region-growing approach, both projection_method and
% opt.init_seg_method must be set to 'corr'
if isempty(opt.resume_CalciSeg) && strcmp(opt.init_seg_method, 'corr')
    if ~strcmp(opt.projection_method, 'corr')
        warning('CalciSeg:forcedProjection', 'When the input "init_seg_method" is set to "corr", the input "projection_method" must be set to "corr" as well. We fixed this.')
        opt.projection_method = 'corr';
    end%if
end%if
% 'both' regmax_method is only applicable with Voronoi tessellation
if strcmp(opt.regmax_method, 'both') && ~strcmp(opt.init_seg_method, 'voronoi')
    warning('CalciSeg:regmaxReset', 'regmax_method "both" is only applicable with init_seg_method "voronoi". Setting to "raw".');
    opt.regmax_method = 'raw';
end%if
% Check explicit options against the effective computation after method resets.
initializing = isempty(opt.resume_CalciSeg);
growing = initializing && strcmp(opt.init_seg_method, 'corr');
fittingICA = initializing && strcmp(opt.init_seg_method, 'rICA');
usingPCA = strcmp(opt.projection_method, 'pca') || (fittingICA && opt.n_rICA == 0);
if ismember('n_PCA', supplied) && ~usingPCA
    warning('CalciSeg:inertParameter', '"n_PCA" is unused without PCA projection or automatic rICA feature selection.');
end
if any(ismember({'n_rICA','rICA_thresh','rica_sd_threshold','seed'}, supplied)) && ~fittingICA
    warning('CalciSeg:inertParameter', 'The supplied rICA options require rICA initialization.');
end
usesThreshold = opt.pool_fragments || (growing && (isempty(opt.seed_thresh) || isempty(opt.grow_thresh)));
if ismember('corr_thresh', supplied) && ~usesThreshold
    warning('CalciSeg:inertParameter', '"corr_thresh" is unused with these growth and pooling settings.');
end
if any(ismember({'seed_thresh','grow_thresh'}, supplied)) && ~growing
    warning('CalciSeg:inertParameter', 'Seed and growth thresholds require correlation initialization.');
end
if any(ismember({'keep_border','seed_on','pixel_size'}, supplied)) && ~(initializing && strcmp(opt.init_seg_method, 'voronoi'))
    warning('CalciSeg:inertParameter', 'Voronoi seed options require Voronoi initialization.');
end
end%FCN:validateInputs

% -------------------------------------------------------------------------

function [reshaped_stack, projection, reshaped_stack_raw] = initialPreprocessing(stack, opt)
% Compute the temporal projection (collapsing the time dimension into a 2D
% image) and reshape the stack into a 2D pixels-by-time matrix.
% reshaped_stack   : data used for segmentation/refinement (may be
%                    transformed, e.g. PC scores when projection_method='pca').
% reshaped_stack_raw: pixels by time before detrending or projection, at the
%                    input resolution of this helper (possibly binned).

% Get the size of the stack
[x, y, t] = size(stack);

% Capture the raw pixels-by-time matrix BEFORE any projection method
% modifies the stack (only 'pca' replaces it with PC scores).
reshaped_stack_raw = reshape(stack, [x*y, t]);
% Optional detrending or baseline normalization of segmentation data.
if ~strcmp(opt.detrend, 'none')
    stack = detrendStack(stack, opt.detrend);
end

% Calculate the projection based on the specified method
switch opt.projection_method
    case 'std'
        projection = std(stack, 0, 3, 'omitnan');
    case 'mean'
        projection = mean(stack, 3, 'omitnan');
    case 'median'
        projection = median(stack, 3, 'omitnan');
    case 'max'
        projection = max(stack, [], 3, 'omitnan');
    case 'min'
        projection = min(stack, [], 3, 'omitnan');
    case 'pca'
        % Project into PC space and keep components explaining n_PCA% variance
        [~, stack, ~, ~, explained] = pca(reshape(stack, [x*y, t]));
        t_pca = find(cumsum(explained) >= opt.n_PCA, 1);
        if isempty(t_pca)
            % Threshold not reached (float imprecision or fewer components
            % than time frames) — keep all available components
            t = size(stack, 2);
        else
            t = t_pca;
        end%if variance threshold reached
        stack = reshape(stack(:,1:t), [x, y, t]);
        % Variance-weighted projection of retained components
        weight_sum = sum(explained(1:t));
        if weight_sum > 0
            weights = explained(1:t) / weight_sum;
        else
            % Degenerate: zero explained variance — use equal weights
            weights = ones(t, 1) / t;
        end
        % Weighted sum of the retained PC score images (vectorized)
        projection = reshape(reshape(stack, [x*y, t]) * weights(:), [x, y]);
        if opt.verbose
            disp(['PCA projection | ', num2str(t), ' components | Explained variance: ', num2str(round(sum(explained(1:t)),2)), '%'])
        end
    case 'corr'
        projection = correlationImage(stack, opt);
        if isempty(opt.resume_CalciSeg) && strcmp(opt.init_seg_method, 'corr') && strcmp(opt.regmax_method, 'filtered')
            projection = imboxfilt(projection, 2*opt.filter_size + 1, 'Padding', 'replicate');
        end
    case 'none'
        projection = stack(:,:,1);
end%switch projection method

% Optional spatial smoothing of the projection to suppress spurious local
% extrema (helps counter over-segmentation on smooth maps).
if opt.spatial_smooth > 0 && ~isempty(projection)
    projection = imgaussfilt(projection, opt.spatial_smooth);
end

% Reshape the stack to a 2D matrix for further processing
reshaped_stack = reshape(stack, [x*y, t]);

end%FCN:initialPreprocessing

% -------------------------------------------------------------------------

function stack = detrendStack(stack, method)
% Temporally detrend each pixel's trace. 'linear' removes a per-pixel linear
% trend; 'dff' converts to Delta F / F using a per-pixel 10th-percentile
% baseline. Operates along the time dimension (3rd).
[x, y, t] = size(stack);
X = reshape(stack, [x*y, t]);
switch method
    case 'linear'
        % detrend works column-wise; transpose so each pixel is a column
        X = detrend(double(X)', 'linear')';
        stack = reshape(single(X), [x, y, t]);
    case 'dff'
        F0 = prctile(X, 10, 2);
        F0(F0 <= 0) = eps('single');   % guard non-positive baselines
        X = (X - F0) ./ F0;
        thresh = reshape(F0, [x y]);
        thresh = thresh > (0.5*min(multithresh(thresh, 20)));
        thresh = imfill(thresh, 'holes');
        stack = reshape(single(X), [x, y, t]) .* thresh;
end%switch

end%FCN:detrendStack

% -------------------------------------------------------------------------

function [granules_labeled, opt] = initialSegmentation(reshaped_stack, projection, opt)
% Perform the initial spatial segmentation using one of three methods:
% Voronoi tessellation, correlation-based region growing, or rICA.
switch opt.init_seg_method
    case 'voronoi'
        % Get local maxima and minima
        [regional_maxima, regional_minima] = identifyLocalExtrema(projection, opt);
        % Assign each pixel to its nearest seed in scaled spatial coordinates.
        granules_labeled = triang_extrema(regional_maxima, regional_minima, opt);
    case 'corr'
        granules_labeled = growing_regions(projection, reshaped_stack, opt);
    case 'rICA'
        opt.n_rICA_requested = opt.n_rICA;
        [granules_labeled, opt.n_rICA] = rICA_segmentation(projection, reshaped_stack, opt);
        % Check for splits
        granules_labeled = checkSplits(granules_labeled, opt.use_parallel);
end%switch opt.init_seg_method
end%FCN:initialSegmentation

% -------------------------------------------------------------------------

function granules_labeled = triang_extrema(regional_maxima, regional_minima, opt)
% Segment the image via Voronoi tessellation seeded by local extrema.
% Each foreground seed pixel from ultimate erosion is a separate site.
% opt.seed_on chooses the extrema and opt.pixel_size scales the distances.
% With keep_border false, add sites on the edge and discard their cells.

% Apply ultimate erosion to the selected extrema masks.
switch opt.seed_on
    case 'both'
        ultimateErosion = bwulterode(regional_maxima>0) + bwulterode(regional_minima>0);
    case 'max'
        ultimateErosion = bwulterode(regional_maxima>0);
    case 'min'
        ultimateErosion = bwulterode(regional_minima>0);
end
ultimateErosion = ultimateErosion > 0;
% Force edge pixels to be seeds only when border regions are to be dropped
if ~opt.keep_border
    ultimateErosion(1,:) = 1; ultimateErosion(:,1) = 1;
    ultimateErosion(end,:) = 1; ultimateErosion(:,end) = 1;
end
% Seed coordinates as (row, col) subscripts, scaled by pixel size (dy, dx).
% Euclidean nearest-seed assignment is symmetric in the axes, so the row/col
% (vs x/y) naming does not affect results for isotropic pixels.
sy = opt.pixel_size(1); sx = opt.pixel_size(2);
linearIndices = find(ultimateErosion);
[Cr, Cc] = ind2sub(size(ultimateErosion), linearIndices);
C = [Cr(:)*sy, Cc(:)*sx]; clear Cr Cc
% Every pixel's (row, col) coordinates, scaled likewise
[allR, allC] = ind2sub(size(regional_maxima), (1:numel(regional_maxima))');
allVoxels = [allR*sy, allC*sx];
clear allR allC
% Use knnsearch to assign each pixel to the nearest seed point in 'C'
granules_labeled = knnsearch(C, allVoxels);
% Reshape the resulting index array to the dimensions of the 2D image
granules_labeled = reshape(granules_labeled, size(regional_maxima));
% Remove border-seeded regions unless the caller asked to keep them
if ~opt.keep_border
    ultimateErosion(1,:) = 0; ultimateErosion(:,1) = 0;
    ultimateErosion(end,:) = 0; ultimateErosion(:,end) = 0;
    linearIndices = find(ultimateErosion);
    mask = ismember(granules_labeled, granules_labeled(linearIndices));
    granules_labeled(~mask) = 0;
end
end%FCN:triang_extrema

% -------------------------------------------------------------------------

function [regional_maxima, regional_minima] = identifyLocalExtrema(projection, opt)
% Identify local maxima and minima in the 2D projection image using one of
% three methods: raw (direct), filtered (morphological pre-smoothing), or
% both (union of raw and filtered).

% The structuring element is only needed by the morphological (non-'raw')
% methods, so build it lazily to avoid the cost for 'raw'.
if ~strcmp(opt.regmax_method, 'raw')
    SE = strel('disk', opt.filter_size);
end
% Identify the local extrema based on the specified method
switch opt.regmax_method
    case 'raw'
        regional_maxima = imregionalmax(projection, 8);
        regional_minima = imregionalmin(projection, 8);
    case 'filtered'
        projection_max_filtered = imdilate(projection, SE);
        projection_min_filtered = imerode(projection, SE);
        regional_maxima = imregionalmax(projection_max_filtered, 8);
        regional_minima = imregionalmin(projection_min_filtered, 8);
    case 'both'
        projection_max_filtered = imdilate(projection, SE);
        projection_min_filtered = imerode(projection, SE);
        regional_maxima = imregionalmax(projection, 8) | imregionalmax(projection_max_filtered, 8);
        regional_minima = imregionalmin(projection, 8) | imregionalmin(projection_min_filtered, 8);
end%switch opt.regmax_method
end%FCN:identifyLocalExtrema

% -------------------------------------------------------------------------

function projection = correlationImage(stack, opt)
% Correlate each pixel trace with the mean trace of its in-bounds neighbors.
[rows, cols, t] = size(stack);
% Get a mask to identify neighbors
nMask = strel('disk', opt.filter_size);
nMask = nMask.Neighborhood;
center = (size(nMask,1) + 1) / 2;
nMask(center, center) = 0;
% Get neighbor offsets from mask
[dr, dc] = find(nMask);
dr = dr - center;
dc = dc - center;
% Compute neighbor mean using shifted stacks
stack = double(stack);   % accumulate before computing centered sums of squares
neighbor_sum = zeros(rows, cols, t);
neighbor_count = zeros(rows, cols);
for iN = 1:length(dr)
    src_r = max(1, 1-dr(iN)):min(rows, rows-dr(iN));
    src_c = max(1, 1-dc(iN)):min(cols, cols-dc(iN));
    dst_r = src_r + dr(iN);
    dst_c = src_c + dc(iN);
    neighbor_sum(dst_r, dst_c, :) = neighbor_sum(dst_r, dst_c, :) + stack(src_r, src_c, :);
    neighbor_count(dst_r, dst_c) = neighbor_count(dst_r, dst_c) + 1;
end%iN
neighbor_mean = neighbor_sum ./ neighbor_count;
clear neighbor_sum neighbor_count
% Compute pixel-wise Pearson correlation between each pixel and its
% neighbor mean: r = sum((X-muX).*(Y-muY)) / sqrt(sum((X-muX)^2)*sum((Y-muY)^2))
X = reshape(stack, [rows*cols, t]);
Y = reshape(neighbor_mean, [rows*cols, t]);
clear neighbor_mean
X = X - mean(X, 2);
Y = Y - mean(Y, 2);
projection = sum(X .* Y, 2) ./ sqrt(sum(X.^2, 2) .* sum(Y.^2, 2));
% Replace NaN (from constant-signal pixels) with 0
projection(isnan(projection)) = 0;
projection = single(reshape(max(-1, min(1, projection)), [rows, cols]));
end%FCN:correlationImage

% -------------------------------------------------------------------------

function granules_labeled = growing_regions(corr_img, reshaped_stack, opt)
% Greedy correlation-based region growing using BFS boundary tracking.
% Starting from the pixel with the highest local correlation, a region is
% grown by adding neighboring pixels whose time-course correlation
% meets grow_thresh against the region mean. Unassigned pixels can seed a
% new region when their local correlation meets seed_thresh.

% Separate thresholds for seeding vs. growing (default to corr_thresh)
seed_thr = opt.corr_thresh;
if ~isempty(opt.seed_thresh); seed_thr = opt.seed_thresh; end
grow_thr = opt.corr_thresh;
if ~isempty(opt.grow_thresh); grow_thr = opt.grow_thresh; end
% Check whether there is an upper limit to the region size
if isnumeric(opt.limitPixCount)
    max_size = opt.limitPixCount(2);
else
    max_size = inf;  % 'auto' — let post-segmentation enforce limits
end
[nRows, nCols] = size(corr_img);
granules_labeled = zeros(nRows, nCols);
cnt_id = 1;
px_cnt = numel(corr_img);
% Pre-sort pixels by correlation (descending) to avoid O(n) max search
% per iteration. This reduces total seeding cost from O(n^2) to O(n log n).
[~, sorted_idx] = sort(corr_img(:), 'descend');
% Progress bar only when verbose and a graphical desktop is available.
% onCleanup guarantees the handle is closed even if the loop errors out.
show_wb = opt.verbose && usejava('desktop');
wb_step = max(1, round(px_cnt / 100));
if show_wb
    hWait = waitbar(0, ['Growing regions ... 0/',num2str(px_cnt)]);
    wbCleanup = onCleanup(@() closeIfValid(hWait));
end
for iPx = 1:px_cnt
    % Throttle waitbar updates to every ~1% to avoid GUI overhead
    if show_wb && mod(iPx, wb_step) == 0
        waitbar(iPx/px_cnt, hWait, ['Growing regions ... ', num2str(iPx),'/',num2str(px_cnt)])
    end
    % Get the next unassigned pixel from the pre-sorted list
    max_id = sorted_idx(iPx);
    if granules_labeled(max_id) ~= 0
        continue;
    end%if already assigned
    if corr_img(max_id) < seed_thr
        break;
    end%if below seed threshold — all subsequent pixels are lower, so stop
    % Seed a new region
    n_region = 1;
    granules_labeled(max_id) = cnt_id;
    tc_sum = double(reshaped_stack(max_id, :));
    % Initialize boundary from the seed's 4-connected neighbors
    [sr, sc] = ind2sub([nRows, nCols], max_id);
    nr = [sr-1, sr+1, sr, sr];
    nc = [sc, sc, sc-1, sc+1];
    valid = nr>=1 & nr<=nRows & nc>=1 & nc<=nCols;
    boundary = sub2ind([nRows, nCols], nr(valid), nc(valid))';
    boundary = boundary(granules_labeled(boundary) == 0);
    % Grow region by testing boundary pixels against the mean TC
    while ~isempty(boundary) && n_region < max_size
        curr_TC = tc_sum / n_region;
        % Compare each candidate with the same region mean in double precision.
        B_tc = double(reshaped_stack(boundary, :));
        Bc = B_tc - mean(B_tc, 2);
        vc = curr_TC - mean(curr_TC);
        cc = (Bc * vc(:)) ./ (sqrt(sum(Bc.^2, 2)) * sqrt(sum(vc.^2)));
        % Zero-variance pixels give NaN, which fails the threshold (excluded).
        % Apply the size cap only after testing the entire boundary.
        pass = find(cc >= grow_thr);
        pass = pass(1:min(numel(pass), max_size - n_region));
        added_idx = boundary(pass);
        n_added = numel(added_idx);
        if n_added == 0
            break;
        end%if no pixels added
        granules_labeled(added_idx) = cnt_id;
        n_region = n_region + n_added;
        tc_sum = tc_sum + sum(B_tc(pass, :), 1);
        added_idx = added_idx(:);
        % Collect new boundary from added pixels' 4-connected neighbors
        new_boundary = zeros(n_added * 4, 1);
        n_nb = 0;
        for iA = 1:n_added
            [ar, ac] = ind2sub([nRows, nCols], added_idx(iA));
            nr = [ar-1, ar+1, ar, ar];
            nc = [ac, ac, ac-1, ac+1];
            valid = nr>=1 & nr<=nRows & nc>=1 & nc<=nCols;
            nb = sub2ind([nRows, nCols], nr(valid), nc(valid))';
            n_new = length(nb);
            new_boundary(n_nb+1:n_nb+n_new) = nb;
            n_nb = n_nb + n_new;
        end%iA
        new_boundary = unique(new_boundary(1:n_nb));
        boundary = new_boundary(granules_labeled(new_boundary) == 0);
    end%while growing
    cnt_id = cnt_id + 1;
end%iPx
% waitbar (if any) is closed by its onCleanup guard on function exit
end%FCN:growing_regions

% -------------------------------------------------------------------------

function [granules_labeled, nFeatures] = rICA_segmentation(projection, reshaped_stack, opt)
% Segment using reconstruction independent component analysis (rICA).
% The data is mean-subtracted, PCA-transformed, reduced, prewhitened, and
% then decomposed into independent components via rica(). Each pixel is
% assigned to the component with the strongest standardized loading.
% Get size of the dataset
[x, y] = size(projection);
% Center each input feature across pixels before PCA.
reshaped_stack = double(reshaped_stack);
reshaped_stack = reshaped_stack-mean(reshaped_stack,1,'omitnan');
% Apply PCA
[~,reshaped_stack,latent,~,explained] = pca(reshaped_stack, ...
    'Centered', false);
if isempty(latent) || latent(1) <= 0
    error('CalciSeg:degenerateData', 'rICA requires at least one varying PCA direction.');
end
tol = max(size(reshaped_stack)) * eps(latent(1));
nAvailable = nnz(latent > tol);
if nAvailable == 0
    error('CalciSeg:degenerateData', 'rICA requires at least one varying PCA direction.');
end
nFeatures = opt.n_rICA;
if nFeatures == 0
    nFeatures = find(cumsum(explained(1:nAvailable)) >= opt.n_PCA, 1);
    if isempty(nFeatures)
        nFeatures = nAvailable;
    end
end
% Overcomplete fits learn more features than there are input directions.
reshaped_stack = reshaped_stack(:, 1:min(nFeatures, nAvailable));
% Prewhiten the signal so that it has zero mean and identity covariance.
reshaped_stack = prewhiten(reshaped_stack);
% Apply reconstruction ICA (rICA)
rng_state = rng;
rng_cleanup = onCleanup(@() rng(rng_state));% Restore RNG even if rica throws
rng(opt.seed, 'twister'); % Reproducible independently of the caller's generator.
ricaTransform = rica(single(reshaped_stack), nFeatures, 'VerbosityLevel', double(opt.verbose), 'IterationLimit', 5000);
clear rng_cleanup;% Normal exit: restore now rather than at scope end
% Transform the data using the fitted model
transformedData = transform(ricaTransform, reshaped_stack);
% Standardize each component across pixels before comparing spatial loadings.
s = std(transformedData, 0, 1);
s(s == 0) = 1;% Guard zero-variance components against division by zero
TDv = (transformedData - mean(transformedData,1)) ./ s;
% Use skewness by default or compare excursions beyond the requested SD.
if isempty(opt.rica_sd_threshold)
    sign_score = mean(TDv.^3, 1);
else
    threshold = opt.rica_sd_threshold;
    positiveTail = mean(max(TDv - threshold, 0), 1);
    negativeTail = mean(max(-TDv - threshold, 0), 1);
    sign_score = positiveTail - negativeTail;
    sign_score(abs(sign_score) <= 8*eps(max(positiveTail, negativeTail))) = 0;
end
flip_sign = ones(1, size(TDv,2));
flip_sign(sign_score < 0) = -1;
transformedData = TDv .* flip_sign;
% Reshape the independent components back into a 3-D format
transformedData = reshape(transformedData, [x, y, size(transformedData,2)]);
if strcmp(opt.regmax_method, 'filtered')
    for component = 1:nFeatures
        transformedData(:,:,component) = imboxfilt(transformedData(:,:,component), ...
            2*opt.filter_size + 1, 'Padding', 'replicate');
    end
end
% Get the component with the strongest effect on a given pixel
[label_load, granules_labeled] = max(transformedData, [], 3);
% Threshold the median-filtered winning load at the lowest Otsu level.
if opt.rICA_thresh > 0
    min_val = quantile(label_load(:), 1-0.995);
    max_val = quantile(label_load(:), 0.995);
    label_load(label_load>max_val) = max_val;
    label_load(label_load<min_val) = min_val;
    threshold = min(multithresh(label_load, opt.rICA_thresh));
    granules_labeled(label_load < threshold) = 0;
end
end%FCN:rICA_segmentation

% -------------------------------------------------------------------------

function Z = prewhiten(X)
% From https://de.mathworks.com/help/stats/extract-mixed-signals.html
% Decompose the feature covariance, retaining directions with positive variance.
[U,Sig] = svd(cov(X));
Sig     = diag(Sig);
Sig     = Sig(:)';
% Keep singular values above a tolerance relative to the largest value.
tol = eps(class(X));
idx = (Sig > max(Sig)*tol);
if all(idx == 0)
    % All singular values are ~0: the data has no variance to whiten
    % (e.g. a constant stack). ICA is undefined here, so fail clearly
    % rather than with a bare assertion deep in the call stack.
    error('CalciSeg:degenerateData', ...
        'rICA cannot run on data with no variance (constant or rank-deficient stack).');
end
% Select the retained variances and their covariance directions.
Sig = Sig(idx);
U   = U(:,idx);
% Compute prewhitened data.
mu = mean(X,1);
Z = X-mu;
Z = (Z*U)./sqrt(Sig);
end%FCN:prewhiten

% -------------------------------------------------------------------------

function [granules_labeled, n_used] = refineSegmentation(reshaped_stack, granules_labeled, opt)
% Iteratively refine granule boundaries. In each iteration, border pixels
% are reassigned to the neighboring granule with the best fit (lowest root-median-square error
% or highest correlation). After each re-assignment pass, split granules
% are detected and size constraints are enforced. Stops early if the label
% map converges or begins to oscillate between two states (2-cycle).
% n_used returns the number of iterations actually performed.
previous_granules_labeled = -ones(size(granules_labeled));
state_two_ago = -ones(size(granules_labeled));   % state from two iters ago
n_used = opt.n_rep;
% Progress bar only when verbose and a graphical desktop is available.
% onCleanup guarantees the handle is closed even if an iteration errors out.
show_wb = opt.verbose && usejava('desktop');
if show_wb
    hWait = waitbar(0, ['Refine regions ... 0/',num2str(opt.n_rep)]);
    wbCleanup = onCleanup(@() closeIfValid(hWait));
end
for iRep = 1:opt.n_rep
    if show_wb
        waitbar(iRep/opt.n_rep, hWait, ['Refine regions ... ', num2str(iRep),'/',num2str(opt.n_rep)])
    end
    % Get avg signal trace per granule via grouped indexing
    [allIDs, ~, gIdx] = unique(granules_labeled(:));
    bg = allIDs == 0;
    granuleList = allIDs(~bg);
    pixel_groups = accumarray(gIdx, (1:numel(granules_labeled))', [], @(x){x});
    pixel_groups(bg) = [];
    granuleTC = nan(length(granuleList), size(reshaped_stack,2));
    for iGranule = 1:length(granuleList)
        granuleTC(iGranule, :) = mean(reshaped_stack(pixel_groups{iGranule},:),1,'omitnan');
    end%iGranule
    % Identify border pixels that need re-evaluation. Only consider
    % regions that changed since the last iteration (for efficiency).
    ind = find(granules_labeled~=previous_granules_labeled);
    keep_changed_IDs = unique([granules_labeled(ind(:)); previous_granules_labeled(ind(:))]);
    mask = ismember(granules_labeled, keep_changed_IDs);
    candidates = imgradient(granules_labeled, "intermediate")>0;
    candidates(~mask) = 0;
    idx_candidates = find(candidates);
    % Iterate over all pixels in border regions and check
    % whether they have to be re-assigned.
    previous_granules_labeled = granules_labeled;
    idx_candidates_newIdentity = refinement_parfor(idx_candidates, granules_labeled, reshaped_stack, granuleTC, granuleList, opt.refinement_method, opt.use_parallel);
    % Update changed identities
    granules_labeled(idx_candidates) = idx_candidates_newIdentity;
    % Account for splitting of granules. For this, iterate over
    % all granules and check whether the corresponding identity
    % can be found in more than one coherent region
    granules_labeled = checkSplits(granules_labeled, opt.use_parallel);
    % Enforce size constraints after re-assignment
    granules_labeled = removeHugeRegions(granules_labeled, opt.limitPixCount(2));
    granules_labeled = removeTinyRegions(reshaped_stack, granules_labeled, opt);
    % Stop early if the label map has converged, or if it oscillates between
    % two states (a 2-cycle), which exact-equality alone would never catch.
    if isequal(previous_granules_labeled, granules_labeled)
        n_used = iRep; break
    elseif isequal(granules_labeled, state_two_ago)
        if opt.verbose
            fprintf('Refinement stopped after %d iterations (2-cycle oscillation).\n', iRep);
        end
        n_used = iRep; break
    end%if converged/oscillating
    state_two_ago = previous_granules_labeled;
end%iRep
% waitbar (if any) is closed by its onCleanup guard on function exit
% Make sure the highest ID is equal to the number of granules
granules_labeled = finalRefinement(granules_labeled);
end%FCN:refineSegmentation

% -------------------------------------------------------------------------

function idx_candidates_newIdentity = refinement_parfor(idx_candidates, granules_labeled, reshaped_stack, granuleTC, granuleList, method, use_parallel)
% Re-assign border pixels to the best-fitting neighboring granule using one
% of two distance measures:
%   'corr' : Pearson correlation distance |1 - r|
%   'rmse' : root-median-square error sqrt(median((cluster - px).^2))
% Each candidate pixel is compared against its 4-connected neighboring
% granules and assigned to the one with the smallest distance. Per-iteration
% data is pre-extracted so that the full reshaped_stack and granules_labeled
% matrices are not broadcast to parallel workers.

% Prepare output of new labels
idx_candidates_newIdentity = granules_labeled(idx_candidates);
if isempty(idx_candidates) || isempty(granuleList)
    return
end
[nRows, nCols] = size(granules_labeled);
nCandidates = length(idx_candidates);

% Pre-extract per-candidate data to avoid broadcasting the full
% reshaped_stack (pixels x time) and granules_labeled (rows x cols)
candidate_TCs = reshaped_stack(idx_candidates, :);
neighbor_labels = cell(nCandidates, 1);
for iCan = 1:nCandidates
    [pr, pc] = ind2sub([nRows, nCols], idx_candidates(iCan));
    nr = [pr-1, pr+1, pr, pr];
    nc = [pc, pc, pc-1, pc+1];
    valid = nr>=1 & nr<=nRows & nc>=1 & nc<=nCols;
    neighbor_idx = sub2ind([nRows, nCols], nr(valid), nc(valid));
    nC = unique(granules_labeled(neighbor_idx));
    nC(nC == 0) = [];
    neighbor_labels{iCan} = nC;
end%iCan

% Map each granule ID to its row in the matrix of mean traces.
id_to_row = zeros(max(granuleList), 1);
id_to_row(granuleList) = 1:length(granuleList);

% Prepare metric-specific per-row data. For 'corr', pre-normalize (center +
% unit L2 norm, in double) so the per-neighbor Pearson correlation reduces to
% a dot product; for 'rmse', use the raw traces. Zero-variance rows become
% NaN, which is treated as a non-finite distance below.
use_corr = strcmp(method, 'corr');
if use_corr
    cand_data = normalizeRows(candidate_TCs);
    tc_data = normalizeRows(granuleTC);
else
    cand_data = candidate_TCs;
    tc_data = granuleTC;
end

% Parallel reassignment — only sliced (cand_data, neighbor_labels) and small
% broadcast (tc_data, id_to_row, use_corr) variables enter the loop. Gate
% parallelism on the request and a work-size threshold.
parMax = 0;
if use_parallel && nCandidates >= 256; parMax = 1e6; end
parfor (iCan = 1:nCandidates, parMax)
    nC = neighbor_labels{iCan};
    if isempty(nC); continue; end
    % Get the pixel's TC
    activity_px = cand_data(iCan, :);
    % Check which neighborhood is a better fit
    voronoi_R = zeros(length(nC),1);
    for iN = 1:length(nC)
        activity_cluster = tc_data(id_to_row(nC(iN)), :);
        if use_corr
            voronoi_R(iN) = abs(1 - (activity_px * activity_cluster'));
        else
            voronoi_R(iN) = sqrt(median((activity_cluster - activity_px).^2));
        end
    end%iN
    % Assign new identity (keep current if all distances are NaN)
    voronoi_R(~isfinite(voronoi_R)) = NaN;
    if all(isnan(voronoi_R)); continue; end
    [~, iBest] = min(voronoi_R);
    idx_candidates_newIdentity(iCan) = nC(iBest);
end%iCan
end%FCN:refinement_parfor

% -------------------------------------------------------------------------

function granules_labeled = checkSplits(granules_labeled, use_parallel)
% Detect and resolve granules that have been split into multiple
% disconnected components. Each non-contiguous fragment beyond the largest
% receives a new unique ID when it contains multiple pixels. A singleton
% joins the most frequent adjacent positive ID (lowest ID on a tie),
% or becomes background if isolated.

% Identify splitter candidates (granules with >1 connected component).
% Pre-extract per-granule pixel index lists so that the full
% granules_labeled matrix is not broadcast to every parallel worker.
unique_granules = unique(granules_labeled);
unique_granules(unique_granules == 0) = [];
nGranules = length(unique_granules);
pixel_lists = cell(nGranules, 1);
for iP = 1:nGranules
    pixel_lists{iP} = find(granules_labeled == unique_granules(iP));
end%iP
[nRows, nCols] = size(granules_labeled);
is_splitter = false(nGranules, 1);
% Gate parallelism: use workers only when requested and the workload is
% large enough to outweigh parfor overhead. Component counting operates on a
% per-granule bounding box rather than a full-frame mask.
parMax = 0;
if use_parallel && nGranules >= 64; parMax = 1e6; end
parfor (iP = 1:nGranules, parMax)
    px = pixel_lists{iP};
    [rr, cc] = ind2sub([nRows, nCols], px);
    r0 = min(rr); c0 = min(cc);
    sub = false(max(rr)-r0+1, max(cc)-c0+1);
    sub(sub2ind(size(sub), rr-r0+1, cc-c0+1)) = true;
    [~, num] = bwlabel(sub, 4);
    is_splitter(iP) = (num > 1);
end%iP
splitter_candidates = unique_granules(is_splitter);
% Iterate over all splitters. The LARGEST connected component keeps the
% original ID; the others receive new IDs, are merged into a neighbour (single
% pixels), or are dropped to background (isolated single pixels).
for iGranule = 1:length(splitter_candidates)
    id = splitter_candidates(iGranule);
    bw = granules_labeled==id;
    L = bwlabel(bw,4);
    nComp = max(L(:));
    if nComp <= 1
        continue;
    end%if only one connected component
    % Component sizes; the largest fragment retains the original ID
    componentIDs = L(bw);
    compSizes = accumarray(componentIDs(:), 1);
    [~, largest] = max(compSizes);
    for c = 1:nComp
        if c == largest
            continue;
        end%if largest fragment keeps original ID
        idx = find(L == c);
        if isscalar(idx)
            % Assign to an adjacent positive ID, or to background if isolated.
            [pr, pc] = ind2sub(size(granules_labeled), idx);
            nr = [pr-1, pr+1, pr, pr];
            nc = [pc, pc, pc-1, pc+1];
            valid = nr>=1 & nr<=size(granules_labeled,1) & nc>=1 & nc<=size(granules_labeled,2);
            nb_idx = sub2ind(size(granules_labeled), nr(valid), nc(valid));
            nb_idx = nb_idx(~bw(nb_idx));
            neighbor_IDs = granules_labeled(nb_idx);
            neighbor_IDs(neighbor_IDs == 0) = [];
            if isempty(neighbor_IDs)
                granules_labeled(idx) = 0;
            else
                granules_labeled(idx) = mode(neighbor_IDs(:));
            end
        else
            newID = max(unique_granules) + 1;
            granules_labeled(idx) = newID;
            unique_granules = [unique_granules(:); newID];
        end
    end%c
end%iGranule
end%FCN:checkSplits

% -------------------------------------------------------------------------

function granules_labeled = removeTinyRegions(reshaped_stack, granules_labeled, opt)
% Merge granules smaller than limitPixCount(1) into their best-fitting
% neighbor. Uses a batch approach: within each pass, a set of tiny regions
% that share no neighbors is selected so that merging them cannot affect one
% another's best-fit decision; those independent merges are then applied
% sequentially. The outer loop repeats until no granule is below the minimum
% size.
while true
    % Get pixel count for each region via single-pass accumarray
    [allIDs, ~, idx] = unique(granules_labeled(:));
    pixCounts_all = accumarray(idx, 1);
    bg = allIDs == 0;
    granuleList = allIDs(~bg);
    pixCounts = pixCounts_all(~bg);
    if isempty(granuleList); break; end
    C = [pixCounts, granuleList(:)];
    % Sort based on count
    C = sortrows(C,1,"ascend");
    % Find all regions below the minimum size
    small_granules = C(C(:, 1) < opt.limitPixCount(1), 2);
    if isempty(small_granules)
        break;
    end

    % Identify independent tiny regions for batch processing.
    % Two tiny regions are independent if they share no neighbors
    % (neither each other nor any common neighbor), so merging one
    % cannot affect the other's best-fit decision.
    n_small = length(small_granules);
    [nRows_t, nCols_t] = size(granules_labeled);
    neighbor_sets = cell(n_small, 1);
    for iS = 1:n_small
        pix = find(granules_labeled == small_granules(iS));
        [pr, pc] = ind2sub([nRows_t, nCols_t], pix);
        % Collect all 4-connected neighbors of the region's pixels
        all_nr = [pr-1; pr+1; pr; pr];
        all_nc = [pc; pc; pc-1; pc+1];
        valid = all_nr>=1 & all_nr<=nRows_t & all_nc>=1 & all_nc<=nCols_t;
        nb_idx = sub2ind([nRows_t, nCols_t], all_nr(valid), all_nc(valid));
        % Exclude pixels belonging to the region itself
        nb_idx = nb_idx(granules_labeled(nb_idx) ~= small_granules(iS));
        nIDs = unique(granules_labeled(nb_idx));
        nIDs(nIDs == 0) = [];
        neighbor_sets{iS} = nIDs(:)';
    end%iS
    % Greedy selection (smallest first): skip any region that shares
    % a neighbor with an already-selected one
    selected = false(n_small, 1);
    claimed = zeros(1, 0);
    for iS = 1:n_small
        curr_id = small_granules(iS);
        curr_nb = neighbor_sets{iS};
        if ismember(curr_id, claimed) || any(ismember(curr_nb, claimed))
            continue;
        end
        selected(iS) = true;
        claimed = [claimed, curr_id, curr_nb];
    end%iS
    batch_ids = find(selected);

    % Iterate over all tiny regions that do not share neighbors
    for iB = 1:length(batch_ids)
        % Get the ID of the current granule
        curr_ID = small_granules(batch_ids(iB));
        % Get a logical image for the current granule
        bw = granules_labeled==curr_ID;
        % Get the pixels of the current granule
        ind_pix = find(bw);
        % Get the current time course
        curr_TC = mean(reshaped_stack(ind_pix, :), 1, 'omitnan');
        % Selected regions have disjoint neighborhoods throughout this batch.
        neighbor_IDs = neighbor_sets{batch_ids(iB)};
        % Guard: no valid neighbors — assign to background
        if isempty(neighbor_IDs)
            granules_labeled(ind_pix) = 0;
            continue;
        end
        % Get the neighbors' time courses
        TCs = nan(length(neighbor_IDs), size(reshaped_stack,2));
        for iTC = 1:length(neighbor_IDs)
            ind = find(granules_labeled==neighbor_IDs(iTC));
            TCs(iTC,:) = mean(reshaped_stack(ind,:),1,'omitnan');
        end%iTC
        % Check which neighboring granule is the best fit
        id_table = nan(1, length(neighbor_IDs));
        for iN = 1:length(neighbor_IDs)
            % Get the rmse or distance based on correlation
            switch opt.refinement_method
                case 'corr'
                    min_dist = abs(1 - normalizeRows(curr_TC) * normalizeRows(TCs(iN,:))');
                case 'rmse'
                    min_dist = sqrt(median((curr_TC - TCs(iN,:)).^2));
            end
            id_table(iN) = min_dist;
        end%iN
        % Assign all pixels of the current granule to the best fitting neighbor
        % (if all distances are NaN, assign to first neighbor as fallback)
        if all(isnan(id_table))
            best_id = 1;
        else
            [~,best_id] = min(id_table);
        end
        granules_labeled(ind_pix) = neighbor_IDs(best_id(1));
    end%iB
end%while
end%FCN:removeTinyRegions

% -------------------------------------------------------------------------

function granules_labeled = removeHugeRegions(granules_labeled, maxPixCount)
% Split granules larger than maxPixCount. Three splitting strategies are
% tried in order: (1) watershed on the distance transform, (2) cutting
% perpendicular to the major axis, (3) cutting along the major axis.
% Regions that cannot be split by any method are skipped.
unsplittable = [];
while true
    % Get pixel count for each region via single-pass accumarray
    [allIDs, ~, idx] = unique(granules_labeled(:));
    pixCounts_all = accumarray(idx, 1);
    bg = allIDs == 0;
    granuleList = allIDs(~bg);
    pixCounts = pixCounts_all(~bg);
    if isempty(granuleList); break; end
    C = [pixCounts, granuleList(:)];
    % Sort based on count
    C = sortrows(C,1,"descend");
    % Find large granules, excluding those that failed to split
    large_granules = C(C(:, 1) > maxPixCount, 2);
    large_granules = large_granules(~ismember(large_granules, unsplittable));
    if isempty(large_granules)
        break;
    end
    for iGranule = 1:length(large_granules)
        % Get the ID of the current granule
        curr_ID = large_granules(iGranule);
        % Get a logical image for the current granule
        bw = granules_labeled==curr_ID;
        % Strategy 1: Watershed on the distance transform
        D = -bwdist(~bw);
        L = watershed(D);
        L(~bw) = 0;
        bw_cut = L>0;
        bw_cut = bwlabel(bw_cut, 4);   % 4-conn keeps sub-regions 4-connected
        % Strategy 2: If watershed did not split the region (e.g.,
        % circular shapes), try cutting perpendicular to the major axis.
        if max(bw_cut(:)) <= 1
            % Get the fitted ellipse properties. If the region happens to
            % have several connected components, base the cut on the largest
            % one rather than an arbitrary first component.
            stats = regionprops(bw, 'Area', 'Centroid', 'MajorAxisLength', 'Orientation');
            [~, largestComp] = max([stats.Area]);
            ellipse = stats(largestComp);
            longAxis = ellipse.MajorAxisLength;
            % Convert the Cartesian orientation to image coordinates (rows increase downward).
            centerX = ellipse.Centroid(1);
            centerY = ellipse.Centroid(2);
            angle = ellipse.Orientation;
            % Strategy 2: cut perpendicular to the major axis
            x1 = floor(centerX - (longAxis/2) * sind(angle));
            y1 = floor(centerY - (longAxis/2) * cosd(angle));
            x2 = ceil(centerX + (longAxis/2) * sind(angle));
            y2 = ceil(centerY + (longAxis/2) * cosd(angle));
            bw_cut = cutMaskAlongLine(bw, x1, y1, x2, y2);
            % Strategy 3: If the perpendicular cut did not work, try
            % cutting along the major axis instead.
            if max(bw_cut(:)) <= 1
                x1 = floor(centerX - (longAxis/2) * cosd(angle));
                y1 = floor(centerY + (longAxis/2) * sind(angle));
                x2 = ceil(centerX + (longAxis/2) * cosd(angle));
                y2 = ceil(centerY - (longAxis/2) * sind(angle));
                bw_cut = cutMaskAlongLine(bw, x1, y1, x2, y2);
            end%if perpendicular cut did not work
        end%if watershed did not work
        % Check if any splitting method succeeded
        if max(bw_cut(:)) <= 1
            unsplittable = [unsplittable; curr_ID];
        else
            % Fill pixels on the cut line (bw_cut==0 but bw==1) by assigning
            % them to the neighboring sub-region with the highest label.
            bw_fill = bw_cut;
            ind = find(bw_cut==0 & bw==1);
            [nRows_h, nCols_h] = size(bw);
            for iPx = 1:length(ind)
                % Get 4-connected neighbors of the cut-line pixel
                [pr, pc] = ind2sub([nRows_h, nCols_h], ind(iPx));
                nr = [pr-1, pr+1, pr, pr];
                nc = [pc, pc, pc-1, pc+1];
                valid = nr>=1 & nr<=nRows_h & nc>=1 & nc<=nCols_h;
                nb_idx = sub2ind([nRows_h, nCols_h], nr(valid), nc(valid));
                neighbor_IDs = unique(bw_cut(nb_idx));
                neighbor_IDs = neighbor_IDs(neighbor_IDs > 0);
                if isempty(neighbor_IDs); continue; end
                bw_fill(ind(iPx)) = max(neighbor_IDs);
            end%iPx
            bw_cut = bw_fill;
            % Shift all sub-region labels above the current global maximum
            % to avoid ID collisions
            bw_cut = bw_cut+max(granules_labeled(:))+1;
            granules_labeled(bw) = bw_cut(bw);
        end%if split succeeded
    end%iGranule
end%while
end%FCN:removeHugeRegions

function bw_cut = cutMaskAlongLine(bw, x1, y1, x2, y2)
% Zero out a 1-pixel-wide straight line from (x1,y1) to (x2,y2) on the binary
% mask 'bw' and return the 4-connected labeling of the result. Used to split
% a region that watershed could not separate. Each sampled point is clamped
% to the image bounds.
bw_cut = bw;
numPoints = round(sqrt(size(bw,1)^2 + size(bw,2)^2));
x = round(linspace(x1, x2, numPoints)); x(x<1)=1; x(x>size(bw_cut,2))=size(bw_cut,2);
y = round(linspace(y1, y2, numPoints)); y(y<1)=1; y(y>size(bw_cut,1))=size(bw_cut,1);
for iPx = 1:numPoints
    bw_cut(y(iPx), x(iPx)) = 0;
end%iPx
bw_cut = bwlabel(bw_cut, 4);
end%FCN:cutMaskAlongLine

% -------------------------------------------------------------------------

function granules_labeled = poolFragmentedGranules(reshaped_stack, granules_labeled, opt)
% Merge neighboring granules whose mean time courses are highly correlated
% (above corr_thresh). Uses an adjacency-first approach: the spatial
% neighbor graph is built from the label map, correlations are computed
% only for adjacent pairs, and connected components above the threshold
% are merged. Repeats until no pair exceeds the threshold.
while true
    % Get a list of granules via grouped indexing
    [allIDs, ~, gIdx] = unique(granules_labeled(:));
    bg = allIDs == 0;
    granuleList = allIDs(~bg);
    if length(granuleList) < 2; break; end
    pixel_groups = accumarray(gIdx, (1:numel(granules_labeled))', [], @(x){x});
    pixel_groups(bg) = [];
    % Get each granule's time course
    avgTCs = nan(length(granuleList), size(reshaped_stack,2));
    for iGranule = 1:length(granuleList)
        avgTCs(iGranule,:) = mean(reshaped_stack(pixel_groups{iGranule},:),1,'omitnan');
    end%iGranule
    % Only granules sharing a horizontal or vertical edge can be pooled.
    % Horizontal pairs
    left = granules_labeled(:, 1:end-1);
    right = granules_labeled(:, 2:end);
    diff_mask = left ~= right & left > 0 & right > 0;
    left = left(diff_mask); right = right(diff_mask);
    pairs_h = [left(:), right(:)];
    % Vertical pairs
    top = granules_labeled(1:end-1, :);
    bottom = granules_labeled(2:end, :);
    diff_mask = top ~= bottom & top > 0 & bottom > 0;
    top = top(diff_mask); bottom = bottom(diff_mask);
    pairs_v = [top(:), bottom(:)];
    % Get unique adjacent pairs
    all_pairs = unique(sort([pairs_h; pairs_v], 2), 'rows');
    if isempty(all_pairs); break; end
    % Compute correlations only for adjacent pairs (vectorized)
    [~, pair_idx1] = ismember(all_pairs(:, 1), granuleList);
    [~, pair_idx2] = ismember(all_pairs(:, 2), granuleList);
    tc1 = avgTCs(pair_idx1, :);
    tc2 = avgTCs(pair_idx2, :);
    tc1 = tc1 - mean(tc1, 2);
    tc2 = tc2 - mean(tc2, 2);
    r = sum(tc1 .* tc2, 2) ./ sqrt(sum(tc1.^2, 2) .* sum(tc2.^2, 2));
    % Keep above-threshold pairs
    combineIDs = all_pairs(r > opt.corr_thresh, :);
    % Stop when there is nothing left to pool
    if isempty(combineIDs)
        break
    end%if empty
    % Graph node numbers are the labels in the current map.
    G = graph(combineIDs(:,1), combineIDs(:,2));
    bin = conncomp(G);
    % Check all groups
    for iG = 1:max(bin)
        % Get ID list
        ID_list = find(bin == iG);
        if length(ID_list)>1
            % Assign a new ID to those connected
            newID = max(ID_list);
            for iID = 1:length(ID_list)
                granules_labeled(granules_labeled==ID_list(iID)) = newID;
            end%iID
        end%if connected
    end%iG
end%while
end%FCN:poolFragmentedGranules

% -------------------------------------------------------------------------

function granules_RElabeled = finalRefinement(granules_labeled)
% Relabel granules with contiguous integer IDs ordered by DESCENDING pixel
% count (ID 1 = largest granule; ties broken by original ID). This makes the
% labeling deterministic and comparable across runs. Background (0) stays 0.
sz = size(granules_labeled);
[uv, ~, ic] = unique(granules_labeled(:));   % uv sorted; ic maps pixels -> uv index
counts = accumarray(ic, 1);                   % pixel count per unique value
newId = zeros(numel(uv), 1);                  % new id per unique value (bg -> 0)
fgIdx = find(uv ~= 0);                         % foreground entries in uv
% Rank foreground values by descending count, tie-break by ascending value
[~, ord] = sortrows([-counts(fgIdx), double(uv(fgIdx))]);
newId(fgIdx(ord)) = 1:numel(fgIdx);
granules_RElabeled = reshape(newId(ic), sz);
end%FCN:finalRefinement

% -------------------------------------------------------------------------

function summary_stats = regionStats(reshaped_stack, granules_labeled, opt)
% Compute summary statistics for each granule: mean time course, within-
% granule correlation, per-granule projection images (mean, std, min, max,
% correlation), and quality metrics (SNR, compactness). Also estimates an
% active-region binary map using Otsu's method on the metric matching the
% projection method.

% Get list of granules via grouped indexing
[allIDs, ~, gIdx] = unique(granules_labeled(:));
bg = allIDs == 0;
granuleList = allIDs(~bg);
pixel_groups = accumarray(gIdx, (1:numel(granules_labeled))', [], @(x){x});
pixel_groups(bg) = [];
% Preallocation
nG = length(granuleList);
nT = size(reshaped_stack,2);
[nRows, nCols] = size(granules_labeled);
summary_stats.avgTCs = nan(nG, nT);
summary_stats.granule_Corr = nan(nG, 1);
summary_stats.granule_Avg_img = nan(size(granules_labeled));
summary_stats.granule_Min_img = nan(size(granules_labeled));
summary_stats.granule_Max_img = nan(size(granules_labeled));
summary_stats.granule_Std_img = nan(size(granules_labeled));
summary_stats.granule_Median_img = nan(size(granules_labeled));
summary_stats.granule_Corr_img = nan(size(granules_labeled));
summary_stats.area = nan(nG, 1);              % pixel count per granule
summary_stats.centroid = nan(nG, 2);          % [x(col) y(row)]
summary_stats.bbox = nan(nG, 4);              % [x y width height]
summary_stats.pixel_idx = cell(nG, 1);        % linear pixel indices per granule
summary_stats.quality.SNR = nan(nG, 1);
summary_stats.quality.compactness = nan(nG, 1);
% Iterate over all granules, computing time-course statistics, geometry and
% quality metrics in a single pass.
for iGranule = 1:nG
    % Get their index positions
    idx_granule = pixel_groups{iGranule};
    % Geometry (each granule is a single 4-connected component by invariant)
    [rr, cc] = ind2sub([nRows, nCols], idx_granule);
    summary_stats.area(iGranule)      = numel(idx_granule);
    summary_stats.centroid(iGranule,:) = [mean(cc), mean(rr)];
    summary_stats.bbox(iGranule,:)     = [min(cc), min(rr), max(cc)-min(cc)+1, max(rr)-min(rr)+1];
    summary_stats.pixel_idx{iGranule}  = idx_granule;
    % Get the avg time course of the current granule
    summary_stats.avgTCs(iGranule,:) = mean(reshaped_stack(idx_granule,:),1,'omitnan');
    % Get within-granule correlation
    summary_stats.granule_Corr(iGranule,1) = calculateAverageCorrelation(summary_stats.avgTCs(iGranule,:), reshaped_stack(idx_granule,:));
    % Project different metrics into 2D (one scalar per granule)
    summary_stats.granule_Avg_img(idx_granule)  = mean(summary_stats.avgTCs(iGranule,:), 'omitnan');
    summary_stats.granule_Min_img(idx_granule)  = min(summary_stats.avgTCs(iGranule,:), [], 'omitnan');
    summary_stats.granule_Max_img(idx_granule)  = max(summary_stats.avgTCs(iGranule,:), [], 'omitnan');
    summary_stats.granule_Std_img(idx_granule)  = std(summary_stats.avgTCs(iGranule,:), 0, 'omitnan');
    summary_stats.granule_Median_img(idx_granule) = median(summary_stats.avgTCs(iGranule,:), 'omitnan');
    summary_stats.granule_Corr_img(idx_granule) = summary_stats.granule_Corr(iGranule,1);
    % --- Quality metrics ---
    % Estimate noise from successive differences of the granule mean trace.
    % For independent Gaussian noise, median(abs(diff)) is about 0.9539*sigma.
    tc = summary_stats.avgTCs(iGranule, :);
    baseline = prctile(tc, 10);
    noise = median(abs(diff(tc)), 'omitnan') / 0.9539;
    if isfinite(noise) && noise > 0
        summary_stats.quality.SNR(iGranule) = (max(tc) - baseline) / noise;
    end
    % Compactness: 4*pi*area / perimeter^2 (1 = perfect circle; can slightly
    % exceed 1 for small rasterised shapes due to perimeter discretisation).
    roi_mask = granules_labeled == granuleList(iGranule);
    stats_roi = regionprops(roi_mask, 'Area', 'Perimeter');
    if ~isempty(stats_roi)
        total_area = sum([stats_roi.Area]);
        total_perimeter = sum([stats_roi.Perimeter]);
        if total_perimeter > 0
            summary_stats.quality.compactness(iGranule) = ...
                4 * pi * total_area / total_perimeter^2;
        end
    end
end%iGranule
% Report any granules that still exceed the maximum size (max-size enforcement
% is best-effort: some shapes cannot be split by any strategy).
if isnumeric(opt.limitPixCount) && isfinite(opt.limitPixCount(2))
    summary_stats.oversized_ids = granuleList(summary_stats.area > opt.limitPixCount(2));
else
    summary_stats.oversized_ids = [];
end
% PCA uses temporal standard deviation for the activity mask. 'none' skips it.
metric = opt.projection_method;
if strcmp(metric, 'pca')
    metric = 'std';
end
switch metric
    case 'std',    active_img = summary_stats.granule_Std_img;
    case 'mean',   active_img = summary_stats.granule_Avg_img;
    case 'median', active_img = summary_stats.granule_Median_img;
    case 'max',    active_img = summary_stats.granule_Max_img;
    case 'min',    active_img = summary_stats.granule_Min_img;
    case 'corr',   active_img = summary_stats.granule_Corr_img;
    case 'none',   active_img = [];
end%switch projection method
summary_stats.active_region.method = metric;
if ~isempty(active_img)
    % Otsu threshold computed over FOREGROUND (granule) pixels only, so the
    % split is not biased by the amount of background in the field of view.
    fgmask = ~isnan(active_img);
    vals = active_img(fgmask);
    map = false(size(active_img));
    if ~isempty(vals) && max(vals) > min(vals)
        vmin = min(vals); vmax = max(vals);
        level = graythresh((vals - vmin) / (vmax - vmin));   % Otsu on foreground
        map(fgmask) = vals > (vmin + level * (vmax - vmin));
    elseif ~isempty(vals)
        map(fgmask) = true;   % all foreground values equal -> all active
    end
    summary_stats.active_region.map = map;
else
    summary_stats.active_region.map = [];
end%if active image
end%FCN:regionStats

% -------------------------------------------------------------------------

function overallAvg = calculateAverageCorrelation(avg_tc, ind_tc)
% Compute the mean Pearson correlation between each pixel's time course
% (rows of ind_tc) and the granule's mean time course (avg_tc).
% Uses vectorized correlation for efficiency.
if size(avg_tc, 2)>1
    corr_values = normalizeRows(ind_tc) * normalizeRows(avg_tc)';
    % Sanitize non-finite values (from zero-variance traces) before averaging
    corr_values(~isfinite(corr_values)) = NaN;
    overallAvg = mean(corr_values, 'omitnan');
else
    overallAvg = NaN;
end%if more dimensions
end%FCN:calculateAverageCorrelation

% -------------------------------------------------------------------------

function Xn = normalizeRows(X)
% Center each row and scale it to unit L2 norm (in double precision) so that
% the dot product of two normalized rows equals their Pearson correlation.
% Rows with zero variance become NaN (their correlation is undefined), which
% callers treat as a non-finite distance.
X = double(X);
Xc = X - mean(X, 2);
Xn = Xc ./ sqrt(sum(Xc.^2, 2));
end%FCN:normalizeRows