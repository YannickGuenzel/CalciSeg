function [granules_labeled, summary_stats] = CalciSeg(stack, varargin)
% [granules_labeled, summary_stats] = CalciSeg(stack, Name, Value)
% segments the spatial aspects of a 3-D stack (x*y*time) into individual
% regions ('granules') based on a refined 2D Delaunay triangulation.
% Note: CalciSeg uses parallel for-loops for the local growing
% algorithm and the refinement process. It is recommended to start a
% parallel pool of workers before calling the function.
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
%                       - 'corr'   : Correlation space as local correlation
%                                    between neighboring pixels
%                       - 'none'   : Applicable when "stack" has no third 
%                                    dimension or when to take the first 
%                                    slice of the 3-D stack
%
% 'init_seg_method'   : Method for the initial segmentation.
% (string)              Default: 'voronoi'
%                       - 'voronoi' : Delaunay triangulation
%                       - 'corr'    : Local growing based on correlation.
%                                     This will set 'projection_method' to
%                                     be 'corr'
%                       - 'rICA'    : Reconstruction independent component 
%                                     analysis.
%
% 'regmax_method'     : Method for determining how to identify local
% (string)              extrema in the intensity projection.
%                       Default: 'raw'
%                       - 'raw'      : Simply apply imregionalmax/-min on
%                                      the projection. 
%                       - 'filtered' : Dilate/erode image to apply moving
%                                      max/min filter before applying
%                                      imregionalmax/-min. Note, this
%                                      will also smooth the correlation 
%                                      space when init_seg_method is set to
%                                      be 'corr', or the independent
%                                      components accordingly when it is 
%                                      set to be 'rICA'.
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
%                       - 'rmse' : Root median square error
%                       - 'corr' : Pearson correlation coefficient
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
%                                  refinement. Here, a is set to be the
%                                  5th quantile and b to be the 95th
%                                  quantile. The filter size for    
%                                  regmax_method is set to 1.
%
% 'corr_thresh'       : Threshold for the Pearson correlation coefficient
% (number)              when refinement_method is set to be 'corr'.
%                       Default: 0.85
%
% 'n_rICA'            : Number of features to extract during reconstruction 
% (integer)             independent component analysis.
%                       Default: 0 to use the full number of components.
%                       Otherwise provide a number larger than zero for
%                       either an undercomplete (n_rICA < number of frames)
%                       or overcomplete (n_rICA > number of frames) feature 
%                       representations.
%
% 'n_PCA'             : Percentage of the total variance of the data that
% (number)              should be kept if the projection_method is set to
%                       'pca'. The value must be larger than zero and not
%                       exceed 100.
%                       Default: 100
%
% 'fillmissing'       : Replace NaNs or Infs by linear interpolation of 
% (logical)             neighboring, nonmissing values.
%                       Default: false
%
% 'resume_CalciSeg'    : Map with granule IDs from a previous segmentation
% (2-D matrix)          session. Providing this input will skip all initial
%                       segmentation step and directly jump to refinement.
%                       Default: []
%
% 'pool_fragments'    : Sometimes, granules belonging to the same bio-
% (logical)             logical structure are fragmented into one or
%                       multiple pieces. To account for this, we
%                       calculate the Pearson correlation coefficient
%                       between neighboring granules and pool them if the
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
%                       every N consecutive frames before segmentation.
%                       Full-resolution time courses are used for the
%                       final summary statistics.
%                       Default: 1 (no binning)
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
%                      - projection           : Projection image used for
%                                               initial segmentation.
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
%                      - active_region.method : Method used for binarizing
%                      - quality.SNR          : Per-granule signal-to-noise
%                                               ratio (peak value over
%                                               baseline noise std)
%                      - quality.compactness  : Per-granule spatial
%                                               compactness (4*pi*area /
%                                               perimeter^2; 1 = circle)
%
%
% Version 1.2 | 13-Apr-26 (R2024a)
%
% =========================================================================

% Check whether this is the latest version
try
    checkVersion(...
        fileread([mfilename('fullpath'),'.m']),...
        webread('https://raw.githubusercontent.com/yannickguenzel/CalciSeg/main/code/CalciSeg.m'))
catch
    % Version check failed (e.g., no internet connection) - continue
end

% Validate inputs (only squeeze 4D+ inputs to collapse trailing singletons;
% squeeze on 3D data can destroy a valid singleton spatial dimension)
if ndims(stack) > 3
    stack = squeeze(stack);
end
opt = validateInputs(stack, varargin);

% Reduce precision of stack to be single
stack_size = whos('stack');
if ~strcmp(stack_size.class, 'single')
    stack = single(stack);
end%if not single precision

% Check for missing data
if any(isnan(stack(:))) || any(isinf(stack(:)))
    if opt.fillmissing
        stack(isinf(stack)) = NaN;
        % Fill missing data      
        stack = fillmissing(stack, "linear", 1); stack = fillmissing(stack, "nearest", 1);
        stack = fillmissing(stack, "linear", 2); stack = fillmissing(stack, "nearest", 2);
        if ndims(stack) == 3
            stack = fillmissing(stack, "linear", 3);
            stack = fillmissing(stack, "nearest", 3);
        end
    else
        error('The stack contains NaN or inf values. Please check your data or set the input "fillmissing" to true.')
    end%if fill
end%if missing data

% Temporal downsampling for faster segmentation
reshaped_stack_fullres = [];
if opt.temporal_bin > 1
    [x_s, y_s, t_s] = size(stack);
    t_binned = floor(t_s / opt.temporal_bin);
    % Guard: binning must leave at least 2 frames
    if t_binned < 2
        warning('CalciSeg:temporalBin', ...
            'Temporal binning (bin=%d, frames=%d) would leave %d frame(s). Disabling binning.', ...
            opt.temporal_bin, t_s, t_binned);
        opt.temporal_bin = 1;
    else
        reshaped_stack_fullres = reshape(stack, [x_s*y_s, t_s]);
        stack = reshape(stack(:,:,1:t_binned*opt.temporal_bin), [x_s, y_s, opt.temporal_bin, t_binned]);
        % Use mean then reshape instead of squeeze to prevent collapsing
        % the time dimension when t_binned happens to equal 1
        stack = reshape(mean(stack, 3), [x_s, y_s, t_binned]);
        disp(['Temporal binning: ', num2str(t_s), ' -> ', num2str(t_binned), ' frames (bin factor: ', num2str(opt.temporal_bin), ')'])
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
[stack, reshaped_stack, projection] = initialPreprocessing(stack, opt);

% Perform initial segmentation (or use a previous segmentation result)
if ~isempty(opt.resume_CalciSeg)
    granules_labeled = opt.resume_CalciSeg;
else
    granules_labeled = initialSegmentation(stack, reshaped_stack, projection, opt);
end

% Automatically estimate min/max granule size from the initial segmentation
if (isstring(opt.limitPixCount) || ischar(opt.limitPixCount)) && strcmp(opt.limitPixCount, 'auto')
    % Get pixel count for each region (exclude background ID 0)
    granuleList = unique(granules_labeled);
    [~,~,C] = unique(granules_labeled);
    C = histcounts(C, 1:length(granuleList)+1);
    C = [C(:), granuleList(:)];
    C = C(C(:,2) ~= 0, :);
    opt.limitPixCount = [floor(quantile(C(:,1), 0.05)), floor(quantile(C(:,1), 0.95))];
    disp(['automatic assessment of the min./max. granule size: [', num2str(opt.limitPixCount),']'])
end%if auto min size

% Enforce granule size constraints (split oversized, merge undersized).
% Note: removeHugeRegions runs first because splitting may create tiny
% fragments that removeTinyRegions then merges.
granules_labeled = removeHugeRegions(granules_labeled, opt.limitPixCount(2));
granules_labeled = removeTinyRegions(reshaped_stack, granules_labeled, opt);
if opt.n_rep>0
    granules_labeled = refineSegmentation(reshaped_stack, granules_labeled, opt);
end%if refine

% Optionally pool fragmented granules based on neighbor correlation
if opt.pool_fragments
    granules_labeled = poolFragmentedGranules(reshaped_stack, granules_labeled, opt);
end%if pool fragments

% Relabel granules with contiguous IDs (1, 2, ..., N)
granules_labeled = finalRefinement(granules_labeled);

% Compute summary statistics if requested
if nargout>1
    if opt.temporal_bin > 1
        [granules_labeled, summary_stats] = regionStats(reshaped_stack_fullres, granules_labeled, opt);
    else
        [granules_labeled, summary_stats] = regionStats(reshaped_stack, granules_labeled, opt);
    end
    summary_stats.projection = projection;
end%if summary stats

end%FCN:CalciSeg

% -------------------------------------------------------------------------

function checkVersion(localVersion, GitHubVersion)
% Compare the local file content against the latest GitHub version and warn
% the user if they differ.
if ~strcmp(localVersion, GitHubVersion)
    % Display a warning message
    msg = sprintf(['The local version of the CalciSeg does not match the latest version on GitHub. \n', ...
        'Please update your script to ensure you are using the most recent version. \n', ...
        'Latest version can be found at: https://github.com/YannickGuenzel/CalciSeg']);
    warning('CalciSeg:version','%s',msg);
end
end%FCN:checkVersion

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
opt.limitPixCount = [1 inf];
opt.corr_thresh = 0.85;
opt.n_rICA = 0;
opt.n_PCA = 100;
opt.fillmissing = false;
opt.resume_CalciSeg = [];
opt.pool_fragments = false;
opt.temporal_bin = 1;

% Check optional input arguments (Name-Value pairs)
varargin = varargin{1};
if ~isempty(varargin)
    if mod(length(varargin),2)>0
        error('Wrong number of arguments.')
    else
        % Iterate over all input arguments
        for iArg = 1:2:length(varargin)
            switch varargin{iArg}
                % Projection method -------------------------------------------
                case 'projection_method'
                    valid_projection_methods = {'std', 'mean', 'median', 'max', 'min', 'pca', 'corr', 'none'};
                    if ~ismember(varargin{iArg+1}, valid_projection_methods)
                        error('Invalid "projection_method". Valid options are: %s.', strjoin(valid_projection_methods, ', '));
                    else
                        opt.projection_method = varargin{iArg+1};
                    end
                    % Initial segmentation method ---------------------------------
                case 'init_seg_method'
                    valid_init_seg_methods = {'voronoi', 'corr', 'rICA'};
                    if ~ismember(varargin{iArg+1}, valid_init_seg_methods)
                        error(sprintf('Invalid "init_seg_method". Valid options are: %s.', strjoin(valid_init_seg_methods, ', ')));
                    else
                        opt.init_seg_method = varargin{iArg+1};
                    end
                    % Method for smoothing regional maxima ------------------------
                case 'regmax_method'
                    valid_regmax_methods = {'raw', 'filtered', 'both'};
                    if ~ismember(varargin{iArg+1}, valid_regmax_methods)
                        error('Invalid "regmax_method". Valid options are: %s.', strjoin(valid_regmax_methods, ', '));
                    else
                        opt.regmax_method = varargin{iArg+1};
                    end
                    % Number of refinement iterations ----------------------------
                case 'n_rep'
                    if ~isnumeric(varargin{iArg+1}) || ~isscalar(varargin{iArg+1}) || varargin{iArg+1} < 0 || mod(varargin{iArg+1},1) ~= 0
                        error('Input "n_rep" must be a non-negative scalar integer.');
                    else
                        opt.n_rep = varargin{iArg+1};
                    end
                    % Distance metric during refinement ---------------------------
                case 'refinement_method'
                    valid_refinement_methods = {'corr', 'rmse'};
                    if ~ismember(varargin{iArg+1}, valid_refinement_methods)
                        error('Invalid "refinement_method". Valid options are: %s.', strjoin(valid_refinement_methods, ', '));
                    else
                        opt.refinement_method = varargin{iArg+1};
                    end
                    % Minimum and maximum number of pixels per region -------------
                case 'limitPixCount'
                    val = varargin{iArg+1};
                    if ~isnumeric(val) && ~strcmp(val, 'auto')
                        error('Input "limitPixCount" must be a non-negative pair [min max] or "auto".');
                    elseif isnumeric(val)
                        if numel(val) ~= 2 || any(val < 0) || val(1) > val(2) || ~isfinite(val(1)) || any(mod(val(isfinite(val)),1) ~= 0)
                            error('Input "limitPixCount" must be [min max] with non-negative integers and finite min (max may be inf).');
                        end
                        opt.limitPixCount = val;
                    else
                        opt.limitPixCount = val;
                    end
                    % Correlation threshold ---------------------------------------
                case 'corr_thresh'
                    if ~isnumeric(varargin{iArg+1}) || ~isscalar(varargin{iArg+1}) || varargin{iArg+1} < 0 || varargin{iArg+1} > 1
                        error('Input "corr_thresh" must be a scalar number between 0 and 1.');
                    else
                        opt.corr_thresh = varargin{iArg+1};
                    end
                    % Number of independent components ----------------------------
                case 'n_rICA'
                    if ~isnumeric(varargin{iArg+1}) || ~isscalar(varargin{iArg+1}) || varargin{iArg+1} < 0 || mod(varargin{iArg+1},1) ~= 0
                        error('Input "n_rICA" must be a non-negative scalar integer.');
                    else
                        opt.n_rICA = varargin{iArg+1};
                    end
                case 'n_PCA'
                    if ~isnumeric(varargin{iArg+1}) || ~isscalar(varargin{iArg+1}) || varargin{iArg+1} <= 0 || varargin{iArg+1} > 100
                        error('Input "n_PCA" must be a scalar number larger than 0 and not exceeding 100.');
                    else
                        opt.n_PCA = varargin{iArg+1};
                    end
                case 'fillmissing'
                    if ~islogical(varargin{iArg+1}) || ~isscalar(varargin{iArg+1})
                        error('Input "fillmissing" must be a scalar logical (true or false).');
                    else
                        opt.fillmissing = varargin{iArg+1};
                    end
                case 'resume_CalciSeg'
                    val = varargin{iArg+1};
                    if ~isnumeric(val) || ~ismatrix(val) || size(val,1) ~= size(stack,1) || size(val,2) ~= size(stack,2)
                        error('Input "resume_CalciSeg" must be a 2-D matrix with the same spatial size as "stack".');
                    elseif any(val(:) < 0) || any(~isfinite(val(:))) || any(mod(val(:),1) ~= 0)
                        error('Input "resume_CalciSeg" must contain only finite, non-negative integers.');
                    else
                        opt.resume_CalciSeg = single(val);
                    end
                case 'pool_fragments'
                    if ~islogical(varargin{iArg+1}) || ~isscalar(varargin{iArg+1})
                        error('Input "pool_fragments" must be a scalar logical (true or false).');
                    else
                        opt.pool_fragments = varargin{iArg+1};
                    end
                case 'temporal_bin'
                    if ~isnumeric(varargin{iArg+1}) || ~isscalar(varargin{iArg+1}) || varargin{iArg+1} < 1 || mod(varargin{iArg+1},1) ~= 0
                        error('Input "temporal_bin" must be a positive integer.');
                    else
                        opt.temporal_bin = varargin{iArg+1};
                    end
                otherwise
                    error(['Unknown argument "',varargin{iArg},'"'])
            end%switch
        end%iArg
    end%if mod
end%if custom arguments

% Validate 'stack' dimensions and enforce parameter dependencies
if ndims(stack) > 3
    error('Number of dimensions of the input "stack" cannot exceed 3.')
elseif ismatrix(stack)
    warning('Input "stack" has no temporal component. Inputs are adjusted accordingly (projection_method="none"; init_segment_method="voronoi"; refinement_method="rmse"; temporal_bin=1).');
    opt.projection_method = 'none';
    opt.init_seg_method = 'voronoi';
    opt.refinement_method = 'rmse';
    opt.temporal_bin = 1;
end%if invalid 'stack' input
% For the region-growing approach, both projection_method and
% opt.init_seg_method must be set to 'corr'
if strcmp(opt.init_seg_method, 'corr')
    if ~strcmp(opt.projection_method, 'corr')
        warning('When the input "init_seg_method" is set to "corr", the input "projection_method" must be set to "corr" as well. We fixed this.')
        opt.projection_method = 'corr';
    end%if
end%if
% 'both' regmax_method is only applicable with Voronoi tessellation
if strcmp(opt.regmax_method, 'both') && ~strcmp(opt.init_seg_method, 'voronoi')
    warning('regmax_method "both" is only applicable with init_seg_method "voronoi". Setting to "raw".');
    opt.regmax_method = 'raw';
end%if
end%FCN:validateInputs

% -------------------------------------------------------------------------

function [stack, reshaped_stack, projection] = initialPreprocessing(stack, opt)
% Compute the temporal projection (collapsing the time dimension into a 2D
% image) and reshape the stack into a 2D pixels-by-time matrix.

% Get the size of the stack
[x, y, t] = size(stack);

% Initialize the projection variable
projection = [];

% Calculate the projection based on the specified method
switch opt.projection_method
    case 'std'
        projection = nanstd(stack, [], 3);
    case 'mean'
        projection = nanmean(stack, 3);
    case 'median'
        projection = nanmedian(stack, 3);
    case 'max'
        projection = nanmax(stack, [], 3);
    case 'min'
        projection = nanmin(stack, [], 3);
    case 'pca'
        % Project into PC space and keep components explaining n_PCA% variance
        [~, stack, ~, ~, explained] = pca(reshape(stack, [x*y, t]));
        if ~isempty(find(cumsum(explained)>=opt.n_PCA,1))
            t = find(cumsum(explained)>=opt.n_PCA,1);
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
        projection = zeros(x, y, 'single');
        for iPC = 1:t
            projection = projection + weights(iPC) * squeeze(stack(:,:,iPC));
        end%iPC
        disp(['PCA projection | ', num2str(t), ' components | Explained variance: ', num2str(round(sum(explained(1:t)),2)), '%'])       
    case 'corr'
        projection = correlationImage(stack, opt);
    case 'none'
        projection = squeeze(stack(:,:,1));
end%switch projection method

% Reshape the stack to a 2D matrix for further processing
reshaped_stack = reshape(stack, [x*y, t]);

end%FCN:initialPreprocessing

% -------------------------------------------------------------------------

function granules_labeled = initialSegmentation(stack, reshaped_stack, projection, opt)
% Perform the initial spatial segmentation using one of three methods:
% Voronoi tessellation, correlation-based region growing, or rICA.
switch opt.init_seg_method
    case 'voronoi'
        % Get local maxima and minima
        [regional_maxima, regional_minima] = identifyLocalExtrema(projection, opt);
        % 2-D Delaunay triangulation
        granules_labeled = triang_extrema(regional_maxima, regional_minima);
    case 'corr'
        granules_labeled = growing_regions(stack, projection, reshaped_stack, opt);
    case 'rICA'
        granules_labeled = rICA_segmentation(projection, reshaped_stack, opt);
        % Check for splits
        granules_labeled = checkSplits(granules_labeled);
end%switch opt.init_seg_method
end%FCN:initialSegmentation

% -------------------------------------------------------------------------

function granules_labeled = triang_extrema(regional_maxima, regional_minima)
% Segment the image via Voronoi tessellation seeded by local extrema.
% Each pixel is assigned to the nearest extremum centroid using knnsearch
% (equivalent to a 2-D Delaunay triangulation). Border regions are
% excluded by setting edge pixels to zero.

% Erode the extrema locations to keep the centroids only
ultimateErosion = bwulterode(regional_maxima>0) + bwulterode(regional_minima>0);
% Include points that are on the edge of the image. Later, we will set all
% regions that are at the edge to zero
ultimateErosion(1,:) = 1; ultimateErosion(:,1) = 1;
ultimateErosion(end,:) = 1; ultimateErosion(:,end) = 1;
% Segment image based on these maxima using a 2-D Delaunay
% triangulation
% Get linear indices of local maxima
linearIndices = find(ultimateErosion);
% Convert linear indices to (x, y) coordinates
[Cx, Cy] = ind2sub(size(ultimateErosion), linearIndices);
% Combine coordinates into Nx2 matrix
C = [Cx, Cy]; clear Cx Cy
% Segment image based on these maxima using a 2-D Delaunay
% triangulation (i.e. assign pixels to the closest extreme point
% Convert linear indices to (x, y) coordinates
linearIndices = 1:numel(regional_maxima);
[x, y] = ind2sub(size(regional_maxima), linearIndices);
% Combine all and adjust aspect ratio;
allVoxels = [x(:), y(:)];
clear x y linearIndices
% Use knnsearch to find the index of the closest point in 'C' for each voxel
granules_labeled = knnsearch(C, allVoxels);
% Reshape the resulting index array to the dimensions of the 2D image
granules_labeled = reshape(granules_labeled, size(regional_maxima));
% Remove border
ultimateErosion(1,:) = 0; ultimateErosion(:,1) = 0;
ultimateErosion(end,:) = 0; ultimateErosion(:,end) = 0;
linearIndices = find(ultimateErosion);
% Create a logical mask
mask = ismember(granules_labeled, granules_labeled(linearIndices));
granules_labeled(~mask) = 0;
end%FCN:triang_extrema

% -------------------------------------------------------------------------

function [regional_maxima, regional_minima] = identifyLocalExtrema(projection, opt)
% Identify local maxima and minima in the 2D projection image using one of
% three methods: raw (direct), filtered (morphological pre-smoothing), or
% both (union of raw and filtered).

% Initialize output variables
regional_maxima = [];
regional_minima = [];
SE = strel('disk', opt.filter_size);
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
% Get a map of local correlations of a pixel with each of its neighbors
% Uses a vectorized shifted-stack approach for efficiency.
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
neighbor_sum = zeros(rows, cols, t, 'single');
neighbor_count = zeros(rows, cols, 'single');
for iN = 1:length(dr)
    src_r = max(1, 1-dr(iN)):min(rows, rows-dr(iN));
    src_c = max(1, 1-dc(iN)):min(cols, cols-dc(iN));
    dst_r = src_r + dr(iN);
    dst_c = src_c + dc(iN);
    neighbor_sum(dst_r, dst_c, :) = neighbor_sum(dst_r, dst_c, :) + stack(src_r, src_c, :);
    neighbor_count(dst_r, dst_c) = neighbor_count(dst_r, dst_c) + 1;
end%iN
neighbor_mean = neighbor_sum ./ neighbor_count;
% Compute pixel-wise Pearson correlation between each pixel and its
% neighbor mean: r = sum((X-muX).*(Y-muY)) / sqrt(sum((X-muX)^2)*sum((Y-muY)^2))
X = reshape(stack, [rows*cols, t]);
Y = reshape(neighbor_mean, [rows*cols, t]);
X = X - mean(X, 2);
Y = Y - mean(Y, 2);
projection = sum(X .* Y, 2) ./ sqrt(sum(X.^2, 2) .* sum(Y.^2, 2));
% Replace NaN (from constant-signal pixels) with 0
projection(isnan(projection)) = 0;
projection = reshape(projection, [rows, cols]);
end%FCN:correlationImage

% -------------------------------------------------------------------------

function granules_labeled = growing_regions(stack, corr_img, reshaped_stack, opt)
% Greedy correlation-based region growing using BFS boundary tracking.
% Starting from the pixel with the highest local correlation, a region is
% grown by iteratively adding neighboring pixels whose time course
% correlates above corr_thresh with the region's mean time course. Once a
% region can no longer grow, the next highest-correlation unassigned pixel
% seeds a new region.

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
hWait = waitbar(0, ['Growing regions ... 0/',num2str(px_cnt)]);
for iPx = 1:px_cnt
    waitbar(iPx/px_cnt, hWait, ['Growing regions ... ', num2str(iPx),'/',num2str(px_cnt)])
    % Get the unassigned pixel with the highest local correlation
    [~, max_id] = max(corr_img(:));
    if corr_img(max_id) < opt.corr_thresh
        break;
    end%if below threshold
    % Seed a new region
    granules_labeled(max_id) = cnt_id;
    region_pixels = max_id;
    tc_sum = reshaped_stack(max_id, :);
    % Initialize boundary from the seed's 4-connected neighbors
    [sr, sc] = ind2sub([nRows, nCols], max_id);
    nr = [sr-1, sr+1, sr, sr];
    nc = [sc, sc, sc-1, sc+1];
    valid = nr>=1 & nr<=nRows & nc>=1 & nc<=nCols;
    boundary = sub2ind([nRows, nCols], nr(valid), nc(valid))';
    boundary = boundary(granules_labeled(boundary) == 0);
    % Grow region by testing boundary pixels against the mean TC
    while ~isempty(boundary) && length(region_pixels) < max_size
        curr_TC = tc_sum / length(region_pixels);
        % Limit candidates to respect max_size (skip when max_size is inf)
        if isfinite(max_size)
            n_available = max_size - length(region_pixels);
            if length(boundary) > n_available
                boundary = boundary(1:n_available);
            end
        end
        % Test each boundary pixel
        added_idx = [];
        for iB = 1:length(boundary)
            if granules_labeled(boundary(iB)) ~= 0
                continue;
            end%if already assigned
            px_tc = reshaped_stack(boundary(iB), :);
            if corr(curr_TC(:), px_tc(:)) >= opt.corr_thresh
                granules_labeled(boundary(iB)) = cnt_id;
                region_pixels = [region_pixels; boundary(iB)];
                tc_sum = tc_sum + px_tc;
                added_idx = [added_idx; boundary(iB)];
                if length(region_pixels) >= max_size
                    break;
                end%if max size reached
            end%if above threshold
        end%iB
        if isempty(added_idx)
            break;
        end%if no pixels added
        % Collect new boundary from added pixels' 4-connected neighbors
        new_boundary = [];
        for iA = 1:length(added_idx)
            [ar, ac] = ind2sub([nRows, nCols], added_idx(iA));
            nr = [ar-1, ar+1, ar, ar];
            nc = [ac, ac, ac-1, ac+1];
            valid = nr>=1 & nr<=nRows & nc>=1 & nc<=nCols;
            nb = sub2ind([nRows, nCols], nr(valid), nc(valid))';
            new_boundary = [new_boundary; nb(:)];
        end%iA
        new_boundary = unique(new_boundary);
        boundary = new_boundary(granules_labeled(new_boundary) == 0);
    end%while growing
    % Mark current region's pixels so they are not revisited
    corr_img(granules_labeled == cnt_id) = -inf;
    cnt_id = cnt_id + 1;
end%iPx
close(hWait)
end%FCN:growing_regions

% -------------------------------------------------------------------------

function granules_labeled = rICA_segmentation(projection, reshaped_stack, opt)
% Segment using reconstruction independent component analysis (rICA).
% The data is mean-subtracted, prewhitened, PCA-transformed, and then
% decomposed into independent components via rica(). Each pixel is assigned
% to the component with the strongest loading.
% Get size of the dataset
[x, y] = size(projection);
% Subtract mean since ICA cannot separate sources with a mean signal
% effect.
reshaped_stack = reshaped_stack-nanmean(reshaped_stack,1);
% Prewhiten the signal so that it has zero mean and identity covariance.
reshaped_stack = prewhiten(reshaped_stack);
% Apply PCA
[~,reshaped_stack] = pca(reshaped_stack);
% Check how many components to request
if opt.n_rICA==0
    opt.n_rICA = size(reshaped_stack, 2);
end
% Apply reconstruction ICA (rICA)
rng_state = rng;
rng_cleanup = onCleanup(@() rng(rng_state));% Restore RNG even if rica throws
rng(1);% For reproducibility
ricaTransform = rica(reshaped_stack, opt.n_rICA, 'VerbosityLevel', 1);
clear rng_cleanup;% Normal exit: restore now rather than at scope end
% Transform the data using the fitted model
transformedData = transform(ricaTransform, reshaped_stack);
% Correct sign of each component so that the dominant deflection is
% positive. This is done by z-scoring each component, measuring the mean
% positive excursion beyond +2 SD and the mean negative excursion beyond
% -2 SD, and flipping the sign if the negative tail dominates.
TD = transformedData;
s = std(TD);
s(s == 0) = 1;% Guard zero-variance components against division by zero
TDv = TD ./ repmat(s, size(TD,1), 1);             % z-score
TDzp = TDv - 2;  TDzp(TDzp < 0) = 0;              % positive tail (>+2 SD)
TDzpm = mean(TDzp);
TDzn = TDv + 2;  TDzn(TDzn > 0) = 0;              % negative tail (<-2 SD)
TDznm = mean(TDzn);
flip_sign = sign(TDzpm + TDznm);
flip_sign(flip_sign == 0) = 1;% Ambiguous sign defaults to no flip
transformedData = TD .* repmat(flip_sign, size(TD, 1), 1);
% Reshape the independent components back into a 3-D format
reshaped_stack = reshape(transformedData, [x, y, size(transformedData,2)]);
% Smooth components
if strcmp(opt.regmax_method, 'filtered')
    reshaped_stack = imboxfilt3(reshaped_stack, [opt.filter_size, opt.filter_size, 1]);
end%if filter
% Get the component with the strongest effect on a given pixel
[~, granules_labeled] = max(reshaped_stack, [], 3);
end%FCN:rICA_segmentation

% -------------------------------------------------------------------------

function Z = prewhiten(X)
    % From https://de.mathworks.com/help/stats/extract-mixed-signals.html
    % 1. Size of X.
    [N,P] = size(X); %#ok<ASGLU> — N not used; SVD rank handles wide matrices
    % 2. SVD of covariance of X. We could also use svd(X) to proceed but N
    % can be large and so we sacrifice some accuracy for speed.
    [U,Sig] = svd(cov(X));
    Sig     = diag(Sig);
    Sig     = Sig(:)';
    % 3. Figure out which values of Sig are non-zero.
    tol = eps(class(X));
    idx = (Sig > max(Sig)*tol);
    assert(~all(idx == 0));
    % 4. Get the non-zero elements of Sig and corresponding columns of U.
    Sig = Sig(idx);
    U   = U(:,idx);
    % 5. Compute prewhitened data.
    mu = mean(X,1);
    Z = X-mu;
    Z = (Z*U)./sqrt(Sig);
end%FCN:prewhiten

% -------------------------------------------------------------------------

function granules_labeled = refineSegmentation(reshaped_stack, granules_labeled, opt)
% Iteratively refine granule boundaries. In each iteration, border pixels
% are re-assigned to the neighboring granule with the best fit (lowest RMSE
% or highest correlation). After each re-assignment pass, split granules
% are detected and size constraints are enforced. Stops early if the label
% map converges.
previous_granules_labeled = -ones(size(granules_labeled));
hWait = waitbar(0, ['Refine regions ... 0/',num2str(opt.n_rep)]);
for iRep = 1:opt.n_rep
    waitbar(iRep/opt.n_rep, hWait, ['Refine regions ... ', num2str(iRep),'/',num2str(opt.n_rep)])
    % Get avg signal trace per granule (exclude background ID 0)
    granuleList = unique(granules_labeled);
    granuleList(granuleList == 0) = [];
    granuleTC = nan(length(granuleList), size(reshaped_stack,2));
    for iGranule = 1:length(granuleList)
        idx = find(granules_labeled==granuleList(iGranule));
        granuleTC(iGranule, :) = nanmean(reshaped_stack(idx,:),1);
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
    switch opt.refinement_method
        case 'corr'
            idx_candidates_newIdentiy = refinement_parfor_corr(idx_candidates, granules_labeled, reshaped_stack, granuleTC, granuleList);
        case 'rmse'
            idx_candidates_newIdentiy = refinement_parfor_rmse(idx_candidates, granules_labeled, reshaped_stack, granuleTC, granuleList);
    end
    % Update changed identities
    granules_labeled(idx_candidates) = idx_candidates_newIdentiy;
    % Account for splitting of granules. For this, iterate over
    % all granules and check whether the corresponding identity
    % can be found in more than one coherent region
    granules_labeled = checkSplits(granules_labeled);
    % Enforce size constraints after re-assignment
    granules_labeled = removeHugeRegions(granules_labeled, opt.limitPixCount(2));
    granules_labeled = removeTinyRegions(reshaped_stack, granules_labeled, opt);
    % Stop early if the label map has converged
    if sum(previous_granules_labeled(:) == granules_labeled(:)) == numel(granules_labeled(:))
        break
    end%if converged
end%iRep
close(hWait)
% Make sure the highest ID is equal to the number of granules
has_background = any(granules_labeled(:) == 0);
[~, ~, newIDs] = unique(granules_labeled);
if has_background
    newIDs = newIDs - 1;
end
granules_labeled = reshape(newIDs, size(granules_labeled));
end%FCN:refineSegmentation

% -------------------------------------------------------------------------

function idx_candidates_newIdentiy = refinement_parfor_corr(idx_candidates, granules_labeled, reshaped_stack, granuleTC, granuleList)
% Re-assign border pixels using Pearson correlation distance (|1 - r|).
% Each candidate pixel is compared against its 4-connected neighboring
% granules and assigned to the one with the smallest correlation distance.

% Prepare output of new labels
idx_candidates_newIdentiy = granules_labeled(idx_candidates);
[nRows, nCols] = size(granules_labeled);
% Use a parallel for loop to speed things up
parfor iCan = 1:length(idx_candidates)
    % Get 4-connected neighbor indices directly
    [pr, pc] = ind2sub([nRows, nCols], idx_candidates(iCan));
    nr = [pr-1, pr+1, pr, pr];
    nc = [pc, pc, pc-1, pc+1];
    valid = nr>=1 & nr<=nRows & nc>=1 & nc<=nCols;
    neighbor_idx = sub2ind([nRows, nCols], nr(valid), nc(valid));
    % Get neighborhood clusters (exclude background ID 0)
    nC = unique(granules_labeled(neighbor_idx));
    nC(nC == 0) = [];
    if isempty(nC); continue; end
    % Get the pixel's TC
    activity_px = reshaped_stack(idx_candidates(iCan), :);
    % Check which neighborhood is a better fit
    voronoi_R = zeros(length(nC),1);
    for iN = 1:length(nC)
        activity_cluster = granuleTC(granuleList == nC(iN), :);
        rr = corrcoef(activity_px, activity_cluster); voronoi_R(iN) = abs(diff([1, rr(2)]));
    end%iN
    % Assign new identity (keep current if all distances are NaN)
    voronoi_R(~isfinite(voronoi_R)) = NaN;
    if all(isnan(voronoi_R)); continue; end
    [~, iBest] = min(voronoi_R);
    idx_candidates_newIdentiy(iCan) = nC(iBest);
end%iCan
end%FCN:refinement_parfor_corr

% -------------------------------------------------------------------------

function idx_candidates_newIdentiy = refinement_parfor_rmse(idx_candidates, granules_labeled, reshaped_stack, granuleTC, granuleList)
% Re-assign border pixels using root-median-square error (RMSE).
% Each candidate pixel is compared against its 4-connected neighboring
% granules and assigned to the one with the smallest RMSE.

% Prepare output of new labels
idx_candidates_newIdentiy = granules_labeled(idx_candidates);
[nRows, nCols] = size(granules_labeled);
% Use a parallel for loop to speed things up
parfor iCan = 1:length(idx_candidates)
    % Get 4-connected neighbor indices directly
    [pr, pc] = ind2sub([nRows, nCols], idx_candidates(iCan));
    nr = [pr-1, pr+1, pr, pr];
    nc = [pc, pc, pc-1, pc+1];
    valid = nr>=1 & nr<=nRows & nc>=1 & nc<=nCols;
    neighbor_idx = sub2ind([nRows, nCols], nr(valid), nc(valid));
    % Get neighborhood clusters (exclude background ID 0)
    nC = unique(granules_labeled(neighbor_idx));
    nC(nC == 0) = [];
    if isempty(nC); continue; end
    % Get the pixel's TC
    activity_px = reshaped_stack(idx_candidates(iCan), :);
    % Check which neighborhood is a better fit
    voronoi_R = zeros(length(nC),1);
    for iN = 1:length(nC)
        activity_cluster = granuleTC(granuleList == nC(iN), :);
        voronoi_R(iN) = sqrt(median((activity_cluster-activity_px).^2));
    end%iN
    % Assign new identity (keep current if all distances are NaN)
    voronoi_R(~isfinite(voronoi_R)) = NaN;
    if all(isnan(voronoi_R)); continue; end
    [~, iBest] = min(voronoi_R);
    idx_candidates_newIdentiy(iCan) = nC(iBest);
end%iCan
end%FCN:refinement_parfor_rmse

% -------------------------------------------------------------------------

function granules_labeled = checkSplits(granules_labeled)
% Detect and resolve granules that have been split into multiple
% disconnected components. Each non-contiguous fragment beyond the first
% receives a new unique ID. Single-pixel fragments are merged into the
% nearest neighboring granule instead.

% Identify splitter candidates (granules with >1 connected component)
unique_granules = unique(granules_labeled);
unique_granules(unique_granules == 0) = [];
is_splitter = false(length(unique_granules), 1);
parfor iP = 1:length(unique_granules)
    % Get region mask and count connected components
    region_mask = (granules_labeled == unique_granules(iP));
    [~, num] = bwlabel(region_mask, 4);
    is_splitter(iP) = (num > 1);
end%iP
splitter_candidates = unique_granules(is_splitter);
% Iterate over all splitters and account for them
for iGranule = 1:length(splitter_candidates)
    % Get b/w image and identify number of regions
    bw = granules_labeled==splitter_candidates(iGranule);
    L = bwlabel(bw,4);
    % If there is more than one region (i.e. more than
    % the labels 0 and 1, assign new identities
    unique_L = unique(L(:));
    if length(unique_L)>2
        % Iterate over all splitter
        for iSplitter = 3:length(unique_L)
            idx = find(L == unique_L(iSplitter));
            % If the new region is just a single pixel, do not make it
            % a new region
            if length(idx)==1
                % Get pixels on the edge
                neighbors = imgradient(L == unique_L(iSplitter), 'central')>0;
                neighbors(bw) = 0;
                neighbor_IDs = unique(granules_labeled(neighbors));
                neighbor_IDs(neighbor_IDs == 0) = [];
                if isempty(neighbor_IDs)
                    granules_labeled(idx) = max(unique_granules) + 1;
                    unique_granules = [unique_granules(:); max(unique_granules) + 1];
                else
                    granules_labeled(idx) = mode(neighbor_IDs(:));
                end
            else
                granules_labeled(idx) = max(unique_granules)+1;
                unique_granules = [unique_granules(:); max(unique_granules)+1];
            end
        end%iSplitter
    end%if more than one connected component
end%iGranule
end%FCN:checkSplits

% -------------------------------------------------------------------------

function granules_labeled = removeTinyRegions(reshaped_stack, granules_labeled, opt)
% Merge granules smaller than limitPixCount(1) into their best-fitting
% neighbor. Uses a batch approach: tiny regions that share no neighbors
% are merged in parallel within a single pass. The outer loop repeats
% until no granule is below the minimum size.
while true
    % Get pixel count for each region (exclude background ID 0)
    granuleList = unique(granules_labeled);
    granuleList(granuleList == 0) = [];
    if isempty(granuleList); break; end
    pixCounts = zeros(length(granuleList), 1);
    for iG = 1:length(granuleList)
        pixCounts(iG) = sum(granules_labeled(:) == granuleList(iG));
    end%iG
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
    neighbor_sets = cell(n_small, 1);
    for iS = 1:n_small
        bw_s = granules_labeled == small_granules(iS);
        border_s = imgradient(bw_s, 'central') > 0;
        border_s(bw_s) = 0;
        nIDs = unique(granules_labeled(border_s));
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
        curr_TC = nanmean(reshaped_stack(ind_pix, :), 1);
        % Get pixels on the edge (exclude background ID 0)
        neighbors = imgradient(bw, 'central')>0;
        neighbors(bw) = 0;
        neighbor_IDs = unique(granules_labeled(neighbors));
        neighbor_IDs(neighbor_IDs == 0) = [];
        % Guard: no valid neighbors — assign to background
        if isempty(neighbor_IDs)
            granules_labeled(ind_pix) = 0;
            continue;
        end
        % Get the neighbors' time courses
        TCs = nan(length(neighbor_IDs), size(reshaped_stack,2));
        for iTC = 1:length(neighbor_IDs)
            ind = find(granules_labeled==neighbor_IDs(iTC));
            TCs(iTC,:) = nanmean(reshaped_stack(ind,:),1);
        end%iTC
        % Check which neighboring granule is the best fit
        id_table = nan(1, length(neighbor_IDs));
        for iN = 1:length(neighbor_IDs)
            % Get the rmse or distance based on correlation
            switch opt.refinement_method
                case 'corr'
                    min_dist = corrcoef(curr_TC, TCs(iN,:)); min_dist = abs(diff([1, min_dist(2)]));
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
    % Get pixel count for each region (exclude background ID 0)
    granuleList = unique(granules_labeled);
    granuleList(granuleList == 0) = [];
    if isempty(granuleList); break; end
    pixCounts = zeros(length(granuleList), 1);
    for iG = 1:length(granuleList)
        pixCounts(iG) = sum(granules_labeled(:) == granuleList(iG));
    end%iG
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
        bw_cut = bwlabel(bw_cut);
        % Strategy 2: If watershed did not split the region (e.g.,
        % circular shapes), try cutting perpendicular to the major axis.
        if max(bw_cut(:)) <= 1
            % Get the fitted ellipse properties
            stats = regionprops(bw, 'Centroid', 'MajorAxisLength', 'MinorAxisLength', 'Orientation');
            ellipse = stats(1);
            longAxis = ellipse.MajorAxisLength;
            % Get the position and angle of the ellipse
            centerX = ellipse.Centroid(1);
            centerY = ellipse.Centroid(2);
            angle = ellipse.Orientation;
            % Define a line perpendicular to the major axis
            x1 = floor(centerX + (longAxis/2) * sind(angle));
            y1 = floor(centerY - (longAxis/2) * cosd(angle));
            x2 = ceil(centerX - (longAxis/2) * sind(angle));
            y2 = ceil(centerY + (longAxis/2) * cosd(angle));
            % Draw a line
            bw_cut = bw;
            numPoints = round(sqrt(size(bw,1)^2+size(bw,2)^2));
            x = round(linspace(x1, x2, numPoints)); x(x<1)=1; x(x>size(bw_cut,2))=size(bw_cut,2);
            y = round(linspace(y1, y2, numPoints)); y(y<1)=1; y(y>size(bw_cut,1))=size(bw_cut,1);
            for iPx = 1:numPoints
                bw_cut(y(iPx), x(iPx)) = 0;
            end%iPx
            bw_cut = bwlabel(bw_cut,4);
            % Strategy 3: If perpendicular cut did not work, try
            % cutting along the major axis
            if max(bw_cut(:)) <= 1
                % Define a line along the major axis
                x1 = floor(centerX - (longAxis/2) * cosd(angle));
                y1 = floor(centerY - (longAxis/2) * sind(angle));
                x2 = ceil(centerX + (longAxis/2) * cosd(angle));
                y2 = ceil(centerY + (longAxis/2) * sind(angle));
                % Draw a line
                bw_cut = bw;
                numPoints = round(sqrt(size(bw,1)^2+size(bw,2)^2));
                x = round(linspace(x1, x2, numPoints)); x(x<1)=1; x(x>size(bw_cut,2))=size(bw_cut,2);
                y = round(linspace(y1, y2, numPoints)); y(y<1)=1; y(y>size(bw_cut,1))=size(bw_cut,1);
                for iPx = 1:numPoints
                    bw_cut(y(iPx), x(iPx)) = 0;
                end%iPx
                bw_cut = bwlabel(bw_cut,4);
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
        for iPx = 1:length(ind)
            % Get the current pixel's neighborhood
            neighborhood = zeros(size(bw));
            neighborhood(ind(iPx)) = 1;
            neighbors = imgradient(neighborhood, 'central')>0;
            neighbors(ind(iPx)) = 0;
            neighbor_IDs = unique(bw_cut(neighbors));
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

% -------------------------------------------------------------------------

function granules_labeled = poolFragmentedGranules(reshaped_stack, granules_labeled, opt)
% Merge neighboring granules whose mean time courses are highly correlated
% (above corr_thresh). Uses an adjacency-first approach: the spatial
% neighbor graph is built from the label map, correlations are computed
% only for adjacent pairs, and connected components above the threshold
% are merged. Repeats until no pair exceeds the threshold.
while true
    % Get a list of granules (exclude background ID 0)
    granuleList = unique(granules_labeled);
    granuleList(granuleList == 0) = [];
    if length(granuleList) < 2; break; end
    % Get each granule's time course
    avgTCs = nan(length(granuleList), size(reshaped_stack,2));
    for iGranule = 1:length(granuleList)
        idx_granule = find(granules_labeled==granuleList(iGranule));
        avgTCs(iGranule,:) = nanmean(reshaped_stack(idx_granule,:),1);
    end%iGranule
    % Build adjacency from label map (8-connected neighborhood)
    % Horizontal pairs
    left = granules_labeled(:, 1:end-1);
    right = granules_labeled(:, 2:end);
    diff_mask = left ~= right & left > 0 & right > 0;
    pairs_h = [left(diff_mask), right(diff_mask)];
    % Vertical pairs
    top = granules_labeled(1:end-1, :);
    bottom = granules_labeled(2:end, :);
    diff_mask = top ~= bottom & top > 0 & bottom > 0;
    pairs_v = [top(diff_mask), bottom(diff_mask)];
    % Diagonal pairs (top-left to bottom-right)
    tl = granules_labeled(1:end-1, 1:end-1);
    br = granules_labeled(2:end, 2:end);
    diff_mask = tl ~= br & tl > 0 & br > 0;
    pairs_d1 = [tl(diff_mask), br(diff_mask)];
    % Diagonal pairs (top-right to bottom-left)
    tr = granules_labeled(1:end-1, 2:end);
    bl = granules_labeled(2:end, 1:end-1);
    diff_mask = tr ~= bl & tr > 0 & bl > 0;
    pairs_d2 = [tr(diff_mask), bl(diff_mask)];
    % Get unique adjacent pairs
    all_pairs = unique(sort([pairs_h; pairs_v; pairs_d1; pairs_d2], 2), 'rows');
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
    % Get all connected IDs via graph-based component analysis
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
% Relabel all granules with contiguous integer IDs (1, 2, ..., N) so that
% the highest ID equals the total number of granules. Background (0) is
% preserved as 0.
has_background = any(granules_labeled(:) == 0);
[~, ~, granules_RElabeled] = unique(granules_labeled);
if has_background
    granules_RElabeled = granules_RElabeled - 1;
end
granules_RElabeled = reshape(granules_RElabeled, size(granules_labeled));
end%FCN:finalRefinement

% -------------------------------------------------------------------------

function [granules_labeled, summary_stats] = regionStats(reshaped_stack, granules_labeled, opt)
% Compute summary statistics for each granule: mean time course, within-
% granule correlation, per-granule projection images (mean, std, min, max,
% correlation), and quality metrics (SNR, compactness). Also estimates an
% active-region binary map using Otsu's method on the metric matching the
% projection method.

% Get list of granules (exclude background ID 0)
granuleList = unique(granules_labeled);
granuleList(granuleList == 0) = [];
% Preallocation
summary_stats.avgTCs = nan(length(granuleList), size(reshaped_stack,2));
summary_stats.granule_Corr = nan(length(granuleList), 1);
summary_stats.granule_Avg_img = nan(size(granules_labeled));
summary_stats.granule_Min_img = nan(size(granules_labeled));
summary_stats.granule_Max_img = nan(size(granules_labeled));
summary_stats.granule_Std_img = nan(size(granules_labeled));
summary_stats.granule_Median_img = nan(size(granules_labeled));
summary_stats.granule_Corr_img = nan(size(granules_labeled));
summary_stats.quality.SNR = nan(length(granuleList), 1);
summary_stats.quality.compactness = nan(length(granuleList), 1);
% Iterate over all granules
for iGranule = 1:length(granuleList)
    % Get their index positions
    idx_granule = find(granules_labeled==granuleList(iGranule));
    % Get the avg time course of the current granule
    summary_stats.avgTCs(iGranule,:) = nanmean(reshaped_stack(idx_granule,:),1);
    % Get within-granule correlation
    summary_stats.granule_Corr(iGranule,1) = calculateAverageCorrelation(summary_stats.avgTCs(iGranule,:), reshaped_stack(idx_granule,:));
    % Project different metrics into 2D (one scalar per granule)
    summary_stats.granule_Avg_img(idx_granule)  = nanmean(summary_stats.avgTCs(iGranule,:));
    summary_stats.granule_Min_img(idx_granule)  = nanmin(summary_stats.avgTCs(iGranule,:));
    summary_stats.granule_Max_img(idx_granule)  = nanmax(summary_stats.avgTCs(iGranule,:));
    summary_stats.granule_Std_img(idx_granule)  = nanstd(summary_stats.avgTCs(iGranule,:));
    summary_stats.granule_Median_img(idx_granule) = nanmedian(summary_stats.avgTCs(iGranule,:));
    summary_stats.granule_Corr_img(idx_granule) = summary_stats.granule_Corr(iGranule,1);
end%iGranule
% Compute per-granule quality metrics
for iGranule = 1:length(granuleList)
    % SNR: peak dF/F over baseline noise standard deviation
    tc = summary_stats.avgTCs(iGranule, :);
    baseline = prctile(tc, 20);
    below_median = tc(tc <= prctile(tc, 50));
    if length(below_median) > 1
        noise_std = std(below_median);
    else
        noise_std = NaN;
    end
    if noise_std > 0
        summary_stats.quality.SNR(iGranule) = (max(tc) - baseline) / noise_std;
    end
    % Compactness: 4*pi*area / perimeter^2 (1 = perfect circle)
    % Sum across all connected components to handle pooled/fragmented ROIs
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
end%iGranule quality
% Estimate active regions via Otsu's method on the metric that matches
% the projection method
switch opt.projection_method
    case 'std',    active_img = summary_stats.granule_Std_img;
    case 'mean',   active_img = summary_stats.granule_Avg_img;
    case 'median', active_img = summary_stats.granule_Median_img;
    case 'max',    active_img = summary_stats.granule_Max_img;
    case 'min',    active_img = summary_stats.granule_Min_img;
    case 'pca',    active_img = summary_stats.granule_Std_img;
    case 'corr',   active_img = summary_stats.granule_Corr_img;
    case 'none',   active_img = [];
end%switch projection method
summary_stats.active_region.method = opt.projection_method;
if ~isempty(active_img)
    % Mask NaN pixels before binarization (imbinarize does not handle NaN)
    nan_mask = isnan(active_img);
    active_img(nan_mask) = 0;
    summary_stats.active_region.map = imbinarize(active_img);
    summary_stats.active_region.map(nan_mask) = false;
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
    X = ind_tc - mean(ind_tc, 2);
    Y = avg_tc - mean(avg_tc);
    corr_values = (X * Y') ./ (sqrt(sum(X.^2, 2)) * sqrt(sum(Y.^2)));
    % Sanitize non-finite values (from zero-variance traces) before averaging
    corr_values(~isfinite(corr_values)) = NaN;
    overallAvg = nanmean(corr_values);
else
    overallAvg = NaN;
end%if more dimensions
end%FCN:calculateAverageCorrelation