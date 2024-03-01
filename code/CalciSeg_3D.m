function [granules_labeled, summary_stats] = CalciSeg_3D(stack, varargin)
% [granules_labeled, summary_stats] = CalciSeg_3D(stack, Name, Value)
% segments the spatial aspects of a 4-D stack (x*y*z*time) into individual
% regions ('granules') based on a refined 3D Delaunay triangulation.
% Note: CalciSeg uses parallel for-loops at for the local growing
% algorithm and the refinement process. It is recommended to start a
% parallel pool of workers before calling the function.
%
%
% OPTIONAL INPUT PARAMETERS:
% Parameter name        Values and description
% =========================================================================
%
% 'aspect_ratio'        Aspect ratio of the the data (x, y, z)
% (array)               Default: [1 1 1]
%
% 'projection_method'   Method for calculating the projection
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
%                                    between  neighboring pixels
%                       - 'none'   : Applicable when "stack" has no fourth 
%                                    dimension or when to take the first 
%                                    slice of the 4-D stack
%
% 'init_seg_method'    Method for the initial segmentation.
% (string)             Default: 'voronoi'
%                       - 'voronoi' : Delaunay triangulation
%                       - 'corr'    : Local growing based on correlation.
%                                   : This will set 'projection_method' to
%                                     be 'corr'
%                       - 'rICA'    : Reconstruction independent component 
%                                     analysis.
%
% 'regmax_method'       Method for determining how to identify local
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
% 'n_rep'               Number of iterations for refining the regions.
% (integer)             Default: 0
%
% 'refinement_method'   Measure to fine-tune granule assignment during 
% (string)              refinement iterations.
%                       Default: 'rmse'
%                       - 'rmse' : Root median square error
%                       - 'corr' : Pearson  correlation coefficient
%                       
% 'minPixCount'     :   Limits the pixel area per granule
% (matrix or string)    Default:  1
%                       - [a]    : The minimum (a) number of pixels that 
%                                  can be assigned to a  granule. Note 
%                                  that, it affects the filter size for the 
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
%
% OUTPUT
% Parameter name      Values and description
% =========================================================================
%
% granules_labeled    : Map with granule IDs. The size of granules_labeled
%                       corresponds to the first two dimensions of stack.
%
% summary_stats       : Summary statistics for each granule and additional
%                       inforamtion on the segmentation process saved as
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
%                      - granule_Corr_img     : Within-granule correlation
%                      - active_region.map    : Binary map of active 
%                                               regions based on Otsu's
%                                               method
%                      - active_region.method : Method used for binarizing
%
%
% Version: 29-Feb-24 (R2023a)
% =========================================================================

% Validate inputs
stack = squeeze(stack);
opt = validateInputs(stack, varargin);

% Reduce precision of stack to be single
stack_size = whos('stack');
if ~strcmp(stack_size.class, 'single')
    stack = single(stack);
end%if not single precision

% Define size for local neighborhood
if isnumeric(opt.minPixCount)
    opt.filter_size = ceil(sqrt(opt.minPixCount(1)/pi)) + double(opt.minPixCount(1)==0);
elseif strcmp(opt.minPixCount, 'auto')
    opt.filter_size = 1;
end

% Initial preprocessing
[reshaped_stack, projection] = initialPreprocessing(stack, opt);

% Perform initial segmentation
granules_labeled = initialSegmentation(stack, reshaped_stack, projection, opt);

% Maybe, the user wants to have the min/max size estimated based on the data
if (isstring(opt.minPixCount) || ischar(opt.minPixCount)) && strcmp(opt.minPixCount, 'auto')
    % Get pixel count for each region
    granuleList = unique(granules_labeled);
    [~,~,C] = unique(granules_labeled);
    C = histcounts(C, 1:length(granuleList)+1);
    C = [C(:), granuleList(:)];
    opt.minPixCount = floor(quantile(C(2:end,1), 0.05));
    disp(['automatic assessment of the min. granule size: ', num2str(opt.minPixCount)])
end%if auto min size

% Refine segmentation
granules_labeled = removeTinyRegions(granules_labeled, opt);
if opt.n_rep>0
    granules_labeled = refineSegmentation(reshaped_stack, granules_labeled, opt);
end%if refine

% Final refinement steps and statistics calculation
if nargout>1
    [granules_labeled, summary_stats] = finalRefinementAndStats(reshaped_stack, granules_labeled, opt);
    summary_stats.projection = projection;
end%if summary stats

end%FCN:CalciSeg_3D

% -------------------------------------------------------------------------

function opt = validateInputs(stack, varargin)
% First, set all parameters to their default value
opt.aspect_ratio = [1 1 1];
opt.projection_method = 'std';
opt.init_seg_method = 'voronoi';
opt.regmax_method = 'raw';
opt.n_rep = 50;
opt.refinement_method = 'rmse';
opt.minPixCount = 1;
opt.corr_thresh = 0.85;
opt.n_rICA = 0;

% Now, check optional input variable
varargin = varargin{1};
if ~isempty(varargin)
    if mod(length(varargin),2)>0
        error('Wrong number of arguments.')
    else
        % Iterate over all inuput arguments
        for iArg = 1:2:length(varargin)
            switch varargin{iArg}
                % Aspect ratio ------------------------------------------------
                case 'aspect_ratio'
                    if ~isnumeric(varargin{iArg+1}) || length(varargin{iArg+1}) ~= 3 || isempty(varargin{iArg+1}) || any(varargin{iArg+1}<=0)
                        error('Input "aspect_ratio" has to be a vector of positive values with a length of 3.');
                    else
                        opt.aspect_ratio = varargin{iArg+1};
                    end
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
                        error(sprintf('Invalid "projection_method". Valid options are: %s.', strjoin(valid_projection_methods, ', ')));
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
                    % Number of refinement interations ----------------------------
                case 'n_rep'
                    if ~isnumeric(varargin{iArg+1}) || varargin{iArg+1} < 0
                        error('Input "n_rep" must be a non-negative integer value.');
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
                case 'minPixCount'
                    if ~isnumeric(varargin{iArg+1}) && ~strcmp(varargin{iArg+1}, 'auto')
                        error( 'Input "minPixCount" must be a non-negative pair of two integer values or "auto".');
                    elseif isnumeric(varargin{iArg+1}) && length(varargin{iArg+1}) ~= 2
                        error('Input "minPixCount" must be a non-negative pair of two integer values or "auto".');
                    else
                        opt.minPixCount = varargin{iArg+1};
                    end
                    % Correlation threshold ---------------------------------------
                case 'corr_thresh'
                    if ~isnumeric(varargin{iArg+1}) || varargin{iArg+1} < 0 || varargin{iArg+1} > 1
                        error('Input "corr_thresh" must be a number between 0 and 1.');
                    else
                        opt.corr_thresh = varargin{iArg+1};
                    end
                    % Number of independent components ----------------------------
                case 'n_rICA'
                    if ~isnumeric(varargin{iArg+1}) || varargin{iArg+1} < 0
                        error('Input "n_rICA" must be an integer larger than 0.');
                    else
                        opt.n_rICA = round(varargin{iArg+1});
                    end
                otherwise
                    error(['Unknown argument "',varargin{iArg},'"'])
            end%switch
        end%iArg
    end%if mod
end%if custom arguments

% Last, validate 'stack' input and potentially adjust some input parameters
if ndims(stack) < 2
    error('Input "stack" has to be at least a 2-D (x,y) matrix.')
elseif ndims(stack) > 3
    error('Number of dimension of the input "stack" cannot exceed 4.')
elseif ismatrix(stack)
    warning('Input "stack" has no temporal component. Inputs are adjusted accordingly (projection_method="none"; init_segment_method="voronoi"; refinement_method="rmse").');
    opt.projection_method = 'none';
    opt.init_seg_method = 'voronoi';
    opt.refinement_method = 'rmse';
end%if invalid 'stack' input
% For the region-growing approach, both projection_method and
% opt.init_seg_method must be set to 'corr'
if strcmp(opt.init_seg_method, 'corr')
    if ~strcmp(opt.projection_method, 'corr')
        warning('When the input "init_seg_method" is set to "corr", the input "projection_method" must be set to "corr" as well. We fixed this.')
        opt.projection_method = 'corr';
    end%if
end%if
end%FCN:validateInputs

% -------------------------------------------------------------------------

function [reshaped_stack, projection] = initialPreprocessing(stack, opt)
% Get the size of the stack
[x, y, z, t] = size(stack);

% Initialize the projection variable
projection = [];

% Calculate the projection based on the specified method
switch opt.projection_method
    case 'std'
        projection = nanstd(stack, [], 4);
    case 'mean'
        projection = nanmean(stack, 4);
    case 'median'
        projection = nanmedian(stack, 4);
    case 'max'
        projection = nanmax(stack, [], 4);
    case 'min'
        projection = nanmin(stack, [], 4);
    case 'pca'
        [~, stack, ~, ~, explained] = pca(reshape(stack, [x*y*z, t]));
        if ~isempty(find(cumsum(explained)>=99.99,1))
            t = find(cumsum(explained)>=99.99,1);
        end%if
        stack = reshape(stack(:,1:t), [x, y, z, t]);
        projection = squeeze(stack(:,:,:,1));
        disp(['PC1 | Explained variance: ', num2str(round(explained(1),2)), '%'])       
    case 'corr'
        projection = correlationImage(stack, opt.filter_size);
    case 'none'
        projection = squeeze(stack(:,:,:,1));
end%switch projection method

% Reshape the stack to a 2D matrix for further processing
reshaped_stack = reshape(stack, [x*y, t]);

end%FCN:initialPreprocessing

% -------------------------------------------------------------------------

function granules_labeled = initialSegmentation(stack, reshaped_stack, projection, opt)
switch opt.init_segment_method
    case 'voronoi'
        % Get local maxima and minima
        [regional_maxima, regional_minima] = identifyLocalExtrema(projection, opt);
        % 2-D Delaunay triangulation
        granules_labeled = triang_extrema(regional_maxima, regional_minima, opt);
    case 'corr'
        granules_labeled = growing_regions(stack, projection, reshaped_stack, opt);
    case 'rICA'
        granules_labeled = rICA_segmentation(projection, reshaped_stack, opt);
end%switch opt.init_segment_method
end%FCN:initialSegmentation

% -------------------------------------------------------------------------

function granules_labeled = triang_extrema(regional_maxima, regional_minima, opt)
% Erode the extrema location to keep the controids only
ultimateErosion = bwulterode((double(regional_maxima) + double(regional_minima)), 26);
% Reduce to single voxels
% --- First, label all centroids
L = bwlabeln(ultimateErosion, 26);
%  --- Then, get their properties
stats = regionprops3(L, 'Volume', 'Centroid', "VoxelIdxList", "VoxelList");
% --- Iterate over all that are comprised of more than one voxel
V = stats.Volume;
large_Reg = find(V>1);
% Iterate over all regions
for iVox = 2057:length(large_Reg)
    % Get the voxels and centroid
    ind = stats.VoxelIdxList{large_Reg(iVox)};
    C = round(stats.Centroid(large_Reg(iVox),:));
    % Kick out all pixels keeping only the centroid
    ultimateErosion(ind) = 0;
    ultimateErosion(C(2),C(1),C(3)) = 1; % Note to switch x and y
end%iVox
clear L stats V large_Reg ind C 
% Include points that are on the edge of the image. Later, we will set all
% regions that are at the edge to zero
ultimateErosion(1,:,:) = 1; ultimateErosion(:,1,:) = 1; ultimateErosion(:,:,1) = 1;
ultimateErosion(end,:,:) = 1; ultimateErosion(:,end,:) = 1; ultimateErosion(:,:,end) = 1;
% Segment image based on these maxima using a 2-D Delaunay
% triangulation
% Get linear indices of local maxima
linearIndices = find(ultimateErosion);
% Convert linear indices to (x, y, z) coordinates
[Cx, Cy, Cz] = ind2sub(size(ultimateErosion), linearIndices);
% Combine coordinates into Nx3 matrix
C = [Cx, Cy, Cz]; clear Cx Cy Cz
% Segment image based on these maxima using a 2-D Delaunay
% triangulation (i.e. assign pixels to the closest extrem point
% Convert linear indices to (x, y, z) coordinates
linearIndices = 1:numel(regional_maxima);
[x, y, z] = ind2sub(size(regional_maxima), linearIndices);
% Combine all and adjust aspect ratio;
allVoxels = [x(:), y(:), z(:)];
clear x y z linearIndices
% Use knnsearch to find the index of the closest point in 'C' for each voxel
granules_labeled = knnsearch(C.*opt.aspect_ratio, allVoxels.*opt.aspect_ratio);
% Reshape the resulting index array to the dimensions of the 3D volume
granules_labeled = reshape(granules_labeled, size(regional_maxima));
% Remove border
ultimateErosion(1,:,:) = 0; ultimateErosion(:,1,:) = 0; ultimateErosion(:,:,1) = 0;
ultimateErosion(end,:,:) = 0; ultimateErosion(:,end,:) = 0; ultimateErosion(:,:,end) = 0;
linearIndices = find(ultimateErosion);
% Create a logical mask
mask = ismember(granules_labeled, granules_labeled(linearIndices));
granules_labeled(~mask) = 0;
end%FCN:triang_extrema

% -------------------------------------------------------------------------

function [regional_maxima, regional_minima] = identifyLocalExtrema(projection, opt)
% Initialize the output variables
regional_maxima = [];
regional_minima = [];
% Identify the local extrema based on the specified method
SE = strel('sphere', opt.filter_size);
switch opt.regmax_method
    case 'raw'
        regional_maxima = imregionalmax(projection);
        regional_minima = imregionalmin(projection);
    case 'filtered'
        projection_max_filtered = imdilate(projection, SE);
        projection_min_filtered = imerode(projection, SE);
        regional_maxima = imregionalmax(projection_max_filtered, 26);
        regional_minima = imregionalmin(projection_min_filtered, 26);
    case 'both'
        projection_max_filtered = imdilate(projection, SE);
        projection_min_filtered = imerode(projection, SE);
        regional_maxima = imregionalmax(projection, 26) | imregionalmax(projection_max_filtered, 26);
        regional_minima = imregionalmin(projection, 26) | imregionalmin(projection_min_filtered, 26);
end%switch opt.regmax_method
end%FCN:identifyLocalExtrema

% -------------------------------------------------------------------------

function projection = correlationImage(stack, opt)
% --- Preallocation
corr_img = zeros(size(stack,1), size(stack,2), size(stack,3));
projection = corr_img;
% --- Get number of pixels
px_cnt = numel(corr_img);
% --- Get a mask to identify neighbors
SE = strel('sphere', opt.filter_size);
% Iterate over all pixel
parfor iPx = 1:px_cnt
    % Create mask  including all neighbors
    neighbors = zeros(size(projection));
    neighbors(iPx) = 1;
    neighbors = imdilate(neighbors, SE) - neighbors;
    % Get neighbors' indices
    [neighbors_ind_X, neighbors_ind_Y, neighbors_ind_Z] = ind2sub(size(neighbors), find(neighbors));
    [px_ind_X, px_ind_Y, px_ind_Z] = ind2sub(size(neighbors), iPx);
    % Get the correlation
    neighbors_TC = squeeze(mean(mean(mean(stack(neighbors_ind_X, neighbors_ind_Y, neighbors_ind_Z,:)))));
    corr_img(iPx) = corr(neighbors_TC, squeeze(stack(px_ind_X, px_ind_Y, px_ind_Z, :)));
end%iPx
projection = corr_img;
end%FCN:correlationImage

% -------------------------------------------------------------------------

function granules_labeled = growing_regions(stack, corr_img, reshaped_stack, opt)
% Check whether there is an upper limit to the region size
if isnumeric(opt.minPixCount)
    max_size = opt.minPixCount(2);
else
    max_size = 1;
end
% Now, start with the best local correlation and gradually increase the
% region until we reach a threshold of a correlation that got too worse
% --- Preallocation
granules_labeled = zeros(size(corr_img));
% --- Counter for regions
cnt_id = 1;
% Iterate over all pixels
hWait = waitbar(0, ['Growing regions. Please wait ... (0/', num2str(px_cnt),')']);
for iPx = 1:px_cnt
    waitbar(iPx/px_cnt, hWait, ['Growing regions. Please wait ... (', num2str(iPx), '/', num2str(px_cnt), ')']);
    % Get the current correlation value
    [~, max_id] = max(corr_img(:));
    if corr_img(max_id) >= sqrt(0.85)
        % This will be the origin of a new region
        granules_labeled(max_id) = cnt_id;
        while true
            % Get the current region's TC
            curr_TC = mean(reshaped_stack(find(granules_labeled==cnt_id),:),1);
            % Create mask  including all neighbors
            neighbors = granules_labeled == cnt_id;
            neighbors = imdilate(neighbors, SE) - neighbors;
            % Get neighbors' indices
            [neighbors_ind_X, neighbors_ind_Y, neighbors_ind_Z] = ind2sub(size(neighbors), find(neighbors));
            % Keep track of whether pixels got added
            added_px = zeros(1, length(neighbors_ind_X));
            % Check neighbors
            for iN = 1:length(neighbors_ind_X)
                if (granules_labeled(neighbors_ind_X(iN), neighbors_ind_Y(iN), neighbors_ind_Z(iN)) == 0) && (corr(curr_TC(:), squeeze(stack(neighbors_ind_X(iN), neighbors_ind_Y(iN), neighbors_ind_Z(iN),:))) >= sqrt(0.85))
                    granules_labeled(neighbors_ind_X(iN), neighbors_ind_Y(iN), neighbors_ind_Z(iN)) = cnt_id;
                    added_px(iN) = 1;
                end%if valid addition
            end%iN
            % Check break criterion
            if sum(added_px) == 0 || sum(granules_labeled(:)==cnt_id) >= max_size
                break
            end% if no more
        end%while
        % Dont check the current pixels anymore
        corr_img(granules_labeled == cnt_id) = -inf;
        cnt_id = cnt_id+1;
    end%if good enough
end%iPx
close(hWait)
end%FCN:growing_regions

% -------------------------------------------------------------------------

function granules_labeled = rICA_segmentation(projection, reshaped_stack, opt)
% Employ a reconstruction independent component analysis (rICA) to extract
% features that can be turned into individual regions. For this, determine
% the component with the strongest effec on a given pixel.
% Get size of the dataset
[x, y, z] = size(projection);
% Subtract mean since ICA cannot separate sources with a mean signal
% effect.
reshaped_stack = reshaped_stack-nanmean(reshaped_stack,1);
% Prewhiten the signal so that so that it has zero mean and identity
% covariance.
reshaped_stack = prewhiten(reshaped_stack);
% Apply PCA
[~,reshaped_stack] = pca(reshaped_stack);
% Check how many components to request
if opt.n_rICA==0
    opt.n_rICA = size(reshaped_stack, 2);
end
% Apply reconstruction ICA (rICA)
rng(1)% For reproducibility
ricaTransform = rica(reshaped_stack, opt.n_rICA, 'VerbosityLevel', 1);
% Transform the data using the fitted model
transformedData = transform(ricaTransform, reshaped_stack);
% Correct sign
TD = transformedData;
TDv = TD ./ repmat(std(TD), size(TD,1), 1);
TDzp = TDv - 2;
TDzp(TDzp < 0) = 0;
TDzpm = mean(TDzp);
TDzn = TDv + 2;
TDzn(TDzn > 0) = 0;
TDznm = mean(TDzn);
transformedData = TD .* repmat(sign(TDzpm + TDznm), size(TD, 1), 1);
% Reshape the independent components back into a 3-D format
reshaped_stack = reshape(transformedData, [x, y, z, size(transformedData,2)]);
% Smooth components
if strcmp(opt.regmax_method, 'filtered')
    reshaped_stack = imboxfilt3(reshaped_stack, [opt.filter_size, opt.filter_size, 1]);
end%if filter
% Get the component with the strongest effect on a given pixel
[~, granules_labeled] = max(reshaped_stack, [], 4);
% Check for splits
granules_labeled = checkSplits(granules_labeled);
end%FCN:rICA_segmentation

% -------------------------------------------------------------------------

function granules_labeled = refineSegmentation(reshaped_stack, granules_labeled, opt)
% Get mask for neighbors
nMask = strel('sphere', 1);
hWait = waitbar(0, ['Refinig regions. Please wait ... (0/', num2str(opt.n_rep),')']);
for iRep = 1:opt.n_rep
    waitbar(iRep/opt.n_rep, hWait, ['Refinig regions. Please wait ... (', num2str(iRep), '/', num2str(opt.n_rep), ')']);
    % Keep track of changes
    previous_granules_labeled = granules_labeled;
    % Get avg signal trace per granule
    % Get grouping factors
    groupID = granules_labeled(:);
    [~, ~, groupIndices] = unique(groupID);
    nGroups = max(groupIndices);
    % Preallocation
    granuleTC = zeros(nGroups, size(reshaped_stack, 2));
    % Applying accumarray for each column
    for iCol = 1:size(reshaped_stack, 2)
        granuleTC(:, iCol) = accumarray(groupIndices, reshaped_stack(:, iCol), [], @mean);
    end
    clear iCol groupID groupIndices nGroups
    % Get pixels in border regions
    candidates = imgradient3(granules_labeled, 'intermediate')>0;
    idx_candidates = find(candidates);
    % Iterate over all pixels in border regions and check
    % whether they have to be re-assigned.
    switch opt.refinement_method
        case 'corr'
            idx_candidates_newIdentiy = refinement_parfor_corr(idx_candidates, granules_labeled, reshaped_stack, granuleTC, granuleList, nMask);
        case 'rmse'
            idx_candidates_newIdentiy = refinement_parfor_rmse(idx_candidates, granules_labeled, reshaped_stack, granuleTC, granuleList, nMask);
    end
    % Update changed identities
    granules_labeled(idx_candidates) = idx_candidates_newIdentiy;
    % Account for splitting of granules. For this, iterate over
    % all granule and check whether the corresponding itentity
    % can be found in more then one coherent regions
    granules_labeled = checkSplits(granules_labeled);
    % Also account for regions that are too small
    granules_labeled = removeTinyRegions(granules_labeled, opt.minPixCount);
    % Stop if no changes
    if sum(previous_granules_labeled(:) == granules_labeled(:)) == numel(granules_labeled(:))
        break
    end% if change
end%iRep
% Make sure the highest ID is equal to the number of granules
[~, ~, newIDs] = unique(granules_labeled);
granules_labeled = reshape(newIDs, size(granules_labeled));
% Close waitbar
close(hWait)
end%FCN:refineSegmentation

% -------------------------------------------------------------------------

function idx_candidates_newIdentiy = refinement_parfor_corr(idx_candidates, granules_labeled, reshaped_stack, granuleTC, granuleList, nMask)
% Get current labels and prepare output of new labels
idx_candidates_Identiy = granules_labeled(idx_candidates);
idx_candidates_newIdentiy = idx_candidates_Identiy;
% Use a parallel for loop to speed things up
parfor iCan = 1:length(idx_candidates)
    % Create mask  including all neighbors
    neighbors = zeros(size(granules_labeled));
    neighbors(idx_candidates(iCan)) = 1;
    neighbors = imdilate(neighbors, nMask) - neighbors;
    % Get neighbors' indices
    [neighbors_ind_X, neighbors_ind_Y, neighbors_ind_Z] = ind2sub(size(neighbors), find(neighbors));
    % Get neighborhood clusters
    nC = unique(granules_labeled(neighbors_ind_X, neighbors_ind_Y, neighbors_ind_Z));
    % Get the pixel's TC
    activity_px = reshaped_stack(idx_candidates(iCan), :);
    % Check which neighborhood is a better fit
    voronoi_R = zeros(length(nC),1);
    for iN = 1:length(nC)
        activity_cluster = granuleTC(granuleList == nC(iN), :);
        r = corrcoef(activity_px, activity_cluster); voronoi_R(iN) = 1/r(2);
    end%iN
    % Assign new identity
    [~, iBest] = min(voronoi_R);
    idx_candidates_newIdentiy(iCan) = nC(iBest);
end%iCan
end%FCN:refinement_parfor_corr

% -------------------------------------------------------------------------

function idx_candidates_newIdentiy = refinement_parfor_rmse(idx_candidates, granules_labeled, reshaped_stack, granuleTC, granuleList, nMask)
% Get current labels and prepare output of new labels
idx_candidates_Identiy = granules_labeled(idx_candidates);
idx_candidates_newIdentiy = idx_candidates_Identiy;
% Use a parallel for loop to speed things up
parfor iCan = 1:length(idx_candidates)
    % Create mask  including all neighbors
    neighbors = zeros(size(granules_labeled));
    neighbors(idx_candidates(iCan)) = 1;
    neighbors = imdilate(neighbors, nMask) - neighbors;
    % Get neighbors' indices
    [neighbors_ind_X, neighbors_ind_Y, neighbors_ind_Z] = ind2sub(size(neighbors), find(neighbors));
    % Get neighborhood clusters
    nC = unique(granules_labeled(neighbors_ind_X, neighbors_ind_Y, neighbors_ind_Z));
    % Get the pixel's TC
    activity_px = reshaped_stack(idx_candidates(iCan), :);
    % Check which neighborhood is a better fit
    voronoi_R = zeros(length(nC),1);
    for iN = 1:length(nC)
        activity_cluster = granuleTC(granuleList == nC(iN), :);
        voronoi_R(iN) = sqrt(median((activity_cluster-activity_px).^2));
    end%iN
    % Assign new identity
    [~, iBest] = min(voronoi_R);
    idx_candidates_newIdentiy(iCan) = nC(iBest);
end%iCan
end%Fcn:refinement_parfor_rmse

% -------------------------------------------------------------------------

function granules_labeled = checkSplits(granules_labeled)
% Account for splitting of granules. For this, iterate over
% all granule and check whether the corresponding itentity
% can be found in more then one coherent regions
% Identifying Splitter Candidates
unique_granules = unique(granules_labeled);
unique_granules(unique_granules==0) = [];
splitter_candidates = [];
for iP = 1:length(unique_granules)
    % Get region mask
    region_id = unique_granules(iP);
    region_mask = (granules_labeled == region_id);
    % Label connected components within the region
    [~, num] = bwlabeln(region_mask, 26);
    % If more than one connected component, then it's a splitter
    if num > 1
        splitter_candidates = [splitter_candidates, region_id];
    end%if splitter
end%iP
% Get a mask to identify neighbors
SE = strel('sphere', 1);
for iGranule = 1:length(splitter_candidates)
    % Neglect identity of zero
    % Get b/w image and identify number of regions
    bw = granules_labeled==splitter_candidates(iGranule);
    L = bwlabeln(bw, 26);
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
                % Create mask  including all neighbors
                neighbors = L == unique_L(iSplitter);
                neighbors = imdilate(neighbors, SE) - neighbors;
                % Get IDs on the edge
                neighbor_IDs = unique(granules_labeled(find(neighbors)));
                granules_labeled(idx) = mode(neighbor_IDs(:));
            else
                granules_labeled(idx) = max(unique_granules)+1;
                unique_granules = [unique_granules(:); max(unique_granules)+1];
            end
        end%iSplitter
    end% if more
end%iGranule
end%FCN:checkSplits

% -------------------------------------------------------------------------

function granules_labeled = removeTinyRegions(granules_labeled, opt)
% Kick out granules that are too small. For this, iterate
% over all granules, and check whether it is large enough. If
% not, assign each pixel to the best neigbour.
% Repeat as long as all small granules are gone.
% --- Get a mask to identify neighbors
SE = strel('sphere', 1);
small_granules = inf;
while ~isempty(small_granules)
    % Get pixel count for each region
    granuleList = unique(granules_labeled);
    [~,~,C] = unique(granules_labeled);
    C = histcounts(C, 1:length(granuleList)+1);
    C = [C(:), granuleList(:)];
    % Sort base don count
    C = sortrows(C,1,"ascend");
    % Iterate over all regions and take care of those that are too small
    % Note: only take the first and the check again. Otherwise we may get
    % stuck by combining all pixels!
    small_granules = C(C(:, 1) < opt.minPixCount, 2);
    % Get the ID of the current granule
    curr_ID = C(1,2);
    % Get a logical image for the current granule
    bw = granules_labeled==curr_ID;
    % If it is a single pixel, use the mode. Otherwise the closest
    % value that is not the current ID
    [x, y, z] = ind2sub(size(bw), find(bw));
    % Create mask  including all neighbors
    neighbors = imdilate(bw, SE) - bw;
    % Get IDs on the edge
    neighbor_IDs = unique(granules_labeled(find(neighbors)));
    if C(1,1) == 1
        granules_labeled(x, y, z) = mode(neighbor_IDs(:));
    else
        % Iterate over all pixels of the region that is too small and
        % assign them to the closest neighbor region
        for iPx = 1:length(x)
            id_table = nan(1, length(neighbor_IDs));
            for iN = 1:length(neighbor_IDs)
                [Nx, Ny, Nz] = ind2sub(size(bw), find(granules_labeled==neighbor_IDs(iN)));
                dist = [Nx(:), Ny(:), Nz(:)] - [x(iPx), y(iPx), z(iPx)];
                min_dist = min(sqrt(sum(dist'.*dist'))');
                id_table(iN) = min_dist;
            end%iN
            [~,best_id] = min(id_table);
            granules_labeled(x(iPx), y(iPx), z(iPx)) = neighbor_IDs(best_id);
        end%iPx
    end% if single pixel
end%while
end%FCN:removeTinyRegions

% -------------------------------------------------------------------------

function [granules_labeled, summary_stats] = finalRefinementAndStats(reshaped_stack, granules_labeled, opt)
% Pool pixels belonging to the same granule
% --- Get list of granules
granuleList = unique(granules_labeled);
% --- Preallocation
summary_stats.avgTCs = nan(length(granuleList), size(reshaped_stack,2));
summary_stats.granule_Corr = nan(length(granuleList), 1);
summary_stats.granule_Avg_img = nan(size(granules_labeled));
summary_stats.granule_Min_img = nan(size(granules_labeled));
summary_stats.granule_Max_img = nan(size(granules_labeled));
summary_stats.granule_Std_img = nan(size(granules_labeled));
summary_stats.granule_Max_img = nan(size(granules_labeled));
summary_stats.granule_Corr_img = nan(size(granules_labeled));
% Iterate over all granules
for iGranule = 1:length(granuleList)
    % Get their index positions
    idx_granule = find(granules_labeled==granuleList(iGranule));
    % Get the avg time course of the current granule
    summary_stats.avgTCs(iGranule,:) = nanmean(reshaped_stack(idx_granule,:),1);
    % Get within-granule correlation
    summary_stats.granule_Corr(iGranule,1) = calculateAverageCorrelation(summary_stats.avgTCs(iGranule,:), reshaped_stack(idx_granule,:));
    % Project different metrics like the avg or std into 2D
    % --- mean
    summary_stats.granule_Avg_img(idx_granule) = nanmean(summary_stats.avgTCs(iGranule,:));
    % --- min
    summary_stats.granule_Min_img(idx_granule) = nanmin(summary_stats.avgTCs(iGranule,:));
    % --- max
    summary_stats.granule_Max_img(idx_granule) = nanmax(summary_stats.avgTCs(iGranule,:));
    % --- std
    summary_stats.granule_Std_img(idx_granule) = nanstd(summary_stats.avgTCs(iGranule,:));
    % --- corr
    summary_stats.granule_Corr_img(idx_granule) = summary_stats.granule_Corr(iGranule,1);
end%iGranule
% Based on the projection method, estimate which regions are active
switch opt.projection_method
    case 'std'
        summary_stats.active_region.map = imbinarize(summary_stats.granule_Std_img);
        summary_stats.active_region.method = 'std';
    case 'mean'
        summary_stats.active_region.map = imbinarize(summary_stats.granule_Avg_img);
        summary_stats.active_region.method = 'mean';
    case 'max'
        summary_stats.active_region.map = imbinarize(summary_stats.granule_Max_img);
        summary_stats.active_region.method = 'max';
    case 'none'
        summary_stats.active_region.map = [];
        summary_stats.active_region.method = 'none';
end%switch projection method
end%FCN:finalRefinementAndStats

% -------------------------------------------------------------------------

function overallAvg = calculateAverageCorrelation(avg_tc, ind_tc)
% Calculate the average correlation with the centroid
if size(avg_tc, 2)>1
    % --- Preallocation
    corr_values = nan(1, size(ind_tc,1));
    for iPx = 1:size(ind_tc,1)
        r = corrcoef(ind_tc(iPx,:), avg_tc);
        corr_values(iPx) = r(2);
    end%iPx
    overallAvg = nanmean(corr_values);
else
    overallAvg = NaN;
end%if more dimensions
end%FCN:calculateAverageCorrelation