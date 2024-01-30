function [granules_labeled, summary_stats] = CalciSeg_3D(stack, aspect_ratio, projection_method, init_segment_method, regmax_method, n_rep, refinement_method, minPixCount)
% CalciSeg_3D(stack, projection_method, init_segment_method, regmax_method, n_rep, refinement_method, minPixCount)
% segments the spatial aspects of a stack (x*y*z*time) into individual
% regions ('granules') based on a refined 3D Delaunay triangulation.
% Note: CalciSeg uses parallel for-loops at for the local growing
% algorithm and the refinement process. It is recommended to start a
% parallel pool of workers before calling the function.
%
% Input:
%   stack             : 4-D matrix (x*y*z*time).
%   aspect_ratio      : x-y-z aspect ratio (e.g., [1 1 1])
%   projection_method : method for calculating the projection across time
%                       - 'std'    - standard deviation projection
%                       - 'mean'   - mean intensity projection
%                       - 'median' - median intensity projection
%                       - 'max'    - maximum intensity projection
%                       - 'min'    - minimum intensity projection
%                       - 'pca'    - principal component projection
%                       - 'corr'   - local correlation between neighboring
%                                    pixels
%
%   init_seg_method   : method for initial segmentation
%                       - 'voronoi'  - Delaunay triangulation
%                       - 'corr'     - local growing based on correlation (r_threshold = sqrt(0.85))
%
%   regmax_method     : method for determining how to identify local extrema
%                       - 'raw'      - simply apply imregionalmax/-min on
%                                      the projection
%                       - 'filtered' - dilate/erode image to apply moving
%                                      max/min filter before applying
%                                      imregionalmax/-min
%                       - 'both'     - combine both above-mentioned methods
%
%   n_rep             : Number of iterations for refining the regions.
%   refinement_method : Measure to fine-tune granule assignment
%                       - 'corr'     - correlation
%                       - 'rmse'     - root median square error
%   minPixCount       : Minimum pixel area. Use 'auto' for an
%                       automatic assessment based on the distribution of
%                       region sizes before refined segmentation (min=5th
%                       quantile). Or provide a single number. Note
%                       that, if not set to 'auto', the min region size
%                       affects the filter size for regmax_method. If set
%                       to auto, the size is set to 1.
%
%
% Output:
%   granules_labeled : granule labels
%   summary_stats    : avg, std, corr for each granule
%                      .projection           : projection image used for initial segmentation
%                      .avgTCs               : each granule's mean time course
%                      .granule_Corr         : avg within-granule correlation
%                      .granule_Avg_img      : within-granule average activity
%                      .granule_Std_img      : each granule's std over time
%                      .granule_Max_img      : each granule's max over time
%                      .granule_Corr_img     : within-granule correlation
%                      .active_region.map    : binary map of active regions
%                      .active_region.method : method used for binarizing
%
% Version: 17-Jan-24 (R2023a)

% Validate inputs
stack = squeeze(stack);
stack_size = whos('stack');
if ~strcmp(stack_size.class, 'single')
    stack = single(stack);
end%if not single precision
error_message = '';
[valid, error_message] = validateInputs(error_message, stack, aspect_ratio, projection_method, init_segment_method, regmax_method, n_rep, refinement_method, minPixCount);
if valid == -1
    error(error_message);
elseif valid == 0
    warning(error_message);
end

% Define size for local neighborhood
if isnumeric(minPixCount)
    filter_size = ceil(sqrt(minPixCount(1)/pi)) + double(minPixCount(1)==0);
elseif strcmp(minPixCount, 'auto')
    filter_size = 1;
end

% Initial preprocessing
[reshaped_stack, projection, refinement_method, init_segment_method] = initialPreprocessing(stack, projection_method, refinement_method, init_segment_method, filter_size);

% Perform initial segmentation
granules_labeled = initialSegmentation(stack, aspect_ratio, reshaped_stack, projection, regmax_method, init_segment_method, filter_size, minPixCount);

% Maybe, the user wants to have the min/max size estimated based on the data
if (isstring(minPixCount) || ischar(minPixCount)) && strcmp(minPixCount, 'auto')
    % Get pixel count for each region
    granuleList = unique(granules_labeled);
    [~,~,C] = unique(granules_labeled);
    C = histcounts(C, 1:length(granuleList)+1);
    C = [C(:), granuleList(:)];
    minPixCount = floor(quantile(C(2:end,1), 0.05));
    disp(['automatic assessment of the min. granule size: ', num2str(minPixCount)])
end%if auto min size

% Refine segmentation
if n_rep>0
    granules_labeled = refineSegmentation(reshaped_stack, granules_labeled, n_rep, refinement_method, minPixCount);
end%if refine

% Final refinement steps and statistics calculation
if nargout>1
    [granules_labeled, summary_stats] = finalRefinementAndStats(reshaped_stack, granules_labeled, projection_method);
    summary_stats.projection = projection;
end%if summary stats

end%FCN:CalciSeg_3D

% -------------------------------------------------------------------------

function [valid, error_message] = validateInputs(error_message, stack, aspect_ratio, projection_method, init_segment_method, regmax_method, n_rep, refinement_method, minPixCount)
% Initialize output variables
valid = 1;

% Validate 'stack' input
if ndims(stack) ~= 4
    valid = 0;
    error_message = 'Input "stack" has no temporal component. Inputs are adjusted accordingly (projection_method="none"; init_segment_method="voronoi"; refinement_method="rmse"). ';
end%if invalid 'stack' input

% Validate 'aspect_ratio' input
if ~isnumeric(aspect_ratio) || length(aspect_ratio) ~= 3 || isempty(aspect_ratio) || any(aspect_ratio<=0)
    valid = -1;
    error_message = [error_message, 'Input "aspect_ratio" has to be a vector of positive values with a length of 3.'];
end%if invalid 'stack' input

% Validate 'projection_method' input
valid_projection_methods = {'std', 'mean', 'median', 'max', 'min', 'pca', 'corr', 'none'};
if ~ismember(projection_method, valid_projection_methods)
    valid = -1;
    error_message = [error_message, sprintf('Invalid "projection_method". Valid options are: %s.', strjoin(valid_projection_methods, ', '))];
    return;
end%if invalid 'projection_method' input

% Validate 'projection_method' input
valid_init_segment_methods = {'voronoi', 'corr'};
if ~ismember(init_segment_method, valid_init_segment_methods)
    valid = -1;
    error_message = [error_message, sprintf('Invalid "init_segment_method". Valid options are: %s.', strjoin(valid_init_segment_methods, ', '))];
    return;
end%if invalid 'projection_method' input

% Validate 'regmax_method' input
valid_regmax_methods = {'raw', 'filtered', 'both'};
if ~ismember(regmax_method, valid_regmax_methods)
    valid = -1;
    error_message = [error_message, sprintf('Invalid "regmax_method". Valid options are: %s.', strjoin(valid_regmax_methods, ', '))];
    return;
end%if invalid 'regmax_method' input

% Validate 'n_rep' input
if ~isnumeric(n_rep) || n_rep < 0
    valid = -1;
    error_message = [error_message, 'Input "n_rep" must be a non-negative integer value.'];
    return;
end%if invalid 'n_rep' input

% Validate 'refinement_method' input
valid_refinement_methods = {'corr', 'rmse'};
if ~ismember(refinement_method, valid_refinement_methods)
    valid = -1;
    error_message = [error_message, sprintf('Invalid "refinement_method". Valid options are: %s.', strjoin(valid_refinement_methods, ', '))];
    return;
end%if invalid 'refinement_method' input

% Validate 'minPixCount' input
if ~isnumeric(minPixCount) && ~strcmp(minPixCount, 'auto')
    valid = -1;
    error_message = [error_message, 'Input "minPixCount" must be a non-negative integer value or "auto".'];
    return;
elseif isnumeric(minPixCount) && length(minPixCount) ~= 1
    valid = -1;
    error_message = [error_message, 'Input "minPixCount" must be a non-negative integer value or "auto".'];
    return;
end%if invalid 'minPixCount' input
end%FCN:validateInputs

% -------------------------------------------------------------------------

function [reshaped_stack, projection, refinement_method, init_segment_method] = initialPreprocessing(stack, projection_method, refinement_method, init_segment_method, filter_size)
% Get the size of the stack
[x, y, z, t] = size(stack);

% Check the time dimension and adjust parameters accordingly
if t == 1
    projection_method = 'none';
    refinement_method = 'rmse';
    init_segment_method = 'voronoi';
end%if no temporal component

% Reshape the stack to a 2D matrix for further processing
reshaped_stack = reshape(stack, [x*y*z, t]);

% Initialize the projection variable
projection = [];

% Calculate the projection based on the specified method
switch projection_method
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
        [~, projection, ~, ~, explained] = pca(reshaped_stack, 'NumComponents',1);
        projection = reshape(projection, [x, y, z]);
        disp(['Explained variance: ', num2str(round(explained(1),2)), '%'])
        clear explained
    case 'corr'
        projection = correlationImage(stack, filter_size);
    case 'none'
        projection = stack;
end%switch projection method

end%FCN:initialPreprocessing

% -------------------------------------------------------------------------

function granules_labeled = initialSegmentation(stack, aspect_ratio, reshaped_stack, projection, regmax_method, init_segment_method, filter_size, minPixCount)
switch init_segment_method
    case 'voronoi'
        % Get local maxima and minima
        [regional_maxima, regional_minima] = identifyLocalExtrema(projection, regmax_method, filter_size);
        % 2-D Delaunay triangulation
        granules_labeled = triang_extrema(regional_maxima, regional_minima, aspect_ratio);
    case 'corr'
        granules_labeled = growing_regions(stack, projection, reshaped_stack, minPixCount);
end%switch init_segment_method
end%FCN:initialSegmentation

% -------------------------------------------------------------------------

function granules_labeled = triang_extrema(regional_maxima, regional_minima, aspect_ratio)
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
granules_labeled = knnsearch(C.*aspect_ratio, allVoxels.*aspect_ratio);
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

function [regional_maxima, regional_minima] = identifyLocalExtrema(projection, regmax_method, filter_size)
% Initialize the output variables
regional_maxima = [];
regional_minima = [];
% Identify the local extrema based on the specified method
SE = strel('sphere', filter_size);
switch regmax_method
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
end%switch regmax_method
end%FCN:identifyLocalExtrema

% -------------------------------------------------------------------------

function projection = correlationImage(stack, filter_size)
% --- Preallocation
corr_img = zeros(size(stack,1), size(stack,2), size(stack,3));
projection = corr_img;
% --- Get number of pixels
px_cnt = numel(corr_img);
% --- Get a mask to identify neighbors
SE = strel('sphere', filter_size);
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

function granules_labeled = growing_regions(stack, corr_img, reshaped_stack, minPixCount)
% Check whether there is an upper limit to the region size
if isnumeric(minPixCount)
    max_size = minPixCount(2);
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

function granules_labeled = refineSegmentation(reshaped_stack, granules_labeled, n_rep, refinement_method, minPixCount)
% Account for regions that are too small
granules_labeled = removeTinyRegions(granules_labeled, minPixCount);
% Get mask for neighbors
nMask = strel('sphere', 1);
hWait = waitbar(0, ['Refinig regions. Please wait ... (0/', num2str(n_rep),')']);
for iRep = 1:n_rep
    waitbar(iRep/n_rep, hWait, ['Refinig regions. Please wait ... (', num2str(iRep), '/', num2str(n_rep), ')']);
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
    switch refinement_method
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
    granules_labeled = removeTinyRegions(granules_labeled, minPixCount);
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

function granules_labeled = removeTinyRegions(granules_labeled, minPixCount)
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
    small_granules = C(C(:, 1) < minPixCount, 2);
    hWait = waitbar(0, ['Removing tiny regions ... 0/', num2str(length(small_granules))]);
    for iGranule = 1:length(small_granules)
        waitbar(iGranule/length(small_granules), hWait, ['Removing tiny regions ... ', num2str(iGranule), '/', num2str(length(small_granules))])
        % Get the ID of the current granule
        curr_ID = C(iGranule,2);
        % Get a logical image for the current granule
        bw = granules_labeled==curr_ID;
        % If it is a single pixel, use the mode. Otherwise the closest
        % value that is not the current ID
        [x, y, z] = ind2sub(size(bw), find(bw));
        % Create mask  including all neighbors
        neighbors = imdilate(bw, SE) - bw;
        % Get IDs on the edge
        neighbor_IDs = unique(granules_labeled(find(neighbors)));
        if C(iGranule,1) == 1
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
    end%iGranule
    close(hWait)
end%while
end%FCN:removeTinyRegions

% -------------------------------------------------------------------------

function [granules_labeled, summary_stats] = finalRefinementAndStats(reshaped_stack, granules_labeled, projection_method)
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
switch projection_method
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