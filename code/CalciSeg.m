function [granules_labeled, summary_stats] = CalciSeg(stack, projection_method, init_segment_method, regmax_method, n_rep, refinement_method, limitPixCount)
% CalciSeg(stack, projection_method, init_segment_method, regmax_method, n_rep, refinement_method, limitPixCount)
% segments the spatial aspects of a stack (x*y*time) into individual
% regions ('granules') based on a refined 2D Delaunay triangulation.
% Note: CalciSeg uses parallel for-loops at for the local growing
% algorithm and the refinement process. It is recommended to start a
% parallel pool of workers before calling the function.
%
% Input:
%   stack             : 3-D matrix (x*y*time).
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
%   limitPixCount     : Minimum and maximum pixel area. Use 'auto' for an
%                       automatic assessment based on the distribution of
%                       region sizes before refined segmentation (min=5th
%                       quantile; max=95thquantile). Or provide a
%                       two-element vectorlimitPixCount=[min, max]. Note
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
[valid, error_message] = validateInputs(error_message, stack, projection_method, init_segment_method, regmax_method, n_rep, refinement_method, limitPixCount);
if valid == -1
    error(error_message);
elseif valid == 0
    warning(error_message);
end

% Define size for local neighborhood
if isnumeric(limitPixCount)
    filter_size = ceil(sqrt(limitPixCount(1)/pi)) + double(limitPixCount(1)==0);
elseif strcmp(limitPixCount, 'auto')
    filter_size = 1;
end

% Initial preprocessing
[reshaped_stack, projection, refinement_method, init_segment_method] = initialPreprocessing(stack, projection_method, refinement_method, init_segment_method, filter_size);

% Perform initial segmentation
granules_labeled = initialSegmentation(stack, reshaped_stack, projection, regmax_method, init_segment_method, filter_size, limitPixCount);

% Maybe, the user wants to have the min/max size estimated based on the data
if (isstring(limitPixCount) || ischar(limitPixCount)) && strcmp(limitPixCount, 'auto')
    % Get pixel count for each region
    granuleList = unique(granules_labeled);
    [~,~,C] = unique(granules_labeled);
    C = histcounts(C, 1:length(granuleList)+1);
    C = [C(:), granuleList(:)];
    limitPixCount = [floor(quantile(C(2:end,1), 0.05)), floor(quantile(C(2:end,1), 0.95))];
    disp(['automatic assessment of the min./max. granule size: [', num2str(limitPixCount),']'])
end%if auto min size

% Refine segmentation
if n_rep>0
    granules_labeled = refineSegmentation(reshaped_stack, granules_labeled, n_rep, refinement_method, limitPixCount);
end%if refine

% Final refinement steps and statistics calculation
if nargout>1
    [granules_labeled, summary_stats] = finalRefinementAndStats(reshaped_stack, granules_labeled, projection_method);
    summary_stats.projection = projection;
end%if summary stats

end%FCN:CalciSeg

% -------------------------------------------------------------------------

function [valid, error_message] = validateInputs(error_message, stack, projection_method, init_segment_method, regmax_method, n_rep, refinement_method, limitPixCount)
% Initialize output variables
valid = 1;

% Validate 'stack' input
if ndims(stack) ~= 3
    valid = 0;
    error_message = 'Input "stack" has no temporal component. Inputs are adjusted accordingly (projection_method="none"; init_segment_method="voronoi"; refinement_method="rmse"). ';
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

% Validate 'limitPixCount' input
if ~isnumeric(limitPixCount) && ~strcmp(limitPixCount, 'auto')
    valid = -1;
    error_message = [error_message, 'Input "limitPixCount" must be a non-negative pair of two integer values or "auto".'];
    return;
elseif isnumeric(limitPixCount) && length(limitPixCount) ~= 2
    valid = -1;
    error_message = [error_message, 'Input "limitPixCount" must be a non-negative pair of two integer values or "auto".'];
    return;
end%if invalid 'limitPixCount' input
end%FCN:validateInputs

% -------------------------------------------------------------------------

function [reshaped_stack, projection, refinement_method, init_segment_method] = initialPreprocessing(stack, projection_method, refinement_method, init_segment_method, filter_size)
% Get the size of the stack
[x, y, t] = size(stack);

% Check the time dimension and adjust parameters accordingly
if t == 1
    projection_method = 'none';
    refinement_method = 'rmse';
    init_segment_method = 'voronoi';
end%if no temporal component

% Reshape the stack to a 2D matrix for further processing
reshaped_stack = reshape(stack, [x*y, t]);

% Initialize the projection variable
projection = [];

% Calculate the projection based on the specified method
switch projection_method
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
        [~, projection, ~, ~, explained] = pca(reshaped_stack, 'NumComponents',1);
        projection = reshape(projection, [x, y]);
        disp(['Explained variance: ', num2str(round(explained(1),2)), '%'])
        clear explained
    case 'corr'
        projection = correlationImage(stack, filter_size);
    case 'none'
        projection = stack;
end%switch projection method

end%FCN:initialPreprocessing

% -------------------------------------------------------------------------

function granules_labeled = initialSegmentation(stack, reshaped_stack, projection, regmax_method, init_segment_method, filter_size, limitPixCount)
switch init_segment_method
    case 'voronoi'
        % Get local maxima and minima
        [regional_maxima, regional_minima] = identifyLocalExtrema(projection, regmax_method, filter_size);
        % 2-D Delaunay triangulation
        granules_labeled = triang_extrema(regional_maxima, regional_minima);
    case 'corr'
        granules_labeled = growing_regions(stack, projection, reshaped_stack, limitPixCount);
end%switch init_segment_method
end%FCN:initialSegmentation

% -------------------------------------------------------------------------

function granules_labeled = triang_extrema(regional_maxima, regional_minima)
% Erode the extrema location to keep the centroids only
ultimateErosion = bwulterode(regional_maxima>0) + bwulterode(regional_minima>0);
% Include points that are on the edge of the image. Later, we will set all
% regions that are at the edge to zero
ultimateErosion(1,:) = 1; ultimateErosion(:,1) = 1;
ultimateErosion(end,:) = 1; ultimateErosion(:,end) = 1;
% Segment image based on these maxima using a 2-D Delaunay
% triangulation
% Get linear indices of local maxima
linearIndices = find(ultimateErosion);
% Convert linear indices to (x, y, z) coordinates
[Cx, Cy] = ind2sub(size(ultimateErosion), linearIndices);
% Combine coordinates into Nx3 matrix
C = [Cx, Cy]; clear Cx Cy
% Segment image based on these maxima using a 2-D Delaunay
% triangulation (i.e. assign pixels to the closest extreme point
% Convert linear indices to (x, y) coordinates
linearIndices = 1:numel(regional_maxima);
[x, y] = ind2sub(size(regional_maxima), linearIndices);
% Combine all and adjust aspect ratio;
allVoxels = [x(:), y(:)];
clear x y z linearIndices
% Use knnsearch to find the index of the closest point in 'C' for each voxel
granules_labeled = knnsearch(C, allVoxels);
% Reshape the resulting index array to the dimensions of the 3D volume
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

function [regional_maxima, regional_minima] = identifyLocalExtrema(projection, regmax_method, filter_size)
% Initialize the output variables
regional_maxima = [];
regional_minima = [];
SE = strel('disk', filter_size);
% Identify the local extrema based on the specified method
switch regmax_method
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
end%switch regmax_method
end%FCN:identifyLocalExtrema

% -------------------------------------------------------------------------

function projection = correlationImage(stack, filter_size)
% Get a map of local correlatios of a pixel with each of its neighbors
% --- Preallocation
projection = zeros(size(stack,1), size(stack,2));
corr_img = projection;
% --- Get number of pixels
px_cnt = numel(corr_img);
% --- Get a mask to identify neighbors
nMask = strel('disk', filter_size);
nMask = nMask.Neighborhood;
elementIndex = size(nMask,1); elementIndex = (elementIndex / 2) + 0.5;
nMask(elementIndex, elementIndex) = 0;
% Iterate over all pixel
parfor iPx = 1:px_cnt
    % Create mask  including all neighbors
    neighbors = zeros(size(projection));
    neighbors(iPx) = 1;
    neighbors = conv2(neighbors, nMask, 'same')==1;
    neighbors(iPx) = 0;
    % Get neighbors' indices
    [neighbors_ind_X, neighbors_ind_Y] = ind2sub(size(neighbors), find(neighbors));
    [px_ind_X, px_ind_Y] = ind2sub(size(neighbors), iPx);
    % Get the correlation
    neighbors_TC = squeeze(mean(mean(stack(neighbors_ind_X, neighbors_ind_Y,:))));
    corr_img(iPx) = corr(neighbors_TC, squeeze(stack(px_ind_X, px_ind_Y, :)));
end%iPx
projection = corr_img;
end%FCN:correlationImage

% -------------------------------------------------------------------------

function granules_labeled = growing_regions(stack, corr_img, reshaped_stack, limitPixCount)
% Check whether there is an upper limit to the region size
if isnumeric(limitPixCount)
    max_size = limitPixCount(2);
else
    max_size = 1;
end
% Now, start with the best local correlation and gradually increase the
% region until we reach a threshold of a correlation that got too bad
% --- Preallocation
granules_labeled = zeros(size(corr_img));
% --- Counter for regions
cnt_id = 1;
% Iterate over all pixels
px_cnt = numel(corr_img);
hWait = waitbar(0, ['Growing regions ... ', num2str(0),'/',num2str(px_cnt)]);
for iPx = 1:px_cnt
    waitbar(iPx/px_cnt, hWait, ['Growing regions ... ', num2str(iPx),'/',num2str(px_cnt)])
    % Get the current correlation value
    [~, max_id] = max(corr_img(:));
    if corr_img(max_id) >= sqrt(0.85)
        % This will be the origin of a new region
        granules_labeled(max_id) = cnt_id;
        while true
            % Get the current region's TC
            curr_TC = mean(reshaped_stack(find(granules_labeled==cnt_id),:),1);
            % Get the neighbors
            neighbors = imdilate(granules_labeled == cnt_id, strel('disk', 1)) - (granules_labeled == cnt_id);
            [neighbors_ind_X, neighbors_ind_Y] = ind2sub(size(neighbors), find(neighbors));
            % Keep track of whether pixels got added
            added_px = zeros(1, length(neighbors_ind_X));
            % Check neighbors
            for iN = 1:length(neighbors_ind_X)
                if (granules_labeled(neighbors_ind_X(iN), neighbors_ind_Y(iN)) == 0) && (corr(curr_TC(:), squeeze(stack(neighbors_ind_X(iN), neighbors_ind_Y(iN),:))) >= sqrt(0.85))
                    granules_labeled(neighbors_ind_X(iN), neighbors_ind_Y(iN)) = cnt_id;
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

function granules_labeled = refineSegmentation(reshaped_stack, granules_labeled, n_rep, refinement_method, limitPixCount)
% Account for regions that are too small or too large. For this,
% check each region that has less or more than the required number
% of pixels. Note, first run removeHugeRegions, as this may create tiny
% regions.
granules_labeled = removeHugeRegions(granules_labeled, limitPixCount(2));
granules_labeled = removeTinyRegions(granules_labeled, limitPixCount(1));
% Get mask for neighbors
nMask = [0 1 0; 1 0 1; 0 1 0];
% Refine granules
previous_granules_labeled = -ones(size(granules_labeled));
hWait = waitbar(0, ['Refine regions ... 0/',num2str(n_rep)]);
for iRep = 1:n_rep
    waitbar(iRep/n_rep, hWait, ['Refine regions ... ', num2str(iRep),'/',num2str(n_rep)])
    % Get avg signal trace per granule
    granuleList = unique(granules_labeled);
    granuleTC = nan(length(granuleList), size(reshaped_stack,2));
    for iGranule = 1:length(granuleList)
        idx = find(granules_labeled==granuleList(iGranule));
        granuleTC(iGranule, :) = nanmean(reshaped_stack(idx,:),1);
    end%iGranule
    clear iC iR iGranule
    % Get pixels in border regions. But only take regions into account
    % that changed
    ind = find(granules_labeled~=previous_granules_labeled);
    keep_changed_IDs = unique([granules_labeled(ind(:)); previous_granules_labeled(ind(:))]);
    mask = ismember(granules_labeled, keep_changed_IDs);
    candidates = imgradient(granules_labeled, "intermediate")>0;
    candidates(~mask) = 0;
    idx_candidates = find(candidates);
    idx_table = reshape(1:numel(granules_labeled), size(granules_labeled));
    % Iterate over all pixels in border regions and check
    % whether thez have to be re-assigned.
    previous_granules_labeled = granules_labeled;
    switch refinement_method
        case 'corr'
            idx_candidates_newIdentiy = refinement_parfor_corr(idx_candidates, idx_table, granules_labeled, reshaped_stack, granuleTC, granuleList, nMask);
        case 'rmse'
            idx_candidates_newIdentiy = refinement_parfor_rmse(idx_candidates, idx_table, granules_labeled, reshaped_stack, granuleTC, granuleList, nMask);
    end
    % Update changed identities
    granules_labeled(idx_candidates) = idx_candidates_newIdentiy;
    % Account for splitting of granules. For this, iterate over
    % all granulee and check whether the corresponding itentity
    % can be found in more than one coherent region
    granules_labeled = checkSplits(granules_labeled);
    % Account for regions that are too small or too large. For this,
    % check each region that has less or more than the required number
    % of pixels. Note, first run removeHugeRegions, as this may create tiny
    % regions.
    granules_labeled = removeHugeRegions(granules_labeled, limitPixCount(2));
    granules_labeled = removeTinyRegions(granules_labeled, limitPixCount(1));
    % Stop if no changes
    if sum(previous_granules_labeled(:) == granules_labeled(:)) == numel(granules_labeled(:))
        break
    end% if change
end%iRep
close(hWait)
% Make sure the highest ID is equal to the number of granules
[~, ~, newIDs] = unique(granules_labeled);
granules_labeled = reshape(newIDs, size(granules_labeled));
end%FCN:refineSegmentation

% -------------------------------------------------------------------------

function idx_candidates_newIdentiy = refinement_parfor_corr(idx_candidates, idx_table, granules_labeled, reshaped_stack, granuleTC, granuleList, nMask)
% Get current labels and prepare output of new labels
idx_candidates_Identiy = granules_labeled(idx_candidates);
idx_candidates_newIdentiy = idx_candidates_Identiy;
% Use a parallel for loop to speed things up
parfor iCan = 1:length(idx_candidates)
    % Get indices of neighbors
    idx_neighbor = conv2( double(idx_table==idx_candidates(iCan)), nMask, 'same')==1;
    % Get neighborhood clusters
    nC = unique(granules_labeled(idx_neighbor));
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

function idx_candidates_newIdentiy = refinement_parfor_rmse(idx_candidates, idx_table, granules_labeled, reshaped_stack, granuleTC, granuleList, nMask)
% Get current labels and prepare output of new labels
idx_candidates_Identiy = granules_labeled(idx_candidates);
idx_candidates_newIdentiy = idx_candidates_Identiy;
% Use a parallel for loop to speed things up
parfor iCan = 1:length(idx_candidates)
    % Get indices of neighbors
    idx_neighbor = conv2( double(idx_table==idx_candidates(iCan)), nMask, 'same')==1;
    % Get neighborhood clusters
    nC = unique(granules_labeled(idx_neighbor));
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
splitter_candidates = [];
parfor iP = 1:length(unique_granules)
    % Get region mask
    region_id = unique_granules(iP);
    region_mask = (granules_labeled == region_id);
    % Label connected components within the region
    [~, num] = bwlabel(region_mask, 4);
    % If more than one connected component, then it's a splitter
    if num > 1
        splitter_candidates = [splitter_candidates, region_id];
    end%if splitter
end%iP
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
% not, assign each pixel to the best neighbour.
% Repeat as long as all small granules are gone.
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
    for iGranule = 1:length(small_granules)
        % Get the ID of the current granule
        curr_ID = C(iGranule,2);
        % Get a logical image for the current granule
        bw = granules_labeled==curr_ID;
        % If it is a single pixel, use the mode. Otherwise the closest
        % value that is not the current ID
        [x, y] = ind2sub(size(bw), find(bw));
        % Get pixels on the edge
        neighbors = imgradient(bw, 'central')>0;
        neighbors(bw) = 0;
        neighbor_IDs = unique(granules_labeled(neighbors));
        if C(iGranule,1) == 1
            granules_labeled(x, y) = mode(neighbor_IDs(:));
        else
            % Iterate over all pixels of the region that is too small and
            % assign them to the closest neighbor region
            for iPx = 1:length(x)
                id_table = nan(1, length(neighbor_IDs));
                for iN = 1:length(neighbor_IDs)
                    [Nx, Ny] = ind2sub(size(bw), find(granules_labeled==neighbor_IDs(iN)));
                    dist = [Nx(:), Ny(:)] - [x(iPx), y(iPx)];
                    min_dist = min(sqrt(sum(dist'.*dist'))');
                    id_table(iN) = min_dist;
                end%iN
                [~,best_id] = min(id_table);
                granules_labeled(x(iPx), y(iPx)) = neighbor_IDs(best_id);
            end%iPx
        end% if single pixel
    end%iGranule
end%while
end%FCN:removeTinyRegions

% -------------------------------------------------------------------------

function granules_labeled = removeHugeRegions(granules_labeled, maxPixCount)
% Kick out granules that are too large. For this, iterate
% over all granules, and check whether it is small enough. If
% not, split the region by watershedding it
% Repeat as long as all small granules are gone.
large_granules = inf;
while ~isempty(large_granules)
    % Get pixel count for each region
    granuleList = unique(granules_labeled);
    [~,~,C] = unique(granules_labeled);
    C = histcounts(C, 1:length(granuleList)+1);
    C = [C(:), granuleList(:)];
    % Sort base don count
    C = sortrows(C,1,"descend");
    % Iterate over all regions and take care of those that are too large
    large_granules = C(C(:, 1) > maxPixCount, 2);
    for iGranule = 1:length(large_granules)
        % Get the ID of the current granule
        curr_ID = C(iGranule,2);
        % Get a logical image for the current granule
        bw = granules_labeled==curr_ID;
        % Watershed
        D = bw;
        D = bwdist(~D);
        D = -D;
        L = watershed(D);
        L(~bw) = 0;
        bw_cut = L>0;
        bw_cut = bwlabel(bw_cut);
        % Sometimes, this does not work for circular shapes. For
        % this, cut along the minor axis.
        if ~any(bw(bw) ~= bw_cut(bw))
            % Get the minor axis
            stats = regionprops(bw, 'Centroid', 'MajorAxisLength', 'MinorAxisLength', 'Orientation');
            ellipse = stats(1);
            longAxis = ellipse.MajorAxisLength;
            shortAxis = ellipse.MinorAxisLength;
            % Get the position and anle of the ellipse
            centerX = ellipse.Centroid(1);
            centerY = ellipse.Centroid(2);
            angle = ellipse.Orientation;
            % Define the minor axis
            x1 = floor(centerX - (shortAxis/2) * cosd(angle));
            y1 = floor(centerY - (shortAxis/2) * sind(angle));
            x2 = ceil(centerX + (shortAxis/2) * cosd(angle));
            y2 = ceil(centerY + (shortAxis/2) * sind(angle));
            % Draw a line
            numPoints = sqrt(size(bw,1)^2+size(bw,2)^2);
            x = round(linspace(x1, x2, numPoints)); x(x<1)=1; x(x>size(bw_cut,2))=size(bw_cut,2);
            y = round(linspace(y1, y2, numPoints)); y(y<1)=1; y(y>size(bw_cut,1))=size(bw_cut,1);
            for iPx = 1:numPoints
                bw_cut(y(iPx), x(iPx), :) = 0;
            end%iPx
            bw_cut = bwlabel(bw_cut,4);
            % If this still did not work, try cutting along the long
            % axis
            if ~any(bw(bw) ~= bw_cut(bw))
                % Define the minor axis
                x1 = floor(centerX - (longAxis/2) * cosd(angle));
                y1 = floor(centerY - (longAxis/2) * sind(angle));
                x2 = ceil(centerX + (longAxis/2) * cosd(angle));
                y2 = ceil(centerY + (longAxis/2) * sind(angle));
                % Draw a line
                numPoints = sqrt(size(bw,1)^2+size(bw,2)^2);
                x = round(linspace(x1, x2, numPoints)); x(x<1)=1; x(x>size(bw_cut,2))=size(bw_cut,2);
                y = round(linspace(y1, y2, numPoints)); y(y<1)=1; y(y>size(bw_cut,1))=size(bw_cut,1);
                for iPx = 1:numPoints
                    bw_cut(y(iPx), x(iPx), :) = 0;
                end%iPx
                bw_cut = bwlabel(bw_cut,4);
            end%did not work
        end%did not work
        % Account for separating lines.
        bw_fill = bw_cut;
        ind = find(bw_cut==0 & bw==1);
        for iPx = 1:length(ind)
            % Get the current pixel's neighborhood
            neighborhood = zeros(size(bw));
            neighborhood(ind(iPx)) = 1;
            neighbors = imgradient(neighborhood, 'central')>0;
            neighbors(ind(iPx)) = 0;
            neighbor_IDs = unique(bw_cut(neighbors));
            bw_fill(ind(iPx)) = max(neighbor_IDs);
        end%iPx
        bw_cut = bw_fill;
        % Update granules_labeled
        bw_cut(ind) = max(bw_cut(:))+1 : max(bw_cut(:))+length(ind);
        bw_cut = bw_cut+max(granules_labeled(:))+1;
        granules_labeled(bw) = bw_cut(bw);
    end%iGranule
end%while
end%FCN:removeHugeRegions

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