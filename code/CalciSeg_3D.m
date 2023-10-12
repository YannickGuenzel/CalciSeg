function [pockets_labeled, summary_stats] = CalciSeg_3D(stack, aspect_ratio, projection_method, init_segment_method, regmax_method, n_rep, refinement_method, minPixCount)
% CalciSeg_3D(stack, aspect_ratio, projection_method, init_segment_method, regmax_method, n_rep, refinement_method, minPixCount)
% segments the spatial aspects of a stack (x*y*y*time) into individual
% regions ('granules') based on a refined 3D Delaunay triangulation.
% Note: CalciSeg_3D uses parallel for-loops at for the local growing 
% algorithm and the refinement process. It is reconmended to starts a 
% parallel pool of workers before calling the function.
%
% Input:
%       stack             : 4-D matrix (x*y*z,time)
%       aspect_ratio      : x-y-z aspect ratio. Leave empty for default [1 1 1]
%       projection_method : method for calculating the projection across time
%                             - 'std'    - standard deviation projection
%                             - 'mean'   - mean intensity projection
%                             - 'median' - median intensity projection
%                             - 'max'    - maximum intensity projection
%                             - 'min'    - minimum intensity projection
%                             - 'pca'    - principal component projection
%       init_seg_method   : method for initial segmentation
%                             - 'voronoi' - Delaunay triangulation
%                             - 'corr'    - local growing based on correlation (r_threshold = sqrt(0.7))
%       regmax_method     : method for determining how to identify local
%                           extrema
%                             - 'raw'      - simply apply imregionalmax/-min on the projection
%                             - 'filtered' - dilate/erode image to apply
%                                            moving max/min filter before
%                                            applying imregionalmax/-min 
%                             - 'both'     - combine both above-mentioned methods
%       n_rep             : Number of iterations for refining the regions.
%       refinement_method : Measure to fine-tune granule assignment
%                             - 'corr'     - correlation
%                             - 'rmse'     - root median square error
%       minPixCount       : Minimum pixel area. Use 'auto' for an automatic
%                           assessment. Or provide a number, e.g.
%                           minPixCount=10 for at least 10 pixels per
%                           region.
%
% Output:
%       pockets            : grayscale showing pocket distribution
%       pockets_labeled    : pocket labels
%       summary_stats      : avg, std, corr for each pocker
%                    .pocket_Avg           : within-pocket average activity
%                    .pocket_Corr_img      : average activity image (x*y*z)
%                    .pocket_Std           : each pocket's std over time
%                    .pocket_Std           : each pocket's std over time as image
%                    .pocket_Corr          : within-pocket correlation
%                    .pocket_Corr_img      : correlation image (x*y*z)
%                    .active_region.map    : binary map of active regions
%                    .active_region.method : method used for binarizing
%
% Version: 03-Oct-23 (R2023a)

% Validate inputs
stack = squeeze(stack);
error_message = '';
[valid, error_message] = validateInputs(error_message, stack, aspect_ratio, projection_method, init_segment_method, regmax_method, n_rep, refinement_method, minPixCount);
if valid == -1
    error(error_message);
elseif valid == 0
    warning(error_message);
end

% Initial preprocessing
[reshaped_stack, projection, refinement_method, init_segment_method] = initialPreprocessing(stack, projection_method, refinement_method, init_segment_method);

% Perform initial segmentation
if isnumeric(minPixCount) 
    filter_size = ceil(((3*minPixCount)/(4*pi))^(1/3)) + double(minPixCount==0);
elseif strcmp(minPixCount, 'auto')
    filter_size = 1;
end
pockets_labeled = initialSegmentation(stack, aspect_ratio, reshaped_stack, projection, regmax_method, init_segment_method, filter_size);

% Maybe, the user wants to have the min size estimated based on the data
if (isstring(minPixCount) || ischar(minPixCount)) && strcmp(minPixCount, 'auto')
    % Get pixel count for each region
    pocketList = unique(pockets_labeled);
    [~,~,C] = unique(pockets_labeled);
    C = histcounts(C, 1:length(pocketList)+1);
    C = [C(:), pocketList(:)];
    minPixCount = floor(quantile(C(2:end,1), 0.10));
    disp(['automatic assessment of the min. granule size: ', num2str(minPixCount)])
end%if auto min size

% Refine segmentation
pockets_labeled = refineSegmentation(reshaped_stack, pockets_labeled, n_rep, refinement_method, minPixCount);

% Final refinement steps and statistics calculation
[pockets_labeled, summary_stats] = finalRefinementAndStats(reshaped_stack, pockets_labeled, projection_method);

end%FCN:CalciSeg_optimized



function [valid, error_message] = validateInputs(error_message, stack, aspect_ratio, projection_method, init_segment_method, regmax_method, n_rep, refinement_method, minPixCount)
% Initialize output variables
valid = 1;

% Validate 'stack' input
if ndims(stack) ~= 4
    valid = 0;
    error_message = 'Input "stack" has no temporal component. Inputs are adjusted accordingly (projection_method="none"; init_segment_method="voronoi"; refinement_method="rmse"). ';
end%if invalid 'stack' input

% Validate 'aspect_ratio' input
if ~isnumeric(aspect_ratio) || length(aspect_ratio) ~= 3 || isempty(aspect_ratio)
    valid = 0;
    error_message = [error_message, 'Input "aspect_ratio" was set to be the default of [1 1 1]. '];
end%if invalid 'stack' input

% Validate 'projection_method' input
valid_projection_methods = {'std', 'mean', 'median', 'max', 'min', 'pca', 'none'};
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
end%if invalid 'minPixCount' input
end%FCN:validateInputs



function [reshaped_stack, projection, refinement_method, init_segment_method] = initialPreprocessing(stack, projection_method, refinement_method, init_segment_method)
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
    case 'none'
        projection = stack;
end%switch projection method

end%FCN:initialPreprocessing



function pockets_labeled = initialSegmentation(stack, aspect_ratio, reshaped_stack, projection, regmax_method, init_segment_method, filter_size)
switch init_segment_method
    case 'voronoi'
        % Get local maxima and minima
        [regional_maxima, regional_minima] = identifyLocalExtrema(projection, regmax_method, filter_size);
        % 2-D Delaunay triangulation
        pockets_labeled = triang_extrema(regional_maxima, regional_minima, aspect_ratio);
    case 'corr'
        pockets_labeled = growing_regions(projection, stack, reshaped_stack, regmax_method, filter_size);
end%switch init_segment_method
end%FCN:initialSegmentation



function pockets_labeled = triang_extrema(regional_maxima, regional_minima, aspect_ratio)
% Erode the extrema location to keep the controids only
ultimateErosion = bwulterode((double(regional_maxima) + double(regional_minima)), 26);
% Reduce to single voxels
[L, num] = bwlabeln(ultimateErosion, 26);
centroid_V = zeros(size(ultimateErosion));
% Iterate over all regions
for iVox = 1:num
    % Get the x, y, z coordinates of the current component
    [x,y,z] = ind2sub(size(ultimateErosion), find(L == iVox)); 
    % Calculate and round the centroid
    centroid = round(mean([x(:),y(:),z(:)],1)); 
    % Mark the centroid in the centroid_V matrix
    centroid_V(centroid(1), centroid(2), centroid(3)) = 1;
end%iVox
% Multiply the original map of extrema location in order to account for
% wrong assignments where there was no extrema location before
ultimateErosion = ultimateErosion.*centroid_V;
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
C = C.*aspect_ratio;
% Segment image based on these maxima using a 2-D Delaunay
% triangulation (i.e. assign pixels to the closest extrem point
% Convert linear indices to (x, y, z) coordinates
linearIndices = 1:numel(regional_maxima);
[x, y, z] = ind2sub(size(regional_maxima), linearIndices);
% Combine all and adjust aspect ratio;
allVoxels = [x(:), y(:), z(:)];
clear x y z linearIndices
% Use knnsearch to find the index of the closest point in 'C' for each voxel
pockets_labeled = knnsearch(C, allVoxels);
% Reshape the resulting index array to the dimensions of the 3D volume
pockets_labeled = reshape(pockets_labeled, size(regional_maxima));
% Remove border
ultimateErosion(1,:,:) = 0; ultimateErosion(:,1,:) = 0; ultimateErosion(:,:,1) = 0;
ultimateErosion(end,:,:) = 0; ultimateErosion(:,end,:) = 0; ultimateErosion(:,:,end) = 0;
linearIndices = find(ultimateErosion);
% Create a logical mask
mask = ismember(pockets_labeled, pockets_labeled(linearIndices));
pockets_labeled(~mask) = 0;
end%FCN:triang_extrema



function [regional_maxima, regional_minima] = identifyLocalExtrema(projection, regmax_method)
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



function pockets_labeled = growing_regions(projection, stack, reshaped_stack, regmax_method, filter_size)
% Get a map of local correlatios of a pixel with each of its neighbors
% --- Preallocation
corr_img = zeros(size(projection));
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
% Check whether to filter or not
if strcmp(regmax_method, 'filtered')
    corr_img = imdilate(corr_img, ones(3,3,3));
end
% Now, start with the best local correlation and gradually increase the
% region until we reach a threshold of a correlation that got too worse
% --- Preallocation
pockets_labeled = zeros(size(projection));
% --- Counter for regions
cnt_id = 1;
% Iterate over all pixels
hWait = waitbar(0, ['Growing regions. Please wait ... (0/', num2str(px_cnt),')']);
for iPx = 1:px_cnt
    waitbar(iPx/px_cnt, hWait, ['Growing regions. Please wait ... (', num2str(iPx), '/', num2str(px_cnt), ')']);
    % Get the current correlation value
    [~, max_id] = max(corr_img(:));
    if corr_img(max_id) >= sqrt(0.7)
        % This will be the origin of a new region
        pockets_labeled(max_id) = cnt_id;
        while true
            % Get the current region's TC
            curr_TC = mean(reshaped_stack(find(pockets_labeled==cnt_id),:),1);
            % Create mask  including all neighbors
            neighbors = pockets_labeled == cnt_id;
            neighbors = imdilate(neighbors, SE) - neighbors;
            % Get neighbors' indices
            [neighbors_ind_X, neighbors_ind_Y, neighbors_ind_Z] = ind2sub(size(neighbors), find(neighbors));
            % Keep track of whether pixels got added
            added_px = zeros(1, length(neighbors_ind_X));
            % Check neighbors
            for iN = 1:length(neighbors_ind_X)
                if (pockets_labeled(neighbors_ind_X(iN), neighbors_ind_Y(iN), neighbors_ind_Z(iN)) == 0) && (corr(curr_TC(:), squeeze(stack(neighbors_ind_X(iN), neighbors_ind_Y(iN), neighbors_ind_Z(iN),:))) >= sqrt(0.7))
                    pockets_labeled(neighbors_ind_X(iN), neighbors_ind_Y(iN), neighbors_ind_Z(iN)) = cnt_id;
                    added_px(iN) = 1;
                end%if valid addition
            end%iN
            % Check break criterion
            if sum(added_px) == 0
                break
            end% if no more
        end%while
        % Dont check the current pixels anymore
        corr_img(pockets_labeled == cnt_id) = -inf;
        cnt_id = cnt_id+1;
    end%if good enough
end%iPx
close(hWait)
end%FCN:growing_regions



function pockets_labeled = refineSegmentation(reshaped_stack, pockets_labeled, n_rep, refinement_method, minPixCount)
% Account for splitting of pockets. For this, iterate over
% all pocket and check whether the corresponding itentity
% can be found in more then one coherent regions
pockets_labeled = checkSplits(pockets_labeled);
% Account for regions that are too small. For this, check each region
% that has less than the required number of pixels and assign each
% pixel of the tiny region to the closest neighbor
pockets_labeled = removeTinyRegions(pockets_labeled, minPixCount);
% Get mask for neighbors
nMask = strel('sphere', 1);
% Refine pockets
hWait = waitbar(0, ['Refinig regions. Please wait ... (0/', num2str(n_rep),')']);
for iRep = 1:n_rep
    waitbar(iRep/n_rep, hWait, ['Refinig regions. Please wait ... (', num2str(iRep), '/', num2str(n_rep), ')']);
    % Keep track of changes
    previous_pockets_labeled = pockets_labeled;
    % Get avg signal trace per pocket
    pocketList = unique(pockets_labeled);
    pocketTC = nan(length(pocketList), size(reshaped_stack,2));
    for iPocket = 1:length(pocketList)
        idx = find(pockets_labeled==pocketList(iPocket));
        pocketTC(iPocket, :) = nanmean(reshaped_stack(idx,:),1);
    end%iPocket
    clear iC iR iPocket
    % Get pixels in border regions
    candidates = imgradient3(pockets_labeled, 'central')>0;
    idx_candidates = find(candidates);
    % Iterate over all pixels in border regions and check
    % whether they have to be re-assigned.
    switch refinement_method
        case 'corr'
            idx_candidates_newIdentiy = refinement_parfor_corr(idx_candidates, pockets_labeled, reshaped_stack, pocketTC, pocketList, nMask);
        case 'rmse'
            idx_candidates_newIdentiy = refinement_parfor_rmse(idx_candidates, pockets_labeled, reshaped_stack, pocketTC, pocketList, nMask);
    end
    % Update changed identities
    for iCan = 1:length(idx_candidates)
        pockets_labeled(idx_candidates(iCan)) = idx_candidates_newIdentiy(iCan);
    end%iCan
    % Account for splitting of pockets. For this, iterate over
    % all pocket and check whether the corresponding itentity
    % can be found in more then one coherent regions
    pockets_labeled = checkSplits(pockets_labeled);
    % Account for regions that are too small. For this, check each region
    % that has less than the required number of pixels and assign each
    % pixel of the tiny region to the closest neighbor
    pockets_labeled = removeTinyRegions(pockets_labeled, minPixCount);
    % Stop if no changes
    if sum(previous_pockets_labeled(:) == pockets_labeled(:)) == numel(pockets_labeled(:))
        break
    end% if change
end%iRep
close(hWait)
end%FCN:refineSegmentation



function idx_candidates_newIdentiy = refinement_parfor_corr(idx_candidates, pockets_labeled, reshaped_stack, pocketTC, pocketList, nMask)
% Get current labels and prepare output of new labels
idx_candidates_Identiy = pockets_labeled(idx_candidates);
idx_candidates_newIdentiy = idx_candidates_Identiy;
% Use a parallel for loop to speed things up
parfor iCan = 1:length(idx_candidates)
    % Create mask  including all neighbors
    neighbors = zeros(size(pockets_labeled));
    neighbors(idx_candidates(iCan)) = 1;
    neighbors = imdilate(neighbors, nMask) - neighbors;
    % Get neighbors' indices
    [neighbors_ind_X, neighbors_ind_Y, neighbors_ind_Z] = ind2sub(size(neighbors), find(neighbors));
    % Get neighborhood clusters
    nC = unique(pockets_labeled(neighbors_ind_X, neighbors_ind_Y, neighbors_ind_Z));
    % Get the pixel's TC
    activity_px = reshaped_stack(idx_candidates(iCan), :);
    % Check which neighborhood is a better fit
    voronoi_R = zeros(length(nC),1);
    for iN = 1:length(nC)
        activity_cluster = pocketTC(pocketList == nC(iN), :);
        r = corrcoef(activity_px, activity_cluster); voronoi_R(iN) = 1/r(2);
    end%iN
    % Assign new identity
    [~, iBest] = min(voronoi_R);
    idx_candidates_newIdentiy(iCan) = nC(iBest);
end%iCan
end%FCN:refinement_parfor_corr



function idx_candidates_newIdentiy = refinement_parfor_rmse(idx_candidates, pockets_labeled, reshaped_stack, pocketTC, pocketList, nMask)
% Get current labels and prepare output of new labels
idx_candidates_Identiy = pockets_labeled(idx_candidates);
idx_candidates_newIdentiy = idx_candidates_Identiy;
% Use a parallel for loop to speed things up
parfor iCan = 1:length(idx_candidates)
    % Create mask  including all neighbors
    neighbors = zeros(size(pockets_labeled));
    neighbors(idx_candidates(iCan)) = 1;
    neighbors = imdilate(neighbors, nMask) - neighbors;
    % Get neighbors' indices
    [neighbors_ind_X, neighbors_ind_Y, neighbors_ind_Z] = ind2sub(size(neighbors), find(neighbors));
    % Get neighborhood clusters
    nC = unique(pockets_labeled(neighbors_ind_X, neighbors_ind_Y, neighbors_ind_Z));
    % Get the pixel's TC
    activity_px = reshaped_stack(idx_candidates(iCan), :);
    % Check which neighborhood is a better fit
    voronoi_R = zeros(length(nC),1);
    for iN = 1:length(nC)
        activity_cluster = pocketTC(pocketList == nC(iN), :);
        voronoi_R(iN) = sqrt(median((activity_cluster-activity_px).^2));
    end%iN
    % Assign new identity
    [~, iBest] = min(voronoi_R);
    idx_candidates_newIdentiy(iCan) = nC(iBest);
end%iCan
end%Fcn:refinement_parfor_rmse



function pockets_labeled = checkSplits(pockets_labeled)
% Account for splitting of pockets. For this, iterate over
% all pocket and check whether the corresponding itentity
% can be found in more then one coherent regions
unique_pockets = unique(pockets_labeled(:));
% Get a mask to identify neighbors
SE = strel('sphere', 1);
for iPocket = 1:length(unique_pockets)
    % Neglect identity of zero
    % Get b/w image and identify number of regions
    bw = pockets_labeled==unique_pockets(iPocket);
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
                neighbor_IDs = unique(pockets_labeled(find(neighbors)));
                pockets_labeled(idx) = mode(neighbor_IDs(:));
            else
                pockets_labeled(idx) = max(unique_pockets)+1;
                unique_pockets = [unique_pockets(:); max(unique_pockets)+1];
            end
        end%iSplitter
    end% if more
end%iPocket
end%FCN:checkSplits



function pockets_labeled = removeTinyRegions(pockets_labeled, minPixCount)
% Kick out pockets that are too small. For this, iterate
% over all pockets, and check whether it is large enough. If
% not, assign each pixel to the best neigbour.
% Repeat as long as all small pockets are gone.
% --- Get a mask to identify neighbors
SE = strel('sphere', 1);
while true
    % Get pixel count for each region
    pocketList = unique(pockets_labeled);
    [~,~,C] = unique(pockets_labeled);
    C = histcounts(C, 1:length(pocketList)+1);
    C = [C(:), pocketList(:)];
    % Sort base don count
    C = sortrows(C,1,"ascend");
    % Iterate over all regions and take care of those that are too small
    if any(C(:,1)<minPixCount)
        for iPocket = 1:length(pocketList)
            if C(iPocket,1)<minPixCount
                % Get the ID of the current pocket
                curr_ID = C(iPocket,2);
                % Get a logical image for the current pocket
                bw = pockets_labeled==curr_ID;
                % If it is a single pixel, use the mode. Otherwise the closest
                % value that is not the current ID
                [x, y, z] = ind2sub(size(bw), find(bw));
                % Create mask  including all neighbors
                neighbors = imdilate(bw, SE) - bw;
                % Get IDs on the edge
                neighbor_IDs = unique(pockets_labeled(find(neighbors)));
                if C(iPocket,1) == 1
                    pockets_labeled(x, y, z) = mode(neighbor_IDs(:));
                else
                    % Iterate over all pixels of the region that is too small and
                    % assign them to the closest neighbor region
                    for iPx = 1:length(x)
                        id_table = nan(1, length(neighbor_IDs));
                        for iN = 1:length(neighbor_IDs)
                            [Nx, Ny, Nz] = ind2sub(size(bw), find(pockets_labeled==neighbor_IDs(iN)));
                            dist = [Nx(:), Ny(:), Nz(:)] - [x(iPx), y(iPx), z(iPx)];
                            min_dist = min(sqrt(sum(dist'.*dist'))');
                            id_table(iN) = min_dist;
                        end%iN
                        [~,best_id] = min(id_table);
                        pockets_labeled(x(iPx), y(iPx), z(iPx)) = neighbor_IDs(best_id);
                    end%iPx
                end% if single pixel
            else
                break
            end% if any pixels left
        end%iPocket
    else
        break
    end%if nothing to correct
end%while
end%FCN:removeTinyRegions



function [pockets_labeled, summary_stats] = finalRefinementAndStats(reshaped_stack, pockets_labeled, projection_method)
% First, rename IDs to start with 0 and end with n-1
pocketList = unique(pockets_labeled);
temp = pockets_labeled;
for iPocket = 1:length(pocketList)
    temp(pockets_labeled==pocketList(iPocket)) = iPocket-1;
end%iPocket
pockets_labeled = temp;

% Next, pool pixels belonging to the same pocket
% --- Get list of pockets
pocketList = unique(pockets_labeled);
% --- Preallocation
summary_stats.avgTCs = nan(length(pocketList), size(reshaped_stack,2));
summary_stats.allTCs = nan(size(reshaped_stack));
summary_stats.pocket_Avg = nan(length(pocketList), 1);
summary_stats.pocket_Std = nan(length(pocketList), 1);
summary_stats.pocket_Max = nan(length(pocketList), 1);
summary_stats.pocket_Corr = nan(length(pocketList), 1);
summary_stats.pocket_Avg_img = nan(size(pockets_labeled));
summary_stats.pocket_Std_img = nan(size(pockets_labeled));
summary_stats.pocket_Max_img = nan(size(pockets_labeled));
summary_stats.pocket_Corr_img = nan(size(pockets_labeled));
% Iterate over all pockets
for iPocket = 1:length(pocketList)
    % Get their index positions
    idx_pocket = find(pockets_labeled==pocketList(iPocket));
    % Get the avg time course of the current pocket
    summary_stats.avgTCs(iPocket,:) = nanmean(reshaped_stack(idx_pocket,:),1);
    summary_stats.allTCs(idx_pocket,:) = repmat(nanmean(reshaped_stack(idx_pocket,:),1), [length(idx_pocket),1]);
    % Get avg activity
    avg = mean(summary_stats.avgTCs(iPocket,:));
    summary_stats.pocket_Avg(iPocket,1) = avg;
    summary_stats.pocket_Avg_img(idx_pocket) = avg;
    % Get std time
    sd = nanstd(summary_stats.avgTCs(iPocket,:));
    summary_stats.pocket_Std(iPocket,1) = sd;
    summary_stats.pocket_Std_img(idx_pocket) = sd;
    % Get std time
    max_val = nanmax(summary_stats.avgTCs(iPocket,:));
    summary_stats.pocket_Max(iPocket,1) = max_val;
    summary_stats.pocket_Max_img(idx_pocket) = max_val;
    % Get within-pocket correlation
    r = tril(corrcoef(reshaped_stack(idx_pocket, :)'), -1);
    mask = tril(ones(size(r)), -1);
    r = nanmean(r(mask==1));
    summary_stats.pocket_Corr(iPocket,1) = r;
    summary_stats.pocket_Corr_img(idx_pocket) = r;
end%iPocket
% Based on the projection method, estimate which regions are active
switch projection_method
    case 'std'
        summary_stats.active_region.map = imbinarize(summary_stats.pocket_Std_img);
        summary_stats.active_region.method = 'std';
    case 'mean'
        summary_stats.active_region.map = imbinarize(summary_stats.pocket_Avg_img);
        summary_stats.active_region.method = 'mean';
    case 'max'
        summary_stats.active_region.map = imbinarize(summary_stats.pocket_Max_img);
        summary_stats.active_region.method = 'max';
    case 'none'
        summary_stats.active_region.map = [];
        summary_stats.active_region.method = 'none';
end%switch projection method
end%FCN:finalRefinementAndStats




























