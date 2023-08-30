function [pockets_labeled, summary_stats] = CalciSeg(stack, projection_method, regmax_method, n_rep, minPixCount)
% CalciSeg(stack, projection_method, regmax_method, n_rep, minPixCount) segments the spatial
% aspects of a 3-dimensional stack (x*y*time) into individual regions
% based on a refined 2D Delaunay triangulation.
%
% Input:
%       stack             : 3-D matrix (x*y*time)
%       projection_method : method for calculating the projection across time
%                             - 'std'  - standard deviation projection
%                             - 'mean' - mean projection
%                             - 'max'  - max projection
%       regmax_method     : method for determining how to identify local maxima%
%                             - 'raw'      - simply apply imregionalmax on the projection
%                             - 'filtered' - apply 2-D Gaussian filtering (sd=1) and imregionalmax on movmax-filtered projection
%                             - 'both'     - combine both above-mentioned methods
%       n_rep             : Number of iterations for refining the regions.
%       minPixCount       : Minimum pixel area. Use 'auto' for an automatic
%                           assessment. Or provide a number, e.g. 
%                           minPixCount=10 for at least 10 pixels per 
%                           region.
%
% Output:
%       pockets            : grayscale showing pocket distribution
%       pockets_labeled    : pocket labels
%       summary_stats      : avg, std, corr for each pocker
%                    .pocket_Avg         : within-pocket average activity
%                    .pocket_Corr_img    : average activity image (x*y)
%                    .pocket_Std         : each pocket's std over time
%                    .pocket_Std         : each pocket's std over time as image
%                    .pocket_Corr        : within-pocket correlation
%                    .pocket_Corr_img    : correlation image (x*y)
%
% Version: 26-Aug-23 (R2023a)

% -------------------------------------------------------------------------

% First, check whether the 3rd dimension of the stack (i.e. time) is larger
% than 1
if size(stack,3)==1
    projection_method = 'none';
    n_rep = 0;
end

% Get a reshaped version of the stack. Later, this will be used
% to correlate pixels and pockets
reshaped_stack = reshape(stack, [size(stack,1)*size(stack,2), size(stack,3)]);

% -------------------------------------------------------------------------

% Get projection from stack
switch projection_method
    case 'std'
        BF = nanstd(stack, [], 3);
    case 'mean'
        BF = nanmean(stack, 3);
    case 'max'
        BF = nanmax(stack, [], 3);
    case 'none'
        BF = stack;
end%switch
% Substitute NaN with the minimum value
BF(isnan(BF)) = min(BF(:));

% -------------------------------------------------------------------------

% Identify local extrema
switch regmax_method
    case 'raw'
        regmax_img = imregionalmax(BF,8)>0;
        regmin_img = imregionalmin(BF,8)>0;
    case 'filtered'
        % Filter to combine pockets that belong to the same large scale
        % maximum/minimum
        BF = imgaussfilt(BF, 1);
        BF_std_max_filtered = movmax(movmax(BF',3)',3);
        regmax_img = imregionalmax(BF_std_max_filtered,8)>0;
        BF_std_min_filtered = movmin(movmin(BF',3)',3);
        regmin_img = imregionalmin(BF_std_min_filtered,8)>0;
    case 'both'
        % Filter to combine pockets that belong to the same large scale
        % maximum/minimum
        BF_std_max_filtered = movmax(movmax(BF',3)',3);
        regmax_img = (imregionalmax(BF_std_max_filtered,8) + imregionalmax(BF,8))>0;
        BF_std_min_filtered = movmin(movmin(BF',3)',3);
        regmin_img = (imregionalmin(BF_std_min_filtered,8) + imregionalmin(BF,8))>0;
end%switch
ultimateErosion = bwulterode(regmax_img>0) + bwulterode(regmin_img>0);
% Include points that are on the edge of the image. Later, we will set all
% regions that are at the edge to zero
ultimateErosion(1,:) = 1; ultimateErosion(:,1) = 1;
ultimateErosion(end,:) = 1; ultimateErosion(:,end) = 1;
[locMax.y, locMax.x] = find(ultimateErosion);

% -------------------------------------------------------------------------

% Segment image based on these maxima using a 2-D Delaunay
% triangulation
% Get linear indices of local maxima
linearIndices = find(ultimateErosion);
% Convert linear indices to (x, y, z) coordinates
[Cx, Cy] = ind2sub(size(ultimateErosion), linearIndices);
% Combine coordinates into Nx3 matrix
C = [Cx, Cy]; clear Cx Cy

% Segment image based on these maxima using a 3-D Delaunay
% triangulation (i.e. assign pixels to the closest extrem point
% Convert linear indices to (x, y, z) coordinates
linearIndices = 1:numel(BF);
[x, y] = ind2sub(size(BF), linearIndices);
% Combine all and adjust aspect ratio;
allVoxels = [x(:), y(:)];
clear x y z linearIndices

% Use knnsearch to find the index of the closest point in 'C' for each voxel
pockets_labeled = knnsearch(C, allVoxels);
% Reshape the resulting index array to the dimensions of the 3D volume
pockets_labeled = reshape(pockets_labeled, size(BF));

% Remove border
ultimateErosion(1,:) = 0; ultimateErosion(:,1) = 0;
ultimateErosion(end,:) = 0; ultimateErosion(:,end) = 0;
linearIndices = find(ultimateErosion);
% Create a logical mask
mask = ismember(pockets_labeled, pockets_labeled(linearIndices));
pockets_labeled(~mask) = 0;

% Check whether we should assess the min number of pixels
if (isstring(minPixCount) || ischar(minPixCount)) && strcmp(minPixCount, 'auto')
    pocketList = unique(pockets_labeled);
    [~,~,C] = unique(pockets_labeled);
    C = histcounts(C, 1:length(pocketList)+1);
    % Maybe, the user wants to have the min size estimated based on the data
    % minPixCount = floor(quantile(C(2:end), 0.10));
    minPixCount = floor((numel(BF)/length(pocketList))/2);
    disp('- - - - -')
    disp(['automatic assessment of the min. granule size: ', num2str(minPixCount)])
    disp('- - - - -')
    clear pocketList C
end%if auto min size

% -------------------------------------------------------------------------

% Refine edges. Check first whether the 3rd dimension of the stack (i.e.
% time) is larger than 1
if ~size(stack,3)>1
    n_rep = 0;
end
% Get mask for neighbors
nMask = [0 1 0; 1 1 1; 0 1 0];
% Refine pockets
hWait = waitbar(0, ['Please wait ... (0/', num2str(n_rep)]);
for iRep = 1:n_rep
    waitbar(iRep/n_rep, hWait, ['Please wait ... (', num2str(iRep), '/', num2str(n_rep), ')']);
    % Keep track of changes
    previous_pockets_labeled = pockets_labeled;
    % Get avg signal trace per pocket
    pocketList = unique(pockets_labeled);
    pocketTC = nan(length(pocketList), size(stack,3));
    for iPocket = 1:length(pocketList)
        idx = find(pockets_labeled==pocketList(iPocket));
        pocketTC(iPocket, :) = nanmean(reshaped_stack(idx,:));
    end%iPocket
    clear iC iR iPocket
    % Get pixels in border regions
    candidates = imgradient(pockets_labeled, 'central')>0;
    idx_candidates = find(candidates);
    idx_table = reshape(1:numel(pockets_labeled), size(pockets_labeled));
    % Iterate over all pixels in border regions and check
    % whether thez have to be re-assigned.
    idx_candidates_newIdentiy = refinement_parfor(idx_candidates, idx_table, pockets_labeled, reshaped_stack, pocketTC, pocketList, nMask);
    % Update changed identities
    for iCan = 1:length(idx_candidates)
        pockets_labeled(idx_candidates(iCan)) = idx_candidates_newIdentiy(iCan);
    end%iCan
    % Account for splitting of pockets. For this, iterate over
    % all pocket and check whether the corresponding itentity
    % can be found in more then one coherent regions
    unique_pockets = unique(pockets_labeled(:));
    for iPocket = 1:length(unique_pockets)
        % Neglect identity of zero
        % Get b/w image and identify number of regions
        bw = pockets_labeled==unique_pockets(iPocket);
        L = bwlabel(bw,4);
        % If there is more than one region (i.e. more than
        % the labels 0 and 1, assign new identities
        unique_L = unique(L(:));
        if length(unique_L)>2
            % Iterate over all splitter
            for iSplitter = 3:length(unique_L)
                idx = find(L == unique_L(iSplitter));
                pockets_labeled(idx) = max(unique_pockets)+1;
                unique_pockets = [unique_pockets(:); max(unique_pockets)+1];
            end%iSplitter
        end% if more
    end%iPocket
    % Stop if no changes
    if sum(previous_pockets_labeled(:) == pockets_labeled(:)) == numel(pockets_labeled(:))
        break
    end% if change
end%iRep
close(hWait)
% Apply mask again
pockets_labeled(~mask) = 0;

% -------------------------------------------------------------------------

% Now, kick out pockets that are too small. For this, iterate
% over all pockets, and check whether it is large enough. If
% not, assign each pixel to the best neigbour.
% Repeat as long as all small pockets are gone.
pocketList = unique(pockets_labeled);
[~,~,C] = unique(pockets_labeled);
C = histcounts(C, 1:length(pocketList)+1);
% Maybe, the user wants to have the min size estimated based on the data
if (isstring(minPixCount) || ischar(minPixCount)) && strcmp(minPixCount, 'auto')
    minPixCount = floor(quantile(C(2:end), 0.10));
    disp(['automatic assessment of the min. granule size: ', num2str(minPixCount)])
end%if auto min size
if any(C<=minPixCount)
    didCorrect = 0;
    while true
        % --- Get avg signal trace per pocket
        pocketList = unique(pockets_labeled);
        pocketTC = nan(length(pocketList), size(stack,3));
        for iPocket = 1:length(pocketList)
            idx = find(pockets_labeled==pocketList(iPocket));
            pocketTC(iPocket, :) = nanmean(reshaped_stack(idx,:),1);
        end%iPocket
        clear iC iR iPocket
        % --- Get tabel with correct indices
        idx_table = reshape(1:numel(pockets_labeled), size(pockets_labeled));
        % --- Get unique pockets
        unique_pockets = unique(pockets_labeled(:));
        for iPocket = 1:length(unique_pockets)
            % Neglect identity of zero
            if unique_pockets(iPocket)~=0
                % Get b/w image to identify number of pixels
                bw = pockets_labeled==unique_pockets(iPocket);
                if sum(bw(:))<minPixCount
                    % Get indices
                    idx_candidates = find(bw);
                    % Get current identity
                    currC = unique_pockets(iPocket);
                    % Iterate over all pixels
                    for iPx = 1:length(idx_candidates)
                        % Get indices of neighbors
                        idx_neighbor = conv2( double(idx_table==idx_candidates(iPx)), nMask, 'same')==1;
                        % Get neighborhood clusters
                        nC = unique(pockets_labeled(idx_neighbor));
                        % Kick out current cluster
                        nC(nC==currC) = [];
                        if ~isempty(nC)
                            nC_cut = nC;
                            % Check whether neighbors are large enough
                            for iN = 1:length(nC)
                                bw = pockets_labeled==nC(iN);
                                if sum(bw(:))+1<minPixCount
                                    nC_cut(nC_cut==nC(iN)) = [];
                                end%if neighbor would be too small
                            end%iN
                        end%if isempty(nC)
                        nC = nC_cut;
                        % Assign identity of best fitting neighbor
                        if ~isempty(nC)
                            voronoi_R = zeros(length(nC),1);
                            for iN = 1:length(nC)
                                activity_px = reshaped_stack(idx_candidates(iPx), :);
                                activity_cluster = pocketTC(pocketList == nC(iN), :);
                                r = corrcoef(activity_px, activity_cluster); voronoi_R(iN) = r(2);
                                % Keep track of whether we have performed a
                                % correction
                                didCorrect = 1;
                            end%iN
                            % Assign new identity
                            [~, iBest] = max(voronoi_R);
                            pockets_labeled(idx_candidates(iPx)) = nC(iBest);
                        end%if isempty(nC)
                    end%iSplitter
                end% if more
            end% if not label 0
        end%iPocket
        % Check whether there was something to correct.
        if ~didCorrect
            break
        else
            didCorrect = 0;
        end%
    end%while
end%if any

% -------------------------------------------------------------------------

% Get each pockets summary statistics
if nargout > 1
    pocketList = unique(pockets_labeled);
    pocketList = pocketList(pocketList~=0);
    % Preallocation
    summary_stats.pocket_Avg = nan(length(pocketList), 1);
    summary_stats.pocket_Std = nan(length(pocketList), 1);
    summary_stats.pocket_Corr = nan(length(pocketList), 1);
    summary_stats.pocket_Avg_img = nan(size(pockets_labeled));
    summary_stats.pocket_Std_img = nan(size(pockets_labeled));
    summary_stats.pocket_Corr_img = nan(size(pockets_labeled));
    % Iterate over all pockets
    for iPocket = 1:length(pocketList)
        % Get their index positions
        idx_pocket = find(pockets_labeled==pocketList(iPocket));
        % Get avg activity
        avg = nanmean(nanmean(reshaped_stack(idx_pocket,:)));
        summary_stats.pocket_Avg(iPocket,1) = avg;
        summary_stats.pocket_Avg_img(idx_pocket) = avg;
        % Get std over time
        sd = nanstd(nanmean(reshaped_stack(idx_pocket,:)));
        summary_stats.pocket_Std(iPocket,1) = sd;
        summary_stats.pocket_Std_img(idx_pocket) = sd;
        % Get within-pocket correlation
        r = tril(corrcoef(reshaped_stack(idx_pocket, :)'), -1);
        mask = tril(ones(size(r)), -1);
        r = nanmean(r(mask==1));
        summary_stats.pocket_Corr(iPocket,1) = r;
        summary_stats.pocket_Corr_img(idx_pocket) = r;
    end%iPocket
    % Get active regions base don std
    summary_stats.active_mask = imbinarize(summary_stats.pocket_Std_img);
end% if arg

end%FCN:stackSegmentation

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

function idx_candidates_newIdentiy = refinement_parfor(idx_candidates, idx_table, pockets_labeled, reshaped_stack, pocketTC, pocketList, nMask)
% Get current labels and prepare output of new labels
idx_candidates_Identiy = pockets_labeled(idx_candidates);
idx_candidates_newIdentiy = idx_candidates_Identiy;
% Use a parallel for loop to speed things up
parfor iCan = 1:length(idx_candidates)
    % Get indices of neighbors
    idx_neighbor = conv2( double(idx_table==idx_candidates(iCan)), nMask, 'same')==1;
    % Get neighborhood clusters
    nC = unique(pockets_labeled(idx_neighbor));
    % Check which neighborhood is a better fit
    voronoi_R = zeros(length(nC),1);
    for iN = 1:length(nC)
        activity_px = reshaped_stack(idx_candidates(iCan), :);
        activity_cluster = pocketTC(pocketList == nC(iN), :);
        r = corrcoef(activity_px, activity_cluster); voronoi_R(iN) = r(2);
    end%iN
    % Assign new identity
    [~, iBest] = max(voronoi_R);
    idx_candidates_newIdentiy(iCan) = nC(iBest);
end%iCan
end%Fcn:refinement_parfor