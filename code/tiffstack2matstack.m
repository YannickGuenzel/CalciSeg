function matstack = tiffstack2matstack(tiffstack)
% matstack = tiffstack2matstack(tiffstack) takes the path to a tiff stack
% file as input and returns a 3D matrix representing the stack of images.
%
% Input:
% - filePath: A string representing the path to the tiff stack file.
%
% Output:
% - img_stack: A 3D matrix where each slice (3rd dimension) is an image
%               from the tiff stack.
%
% Usage:
% img_stack = tiffstack2matstack('path/to/your/tiffstack.tif');
%
% Version: 26-Aug-23 (R2023a)

% Check if the file exists
if ~exist(tiffstack, 'file')
    error('File not found');
end

% Get the information of the tiff stack
info = imfinfo(tiffstack);

% Get the dimensions of the tiff stack
num_images = numel(info);
img_width = info(1).Width;
img_height = info(1).Height;

% Get the data format
test_img = imread(tiffstack, 1, 'Info', info);
test_img_info = whos('test_img');

% Initialize an empty 3D matrix to store the tiff stack
matstack = zeros(img_height, img_width, num_images, test_img_info.class);

% Loop through each slice and store it in the 3D matrix
for k = 1:num_images
    matstack(:, :, k) = imread(tiffstack, k, 'Info', info);
end

