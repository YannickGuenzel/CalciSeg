![alt text][logo]

[logo]: CalciSeg_logo.png

CalciSeg - A versatile approach for unsupervised segmentation of calcium imaging data
=====================================================================================


Contents
--------
1.	[About](#about)  
2.	[Requirements](#requirements)  
3.	[Running CalciSeg](#running-calciseg)  
4.	[Cite code as](#cite)


About
-----
Recent advances in calcium imaging, including the development of fast and sensitive genetically encoded indicators, high-resolution camera chips for wide-field imaging, and resonant scanning mirrors in laser scanning microscopy, have notably improved the temporal and spatial resolution of functional imaging analysis. Nonetheless, the variability of imaging approaches and brain structures challenges the development of versatile and reliable segmentation methods. Standard techniques, such as manual selection of regions of interest or machine learning solutions, often fall short due to either user bias, non-transferability among systems, or computational demand. To overcome these issues, we developed **CalciSeg**, a data-driven and reproducible approach for unsupervised functional calcium imaging data segmentation. **CalciSeg** addresses the challenges associated with brain structure variability and user bias by offering a computationally efficient solution for automatic image segmentation based on two parameters: regions' size limits and number of refinement iterations. We evaluated **CalciSeg** efficacy on datasets of varied complexity, different insect species (locusts, bees, and cockroaches), and imaging systems (wide-field, confocal, and multiphoton), showing the robustness and generality of our approach. Finally, the user-friendly nature and the open-source availability of **CalciSeg** facilitate the integration of this algorithm into existing analysis pipelines.

Requirements
------------
In order to run **CalciSeg** you need to install *MATLAB* together with the *Image Processing Toolbox*, the *Statistics and Machine Learning Toolbox*, and the *Parallel Computing Toolbox*.


Running CalciSeg
----------------
Both **CalciSeg** and **CalciSeg_3D** are MATLAB functions that can be called from anywhere in your code. Note that both depend on the *Parallel Computing Toolbox*. It is reconmended to starts a parallel pool of workers before calling the function.
For **CalciSeg**, use the following input:
CalciSeg(stack, *Name*, *Value*)

### Input parameter ###
- stack (matrix): 3-D dataset (x,y,time) that should be segmented.

### Optional input parameter ###
- 'projection_method' (string) Method for calculating the projection across time. Default: 'std'.
	- 'std' : Standard deviation projection.
	- 'mean' : Mean intensity projection.
	- 'median' : Median intensity projection.
	- 'max' : Maximum intensity projection.
	- 'min' : Minimum intensity projection.
	- 'pca' : Principal component projection. For this, "stack" is projected into PC space.
	- 'corr' : Correlation space as local correlation between  neighboring pixels.
	- 'none' : Applicable when "stack" has no third dimension or when to take the first slice of the 3-D stack.
- 'init_seg_method' (string) Method for the initial segmentation. Default: 'voronoi'.
	- 'voronoi' : Delaunay triangulation.
	- 'corr' : Local growing based on correlation. This will set 'projection_method' to be 'corr'.
	- 'rICA' : Reconstruction independent component analysis.
- 'regmax_method' (string) Method for determining how to identify local extrema in the intensity projection. Default: 'raw'.
	- 'raw' : Simply apply imregionalmax/-min on the projection. 
	- 'filtered' : Dilate/erode image to apply moving max/min filter before applying imregionalmax/-min. Note, this will also smooth the correlation space when init_seg_method is set to be 'corr', or the independent components accordingly when it is set to be 'rICA'.
	- 'both' : Combine both above-mentioned methods. Note, this is only applicable when init_seg_method is set to be 'voronoi'. 
- 'n_rep' (integer) Number of iterations for refining the regions. Default: 0.
- 'refinement_method' (string) Measure to fine-tune granule assignment during refinement iterations. Default: 'rmse'.
	- 'rmse' : Root median square error.
	- 'corr' : Pearson  correlation coefficient.                     
- 'limitPixCount' (matrix or string) Limits the pixel area per granule. Default:  [1 inf].
	- [a, b] : The minimum (a) and maximum (b) number of pixels that can be assigned to a granule. Note that, the minimum size affects the filter size for the regmax_method input.
	- 'auto' : An automatic assessment based on the distribution of granule sizes before refinement. Here, a is set to be the 5th quantile and b to be the 95th quantile. The filter size for regmax_method is set to 1.
- 'corr_thresh' (number) Threshold for the Pearson correlation coefficient when refinement_method is set to be 'corr'. Default: 0.85.
- 'n_rICA' (integer) Number of features to extract during reconstruction independent component analysis. Default: 0 to use the full number of components. Otherwise provide a number larger than zero for either an undercomplete (n_rICA < number of frames) or overcomplete (n_rICA > number of frames) feature representations.
- 'n_PCA' (number) Percentage of the total variance of the data that should be kept if the projection_method is set to 'pca'. The value must be larger than zero and not exceed 100. Default: 100.

For **CalciSeg_3D**, the input remains similar. Note the additional parameter *aspect_ratio* to account for differences in x-y-z dimensions. Further, the input *stack* is now expected to a 4-D matrix (x,y,z,time) and there is only an option for a minimum region size. See the function's documentation for more information.


### Both functions return two variables ###
- pockets_labeled: resulting segmentation. An ID was assigned to each pixel.
- summary_stats: summarizing statisitcs.


Cite code as
------------
If you use **CalciSeg**, please cite it using this DOI:
[![DOI](https://zenodo.org/badge/671860103.svg)](https://zenodo.org/doi/10.5281/zenodo.11190097)


Thank you for using **CalciSeg**!