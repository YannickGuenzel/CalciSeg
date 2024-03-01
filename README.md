![alt text][logo]

[logo]: CalciSeg_logo.png

CalciSeg - A versatile approach for unsupervised segmentation of calcium imaging data
=====================================================================================


Contents
--------
1.	[About](#about)  
2.	[Requirements](#requirements)  
3.	[Running CalciSeg](#running-calciseg)  
4.	[Cite](#cite)


About
-----
Functional calcium imaging has emerged as a powerful tool to provide insight into neuronal activity patterns ranging from local networks to whole-brain dynamics. However, the versatility of imaging approaches and brain structures in non-model organisms poses a challenge in how datasets can be reliably and reproducibly segmented. Common techniques include, *e.g.,* the manual selection of the regions of interest (ROIs) and machine learning solutions such as those obtained with *k*-means clustering. Still, the former suffers from user bias and can be rather labor-intensive, while the latter is not always transferable among study systems and can be computationally expensive. Intending to overcome these issues and provide a universal tool for calcium imaging data segmentation and processing, we developed **CalciSeg**, a data-driven and reproducible approach for unsupervised segmentation of calcium imaging data. **CalciSeg** addresses the challenges associated with brain structure variability and user bias by offering a computationally efficient solution for automatically selecting the regions of interest.


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
- 'projection_method' (string): Method for calculating the projection across time. Default: 'std'.
	- 'std': Standard deviation projection
	- 'mean': Mean intensity projection
	- 'median': Median intensity projection
	- 'max': Maximum intensity projection
	- 'min': Minimum intensity projection
	- 'pca': Principal component projection. For this, "stack" is projected into PC space.
	- 'corr': Correlation space as local correlation between neighboring pixels
	- 'none': Applicable when "stack" has no third dimension or when to take the first slice of the 3-D stack

- 'init_seg_method' (string): Method for the initial segmentation. Default: 'voronoi'.
	- 'voronoi': Delaunay triangulation
	- 'corr': Local growing based on correlation. This will set *projection_method* to be 'corr'.
	- 'rICA': Reconstruction independent component analysis.

- 'regmax_method' (string): Method for determining how to identify local extrema in the intensity projection. Default: 'raw'.
	- 'raw': Simply apply imregionalmax/-min on the projection. 
	- 'filtered': Dilate/erode image to apply moving max/min filter before applying imregionalmax/-min. Note, this will also smooth the correlation space when *init_seg_method* is set to be 'corr', or the independent components accordingly when it is set to be 'rICA'.
	- 'both': Combine both above-mentioned methods. Note, this is only applicable when *init_seg_method* is set to be 'voronoi'. 

- 'n_rep' (integer): Number of iterations for refining the regions. Default: 0.

- 'refinement_method' (string) Measure to fine-tune granule assignment during refinement iterations. Default: 'rmse'.
	- 'rmse': Root median square error
	- 'corr': Pearson correlation coefficient
 
- 'limitPixCount' (matrix or string): Limits the pixel area per granule. Default: [1 inf].
	- [a, b] : The minimum (a) and maximum (b) number of pixels that can be assigned to a granule. Note that, the minimum size affects the filter size for the *regmax_method* input.
 	- 'auto' : An automatic assessment based on the distribution of granule sizes before refinement. Here, a is set to be the 5th quantile and b to be the 95th quantile. The filter size for *regmax_method* is set to 1.

 -'corr_thresh' (number): Threshold for the Pearson correlation coefficient when refinement_method is set to be 'corr'.  Default: 0.85.

 - 'n_rICA' (integer): Number of features to extract during reconstruction independent component analysis. Default: 0 to use the full number of components. Otherwise provide a number larger than zero for either an undercomplete (n_rICA < number of frames) or overcomplete (n_rICA > number of frames) feature representations.

For **CalciSeg_3D**, the input remains similar. Note the additional parameter *aspect_ratio* to account for differences in x-y-z dimensions. Further, the input *stack* is now expected to a 4-D matrix (x,y,z,time) and there is only an option for a minimum region size. See the function's documentation for more information.


### Both functions return two variables ###
- pockets_labeled: resulting segmentation. An ID was assigned to each pixel.
- summary_stats: summarizing statisitcs.

Cite
----
If you use **CalciSeg**, please cite it using this DOI:
[![DOI]( *coming soon* )

Thank you for using **CalciSeg**!