![alt text][logo]

[logo]: CalciSeg_logo.png

## CalciSeg - A versatile approach for unsupervised segmentation of calcium imaging data.


Contents
--------
1.	[About](#about)  
2.	[Requirements](#requirements)  
3.	[Running CalciSeg](#running-calciseg)  
4.	[Cite](#cite)


About
------
Functional calcium imaging has emerged as a powerful tool to provide insight into neuronal activity patterns ranging from local networks to whole-brain dynamics. However, the versatility of imaging approaches and brain structures in non-model organisms poses a challenge in how datasets can be reliably and reproducibly segmented. Common techniques include, *e.g.,* the manual selection of the regions of interest (ROIs) and machine learning solutions such as those obtained with *k*-means clustering. Still, the former suffers from user bias and can be rather labor-intensive, while the latter is not always transferable among study systems and can be computationally expensive. Intending to overcome these issues and provide a universal tool for calcium imaging data segmentation and processing, we developed **CalciSeg**, a data-driven and reproducible approach for unsupervised segmentation of calcium imaging data. **CalciSeg** addresses the challenges associated with brain structure variability and user bias by offering a computationally efficient solution for automatically selecting the regions of interest.


Requirements
------------
In order to run **CalciSeg** you need to install *MATLAB* together with the *Image Processing Toolbox*, the *Statistics and Machine Learning Toolbox*, and the *Parallel Computing Toolbox*.


Running CalciSeg
----------------
Both **CalciSeg** and **CalciSeg_3D** are MATLAB functions that can be called from anywhere in your code. Note that both depend on the *Parallel Computing Toolbox*. It is reconmended to starts a parallel pool of workers before calling the function.
For **CalciSeg**, use the following input:

CalciSeg(stack, projection_method, init_segment_method, regmax_method, n_rep, refinement_method, minPixCount)
- stack : 3-D matrix (x*y*time)
- projection_method : method for calculating the projection across time
	- 'std'      - standard deviation projection
	- 'mean'     - mean projection
	- 'max'      - max projection
- init_seg_method : method for initial segmentation
	- 'voronoi'  - Delaunay triangulation
	- 'corr'     - local growing based on correlation (r_threshold = sqrt(0.7))
- regmax_method : method for determining how to identify local extrema                          
	- 'raw'      - simply apply imregionalmax/-min on the projection
	- 'filtered' - dilate/erode image to apply moving max/min filter before applying imregionalmax/-min                                           
	- 'both'     - combine both above-mentioned methods
- n_rep : Number of iterations for refining the regions.
- refinement_method : Measure to fine-tune granule assignment
	- 'corr' - correlation
	- 'rmse'  - root median square error
- minPixCount : Minimum pixel area. Use 'auto' for an automatic assessment. Or provide a number, *e.g.*, minPixCount=10 for at least 10 pixels per region.

For **CalciSeg_3D**, the input remains similar. Note the additional parameter *aspect_ratio* to account for differences in x-y-z dimensions. Further, the input *stack* is now expected to a 4-D matrix (x*y*z-time).

CalciSeg(stack, projection_method, init_segment_method, regmax_method, n_rep, refinement_method, minPixCount)

Both functions return two variables:
- pockets_labeled : resulting segmentation. An ID was assigned to each pixel.
- summary_stats : summarizing statisitcs.

Cite
----
If you use **CalciSeg**, please cite it using this DOI:
[![DOI]( *coming soon* )

Thank you for using **CalciSeg**!