![alt text][logo]

[logo]: CalciSeg_logo.png

CalciSeg - A versatile approach for unsupervised segmentation of calcium imaging data
=====================================================================================


Contents
--------
1. [About](#about)
2. [Requirements](#requirements)
3. [Running CalciSeg](#running-calciseg)
4. [Output](#output)
5. [Recent changes](#recent-changes)
6. [Cite code as](#cite-code-as)
7. [Funding](#funding)


About
-----
**CalciSeg** is a MATLAB-based method for unsupervised segmentation of calcium imaging data into spatial regions ("granules") using projection-based initialization, temporal similarity, region-size constraints, and optional iterative refinement. It is designed for 2-D or 3-D calcium imaging stacks with spatial dimensions _(x, y)_ and an optional temporal dimension _(time)_.

The current **CalciSeg** function supports multiple projection methods, several initialization strategies, optional refinement, resumable segmentation, fragmented-region pooling, temporal downsampling during segmentation, and per-granule quality metrics.


Requirements
------------
To run **CalciSeg**, install MATLAB together with:

- **Image Processing Toolbox**
- **Statistics and Machine Learning Toolbox**
- **Parallel Computing Toolbox**

**CalciSeg** uses _parfor_ in several stages. Starting a parallel pool before calling the function is recommended for better performance.


Running CalciSeg
----------------
Call the function as:

```
[granules_labeled, summary_stats] = CalciSeg(stack, Name, Value);
```

### Required input

- _stack_ (matrix): calcium imaging dataset.
  - 3-D input is expected as _(x, y, time)_.
  - 2-D input is also accepted. In that case, **CalciSeg** automatically adjusts settings to:
    - _projection_method = 'none'_
    - _init_seg_method = 'voronoi'_
    - _refinement_method = 'rmse'_
    - _temporal_bin = 1_

### Optional input parameters

- **_projection_method_** (char)  
  Method for calculating the projection across time.  
  **Default:** _'std'_

  - _'std'_    : standard deviation projection
  - _'mean'_   : mean intensity projection
  - _'median'_ : median intensity projection
  - _'max'_    : maximum intensity projection
  - _'min'_    : minimum intensity projection
  - _'pca'_    : principal component projection; the stack is projected into PC space
  - _'corr'_   : local correlation image based on neighboring pixels
  - _'none'_   : use the first slice of the stack, or use when no temporal dimension is present

- **_init_seg_method_** (char)  
  Method for the initial segmentation.  
  **Default:** _'voronoi'_

  - _'voronoi'_ : Voronoi-style partitioning based on local extrema
  - _'corr'_    : local growing based on temporal correlation
  - _'rICA'_    : reconstruction independent component analysis

  **Note:** If _init_seg_method_ is set to _'corr'_, _projection_method_ is automatically set to _'corr'_.

- **_regmax_method_** (char)  
  Method for identifying local extrema in the projection.  
  **Default:** _'raw'_

  - _'raw'_      : apply _imregionalmax_ / _imregionalmin_ directly
  - _'filtered'_ : prefilter with morphological dilation / erosion before extrema detection
  - _'both'_     : combine raw and filtered extrema detection

  **Note:** _regmax_method = 'both'_ is only applicable with _init_seg_method = 'voronoi'_. Otherwise it is reset to _'raw'_.

- **_n_rep_** (integer)  
  Number of refinement iterations.  
  **Default:** _0_

- **_refinement_method_** (char)  
  Distance measure used during refinement.  
  **Default:** _'rmse'_

  - _'rmse'_ : root median square error
  - _'corr'_ : Pearson correlation coefficient

- **_limitPixCount_** (numeric vector or char)  
  Minimum and maximum pixel area per granule.  
  **Default:** _[1 inf]_

  - _[a, b]_ : lower and upper pixel-count limits
  - _'auto'_ : estimate limits from the initial segmentation using the 5th and 95th quantiles of granule sizes

  **Note:** The minimum size also affects the filter size used by _regmax_method_.

- **_corr_thresh_** (scalar number)  
  Pearson-correlation threshold used when correlation-based decisions are made.  
  **Default:** _0.85_

- **_n_rICA_** (integer)  
  Number of features to extract during reconstruction ICA.  
  **Default:** _0_

  - _0_ uses the full number of components
  - values _> 0_ allow undercomplete or overcomplete feature representations

- **_n_PCA_** (scalar number)  
  Percentage of total variance to keep when _projection_method = 'pca'_.  
  Must be _> 0_ and _<= 100_.  
  **Default:** _100_

- **_fillmissing_** (logical)  
  Replace _NaN_ or _Inf_ values by interpolation of neighboring nonmissing values.  
  **Default:** _false_

- **_resume**CalciSeg**** (2-D matrix)  
  Label map from a previous segmentation session.  
  Providing this skips the initial segmentation step and starts from refinement / post-processing.  
  **Default:** _[]_

- **_pool_fragments_** (logical)  
  Pool neighboring granules that likely belong to the same biological structure if their average time courses correlate above _corr_thresh_.  
  The neighbor comparison is repeated until no adjacent pair exceeds the threshold.  
  **Default:** _false_

  **Note:** This step is applied after refinement and therefore ignores the upper pixel-count limit.

- **_temporal_bin_** (integer)  
  Temporal downsampling factor for faster segmentation. The stack is binned by averaging every _N_ consecutive frames before segmentation. Final summary statistics are computed from the full-resolution time courses.  
  **Default:** _1_

### Example

``` matlab
[granules_labeled, summary_stats] = CalciSeg(...
    stack, ...
    'projection_method', 'std', ...
    'init_seg_method', 'voronoi', ...
    'regmax_method', 'filtered', ...
    'n_rep', 25, ...
    'refinement_method', 'rmse', ...
    'limitPixCount', [5 100], ...
    'corr_thresh', 0.85, ...
    'fillmissing', true, ...
    'pool_fragments', false, ...
    'temporal_bin', 2);
```

### Note on CalciSeg_3D

**CalciSeg_3D** is provided separately. For current 3D-specific inputs and behavior, consult the inline documentation of that function directly.


Output
------
**CalciSeg** returns two outputs:

- **_granules_labeled_**  
  2-D label map of segmented granules. Its size matches the first two dimensions of _stack_.

- **_summary_stats_**  
  Structure containing per-granule summary statistics and segmentation metadata.

### Fields in summary_stats

- _projection_  
  Projection image used during initial segmentation

- _avgTCs_  
  Mean time course of each granule

- _granule_Corr_  
  Average within-granule correlation

- _granule_Avg_img_  
  Mean activity projected back into image space

- _granule_Std_img_  
  Standard deviation projected back into image space

- _granule_Max_img_  
  Maximum projected back into image space

- _granule_Min_img_  
  Minimum projected back into image space

- _granule_Median_img_  
  Median projected back into image space

- _granule_Corr_img_  
  Within-granule correlation projected back into image space

- _active_region.map_  
  Binary map of active regions

- _active_region.method_  
  Projection method used for active-region binarization

- _quality.SNR_  
  Per-granule signal-to-noise ratio

- _quality.compactness_  
  Per-granule compactness  
  (_4*pi*area / perimeter^2_, where _1_ corresponds to a perfect circle)



Recent changes
--------------
VERSION 1.2: Fixed many segmentation and statistics bugs, including issues with handling tiny and huge regions, background labeling, temporal interpolation, projection-specific outputs, split-and-merge logic, NaN and Inf propagation, RNG restoration, zero-variance cases, and parallel execution. Also made input validation stricter for scalar, integer, finite, logical, and ordered parameter constraints. Resilience was improved in temporal binning, prewhitening, and resume or version-check processes. Optimized several key areas by vectorizing correlation-related tasks, rewriting region growth and fragment pooling, simplifying refinement steps, and batching tiny-region merges to reduce unnecessary iterations.

VERSION 1.1: Added options to correct missing data (_fillmissing_), resume a previous segmentation session (_resumeCalciSeg_), and to pool fragmented granules that presumably belong to the same biological structure (_pool_fragments_). Corrected bug in CalciSeg3D.

VERSION 1.0: The initial, NeuroImage version of [CalciSeg](https://doi.org/10.1016/j.neuroimage.2024.120758).


Cite code as
------------
If you use **CalciSeg**, please cite [Günzel et al. (2024)](https://doi.org/10.1016/j.neuroimage.2024.120758)

	@article{gunzel2024calciseg,
	  title={CalciSeg: A versatile approach for unsupervised segmentation of calcium imaging data},
	  author={G{\"u}nzel, Yannick and Couzin-Fuchs, Einat and Paoli, Marco},
	  journal={NeuroImage},
	  volume={298},
	  pages={120758},
	  year={2024},
	  publisher={Elsevier}
	  url = {https://doi.org/10.1016/j.neuroimage.2024.120758}}


Funding
-------
This work was completed with the support of the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) under Germany’s Excellence Strategy – EXC 2117 – 422037984.


Thank you for using **CalciSeg**!
