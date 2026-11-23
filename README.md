![CalciSeg logo][logo]

[logo]: CalciSeg_logo.png

CalciSeg - A versatile approach for unsupervised segmentation of calcium imaging data
=====================================================================================


Contents
--------

1. [About](#about)
2. [Requirements](#requirements)
3. [Running CalciSeg](#running-calciseg)
4. [Choosing a method](#choosing-a-method)
5. [Output](#output)
6. [Preconditions and caveats](#preconditions-and-caveats)
7. [Recent changes](#recent-changes)
8. [Cite code as](#cite-code-as)
9. [Funding](#funding)


About
-----
**CalciSeg** is a MATLAB-based method for unsupervised segmentation of calcium imaging data into spatial regions ("granules") using projection-based initialization, temporal similarity, region-size constraints, and optional iterative refinement. It is designed for 2-D or 3-D calcium imaging stacks with spatial dimensions _(x, y)_ and an optional temporal dimension _(time)_.

The current **CalciSeg** function supports multiple projection methods, several initialization strategies, optional refinement, resumable segmentation, fragmented-region pooling, temporal downsampling and detrending during segmentation, anisotropic pixel geometry, headless/batch operation, and per-granule quality metrics.


Requirements
------------
To run **CalciSeg**, install MATLAB together with:

- **Image Processing Toolbox** (required)
- **Statistics and Machine Learning Toolbox** (required)
- **Parallel Computing Toolbox** (optional)

The two required toolboxes are checked at startup, and **CalciSeg** fails immediately with an actionable message if either is missing, rather than erroring deep inside a long run.

The Parallel Computing Toolbox is optional. **CalciSeg** uses _parfor_ for the split check and the refinement loop; without the toolbox these run serially with identical results. Parallel execution is also gated on workload size, so small problems run serially regardless. Set _use_parallel_ to _false_ to disable it entirely. Starting a parallel pool before calling the function is recommended when working with large stacks.


Running CalciSeg
----------------
Call the function as:

```matlab
[granules_labeled, summary_stats] = CalciSeg(stack, Name, Value);
```

### Required input

- _stack_ (matrix): calcium imaging dataset.
  - 3-D input is expected as _(rows, columns, time)_.
  - 2-D input is also accepted. In that case, **CalciSeg** automatically adjusts settings to:
    - _projection_method = 'none'_
    - _init_seg_method = 'voronoi'_
    - _refinement_method = 'rmse'_
    - _temporal_bin = 1_
    - _detrend = 'none'_

Inputs with more than three dimensions are squeezed to collapse singleton dimensions; if more than three non-singleton dimensions remain, the call errors.

### Optional input parameters

CalciSeg warns about unused method-specific options after resolving the input dimensions and initialization mode. For example, `n_PCA` applies to PCA projection and to automatic rICA feature selection.

#### Projection and preprocessing

- **_projection_method_** (char)  
  Method for calculating the projection across time.  
  **Default:** _'std'_

  - _'std'_    : standard deviation projection
  - _'mean'_   : mean intensity projection
  - _'median'_ : median intensity projection
  - _'max'_    : maximum intensity projection
  - _'min'_    : minimum intensity projection
  - _'pca'_    : principal component projection; the stack is projected into PC space and a variance-weighted sum of the retained components forms the projection
  - _'corr'_   : Pearson correlation between each pixel’s trace and the mean trace of its in-bounds spatial neighbors
  - _'none'_   : use the first slice of the stack, or use when no temporal dimension is present

- **_detrend_** (char)  
  Temporal preprocessing for **segmentation only**. Summary traces retain the original time resolution and fluorescence scale.  
  **Default:** _'none'_

  - _'none'_   : no detrending
  - _'linear'_ : remove a per-pixel linear trend
  - _'dff'_    : compute `(F − F0) / F0` using each pixel’s 10th percentile as its baseline. Nonpositive baselines are set to `eps('single')`. Pixels above half the lowest of 20 Otsu thresholds on the baseline image form a mask. Holes are filled and traces outside the mask are set to zero.

- **_spatial_smooth_** (scalar number)  
  Gaussian sigma, in pixels, applied to the projection before extrema detection. Larger values suppress spurious local maxima and reduce over-segmentation.  
  **Default:** _0_ (no smoothing)

- **_temporal_bin_** (integer)  
  Temporal downsampling factor for faster segmentation. The stack is binned by averaging complete groups of _N_ consecutive frames. A final partial group is omitted from segmentation. Final summary statistics are computed from the full-resolution time courses. Binning is disabled automatically, with a warning, if it would leave fewer than two frames.  
  **Default:** _1_

- **_fillmissing_** (logical)  
  Fill missing samples along rows, then columns, then time, using linear filling followed by nearest-value filling on each axis. The call errors if any nonfinite samples remain.  
  **Default:** _false_

- **_n_PCA_** (scalar number)  
  Percentage of variance to retain for PCA projection. Automatic rICA also uses this target to select its feature count.  
  Must be _> 0_ and _<= 100_.  
  **Default:** _100_

#### Initial segmentation

- **_init_seg_method_** (char)  
  Method for the initial segmentation.  
  **Default:** _'voronoi'_

  - _'voronoi'_ : Voronoi-style partitioning based on local extrema
  - _'corr'_    : rank unassigned pixels by local correlation and start at the highest eligible seed. Grow through edge-sharing neighbors whose traces correlate with the current region mean at or above the growth threshold
  - _'rICA'_    : center the data, run PCA, whiten the retained directions, and fit reconstruction ICA. Assign each pixel to its strongest sign-corrected, standardized component loading

  Fresh correlation-based initialization uses `projection_method='corr'`. The seed score uses a 2-D disk excluding the center pixel. At the default radius of one pixel, an interior pixel has four neighbors. A larger disk changes the seed score neighborhood, while growth always uses four-connected neighbors. Time supplies the samples for each correlation.

- **_regmax_method_** (char)  
  Method for identifying local extrema in the projection.  
  **Default:** _'raw'_

  - _'raw'_      : apply _imregionalmax_ / _imregionalmin_ directly
  - _'filtered'_ : in Voronoi mode, dilate before detecting maxima and erode before detecting minima. In correlation mode, box-filter the seed map. In rICA mode, box-filter each component loading map. Both box filters use width `2*filter_size+1` with replicated edges
  - _'both'_     : combine raw and filtered extrema detection

  **Note:** _regmax_method = 'both'_ is only applicable with _init_seg_method = 'voronoi'_. Otherwise it is reset to _'raw'_ with a warning.

- **_seed_on_** (char)  
  Which extrema seed the Voronoi tessellation.  
  **Default:** _'both'_

  - _'both'_ : seed on both maxima and minima
  - _'max'_  : seed on maxima only
  - _'min'_  : seed on minima only

- **_keep_border_** (logical)  
  When _false_, the outermost ring of pixels is forced to act as seeds and the resulting edge-seeded granules are then discarded, so border pixels are assigned to the background (ID 0). Set to _true_ to use extrema seeds alone and retain their cells, including cells touching the edge.  
  **Default:** _false_

- **_pixel_size_** (2-element positive vector)  
  Pixel size _[dy dx]_ used to scale coordinates for anisotropic Voronoi tessellation.  
  **Default:** _[1 1]_

- **_seed_thresh_** / **_grow_thresh_** (scalar number in _[0, 1]_, or _[]_)  
  Separate correlation thresholds controlling where a region may start versus how it grows during correlation-based growing. When empty, both fall back to _corr_thresh_.  
  **Default:** _[]_

- **_n_rICA_** (integer)  
  Number of features to extract during reconstruction ICA.  
  **Default:** _0_

  - _0_ chooses a feature count using _n_PCA_ and the number of PCA directions with numerically positive variance
  - values _> 0_ allow undercomplete or overcomplete feature representations

- **_rica_sd_threshold_** (nonnegative scalar or `[]`)  
  The default `[]` orients each component by spatial skewness. A numeric value compares the mean positive and negative excursions beyond that many standard deviations. A tie preserves the fitted sign.

- **_rICA_thresh_** (integer from 0 to 20)  
  The default `0` keeps all winning loadings. A positive value requests that many Otsu thresholds on the 3×3 median-filtered winning-load image. Pixels below the lowest threshold become background.

- **_seed_** (integer)  
  Twister seed for rICA, from 0 through `2^32-1`. The caller's RNG state is saved and restored around the call, so the global stream is left unchanged.  
  **Default:** _1_

#### Refinement, size constraints, and pooling

- **_n_rep_** (integer)  
  Number of refinement iterations. Refinement stops early if the label map converges or begins oscillating between two states.  
  **Default:** _0_

- **_refinement_method_** (char)  
  Distance measure used during refinement.  
  **Default:** _'rmse'_

  - _'rmse'_   : root **median** square error, _sqrt(median((a-b).^2))_ — median-based despite the name
  - _'rmedse'_ : alias for _'rmse'_, accepted as a clearer name and normalized internally
  - _'corr'_   : Pearson correlation coefficient

- **_limitPixCount_** (numeric vector or char)  
  Minimum and maximum pixel area per granule.  
  **Default:** _[1 inf]_

  - _[a, b]_ : lower and upper pixel-count limits
  - _'auto'_ : estimate limits from the initial segmentation using the quantiles given by _auto_quantiles_

  A region smaller than _a_ merges into a neighboring granule or becomes background if isolated. The minimum is enforced again after the final connectivity check. A region larger than _b_ is a candidate for splitting.

  **Note:** The minimum size also affects the filter size used by _regmax_method_. Under _'auto'_ that filter size is fixed at 1.

  **Note:** Maximum-size enforcement is best-effort. Splits can fail, and later merges can produce oversized granules. All final violations appear in _summary_stats.oversized_ids_.

- **_auto_quantiles_** (2-element vector)  
  Quantiles _[lo hi]_ used when _limitPixCount = 'auto'_, with _0 <= lo < hi <= 1_.  
  **Default:** _[0.05 0.95]_

- **_corr_thresh_** (scalar number)  
  Default seed and growth threshold for correlation initialization (both comparisons include equality). Fragment pooling requires correlation strictly above this value. Correlation refinement ranks neighboring granules by fit without a cutoff.  
  **Default:** _0.85_

- **_pool_fragments_** (logical)  
  Pool neighboring granules that likely belong to the same biological structure if their average time courses correlate above _corr_thresh_.  
  The neighbor comparison is repeated until no adjacent pair exceeds the threshold.  
  **Default:** _false_

  **Note:** This step is applied after refinement and therefore ignores the upper pixel-count limit.

- **_resume_CalciSeg_** (2-D matrix)  
  Label map from a previous segmentation session, matching the spatial size of _stack_ and containing only finite, non-negative integers.  
  A nonempty map skips initialization. Preprocessing and size enforcement still run. Labels are compacted before processing to preserve distinct large integer IDs. Passing `[]` starts fresh.  
  **Default:** _[]_

#### Execution control

- **_use_parallel_** (logical)  
  Allow _parfor_ for the split check and refinement. When _false_, or when the workload is too small to outweigh parallel overhead, these run serially with identical results.  
  **Default:** _true_

- **_verbose_** (logical)  
  Show progress bars and print status messages. Set to _false_ to suppress progress output. Warnings remain enabled. Progress bars are additionally suppressed automatically when no graphical desktop is available.  
  **Default:** _true_

- **_check_update_** (logical)  
  When _true_, check GitHub once per MATLAB session for a newer **CalciSeg** version and warn if one is available. Requires an internet connection and fails silently when offline. No network call is made when _false_.  
  **Default:** _false_

### Examples

Default Voronoi segmentation with summary statistics:

```matlab
[granules_labeled, summary_stats] = CalciSeg(stack);
```

Correlation-based growing with baseline normalization and progress output disabled:

```matlab
granules_labeled = CalciSeg(stack, ...
    'init_seg_method', 'corr', ...
    'detrend', 'dff', ...
    'corr_thresh', 0.8, ...
    'verbose', false);
```

Fewer, larger granules via smoothing and a size floor:

```matlab
granules_labeled = CalciSeg(stack, ...
    'spatial_smooth', 2, ...
    'limitPixCount', [25 inf], ...
    'n_rep', 5);
```

Voronoi initialization with filtering, size limits, and temporal binning:

```matlab
[granules_labeled, summary_stats] = CalciSeg(...
    stack, ...
    'projection_method', 'std', ...
    'init_seg_method', 'voronoi', ...
    'regmax_method', 'filtered', ...
    'n_rep', 25, ...
    'refinement_method', 'rmse', ...
    'limitPixCount', [5 100], ...
    'fillmissing', true, ...
    'pool_fragments', false, ...
    'temporal_bin', 2);
```

Automatic rICA feature selection and optional background masking:

```matlab
[L, stats] = CalciSeg(stack, 'init_seg_method', 'rICA', ...
    'n_rICA', 0, 'n_PCA', 95, 'rICA_thresh', 1);
nFeatures = stats.parameters.n_rICA;
```

### Note on CalciSeg_3D

**CalciSeg_3D** is a separate work in progress for volumetric data.


Choosing a method
-----------------
A rough guide to _init_seg_method_:

| Data characteristics | Suggested method |
| --- | --- |
| Sparse, blob-like sources | _'corr'_ (region growing) |
| Dense fields, tiling | _'voronoi'_, optionally with _spatial_smooth_ |
| Overlapping temporal sources | _'rICA'_ |

The _'voronoi'_ default over-segments by design. Increase _spatial_smooth_ or raise the minimum in _limitPixCount_ to obtain fewer, larger granules. When _verbose_ is enabled, **CalciSeg** warns if the result exceeds roughly one granule per four pixels, which usually indicates over-segmentation.


Output
------
**CalciSeg** returns two outputs:

- **_granules_labeled_**  
  2-D label map of segmented granules. Its size matches the first two dimensions of _stack_. Background is 0. Every granule is a single 4-connected component, and IDs are contiguous integers _1..N_ ordered by descending pixel count, so ID 1 is the largest granule.

- **_summary_stats_**  
  Structure containing per-granule summary statistics and segmentation metadata. Computed only when the second output is requested.

Per-granule time courses and their statistics use full-resolution fluorescence, including any imputed samples. Temporal binning, detrending, and PCA affect segmentation. Each intensity summary image assigns one statistic of the granule’s mean trace to all its pixels.

### Fields in summary_stats

#### Time courses and projections

- _projection_  
  Preprocessed projection, potentially binned, detrended, or smoothed. Voronoi and correlation initialization use its values. rICA uses its spatial dimensions. A resumed run also returns the projection.

- _avgTCs_  
  Mean time course of each granule (_nGranules_ x _nFrames_)

- _granule_Corr_  
  Mean Pearson correlation between individual pixel traces and the granule’s mean trace, excluding undefined correlations.

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

#### Geometry

- _area_  
  Per-granule pixel count

- _centroid_  
  Per-granule centroid as _[x y]_

- _bbox_  
  Per-granule bounding box as _[x y width height]_

- _pixel_idx_  
  Cell array of each granule's linear pixel indices

#### Quality and diagnostics

- _quality.SNR_  
  Per-granule signal-to-noise ratio: peak above a 10th-percentile baseline, divided by a robust noise estimate derived from successive differences of the mean trace

- _quality.compactness_  
  Per-granule compactness  
  (_4*pi*area / perimeter^2_, where _1_ corresponds to a perfect circle; values may slightly exceed 1 for small rasterised shapes because of perimeter discretisation)

- _active_region.map_  
  Binary map of active regions, obtained by Otsu thresholding over foreground pixels only, so the split is not biased by how much background the field of view contains

- _active_region.method_  
  Metric used for the activity mask. A PCA projection uses `'std'` here. With `'none'`, the mask is empty.

- _oversized_ids_  
  IDs still exceeding the maximum size after best-effort splitting. Empty when none remain, or when no finite maximum was requested.

- _parameters_  
  Struct of all options actually used for the run, including auto-resolved values such as _limitPixCount_ and the derived filter size. For rICA, `n_rICA_requested` records the input and `n_rICA` the fitted count. When refinement runs, `refine_iterations` records completed passes.


Preconditions and caveats
-------------------------
These are the caller's responsibility; **CalciSeg** cannot detect violations.

- **Motion correction.** The stack should be motion-corrected before segmentation. Residual motion inflates variance and correlations and corrupts the result.
- **Drift and photobleaching.** Correlation-based steps — _'corr'_ projection and initialization, _'corr'_ refinement, and _pool_fragments_ — assume drift has been removed. Use `detrend='linear'` for a linear trend. ΔF/F rescales each trace around a fixed baseline.
- **Correlation thresholds under PCA.** When _projection_method = 'pca'_, refinement, size-based merging, and fragment pooling compare PC scores. Fresh correlation growing selects the correlation projection automatically. Thresholds and RMSE values calibrated on raw fluorescence do not transfer directly to that space.
- **Trace length.** Aggressive _temporal_bin_ or restrictive _n_PCA_ settings can leave very few frames, at which point Pearson correlations become unstable and _corr_thresh_ loses meaning.
- **Memory.** Building the local correlation image holds several full-size copies of the stack simultaneously. Large stacks may need to be processed in blocks or temporally binned first.

### Reproducibility

The rICA path is seeded via the _seed_ input and restores the caller's RNG state afterwards. All other steps are deterministic for a given input, though results may vary across MATLAB and toolbox versions because of tie-breaking in _knnsearch_ and _watershed_.


Recent changes
--------------

**VERSION 1.4**

- rICA now uses PCA before whitening, supports overcomplete fits, and compares standardized loadings. The default sign rule uses skewness. Optional Otsu masking is available through `rICA_thresh`. These changes can alter existing rICA segmentations.
- ΔF/F masks low-baseline pixels. Single-frame inputs skip temporal detrending.
- Correlation calculations use double-precision accumulation. Region growth tests the full boundary before applying a size cap. Filtered correlation and rICA initialization use spatial box filters.
- Resumed maps preserve large IDs. Singleton fragments join the most common adjacent positive ID. Oblique splitting, thin-image pooling, and minimum-size enforcement after connectivity checks have been corrected.
- Input validation and warnings follow the effective settings. Summary metadata records fitted rICA feature counts and the metric used for the activity mask.

VERSION 1.3: Added temporal detrending (_detrend_), spatial smoothing of the projection (_spatial_smooth_), configurable Voronoi seeding (_seed_on_, _keep_border_, _pixel_size_), decoupled seed and growth thresholds for region growing (_seed_thresh_, _grow_thresh_), configurable rICA sign correction and seeding (_rica_sd_threshold_, _seed_), configurable automatic size quantiles (_auto_quantiles_), and explicit execution control (_verbose_, _use_parallel_, _check_update_). Summary statistics now include per-granule geometry (_area_, _centroid_, _bbox_, _pixel_idx_), oversized-granule reporting, and a full record of the parameters used. Toolbox dependencies are verified up front. Refinement now detects two-state oscillation in addition to convergence, Otsu thresholding for active regions is computed over foreground pixels only, and output IDs are ordered by descending granule size.

VERSION 1.2: Fixed many segmentation and statistics bugs, including issues with handling tiny and huge regions, background labeling, temporal interpolation, projection-specific outputs, split-and-merge logic, NaN and Inf propagation, RNG restoration, zero-variance cases, and parallel execution. Also made input validation stricter for scalar, integer, finite, logical, and ordered parameter constraints. Resilience was improved in temporal binning, prewhitening, and resume or version-check processes. Optimized several key areas by vectorizing correlation-related tasks, rewriting region growth and fragment pooling, simplifying refinement steps, and batching tiny-region merges to reduce unnecessary iterations.

VERSION 1.1: Added options to correct missing data (_fillmissing_), resume a previous segmentation session (_resume_CalciSeg_), and to pool fragmented granules that presumably belong to the same biological structure (_pool_fragments_). Corrected bug in CalciSeg3D.

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
	  publisher={Elsevier},
	  url = {https://doi.org/10.1016/j.neuroimage.2024.120758}}


Funding
-------
This work was completed with the support of the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) under Germany's Excellence Strategy – EXC 2117 – 422037984.


Thank you for using **CalciSeg**!
