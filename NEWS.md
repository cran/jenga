# jenga 2.0.0

* Added approximate nearest-neighbor search via custom projection candidates.
* Added optional torch/CUDA distance acceleration with automatic CPU fallback.
* Added probabilistic forecast distributions through configurable simulation draws.
* Added split-conformal prediction intervals calibrated on validation windows.
* Added multiscale sequence representations for neighbor matching.
* Replaced philentropy binary-distance calls with custom vectorized distance metrics and removed the dependency.
* Replaced utility dependencies for iteration, dummy encoding, imputation, smoothing, timing, descriptive moments, entropy, numeric error metrics, and KNN distances with local implementations.
* Expanded distance methods and kernels with additional dependency-free implementations.
* Added sign-respecting filtering on final reconstructed draws before forecast quantiles are computed.
* Enforced historical sign bounds after conformal calibration and in forecast plots.

# jenga 1.0.0

* Added a `NEWS.md` file to track changes to the package.

# jenga 1.1.0

* common sense filter
* added stats to the prediction table (and removed pred_stats)
* added new criterion for model selection: prediction score
* added new metric measures for test error from greybox
* added automatic differentiation module (and removed deriv from hyper-parameters)

# jenga 1.2.0

* prediction score is now managed as a distinct metric and it is not used for model selection
* streamlined code and processing
* matrix distance calculation: improved performance shifting from philentropy to Rfast
* other minutiae

# jenga 1.3.0

* Rationalization and streamlining
* Now you can deal also with categorical variables
