# gcap 0.8.1

- Cleaned logic.

# gcap 0.8.0

- Added options `tightness` and `gap_cn`.
- Renamed `gcap-wes.R` script to `gcap-bam.R`.

# gcap 0.7.0

- Handled void result.
- Updated the background copy number reference and criterion judging a amplicon (#22).

# gcap 0.6.0

- Designed and implemented a class `fCNA` used for storing the workflow key outputs
and downstream analysis and visualization.

# gcap 0.5.0

- Provided a function `convertID()` to convert gene IDs.
- Optimized the output summary.
- Re-constructed scoring and workflow output.

# gcap 0.4.1

- Added `use_best_ntreelimit` in `gcap.runPrediction()` to control the ntree setting.
When it is `FALSE`, we use a custom processing to obtain a more conservative tree number.
- Added `deploy()` to auto-deploy the CLI to `/usr/loca/bin`.
- Added easy-to-use CLI in `inst` directory.
- Filled `NA`s to input when age and gender are not available.
- Automatically appended logs to specific directory with `rappdirs::app_dir("gcap", "ShixiangWang")`.
Users can obtain log path and cat log info with `gcap:::get_log_file()` and `gcap:::cat_log_file()`
for debugging. (#14)
- Supported XGB54 model in workflows. (#13)

# gcap 0.4.0

- After exploration, we found our stepwise model outperform MBO tuned model.
So the models for predicting circular target have been limited to 3.
- `gcap.ASCNworkflow()` now supports input with only total integer copy number,
like the result from ABSOLUTE software (also [DoAbsolute](https://github.com/ShixiangWang/DoAbsolute)).
- Added stepwised model for circle target.
- Changed the way how to select model and run prediction.
- `custom_model` in `gcap.runPrediction()` has been changed to `model`. This is
inconsistent with version below v0.4.

# gcap 0.3.3

* Updated scoring for supporting different thresholds.
* Wrapped ASCAT workflow in `tryCatch()` to avoid abnormal failure.

# gcap 0.3.2

* Fixed input `extra_info` sample (order) issues.

# gcap 0.3.1

* Fixed the issue about rendering wrong data files when skipping existing ASCAT
calling.

# gcap 0.3.0

* Updated models for prediction.

# gcap 0.2.1

* Removed `prob_cutoff` setting in workflows. Directly use prob
0.1, 0.5 and 0.9 for cutting low, medium and high quality amplicon.

# gcap 0.2.0

* Implemented an alternative workflow from allele specific copy number data
to final result files. (#5)

# gcap 0.1.0

* Implemented basic workflow from BAM files to result files.
* Added a `NEWS.md` file to track changes to the package.
