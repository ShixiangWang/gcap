# gcap 1.1.4

- Removed the XGBOOST version limits, instead, a warning is posted.
- Used GitHub action to render pkgdown site.

# gcap 1.1.3

- Updated the logic of using `only_oncogenes` for filtering.

# gcap 1.1.2

- Fixed the data loading for oncogene of mouse.

# gcap 1.1.1

- Updated Sequenza workflow.
- Fixed the re-handling of errored seqz and facets runs.

# gcap 1.1.0

- Supported analysis workflow from FACETS and Sequenza.
- Mouse genome is enabled based on the two implemented workflows above.

# gcap 1.0.0

- The first public and stable release.

# gcap 0.21.0

- Corrected blood label to somatic (https://github.com/BioinformaticsFMRP/TCGAbiolinks/issues/523).

# gcap 0.20.2

- Limited `xgboost` version lower than `1.6` as it will not keep some key info in `.rds` file.

# gcap 0.20.1

- Enhance the `getGeneSummary()` and `getCytobandSummary()` methods to return mutation matrix.

# gcap 0.20.0

- This version is not compatible with previous versions, as the analysis and
visualization functions are moved to an independent package 'gcaputils'.

# gcap 0.14

- Added some utils functions and visualization functions.

# gcap 0.13

- Updated initial setting and CLI.

# gcap 0.12

This version has been discarded from git history.

# gcap 0.11.1

- Supported a `NA` passing as tightness to remove the use of TCGA blood summary data
as a more strict threshold for circular amplicon.

# gcap 0.11.0

- Added `fCNA$subset()` method.
- Added `gcap.plotDistribution()` function.

# gcap 0.10.0

- Added `gcap.plotForest()`.
- Added `gcap.plotKMcurve()`.
- Added `gcap.plotProfile()`.
- Added method `convertGeneID()` to `fCNA` class.
- Set a default value for `pdata` option in `fCNA$new()`.

# gcap 0.9.0

- Supported `gcap` as main command, and previous two commands as subcommands if `GetoptLong` version `>=1.1.0`.
Note: not test yet.

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
