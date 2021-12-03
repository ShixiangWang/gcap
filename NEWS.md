# gcap 0.4.0

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
