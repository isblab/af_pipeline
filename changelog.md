# Changelog

All notable changes to this project will be documented here.

## [1.1.0] - 2026/04/20

### Major changes:

- The concept of job_cycle is artificial and not necessary, hence removed.
- Ranking AlphaFold2 predictions is now supported (from ranking_debug.json).
- Ranking ColabFold predictions is now supported (from file paths).

### Minor changes + fixes:

- Mapping of entities to chains is added to best_af_predictions.json. This can be used subsequently in `save_ppair_interaction` in Interaction class and `save_rigid_bodies`, `assess_rigid_bodies` in RigidBodies class.
- Included ColabFold example structure predictions.
- Example structure prediction output is provided as zippped file.
- Updated tests and test data is provided as zipped file.

## [1.0.1] - 2026/04/09

### Added

- Changelog, contributing guidelines are added.

### Updated

- Added usage on home page.
- Network diagram is only shown for the entire AF-Pipeline at once and not separately for the submodules.

### Fixed

- Bug fixes in `rigid_body_assessment.py`.

## [1.0.0] - 2026/04/07

### Added

- First [release](https://github.com/isblab/af_pipeline/releases/tag/1.0.0) of AF-Pipeline.