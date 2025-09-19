# Changelog

## v0.2.0

## Added
- `CalibrationEditorUI`: user interface for creating, saving, loading, and applying calibrations based on standard data to unknowns.
    - calibrations are saved as json files and include metadata regarding the standards used and standard measurement filter protocols, along with other calibration settings
- `CalibrationApplyUI`: user interface for selecting and applying a calibration file to measurements.

## v0.1.1

### Changed
- Analyses are now "date" + "reading no" due to non-unique "reading no" values

## v0.1.0

### Added
- `StandardUI`: user interface for adding/updating standard measurements in SQLite database