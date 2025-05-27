# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [1.2.2] - 2025-05-27

### Added
- Support for running deterministic simulations using a `seed_file` via the `SigProfilerSimulator(seed_file=...)` parameter.
- Pytest-based test suite to validate reproducibility of simulations with and without seed files.
- `pyproject.toml` added to support modern build tools and packaging standards.

## [1.2.1] - 2025-04-03

### Added
- Support for `mm39` genome assembly.

## [1.2.0] - 2025-02-24

### Changed
- Updated dependencies: Now requires **Pandas >= 2.0.0**, **NumPy >= 2.0.0**, and **Python >= 3.9**.
- Dropped support for **Python 3.8**