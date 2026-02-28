# AI4DITE_GPR_data
Low-Frequency Ground Penetrating Radar (GPR) Dataset Acquired via Ground-Based and UAV Platforms

**Repository Description**

This GitHub repository provides the official Python toolkit accompanying the *Low-Frequency GPR Ground/UAV Dataset* developed by **INGV** and **UNISA**. The repository ensures reproducibility of the dataset organization and supports rigorous comparative analysis across the two acquisition modalities (ground-based cart and UAV-mounted GPR).

The dataset associated with this repository is publicly available on Zenodo at:
https://doi.org/10.5281/zenodo.18769571

The codebase is designed as a **modular and extensible framework**, enabling the complete data handling workflow—from raw file ingestion to cross-modality alignment and export toward other scientific environments. While currently tailored to ground and UAV GPR acquisitions, the architecture has been intentionally designed to support **future extensions to additional sensing modalities and acquisition platforms**, such as vehicle-mounted systems, robotic platforms, or alternative airborne configurations.

### Main Features

* **Data Reading**

  * Native support for both **GeoTable** and **SEG-Y (SEGY)** formats
  * Automatic parsing of metadata and acquisition parameters
  * Consistent internal data structures across formats

* **Data Organization**

  * Unified representation of ground-based and UAV-based acquisitions
  * Trace indexing and GPS-based spatial alignment
  * Metadata harmonization for cross-site and cross-platform experiments

* **Visualization Utilities**

  * Radargram plotting
  * Trace-level inspection
  * Comparative visualization between different acquisition modalities
  * Quick-look tools for quality assessment and artifact inspection

* **Cross-Modality Alignment**

  * Spatial alignment of traces using GPS coordinates
  * Resampling and synchronization utilities
  * Support for domain adaptation, signal comparison, and data fusion studies

* **MATLAB Export**

  * Conversion of processed datasets into `.mat` format
  * Structured export of traces, coordinates, and metadata
  * Compatibility with MATLAB-based signal processing workflows

Thanks to its modular design, the repository can be readily extended to incorporate **new acquisition geometries, additional sensor metadata, or alternative subsurface sensing modalities**, making it suitable as a general-purpose framework for multi-platform GPR and related geophysical datasets.

The code is documented and structured to serve both as a ready-to-use toolkit and as a foundation for further methodological and experimental developments.

