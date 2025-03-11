# Spatially Mapping Locoregional Cancer Recurrence 

A Python-based framework for patient data classification and treatment planning analysis using Scikit-RT.

## Overview
This project provides a set of tools for processing, classifying, and analysing patient imaging and treatment data using Scikit-RT. It includes modules for:
- **Patient data import and management**: Loads patient imaging and associated metadata.
- **Image registration**: Aligns medical images from different scans to a common reference frame.
- **Structure transformation**: Adjusts structures to match registered images.
- **Dose prescription analysis**: Compares treatment plans and radiation doses to expected values.
- **Classification of recurrence structures**: Identifies recurrence patterns based on dose response.
- **Visualisation of treatment planning data**: Generates plots and overlays structures onto images.
  
![FullWorkflow](https://github.com/user-attachments/assets/f072ddcd-9b6e-445f-8b44-d8113377a593)

## Features
- Uses Python package **scikit-rt** for performing image registration of planning CT and relapse CT scans for radiotherapy patients. 
- Determines the optimal spatial mapping transform, and stores the results of the registration in a results directory defined by the user. 
- Spatially maps structures outlined by clinicians from the relapse CT to the planning CT for comparison with radiotherapy dose field data.
- Applies uncertainty quantification methods to determine error in spatial mapping transform and applies to the centroid of the mapped recurrence volume.
- Generates synthetic spherical structure around the recurrence centroid with radius of the quantified uncertainty for comparison to planned radiotherapy dose field. 
- Classifies recurrence based on spatial location and dose delivered to the recurrence sphere. 
- Outputs patient classification data for further analysis. 
- Provides visualisation tools for treatment planning evaluation and comparison with other treatment covariates. 

## Installation
### Prerequisites
Ensure you have the following installed:
- Python 3.10.13
- scikit-rt 0.7.2
 
For further details on the installation of scikit-rt please reference the scikit-rt repository and documentation. 

## Usage
### Running the Classification Pipeline
1. Define the trial and dose prescription parameters.
2. Run the script to process patient data; example code for the pipeline can be found in the classification_workflow.ipynb
3. The results, including classification outputs and visualisation images, will be automatically saved in the `results` directory.

### Code Structure
- **import_analysis.py**: for additional analysis on the import breast cancer trial, additional functionality was incorporated which can be found in this python script. 
- **classification_classes.py**: Defines classification classes and methods.
- **results/**: Directory tores output data.

### Detailed Steps Performed by the Code
1. **Load Patient Data**:
   - Reads patient images and structure sets.
   - Determines the appropriate dose prescription based on trial arms.
     
2. **Register Images**:
   - Aligns the patient’s imaging data for consistency.
   - Uses `Registration_Step` to match fixed and moving scans.
  ![HNC_Registration_Performance](https://github.com/user-attachments/assets/338827d3-6039-4455-8ddc-d0f42027ab8c)

3. **Transform and Process Structures**:
   - Identifies and isolates the **Clinical Target Volume (CTV)**.
   - Transforms recurrence structures based on registration parameters.
   - Creates synthetic structures (e.g., **relapse sphere**) for analysis.
  ![Figure5_PatientYComparisonStructures](https://github.com/user-attachments/assets/c8374513-cee0-44ce-86ae-aa99f972a128)

4. **Classify Treatment Plans**:
   - Determines the **treatment arm** for each patient.
   - Classifies recurrence structures based on dose exposure.
   - Stores results in CSV format for further analysis.
  
    ![redone_classification_schemeHNC](https://github.com/user-attachments/assets/7e766a43-5d4f-4e78-bcbf-23e048e5fb13)

5. **Visualisation**:
   - Generates visual representations of dose fields and structures.
   - Uses `BetterViewer` from the scikit-rt python package to overlay relevant structures onto medical images.

## Configuration
Modify the following parameters in the jupyter notebook example to customise the analysis for the relevant trial data. 

## Output
The pipeline generates:
- **Classifications CSV files** per patient and for all patients.
- **Transformed imaging structures**.
- **Plots and visualisations** for treatment evaluation.

## Licence
This project is licensed under the MIT Licence.


## Future Development

Plans for future development include applications to cancer sites outside of breast cancer and head-and-neck cancer trials, refining image registration techniques, and expanding compatibility for alternative imaging modalities such as MRI.

---
This pipeline is designed to streamline and standardise radiomics analysis, from data import through image registration and recurrence tumour volume classification, in a reproducible and scalable framework for clinical research applications.

This work was supported by Cancer Research UK RadNet Cambridge [C17918/A28870].
![alt text](https://github.com/CeilidhWelsh/radiomics_analysis/blob/main/Radnet%20Cambridge%20%20logo.jpg)
