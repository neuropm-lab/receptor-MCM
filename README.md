# receptor-MCM

# Receptor-enriched Multi-factorial Causal Model (re-MCM)
Scripts for spatiotemporal causal model of multi-modal neuroimaging variables, based on Iturria-Medina et al., 2017 (https://doi.org/10.1016/j.neuroimage.2017.02.05), and extended with molecular template data (gene expression and neurotransmitter receptor maps) in Adewale et al., 2021 (https://doi.org/10.7554/eLife.62589),
Khan et al., 2022 (https://doi.org/10.1093/brain/awab375), and Khan et al., 2023 (https://doi.org/10.1038/s41467-023-41677-w).

## Workflow
Receptor enriched-MCM can be used with multi-modal (MRI, PET, SPECT, etc.) and longitudinal data.

### Step 0: Preprocessing and harmonization
Preprocess raw data as appropriate using tools of the user's choice:
- Motion correction, slice timing correction, artefact removal, etc.
- Alignment to T1 image
- Spatial normalization to MNI template
- Harmonization for site and scanner variations (e.g., using ComBat, https://github.com/Jfortin1/ComBatHarmonization)

### Step 1: Compile the MCM data structure using NeuroPM-box
Organize the harmonized NIFTI files for all subjects and neuroimaging modalities as described in the NeuroPM-box instructions.

NeuroPM-box can be downloaded from 
https://www.neuropm-lab.com/neuropm-box-download.html

### Step 2: Load MCM data structure, molecular templates, and anatomical connectivity matrix

### Step 3: Fit receptor-MCM models

### Step 4: Downstream analysis
Use molecular-MCM outputs (model parameters) for inter-subject comparisons (e.g., with clinical variables).

## Citation

Please cite NeuroPM-box and the relevant molecular-MCM paper from the following:
- **NeuroPM-box**: Iturria-Medina, Y., Carbonell, F., Assadi, A. et al. Integrating molecular, histopathological, neuroimaging and clinical neuroscience data with NeuroPM-box. Commun Biol 4, 614 (2021). https://doi.org/10.1038/s42003-021-02133-x
- **Molecular-MCM (neurotransmitter receptors)**: Khan, A.F., Adewale, Q., Lin, SJ. et al. Patient-specific models link neurotransmitter receptor mechanisms with motor and visuospatial axes of Parkinson’s disease. Nat Commun 14, 6009 (2023). https://doi.org/10.1038/s41467-023-41677-w
- **Molecular-MCM (gene expression)**: Adewale, Q., Khan, A.F., Carbonell, F., & Iturria-Medina, Y., Integrated transcriptomic and neuroimaging brain model decodes biological mechanisms in aging and Alzheimer’s disease. eLife 10:e62589 (2021). https://doi.org/10.7554/eLife.62589 
