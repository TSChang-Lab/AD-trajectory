# 🧠 Identifying Common Disease Trajectories of Alzheimer’s Disease with Electronic Health Records

This repository contains code and documentation for the study:

"Identifying common disease trajectories of Alzheimer’s disease with electronic health records"
by Mingzhou Fu, MPH; Sriram Sankararaman, PhD; Bogdan Pasaniuc, PhD; Keith Vossel, MD, MSc; Timothy S. Chang, MD, PhD.
Published in eBioMedicine, 2025. (https://www.sciencedirect.com/science/article/pii/S2352396425002750#sec2)

📄 Abstract

Background:
Alzheimer’s disease (AD) is a growing public health concern and a major cause of dementia. While many studies have identified AD risk factors, few examine the sequential comorbid trajectories leading to AD.

Methods:
Using longitudinal data from over 24,000 patients in the University of California Health Data Warehouse (UCHDW), we:
* Modeled diagnosis timing using a Fine-Gray subdistribution hazard model
* Identified sequential multi-step disease trajectories
* Clustered trajectories using dynamic time warping (DTW) and k-means
* Characterized trajectory networks
* Performed causal discovery using the Greedy Equivalence Search (GES) algorithm

Results were validated in the All of Us Research Program, a nationally representative, diverse cohort.

Findings:
We identified 6,794 unique AD trajectories across 5,762 patients, grouped into four major trajectory clusters:
* Mental health
* Encephalopathy
* Mild cognitive impairment
* Vascular disease

Over a quarter of diagnostic transitions had consistent temporal ordering (e.g., hypertension → depressive episode → AD). These trajectories predicted AD risk better than single diagnoses and were reproducible across datasets.

🧪 Data Availability

Due to privacy regulations, no patient-level data is included in this repository.
The code operates on de-identified data derived from the UCHDW and All of Us cohorts. Researchers with appropriate access may adapt the pipeline for their own data.

📊 Methods Summary
Survival modeling: Fine-Gray subdistribution hazard model

Trajectory construction: Sequential filtering of significant associations

Trajectory similarity: Dynamic Time Warping (DTW)

Clustering: k-means on DTW-aligned distances

Network analysis: Edge frequency and transition directionality

Causal discovery: Greedy Equivalence Search (GES)

🔁 Validation
Validation was performed in:

The UCHDW cohort via association testing and control comparison

The All of Us Research Program to confirm reproducibility across a diverse population

✍️ Citation
If you use this code or framework, please cite:

Fu M, Sankararaman S, Pasaniuc B, Vossel K, Chang TS.
Identifying common disease trajectories of Alzheimer’s disease with electronic health records.
eBioMedicine. 2025. doi: [DOI]

📬 Contact
For questions or collaboration inquiries, please contact:
📧 Timothy S. Chang – timothychang@mednet.ucla.edu

🧬 Funding
This study was supported by:
* National Institutes of Health (NIH)
* National Institute on Aging (NIA)
* National Science Foundation (NSF)
* Hillblom Foundation
* Fineberg Foundation
* California Department of Public Health
