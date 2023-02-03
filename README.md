# Project Description

This project attempts to replicate the findings from the paper Baron, Maayan, Adrian Veres, Samuel L. Wolock, Aubrey L. Faust, Renaud Gaujoux, Amedeo Vetere, Jennifer Hyoje Ryu, et al. 2016. “A Single-Cell Transcriptomic Map of the Human and Mouse Pancreas Reveals Inter- and Intra-Cell Population Structure.” Cell Systems 3 (4): 346–60.e4. PMID: 27667365

*Final Report: `BF528_Project4_Tinman (1).pdf`*

# Contributors
  
1.  Data Curator: Rizky Kafrawi -@rkafrawi
2.  Programmer: Allison Choy - @AllisonWChoy
3.  Analyst: Reshma Ramaiah - @ReshmaRamaiah10
4.  Biologist: Aneeq Husain - @aneeqh

# Repository Contents

## Data Curator
This directory is split into three directories, two of which represents one of the two primary components of data preparation for project 4: 
barcode whitelisting and salmon. The last directory contains the primary deliverables for each run (mapping summary and UMI count matrix).

## Programmer
- ``Prog2.R`` - R script that creates :
  - a Seurat object with a UMI matrix 
  - filters out low quality cells, normalize and scales data for PCA analysis and supporting plots
  - creates a UMAP with dimensionality based on PCA analysis

## Analyst

## Biologist

