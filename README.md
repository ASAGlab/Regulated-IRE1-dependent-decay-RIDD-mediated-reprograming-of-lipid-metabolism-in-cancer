# Regulated IRE1α-dependent decay (RIDD) mediated reprograming of lipid metabolism in cancer

This repository contains materials and methods for the paper "Regulated IRE1α-dependent decay (RIDD)-mediated reprograming of lipid metabolism in cancer" by Almanza*, Mnich* et al., (2021). All data supporting the findings of this study are available in this repository.

<!-- ABOUT THE PROJECT -->
## About The Project

In this study, we employed unbiased lipidomic and transcriptomic approaches to investigate the contribution of IRE1α signalling to lipid metabolism in TNBC. This system biology framework enables the systematic exploration of regulation of metabolic pathways by IRE1α and represents a novel resource for the study of non-canonical IRE1α functions.

<!-- TABLE OF CONTENTS -->
<details open="open">
  <summary>Table of Contents</summary>
  <ol>
    <li><a href="#Scripts">Scripts</a></li>
    <li><a href="#Data">Data</a></li>
    <li><a href="#Figures">Figures</a></li>
    <li><a href="#Tables">Tables</a></li>
    <li><a href="#contact">Contact</a></li>
  </ol>
</details>

<!-- SCRIPTS -->
## Scripts

* edgeR_final.R: : performs the differential gene expression (DEG) analysis of the transcriptomic data in the paper. It generates Fig.2b,c & S3b,c,d,e figures

* Metasyx_final.R : performs analysis of Metasyx metabolomics data. It generates Figure 1a-m, Figure S1b-i 

* Beatson_final.R : performs analysis of Beatson lipidomics TAG data. It generates Figure 4a-d


<!-- DATA -->
## Data

The following datasets and files are included in the repository:

* Metasyx_MD.txt contains metabolomics measurement day normalized intensities of detected features (annotated and non-annotated). 30 samples that were exposed to the  IRE1α inhibitor or DMSO. Samples were collected at 24, 48 and 72h after the treatment. Five replicates were measured per each condition.

* Metasyx_Sample_Info.txt contains information about group structure and measurement order.

* raw_counts.txt Raw RNAseq counts from featureCounts. 12 samples that were exposed to the IRE1α inhibitor (MKC8866) or DMSO. Samples were collected at 8 and 24h after the treatment. Three replicates measured per each condition.

* sample_info.txt contains sample information for the RNAseq experiment.

* Beatson_MS.txt contains lipidomics measurements. 12 samples divided into 4 biological groups: DMSO, MKC8866, PF-06424439 (DGAT2 inhibitor) and MKC8866 + PF-06424439 cotreatment. Samples were collected 72h after the treatment. Three replicates per each condition.

* Beatson_sample_info.txt contains sample information for the Beatson lipidomics experiment


<!-- FIGURES -->
## Figures

The scripts and data in this repository can generate the following figures from the main text and supplement:

* Main text: Figure 1a-m, Figure 2b-c, Figure 4a-d

* Supplementary information: Figure S1b-i, Figure S2a-c, Figure S3 b-e


<!-- TABLES -->
## Tables

* Differential Gene Expression analysis. For both time points 8 and 24 hours.

* Functional analysis results: with Bioinfominer, GO BP colection. For both time points 8 and 24 hours.

<!-- CONTACT -->
## Contact




