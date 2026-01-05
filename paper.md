---
title: 'IsoPat 3.0: A Python Package for Isotope Pattern Deconvolution in Mass Spectrometry'
tags:
  - Python
  - mass spectrometry
  - isotope labeling
  - deconvolution
  - chemistry
authors:
  - name: Christian C. Gruber
    orcid: 0000-0003-0214-449X
    corresponding: true
    affiliation: "1, 2"
  - name: Wolfgang Kroutil
    orcid: 0000-0002-2151-6394
    affiliation: "3, 4"
affiliations:
  - name: Innophore GmbH, Graz, Austria
    index: 1
  - name: Institute of Molecular Bioscience, University of Graz, Austria
    index: 2
  - name: Institute of Chemistry, University of Graz, Austria
    index: 3
  - name: Field of Excellence BioHealth, BioTechMed Graz, NAWI Graz, University of Graz, Austria
    index: 4
date: 5 January 2025
bibliography: paper.bib
---

# Summary

Isotope labeling is a fundamental technique in chemistry and biochemistry for investigating molecular structures, reaction mechanisms, and metabolic pathways. The accurate determination of isotope incorporation from mass spectrometry data requires deconvolution of overlapping isotope patterns—a non-trivial task when multiple labeled species are present. **IsoPat** provides a Python implementation of a least-squares deconvolution algorithm that determines the relative amounts of each isotope-labeled species from low-resolution mass spectrometry data. The package features minimal dependencies (NumPy only), a clean documented API, command-line interface for batch processing, and comprehensive tests for the core algorithm.

# Statement of need

Mass spectrometry is the primary analytical tool for characterizing isotope-labeled compounds due to its sensitivity and direct measurement of mass differences [@simon1966; @heumann1982; @johnstone1996; @dehoffmann1996]. However, quantifying the degree of labeling in complex mixtures is challenging because natural isotope distributions (e.g., $^{13}$C at 1.1% abundance) create overlapping patterns, multiple labeled derivatives (d$_0$, d$_1$, d$_2$, ..., d$_n$) contribute to the same m/z region, and standard literature procedures relying on single-peak measurements lead to significant errors [@hesse1997; @habler2024].

The original IsoPat algorithm, published in 2007 [@gruber2007], addressed these challenges using a matrix-based least-squares approach. The algorithm has been applied across diverse research areas, including photoredox-catalyzed deuteration and tritiation of pharmaceuticals [@loh2017], hydrogen isotope exchange using iridium and nickel catalysts [@halpin2019; @yang2022; @martinelli2024; @atzrodt2024], $^{13}$C metabolic flux analysis [@zhao2023], enzymatic mechanisms [@lara2009], carbon isotope labeling of carboxylic acids [@kinney2024], regioisomer identification via partial isotopic labeling [@sojdak2025], and the synthesis of tritium- and $^{14}$C-labeled drug candidates [@bergare2024].

However, the original implementation was provided as Excel macros and Mathematica notebooks—formats that are increasingly incompatible with modern computational workflows and difficult to integrate into automated pipelines.

**IsoPat 3.0** reimplements this validated algorithm in Python with minimal dependencies (NumPy only), a clean documented API for programmatic access, command-line interface for batch processing, and comprehensive tests for the core deconvolution algorithm. The package enables researchers to incorporate isotope pattern deconvolution into modern data analysis pipelines, Jupyter notebooks, and high-throughput screening workflows.

# Algorithm

IsoPat solves an overdetermined linear system using least-squares optimization. Given an unlabeled compound d$_0$ with natural isotope abundance pattern ($a_0$, $a_1$, ..., $a_k$) and an analyte mixture containing derivatives with 0 to n isotope labels, with the measured abundance pattern **b**, then **b** is modeled as:

$$\mathbf{A} \cdot \mathbf{x} = \mathbf{b}$$

where **A** is the pattern matrix built by shifting the unlabeled pattern for each derivative, and **x** = [x(d$_0$), x(d$_1$), ..., x(d$_n$)]$^T$ contains the relative fractions of each labeled species.

Since the number of m/z measurements typically exceeds the number of derivatives (overdetermined system), IsoPat minimizes the squared error using the pseudoinverse:

$$\mathbf{x} = (\mathbf{A}^T \mathbf{A})^{-1} \mathbf{A}^T \mathbf{b}$$

This approach is mathematically superior to iterative subtraction methods because it minimizes the total fitting error rather than propagating errors from individual peaks.

\autoref{fig:overview} illustrates the algorithm workflow using 3-octanone as an example compound. Panel A shows the natural isotope distribution of the unlabeled compound. Panel B displays the pattern matrix **A** where each sequential column represents the expected pattern for each derivative with i isotope labels. Panel C compares the measured analyte pattern with the fitted pattern (R$^2$ = 0.9999). Panel D shows the deconvolved fractions for each labeled species, with a labeled ratio of 89.9%.

![IsoPat algorithm overview. (A) Natural isotope distribution of unlabeled 3-octanone. (B) Pattern matrix A showing shifted patterns for each derivative d$_0$–d$_4$. (C) Comparison of measured and fitted analyte patterns with R$^2$ = 0.9999. (D) Deconvolution results showing relative fractions of each labeled species.\label{fig:overview}](figure_1.png)

# Usage

## Python API

```python
from isopat import deconvolve

# Unlabeled compound pattern (natural isotope distribution)
unlabeled = [100, 8.88, 0.37]  # M, M+1, M+2

# Measured mixture pattern
analyte = [10, 20, 40, 25, 5, 0.9, 0.04]

# Deconvolve to get relative amounts (H/D exchange, mass_shift=1)
result = deconvolve(unlabeled, analyte, n_labels=4)

print(result)
# IsotopePattern(d0=10.1%, d1=19.9%, d2=40.0%,
#                d3=25.2%, d4=4.9%, l.r.=89.9%, R²=0.9999)

# For 18O labeling (mass_shift=2)
result_18O = deconvolve(unlabeled, analyte, n_labels=4, mass_shift=2)
```

## Command Line

```bash
# Single pattern analysis (H/D exchange, default mass_shift=1)
isopat deconvolve -u "100,8.88,0.37" -a "10,20,40,25,5,0.9,0.04" -n 4

# For 18O or tritium labeling (mass_shift=2)
isopat deconvolve -u "100,8.88,0.37" -a "10,20,40,25,5,0.9,0.04" -n 4 --mass-shift 2

# Batch processing from files
isopat batch -u reference.csv -a samples.csv -n 4 -o results.csv
```

# Validation

The algorithm was validated against known mixtures in the original publication [@gruber2007] with regression correlation coefficients (R$^2$) of 0.98–0.99. The Python implementation includes unit tests for the core deconvolution algorithm that verify numerical correctness against synthetic test cases, including regression tests based on the 3-octanone example from the original publication. All tests pass on Python 3.9–3.12 across multiple platforms.

# Comparison with Related Software

Several tools exist for isotope pattern analysis. IsoCor [@millard2019] focuses on $^{13}$C metabolic flux analysis in Python but is metabolomics-specific. The mia package [@jungreuthmayer2016] handles $^{13}$C labeling in R but is limited to carbon isotopes. Commercial software like Xcalibur provides general MS analysis but is vendor-specific and proprietary.

IsoPat is distinguished by its support for any isotope (D, $^{13}$C, $^{15}$N, $^{18}$O, etc.), minimal dependencies (NumPy only), validation against published results, and simple well-documented API.

# Acknowledgements

We acknowledge the original co-authors of the 2007 algorithm: Gustav Oberdorfer, Constance V. Voss, Jennifer M. Kremsner, and C. Oliver Kappe. The original work was supported by the Austrian Science Fund (FWF Project P18537-B03). The University of Graz and the Field of Excellence BioHealth of the University of Graz are acknowledged for financial support.

# References
