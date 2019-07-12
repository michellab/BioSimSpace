---
title: 'BioSimSpace: An interoperable Python framework for biomolecular simulation.'
tags:
  - interoperability
  - reproducibility
  - biomolecular-simulation
  - computational-chemistry
  - computational-physics
  - computational-biology
  - molecular dynamics
authors:
  - name: Lester Hedges
    orcid: 0000-0002-5624-0500
    affiliation: 1
  - name: Antonia Mey
    orcid: 0000-0001-7512-5252
    affiliation: 2
  - name: Christopher Woods
    orcid: 0000-0001-6563-9903
    affiliation: 1
  - name: Julien Michel
    orcid: 0000-0003-0360-1760
    affiliation: 2
affiliations:
 - name: Advanced Computing Research Centre, University of Bristol
   index: 1
 - name: Department of Chemistry, University of Edinburgh
   index: 2
date: \today
bibliography: paper.bib
---

# Summary

In most research communities there is not a single, unified software-framework. Instead, researchers are presented with a collection of competing packages from which they pick and choose the functionality that they desire. Interoperability between packages is often poor, with incompatibilities between file formats, hardware, etc. This leads to brittle workflows, poor reproducibility, and lock in to specific software. For the biomolecular simulation community, our solution has been the introduction of an interoperable framework that collects together the core functionality of many packages and exposes it through a simple Python API. By not choosing to reinvent the wheel, we can take advantage of all the fantastic software that has already been written, and can easily plug into new software packages as they appear. Our software can convert between many common molecular file formats and automatically find packages available within the environment on which it is run. This allows allows the user to write portable workflow components that can be run with different input, on different environments, and in completely different ways, e.g. from the command-line, within a Jupyter notebook running on a cloud server, or as part of a Knime workflow.

# Acknowledgments

This work was funded through an EPSRC Flagship Software grant: EP/P022138/1.

# References
