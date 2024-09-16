# Physical geomorphometry

**Physical geomorphometry** is a specialized subfield of geomorphometry that investigates the land surface using the fundamental principles, practices, and concepts of physics such as dimension, energy, work, force, thermodynamics, and equilibria. Unlike descriptive geomorphometry, which focuses on descriptive statistics, physical geomorphometry emphasizes the physical interpretation of land surface features.

By analyzing land surface features through the lens of geomorphometric energies, it provides deeper insights into the relationships between established geomorphometric variables, indexes, and equations. This approach enables the creation of new indexes that express land surface equilibrium and disequilibrium.

For further details, please refer to the article _Physical geomorphometry for elementary land surface segmentation and digital geomorphological mapping_ by Minár et al. (2024)[^1].

## Project: Physical Geomorphometry for Physical-Geographic Research

The project, financially supported by the Slovak Research and Development Agency (reference number APVV-22-0024), is led by **prof.  Jozef Minár**, and runs from June 2023 to June 2027.

This project focuses on studying the geometric properties of landforms and their relationship to the influence of gravitational energy on surface processes (dynamics) and the long-term evolution (genesis) of land surface. One of the fundamental goals is to develop software packages that implement physically based methods for geomorphometric analysis.

## Geoinformatic Outputs

The project aims to develop geoinformatic tools and models for analyzing interconnected physical-geomorphometric characteristics and indices. All developed tools and data will be publicly available in this repository.

The project delivers a series of independently usable components that collectively form a comprehensive approach to land surface segmentation. These components include:

- Generalization of high-resolution Digital Elevation Models (DEMs) to a highly simplified level.
- Calculation of standard and non-standard Land Surface Parameters (LSPs).
- Physically based segmentation of the land surface.

## Project Structure

Each major topic is addressed through a dedicated software package available in the following subdirectories:

- [generalization](https://github.com/xiceph/physical-geomorphometry/tree/main/generalization) – Tool for generalizing DEMs.
- [lsp-calculator](https://github.com/xiceph/physical-geomorphometry/tree/main/lsp-calculator) – Software for the calculation of LSPs.
- [segmentation](https://github.com/xiceph/physical-geomorphometry/tree/main/segmentation) – Implementation of physically based land surface segmentation.

## Usage

Each package includes detailed installation and usage instructions.

[^1]: Minár, J., Drăguţ, L., Evans, I. S., Feciskanin, R., Gallay, M., Jenčo, M., & Popov, A. (2024). Physical geomorphometry for elementary land surface segmentation and digital geomorphological mapping. Earth-Science Reviews, 248, 104631. https://doi.org/10.1016/J.EARSCIREV.2023.104631
