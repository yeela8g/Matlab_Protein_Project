# Protein-Ligand Interaction Analyzer

The Protein-Ligand Interaction Analyzer is a MATLAB-based program designed to analyze and visualize the interactions between proteins
and ligands within the binding site based on a Protein Data Bank (PDB) file.
This program assists researchers in understanding the spatial relationship and proximity between protein amino acids and ligand molecules within the specific binding site outlined in the PDB file.
By parsing the PDB file, the program identifies relevant protein chains and ligands within the binding site, calculates distances between atoms, and visualizes the interaction using 3D scatter plots and heatmaps.

## Usage and Features overview:
Upon loading the file, users can interactively select the protein chain for investigation and subsequently choose a ligand within the binding site for analysis.
The program provides interactive visualizations including 3D scatter plots, histograms, and heatmaps
enabling users to dynamically explore and analyze protein-ligand interactions with ease.
Key features include dynamic distance calculation, customizable distance cutoffs, and interactive visualization of protein-ligand interactions,
leveraging MATLAB's plotting functionalities for intuitive representation within the binding site.

## Running the Program:
To run the Protein-Ligand Interaction Analyzer, simply execute the main script in MATLAB - "protein.m".
Upon execution, the program will prompt the user to select a PDB file containing the protein-ligand structure of interest.
Follow the on-screen instructions to interactively choose the protein chain to investigate and then select the ligand within the defined binding site for analysis.
The program will then generate visualizations depicting the spatial arrangement and interactions between the protein backbone and ligand atoms within the binding site.

## Visualization example:
The program generates several visualizations to aid in the analysis of protein-ligand interactions within the defined binding site:
1. A 3D scatter plot illustrates the spatial distribution of protein backbone atoms and ligand molecules within the binding site (while also displaying the rest of the chain backbone):

![image](https://github.com/yeela8g/Matlab_Protein_Project/assets/118124478/40f5d78a-192c-4a43-8b0d-5d152d5d76e6)

2. A heatmap visualizes the distances between ligand and chain atoms in a specific chosen range:

![image](https://github.com/yeela8g/Matlab_Protein_Project/assets/118124478/6179a74c-a096-4b94-ba64-a20a6933047e)


Additionally, a histogram provides a distribution of ligand-chain atoms distances, allowing users to identify relevant cutoffs for interaction analysis:

![image](https://github.com/yeela8g/Matlab_Protein_Project/assets/118124478/bb5803d9-113f-4772-99c9-7c6721e557a0)


Interactive features enable users to adjust distance parameters dynamically and explore different aspects of the protein-ligand interaction within the binding site.

Enjoy your research!
