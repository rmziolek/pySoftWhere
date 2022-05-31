# pySoftWhere: A Python Package for Analysing Soft Matter Nanostructure Interfaces

pySoftWhere (pySW) has the following classes:

1. ICSI (Intrinsic Core-Shell Interface): an algorithm for calculating the core-shell interface of polymer (or surfactant, etc.) micelles
2. IPI (Intrinsic Planar Interface): an algorithm for calculating the interface between a monolayer (or bilayer, polymer slab, etc.) and surrounding solvent (or secondary phase)
3. Clustering: algorithms to find the the distribution of clusters and largest cluster in a simulation. Also contains a robust micelle unwrapping code.
4. Network Analysis: algorithms to find solvation networks at interfaces and describe water mediated adsorption of solutes

pySW can be used to calculate the following (with respect to both planar and spherical nanostructures):

1. Intrinsic density profiles
2. Track positions of solutes with repsect to the interface (useful for tracking drug molecules, etc.)
3. Intrinsic profiles of any calculated quantity (hydrogen bonds, instantanous mobility, orientation, hydration, etc.)

To do: setup.py, documents/tutorials, add any aux functions (networks: stacking, interface network analysis, ecc/rg/hydro radius) 


