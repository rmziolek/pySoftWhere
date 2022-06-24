.. image:: https://readthedocs.org/projects/pysoftwhere/badge/?version=latest

pySoftWhere:  Automated Interface Analysis of Soft Matter Nanostructures
=========================================================================

**pySoftWhere** is a Python library for analyzing interfaces of soft matter nanostructures, such as *micelles, nanoparticles, monolayers, and polymer films*.	

**pySoftWhere** has the following classes:

1. ICSI (Intrinsic Core-Shell Interface): an algorithm for calculating the core-shell interface of polymer (or surfactant, etc.) micelles
2. IPI (Intrinsic Planar Interface): an algorithm for calculating the interface between a monolayer (or bilayer, polymer slab, etc.) and surrounding solvent (or secondary phase)
3. Clustering: algorithms to find the the distribution of clusters and largest cluster in a simulation. Also contains a robust micelle unwrapping code.

Using these different algorithms, pySW can be used to calculate the following (with respect to both planar and spherical nanostructures):

1. Intrinsic density profiles
2. Track positions of solutes with repsect to the interface (useful for tracking drug molecules, etc.)
3. Intrinsic profiles of any calculated quantity (hydrogen bonds, instantanous mobility, orientation, hydration, etc.)



