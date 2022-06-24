Intrinsic Planar Interface (IPI)
=========


.. _ipi_density:

Calculating an intrinsic density with IPI
----------------------------------

Firstly, define an MDAnalysis Universe

.. code-block:: python3

   u = mda.Universe('system.tpr','system.xtc')


Now call ``pysoftwhere.IPI.ipi()`` to calculate the intrinsic density for multiple frames

.. code-block:: python3
    
   ''' 
   u:           an MDAnalysis Universe that contains bond information
   core_sel:    an MDAnalysis atom selection of the atoms used to construct the intrinsic interface
   density_sel: an MDAnalysis atom selection of the atoms whose intrinsic positions will be calculated
   interface:   'Lower' or 'Upper' - which interface of the core_sel to calculate with respect to
   no_bins:     number of grid edges to use in each lateral dimension
   recombine:   'True' or 'False' - set 'True' to recombine core_sel over the PBC in the vertical direction
   cluster:     'True' or 'False' - set 'True' to check if any core_sel molecules have diffused away from the main selection and remove them from the analysis
   interpolate_interface: 'True' or 'False' - set "True' to linearly interpolate the interface if it is patchy
   '''

intrinsic_z_ = []
interface_vals_ = []

for sim_frame in range(0,100):

    intrinsic_z,interface_vals  =  pysoftwhere.IPI.ipi(u,
                                                       frame=-1,
                                                       core_sel='resname UNK and prop mass >2',
                                                       density_sel='name OW',
                                                       interface='Lower',                        
                                                       no_bins=51,
                                                       recombine=True,                          
                                                       cluster=True,                            
                                                       interpolate_interface=False)()           
    intrinsic_z_.append(intrinsic_z)
    interface_vals_.append(interface_vals)



.. _ipi_tracking:


Tracking solute molecules with IPI
----------------------------------

Firstly, define an MDAnalysis Universe

.. code-block:: python3

   u = mda.Universe('system.tpr','system.xtc')


Now call ``pysoftwhere.IPI.ipi()`` to calculate the intrinsic positions of solute atoms for multiple frames

.. code-block:: python3
    
   ''' 
   u:           an MDAnalysis Universe that contains bond information
   core_sel:    an MDAnalysis atom selection of the atoms used to construct the intrinsic interface
   density_sel: an MDAnalysis atom selection of the atoms whose intrinsic positions will be calculated
   interface:   'Lower' or 'Upper' - which interface of the core_sel to calculate with respect to
   no_bins:     number of grid edges to use in each lateral dimension
   recombine:   'True' or 'False' - set 'True' to recombine core_sel over the PBC in the vertical direction
   cluster:     'True' or 'False' - set 'True' to check if any core_sel molecules have diffused away from the main selection and remove them from the analysis
   interpolate_interface: 'True' or 'False' - set "True' to linearly interpolate the interface if it is patchy
   '''

intrinsic_z_ = []
interface_vals_ = []

for sim_frame in range(0,100):

    intrinsic_z,interface_vals  =  pysoftwhere.IPI.ipi(u,
                                                       frame=-1,
                                                       core_sel='resname UNK and prop mass >2',
                                                       density_sel='name OW',
                                                       interface='Lower',                        
                                                       no_bins=51,
                                                       recombine=True,                          
                                                       cluster=True,                            
                                                       interpolate_interface=False)()           
    intrinsic_z_.append(intrinsic_z)
    interface_vals_.append(interface_vals)

