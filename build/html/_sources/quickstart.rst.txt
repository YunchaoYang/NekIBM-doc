.. _qstart:

==============
Quickstart
==============

-------------------
Directory structure
-------------------

---------------------
Case files
---------------------

-------------------
Scripts
-------------------

Let’s walk through some useful batch scripts:

- ``makenek <case>`` compiles your case
- ``nek/nekb <case>`` runs a serial job in foreground or background
- ``nekmpi/nekbmpi <case> <number of ranks>`` runs a parallel job
- ``neknek <case1> <cas2> <ranks 1> <ranks 2>`` runs Nek5000 with two overlapping component grids 
- ``visnek <case>`` creates metadata file required by `VisIt <https://wci.llnl.gov/simulation/computer-codes/visit/>`_ and `ParaView <https://www.paraview.org/>`_. 
- ``mvn <old name> <new name>`` renames all case files
- ``cpn <old name> <new name>`` copies all case files


-------------------
Meshing
-------------------

Nek5000 is mainly a solver. However, simple box type meshes can be generated with the ``genbox`` tool. For more complex meshes please consider using ``PRENEK`` and the meshing tools ``nekmerge`` and ``n2to3``. We provide mesh converters like ``exo2nek`` and ``msh2nek`` which are quite handy if you want to use your favorite mesh generator. Also check our 
`Bazaar <https://github.com/Nek5000/NekBazaar>`_ for 3rd party meshing tools.

.. _qstart_vis:

-------------------
Visualization
-------------------
Nek5000 output (``.fld`` or ``0.f%05d``) files can be read by `VisIt <https://wci.llnl.gov/simulation/computer-codes/visit/>`_ or 
`ParaView <https://www.paraview.org/>`_. This requires using ``visnek`` to generate a metadata file.  
There is also an build-in X-Window based postprocessor called ``POSTNEK`` located in tools.

  
  
