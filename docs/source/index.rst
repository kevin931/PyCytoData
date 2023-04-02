##########################################################################
PyCytoData: An Elegant Data Analysis Tool for CyTOF
##########################################################################

Do you ever dream of using a single pipeline for all your CyTOF needs?
Look no further: **PyCytoData** is all you need! We offer a complete
pipeline for you to both acquire data as well as do downstream analyses.
This is a flexible and extensible to support future packages and other
workflows.

To get started, you can reference this quick example:

.. code-block:: python

   >>> from PyCytoData import DataLoader
   >>> exprs = DataLoader.load_dataset(dataset = "levine13", preprocess=True)
   Would you like to download levine13? [y/n]y

   Download in progress...
   This may take quite a while, go grab a coffee or cytomulate it!
   >>>  exprs.run_dr_methods(methods = ["PCA", "UMAP", "ICA"])
   Running PCA
   Running ICA
   Running UMAP

For more detailed documentations and examples, look at the topics and
tutorials below! Enjoy your journey working with CyTOF datasets!

.. toctree::
   :maxdepth: 1
   :caption: Getting Started

   installation
   quickstart

.. toctree::
   :maxdepth: 1
   :caption: Tutorial

   tutorial/pycytodata
   tutorial/fileio
   tutorial/preprocessing
   tutorial/dr

.. toctree::
   :maxdepth: 1
   :caption: Development

   change/integration
   change/contribution
   change/build
   change/development
   change/index
   license

.. toctree::
   :maxdepth: 1
   :caption: Full API Reference

   documentation/index

.. toctree::
   :maxdepth: 1
   :caption: Resources

   references
   Dr. Xinlei (Shery) Wang <https://www.uta.edu/academics/faculty/profile?username=wangx9>
   Dr. Tao Wang <https://qbrc.swmed.edu/labs/wanglab/aboutpi.php>
   DBAI <https://dbai.biohpc.swmed.edu/>
   GitHub <https://github.com/kevin931/PyCytoData/>


   
***********************
Resources
***********************

For more resources on our labs, collaborators, and related projects, please visit the following:

   * `Dr. Xinlei (Shery) Wang faculty page <https://www.uta.edu/academics/faculty/profile?username=wangx9>`_
   * `Dr. Tao Wang Lab <https://qbrc.swmed.edu/labs/wanglab/aboutpi.php>`_
   * `Database for Actionable Immunology (DBAI) for more computational immunology-related tools <https://dbai.biohpc.swmed.edu/>`_
   * `GitHub Page for project development <https://github.com/kevin931/PyCytoData/>`_ 
