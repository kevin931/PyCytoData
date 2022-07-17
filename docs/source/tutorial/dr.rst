#################################
Dimension Reuction with CytofDR
#################################

``CytofDR`` is a **PyCytoData Alliance Plus** member, meaning that you can run DR right inside
the ``PyCytoData`` object. Although the option is somewhat limited with only the ``run_dr_methods``
method, which is a close counterpart to that of the same method in ``CytofDR``, you can run most
of the standard DR with this method. Alternatively, you can run your own DR methods and store the
``Reductions`` method in your ``PyCytoData`` object. This tutorial showcases each of the two
methods.

To run DR, **you will need to install the  ``CytofDR`` package first.** To get started, a quick
guide can be found `here <https://cytofdr.readthedocs.io/en/stable/installation.html>`_. ``CytofDR``
is not a mandatory dependency of ``PyCytoData``.

-----------------------

*********************************************
Standard DR Workflow with ``run_dr_methods``
*********************************************

One of the standard workflows of CyTOF data analysis is DR. And if you wish to use the default
settings, you can use this simple API without having to dive into ``CytofDR``. To get it
started, you can do the following:

.. code-block:: python

    >>> from PyCytoData import DataLoader
    >>> exprs = DataLoader.load_dataset(dataset="levine13", preprocess=True)
    >>> exprs.run_dr_methods(methods=["UMAP", "open_tsne"])
    Running UMAP
    Running open_tsne
    ===> Finding 90 nearest neighbors using Annoy approximate search using euclidean distance...
    --> Time elapsed: 8.48 seconds
    ===> Calculating affinity matrix...
    --> Time elapsed: 2.80 seconds
    ===> Running optimization with exaggeration=12.00, lr=6189.67 for 250 iterations...
    Iteration   50, KL divergence 6.3260, 50 iterations in 1.4977 sec
    Iteration  100, KL divergence 5.1775, 50 iterations in 1.5622 sec
    Iteration  150, KL divergence 4.7426, 50 iterations in 1.5219 sec
    Iteration  200, KL divergence 4.4741, 50 iterations in 1.5252 sec
    Iteration  250, KL divergence 4.2860, 50 iterations in 1.5313 sec
    --> Time elapsed: 7.64 seconds
    ===> Running optimization with exaggeration=1.00, lr=6189.67 for 250 iterations...
    Iteration   50, KL divergence 3.5153, 50 iterations in 1.4306 sec
    Iteration  100, KL divergence 2.8738, 50 iterations in 1.5975 sec
    Iteration  150, KL divergence 2.4706, 50 iterations in 2.3137 sec
    Iteration  200, KL divergence 2.1967, 50 iterations in 3.2659 sec
    Iteration  250, KL divergence 1.9980, 50 iterations in 4.6722 sec
    --> Time elapsed: 13.28 seconds

Now, a ``Reductions`` object is stored in the ``reductions`` attribute. Furthermore,
the ``original_data`` and ``original_cell_types`` attributed are automatically
populated if available from the ``PyCytoData`` object. As with the standard workflow
in ``CytofDR``, you can evaluate your DR embeddings:

.. code-block:: python

    >>> exprs.reductions.evaluate(category = ["global", "local", "downstream"], auto_cluster = True, n_clusters = 20)
    Evaluating global...
    Evaluating local...
    Evaluating downstream...

For more on the standard workflow, you can read the `documentation <https://cytofdr.readthedocs.io/en/stable/quickstart.html>`_
of the ``CytofDR`` package.

--------------------------------

*********************************************
Using the Full ``CytofDR`` API
*********************************************

You can use the full ``CytofDR`` API and then add the ``Reductions`` object to the
``PyCytoData`` manually. Here, we showcase how you can use a custom setting for
DR:

.. code-block:: python

    >> from CytofDR import dr
    >>> embedding = dr.NonLinearMethods.UMAP(data = exprs.expression_matrix, out_dims=2, n_neighbors = 30, min_dist = 0)
    >>> results = dr.Reductions(reductions = {"custom_umap": embedding})

Then you can add the ``Reductions`` object into the ``PyCytoData`` object:


.. code-block:: python

    >>> exprs.reductions = results

Then, you can follow the standard procedures as usual. For more documentation on
this, see `this tutorial <https://cytofdr.readthedocs.io/en/stable/tutorial/methods.html>`_.
