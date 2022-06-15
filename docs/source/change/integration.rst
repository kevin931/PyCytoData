===================================
Integration with PyCytoData
===================================

As a package providing a pipeline for users, we would love to integrate
with other tools! In fact, the more packages that join us in development
efforts, the eaiser it will be for users in the field. Therefore, we
introduce two modes of integration for the community! Explore how you
can contribute to our vision here! 

.. note::

    The tiers listed here are not official endorsements or certifications. Rather,
    we are using the nomenclature to determine what we can support and to
    provide a concise idea of what to expect for users.

-------------------------------

*****************************
PyCytoData Alliance
*****************************

This is the base tier of PyCytoData integration. Namely, the workflow
is compatible with the ``PyCytoData`` object. This can be achieved in
one of the following ways:

1. The output is a ``PyCytoData`` object.
2. The output can be casted into a ``PyCytoData`` object (e.g. using a ``to_pycytodata`` method).
3. The output can become an attribute in a ``PyCytoData`` while preserving the functionality of both.

These requirements are frankly quite easy to meet as ``PyCytoData`` object is
an intutive and flexible interface. The first two criteria are mainly for
data packages, such as **Cytomulate**, which generate CyTOF datasets. The last
is more or less targted towards other downstream packages.

For example, if you have a package that perform trajectory inference and your method
returns a ``Trajectory`` object, it may be possible to store it in the ``PyCytoData`` 
object and make use of its functionalities. 


-------------------------

***************************
PyCytoData Alliance Plus
***************************

This is the higher tier of PyCytoData integration. Any packages or methods
that belong in this tier automatically qualifies for the ``PyCytoData Alliance``
tier. Here, the requirement is as follows:

- PyCytoData provides a built-in wrapper for a given method (e.g. CytofDR).

In essence, this means that PyCytoData supports an additional workflow natively.
This of course requires the coordination with our team (or at least through a pull
request). As usual, we would love to incorporate different packages into ``PyCytoData``
so that users can have as complete of a workflow as possible. If your package provides
a type of analyses that is currently missing from us, then we are particularly
interested in prioritizing the support for this!

.. note::

    The **PyCytoData Alliance Plus** does not require porting an entire
    package to PyCytoData. In fact, it doesn't even need to support all the
    customizations. All that is necessary is a wrapper to access the main
    functionalities of a tool. For more advanced usage, users can fall back
    to the **PyCytoData Alliance** mode and work with the original
    packages instead.