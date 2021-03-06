Welcome to libsbn's documentation!
==================================

This is the documentation for interacting with libsbn from Python.
See the README on the `libsbn GitHub repository <https://github.com/phylovi/libsbn>`_ for installation and general description.

For conceptual background, see

.. toctree::

   concepts


Command-line interface
----------------------

The command-line interface is called ``vip``.
See :doc:`modules/vip.cli` for documentation.


Python API
----------

The design principle here is that *as much as is efficient* to be in Python.
So, we want Python users to be able to do everything needed to fit models, while having all of the tree manipulation and traversal happen in C++.
This does mean a reasonable amount is in C++, though: for example, if we want efficient gradients for phylogenetic model parameters, the C++ needs to have data structures for these models.

From the perspective of a Python user, code is broken into two packages:

* ``libsbn`` is the direct interface to the C++ code; this module is generated by `pybind11 <http://pybind11.readthedocs.io/>`_.
* ``vip`` has all of the Python code; this module contains code to infer model parameters and do complete calculations.

A typical workflow looks like:

1. Load some collection of trees into the ``libsbn`` tree cache that define the SBN support
2. Parse those trees to get the SBN support
3. Given this collection of trees, ``libsbn`` sets up data structures to hold SBN model parameters
4. Initialize a phylogenetic model in ``libsbn`` and set parameters
5. Sample trees given the SBN model parameters (these get stored in the ``libsbn`` tree cache, over-writing the previously-cached trees)
6. Set branch lengths for these trees according to some model (``libsbn`` makes branch length vectors accessible as NumPy arrays)
7. Evaluate model likelihoods and gradients in ``libsbn``
8. Modify parameters and go back to 5.

The ``vip`` module facilitates this loop, both by having code that calls ``libsbn``, but also by having models for the ``libsbn`` parameters.
For example, ``libsbn`` doesn't have a data structure that contains a *model* of branch lengths--- it only has a way of setting the branch lengths.
Thus ``vip`` has branch length models, including models based on split and primary subsplit pairs (PSPs).
On the other hand, ``libsbn`` has a lot of functionality that make these models easy to put together on the Python side, such as the ``PSPIndexer`` that gives the indices needed to do PSP model calculations.

Thus, there are typically two levels of state: the "raw" parameters in the ``libsbn`` module, and the parameters of models of those ``libsbn`` parameters, which is contained in ``vip``.
``vip`` also contains routines for running and benchmarking inference.

Erick's suggestions for learning the Python API
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. Learn a little about the ``libsbn`` module by looking at ``test/test_libsbn.py`` and running it via ``pytest -s test/test_libsbn.py``.
2. Learn about the ``vip`` data structures by reading the :doc:`modules/vip.burrito` code and its dependencies.
   You can find additional usage examples by looking at the tests in the ``vip/test`` directory.
3. Get an idea of how to use the burrito module by reading the :doc:`modules/vip.benchmark` code.



Modules
-------

.. autosummary::
   :toctree: modules

   libsbn
   libsbn.beagle_flags
   vip.benchmark
   vip.branch_model
   vip.burrito
   vip.cli
   vip.optimizers
   vip.priors
   vip.sbn_model
   vip.scalar_model
   vip.sgd_server

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
