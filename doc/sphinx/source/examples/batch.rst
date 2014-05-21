.. currentmodule:: ember.input

**********
Batch Jobs
**********

This example demonstrates some features that are useful when running Ember on a
shared multi-node cluster that uses a batch queueing system such as Grid Engine.

``example_batch.py``:

.. literalinclude:: ../../../../python/ember/examples/example_batch.py

Working with Grid Engine
------------------------

The lines above starting with ``#$`` are directives for the Grid Engine batch
queueing system. Since this script is itself executable, it is not necessary to
write a separate job submission script. The various resource requests need to be
specified based on the configuration of the cluster you are using.

Running in parallel
-------------------

The option ``nThreads=8`` means that Ember will use 8 processors simultaneously.
This mirrors the resource request to the queue manager for an 8 processor node,
``-pe singlenode 8``. The name of the "parallel environment", ``singlenode`` in
this case, will depend on the configuration of your cluster. Ember is written as
a multithreaded parallel code, so all of the processors allocated to it need to
be on the same node.

Running several cases sequentially
----------------------------------

In this script, we define a function that solves a case at a particular
equivalence ratio, then call that function for each equivalence ratio of
interest. Python's string interpolation (using the ``%`` operator) is used to
generate consistent names for the output files. This method can be used to set
up parametric studies of various input variables.

To submit this job to the queue, you might use the command::

    $ qsub example_batch.py
