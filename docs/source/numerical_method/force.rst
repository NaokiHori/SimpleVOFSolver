
.. _force:

.. include:: /references.txt

###################################
Curvature and surface tension force
###################################

Interfacial curvature is often needed to describe the surface tension force.
Although not directly relevant to the volume-of-fluid method, here I briefly explain how it is computed.

Mathematically speaking, the surface tension force (only normal component is considered, and the tangential components are ignored) is given by

.. math::

   \frac{1}{We} \kappa \delta n_i,

which is approximated as

.. math::

   \frac{1}{We} \kappa \der{\phi}{x_i},

where :math:`\phi` is the volume-of-fluid.

The pre-factor (inverse of :math:`We`) is computed in the initialisation process:

.. myliteralinclude:: /../../src/interface/init.c
   :language: c
   :tag: compute surface tension coefficient

:math:`\kappa` is the local mean curvature, whose definition is

.. math::

   \kappa
   \equiv
   -\der{n_i}{x_i}.

In this project, this is directly evaluated using the normal vectors computed at the cell corners, i.e.:

.. math::

   \vat{\kappa}{\pic, \pjc}
   & =
   -\vat{\der{n_x}{x}}{\pic, \pjc}
   -\vat{\der{n_y}{y}}{\pic, \pjc} \\
   & \approx
   -\frac{\vat{n_x}{\pip, \pjc} - \vat{n_x}{\pim, \pjc}}{\vat{\Delta x}{\pic}}
   -\frac{\vat{n_y}{\pic, \pjp} - \vat{n_y}{\pic, \pjm}}{     \Delta y       }.

The right-hand-side terms have already been computed in the :ref:`previous step <normal>`.

.. myliteralinclude:: /../../src/interface/curvature_tensor.c
   :language: c
   :tag: compute mean curvature from corner normals

As a result, the surface tension force in the :math:`x` direction is described as

.. math::

   \vat{f_x^{st}}{\pip, \pjc}
   \approx
   \frac{1}{We}
   \frac{
      \vat{\kappa}{\pipp, \pjc}
      +
      \vat{\kappa}{\pic,  \pjc}
   }{2}
   \frac{
      \vat{\phi}{\pipp, \pjc}
      -
      \vat{\phi}{\pic,  \pjc}
   }{\vat{\Delta x}{\pip}}.

Same applies to the other directions.

:math:`x` direction:

.. myliteralinclude:: /../../src/interface/force.c
   :language: c
   :tag: compute surface tension force in x direction

:math:`y` direction:

.. myliteralinclude:: /../../src/interface/force.c
   :language: c
   :tag: compute surface tension force in y direction

:math:`z` direction:

.. myliteralinclude:: /../../src/interface/force.c
   :language: c
   :tag: compute surface tension force in z direction

.. note::

   Here I adopt a very simple approach, but more accurate evaluation of the surface curvature is a challenging problem and not conclusive at all.
   See |POPINET2018| for instance.

