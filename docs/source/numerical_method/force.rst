
.. _force:

#############################
Curvature and surface tension
#############################

Interfacial curvature is often needed to describe the surface tension force.
Although not directly relevant to the volume-of-fluid method, here I briefly explain how it is computed.

Mathematically speaking, the surface tension force (only normal component is considered, and the tangential components are ignored) is written as

.. math::

   \frac{1}{We} \kappa \delta n_i,

which is approximated as

.. math::

   \frac{1}{We} \kappa \der{\phi}{x_i},

where :math:`\phi` is the volume-of-fluid.

Here :math:`\kappa` is the local mean curvature, whose definition leads

.. math::

   \kappa
   \equiv
   -\der{n_i}{x_i}.

In the current framework, this can be directly evaluated using the normal vectors computed at the cell corners, i.e.:

.. math::

   \vat{\kappa}{\pic, \pjc}
   & =
   -\vat{\der{n_x}{x}}{\pic, \pjc}
   -\vat{\der{n_y}{y}}{\pic, \pjc} \\
   & \approx
   -\frac{\vat{n_x}{\pip, \pjc} - \vat{n_x}{\pim, \pjc}}{\vat{\Delta x}{\pic}}
   -\frac{\vat{n_y}{\pic, \pjp} - \vat{n_y}{\pic, \pjm}}{     \Delta y       }.

The right-hand-side terms are nothing else but what were computed as normal vectors at cell corners.

The implementation can be found here:

.. myliteralinclude:: /../../src/interface/tensor.c
   :language: c
   :tag: compute mean curvature from corner normals

Using these values, surface tension force in :math:`x` direction is described as

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

Same holds for the other directions, and are implemented here.

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

   Here I adopt the simplest approach, but accurate evaluation of surface curvature and the resulting surface tension force are not conclusive at all.
   See `this very nice review <https://www.annualreviews.org/doi/abs/10.1146/annurev-fluid-122316-045034>`_ for instance.

