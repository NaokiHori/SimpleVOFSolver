###################
Governing equations
###################

In addition to the incompressibility constraint

.. math::

   \der{u_i}{x_i} = 0,

momentum balance

.. math::

   \der{u_i}{t}
   +
   u_j \der{u_i}{x_j}
   =
   -\der{p}{x_i}
   +
   \frac{\sqrt{Pr}}{\sqrt{Ra}} \frac{\partial u_i}{\partial x_j \partial x_j}
   +
   T \delta_{x i}
   +
   f_i^{st},

and temperature

.. math::

   \der{T}{t}
   +
   u_j \der{T}{x_j}
   =
   \frac{1}{\sqrt{Pr}\sqrt{Ra}} \frac{\partial T}{\partial x_j \partial x_j},

we solve advection equation with respect to :math:`H`

.. math::

   \der{H}{t}
   +
   u_j \der{H}{x_j}
   =
   0,

where :math:`H` is a flag to distinguish two liquid phases taking :math:`0` (in phase 0) or :math:`1` (in phase 1).

Note that :math:`f_i^{st}` is included in the momentum balance, which taking into account the effect of surface tension force.
For the time being, I assume that the surface is extremely clean and no tangential force is caused:

.. math::

   f_i^{st}
   \equiv
   \frac{1}{We} \kappa \delta n_i,

where :math:`We`, :math:`\kappa`, :math:`\delta`, and :math:`n_i` are Weber number, local interfacial (mean) curvature, Dirac delta, and normal vector to the surface, respectively.

For the other equations, see `SimpleNSSolver <https://naokihori.github.io/SimpleNavierStokesSolver/governing_equations/index.html>`_.

