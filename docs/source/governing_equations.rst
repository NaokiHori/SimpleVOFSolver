###################
Governing equations
###################

In addition to the incompressibility constraint

.. math::

   \der{u_i}{x_i} = 0,

the momentum balance

.. math::

   \der{u_i}{t}
   +
   u_j \der{u_i}{x_j}
   =
   -\der{p}{x_i}
   +
   \frac{\sqrt{Pr}}{\sqrt{Ra}} \der{}{x_j} \der{u_i}{x_j}
   +
   T \delta_{x i}
   +
   f_i^{st},

and the equation of the internal energy (temperature)

.. math::

   \der{T}{t}
   +
   u_j \der{T}{x_j}
   =
   \frac{1}{\sqrt{Pr}\sqrt{Ra}} \der{}{x_j} \der{T}{x_j},

I need to solve the advection equation with respect to the indicator function :math:`H`

.. math::

   \frac{DH}{Dt}
   =
   0

or

.. math::

   \der{H}{t}
   +
   u_j \der{H}{x_j}
   =
   0,

where :math:`H` is a flag to distinguish the two liquid phases, which takes :math:`0` (in phase 0) or :math:`1` (in phase 1).

Note that :math:`f_i^{st}` is included in the momentum balance, describing the effect of the surface tension force:

.. math::

   f_i^{st}
   \equiv
   \frac{1}{We} \kappa \delta n_i,

where :math:`We`, :math:`\kappa`, :math:`\delta`, and :math:`n_i` are the Weber number, local interfacial (mean) curvature, Dirac delta function and surface normal vector, respectively.
Here, I assume that the surface is ideally clean and no tangential force (Marangoni stress) is caused.

.. seealso::

   `SimpleNSSolver <https://naokihori.github.io/SimpleNavierStokesSolver/governing_equations/index.html>`_ for details.

