################
Numerical method
################

It would be very easy if I could use naive second-order central difference scheme to discretise the evolution of :math:`H`, which is, as is well-known, impractical because of the numerical instability.
Thus it is (at least for the time being in the context of finite-difference or volume methods) inevitable to use up-wind schemes.
This treatment, however, suffers from numerical diffusion, which tends to smear out the sharp nature of the surface.

One solution is to control the diffusivity of :math:`H` to keep the thickness of the surface, which is categorised as phase-field method.

The other solution is to reconstruct the surface every (or one in a few) time steps, which is categorised as volume-of-fluid method.

Here I adopt the intermediate solution (very roughly speaking), in which I allow the surface to have a finite thickness, while conduct the surface reconstruction: `THINC <https://onlinelibrary.wiley.com/doi/abs/10.1002/fld.975>`_ scheme.

First of all, instead of considering the above equation as it is, I integrate in a control volume (corresponding to computational cell later, but not limited here):

.. math::

   \int_V \der{H}{t} dV
   +
   \int_V u_j \der{H}{x_j} dV
   =
   0.

By assuming commutability, the first term leads

.. math::

   \der{}{t} \int_V H dV.

By adopting the incompressibility constraint, the second terms leads

.. math::

   \int_V u_j \der{H}{x_j} dV
   =
   \int_V u_j \der{H}{x_j} dV
   +
   \int_V H \der{u_j}{x_j} dV
   =
   \int_V \der{u_j H}{x_j} dV.

Thanks to Gauss theorem, this yields

.. math::

   \int_{\partial V} u_j H n_j dS.

Note that :math:`n_j` is *not* the surface normal of the free surface, but surface normal to the control volume.

Now I consider to divide the equation by the volume of the control volume

.. math::

   \int_{V} dV,

giving

.. math::

   \der{\phi}{t}
   +
   \frac{1}{V} \int_{\partial V} u_j H n_j dS
   =
   0,

where :math:`\phi` is ``volume-of-fluid``, whose definition is explained by the name, while the right-hand-side is the numerical flux telling the amount of increase or decrease of :math:`\phi`.

So the question is how to describe the second term numerically.
In particular, I need to explicitly describe :math:`H` so that I can integrate it numerically.
To do so, here I adopt the THINC scheme.

I think it is easy to follow if I explain the numerics and their implementations at the same time.
Please find the following stuffs for details.

.. toctree::
   :maxdepth: 2

   normal
   intercept
   flux
   force

