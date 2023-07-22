
.. _normal:

.. include:: /references.txt

##############
Surface normal
##############

********
Overview
********

I assume that I know the profile of :math:`\phi` at step :math:`n`, e.g.

.. figure:: images/normal0.png
   :width: 75%

   Discrete phase distribution (volume fraction) in a two-dimensioinal domain.

.. note::

   In this project, uniform grid spacings are adopted in the homogeneous (:math:`y`) direction(s), while non-uniform grids can be used in the wall-normal (:math:`x`) direction.

In the above figure, the left-bottom and the right-top cells are fully occupied by one of the two phases, while two liquids coexist in the other cells.

To start the surface reconstruction, I first compute the local surface normal.
In particular, I would like to compute the local normal vector defined at each cell center :math:`n_{\pic,\pjc}`.

There are many options, depending on the accuracy required and the computational cost.
Here, I adopt one of the most simple approaches using the fact that the surface normal is given by

.. math::

   n_i
   =
   \frac{m_i}{\left| m_i \right|},

where :math:`m_i` is the local gradient of the indicator function and is approximated as

.. math::

   m_i
   \approx
   \der{\phi}{x_i}.

To compute this quantity at the cell center :math:`\left( \pic, \pjc \right)`, one of the easiest ways would be to use the second-order-accurate central-difference scheme:

.. math::

   \der{\phi}{x} \approx \frac{\phi_{\pipp, \pjc} - \phi_{\pimm, \pjc}}{2 \Delta x},

which tends to be unstable since this skips :math:`\phi_{\pic, \pjc}` to evaluate the derivative there (not diagonally-dominant).

Instead, I adopt the Youngs approach (introduced by |YOUNGS1982|, extensively explained by |II2012|), which is discussed in the following section.

***************
Youngs approach
***************

#. Compute gradients at the cell corners

   To begin with, I compute the gradients at the cell corners (not the center).
   I use 4 (or 9 in 3D) surrounding cells to compute the gradients at the cell corners.

   Gradients in the :math:`x` direction:

   .. math::

      \vat{m_x}{\pim, \pjm} & \approx \vat{\der{\phi}{x}}{\pim, \pjm}, \\
      \vat{m_x}{\pip, \pjm} & \approx \vat{\der{\phi}{x}}{\pip, \pjm}, \\
      \vat{m_x}{\pim, \pjp} & \approx \vat{\der{\phi}{x}}{\pim, \pjp}, \\
      \vat{m_x}{\pip, \pjp} & \approx \vat{\der{\phi}{x}}{\pip, \pjp}.

   .. myliteralinclude:: /../../src/interface/curvature_tensor.c
      :language: c
      :tag: x gradient

   Gradients in :math:`y`:

   .. math::

      \vat{m_y}{\pim, \pjm} & \approx \vat{\der{\phi}{y}}{\pim, \pjm}, \\
      \vat{m_y}{\pip, \pjm} & \approx \vat{\der{\phi}{y}}{\pip, \pjm}, \\
      \vat{m_y}{\pim, \pjp} & \approx \vat{\der{\phi}{y}}{\pim, \pjp}, \\
      \vat{m_y}{\pip, \pjp} & \approx \vat{\der{\phi}{y}}{\pip, \pjp}.

   .. myliteralinclude:: /../../src/interface/curvature_tensor.c
      :language: c
      :tag: y gradient

   Gradients in :math:`z`:

   .. myliteralinclude:: /../../src/interface/curvature_tensor.c
      :language: c
      :tag: z gradient

   Note that they are gradients and not the normal vectors.
   To obtain normal, I normalise them and the results are assigned to the resulting array ``DVOF``.

   .. myliteralinclude:: /../../src/interface/curvature_tensor.c
      :language: c
      :tag: normalise and obtain corner normals

   .. figure:: images/normal1.png
      :width: 75%

      Normal vectors at the surrounding cell corners.

   .. note::

      In this image, lengths of the arrows are adjusted for the sake of appearance and thus may be inaccurate.

#. Compute normals at the cell centers

   The normal vector defined at each cell center is computed by averaging the normals defined at the surrounding cell corners:

   .. myliteralinclude:: /../../src/interface/curvature_tensor.c
      :language: c
      :tag: average nx

   .. myliteralinclude:: /../../src/interface/curvature_tensor.c
      :language: c
      :tag: average ny

   .. myliteralinclude:: /../../src/interface/curvature_tensor.c
      :language: c
      :tag: average nz

   .. figure:: images/normal2.png
      :width: 75%

      Normal vector at the cell center.

   .. note::

      In this image, lengths of the arrows are adjusted for the sake of appearance and thus may be inaccurate.
      In general, the norm of the new vector is not unity.

#. Coordinate transformation

   The normal vector at each cell center is defined on the *global* coordinate system (e.g., coordinate system used to describe the domain), whose basis is :math:`\left( \underline{e}_x, \underline{e}_y \right)`.
   For convenience, I define a *local* coordinate system which is attached to each cell :math:`\left( i, j \right)`, whose basis is :math:`\left( \underline{e}_{X}, \underline{e}_{Y} \right)`.
   These two coordinate systems are related with

   .. math::

      \begin{pmatrix}
         \underline{e}_X \\
         \underline{e}_Y
      \end{pmatrix}
      =
      \begin{pmatrix}
         \Delta x_i & 0        \\
                  0 & \Delta y
      \end{pmatrix}
      \begin{pmatrix}
         \underline{e}_x \\
         \underline{e}_y
      \end{pmatrix},

   or equivalently

   .. math::

      \begin{pmatrix}
         \underline{e}_x \\
         \underline{e}_y
      \end{pmatrix}
      =
      \begin{pmatrix}
         \frac{1}{\Delta x_i} & 0                 \\
                            0 & \frac{1}{\Delta y}
      \end{pmatrix}
      \begin{pmatrix}
         \underline{e}_X \\
         \underline{e}_Y
      \end{pmatrix}.

   .. figure:: images/normal3.png
      :width: 75%

      (Left) global coordinate system, and (right) local coordinate system used in this project.
      Notice that the directions of the reddish arrows differ, which is because the grid size is scaled.

   Thus, the above normal vector defined on the globa* coordinate system

   .. math::

      n_i^{global}
      =
      \underline{e}_x n_x
      +
      \underline{e}_y n_y

   is written as

   .. math::

      \underline{e}_X \frac{1}{\Delta x_i} n_x
      +
      \underline{e}_Y \frac{1}{\Delta y  } n_y

   on the local coordinate system.

   Finally, the surface normal defined at the cell center of the local coordinate system is computed as

   .. math::

      n_i^{local}
      \equiv
      \underline{e}_X n_X
      +
      \underline{e}_Y n_Y,

   which is

   .. math::

      n_X &= \frac{1}{k} \frac{n_x}{\Delta x_{i}}, \\
      n_Y &= \frac{1}{k} \frac{n_y}{\Delta y    },

   where

   .. math::

      k
      \equiv
      \sqrt{
         \left( \frac{n_x}{\Delta x_{i}} \right)^2
         +
         \left( \frac{n_y}{\Delta y    } \right)^2
      }.

   .. myliteralinclude:: /../../src/interface/curvature_tensor.c
      :language: c
      :tag: normalise and obtain center normals

.. note::

   On the new coordinate system, each control volume is a unit square.

