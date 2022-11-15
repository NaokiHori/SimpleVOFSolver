
xmin = 0.0
xmax = 5.0

ymin = 0.0
ymax = 3.0

mrgn = 0.25

lx = xmax - xmin + 2. * mrgn
ly = ymax - ymin + 2. * mrgn

unset border
set lmargin 0.
set rmargin 0.
set bmargin 0.
set tmargin 0.

set xrange [xmin-mrgn : xmax+mrgn]
set yrange [ymin-mrgn : ymax+mrgn]

unset xtics
unset ytics

set size ratio -1

array xf[4] = [0./9., 4./3., 9./3., 15./3.]
array yf[4] = [0., 1., 2., 3.]
array xc[3]
array yc[3]

do for [i = 1 : 3 : 1] {
  xc[i] = 0.5*(xf[i] + xf[i+1])
}

do for [j = 1 : 3 : 1] {
  yc[j] = 0.5*(yf[j] + yf[j+1])
}

phi00 = 0.0; phi10 = 0.1; phi20 = 0.3
phi01 = 0.2; phi11 = 0.4; phi21 = 0.6
phi02 = 0.4; phi12 = 0.7; phi22 = 1.0

