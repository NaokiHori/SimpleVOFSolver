# three different cases with the same normal
{
  reset

  xmin = 0.
  xmax = 10.
  ymin = 0.
  ymax = 3.

  mrgn = 0.25

  lx = xmax - xmin
  ly = ymax - ymin

  set xrange [xmin-mrgn : xmax+mrgn]
  set yrange [ymin-mrgn : ymax+mrgn]

  set terminal epslatex standalone color size lx, ly font ',17.28'
  set output 'intercept0.tex'

  unset border
  set lmargin 0.
  set rmargin 0.
  set bmargin 0.
  set tmargin 0.

  unset xtics
  unset ytics

  set size ratio -1

  set style line 1 lc rgb '#000000' lw 3
  set style arrow 1 nohead ls 1

  ln = 3.
  do for [i = 0 : 2 : 1] {
    ox = 0.
    oy = 0.
    set arrow from 3.5*i+ox,    oy    to 3.5*i+ox+ln, oy    as 1
    set arrow from 3.5*i+ox,    oy    to 3.5*i+ox,    oy+ln as 1
    set arrow from 3.5*i+ox+ln, oy    to 3.5*i+ox+ln, oy+ln as 1
    set arrow from 3.5*i+ox,    oy+ln to 3.5*i+ox+ln, oy+ln as 1
  }

  do for [i = 0 : 2 : 1] {
    scl = 0.5
    set object polygon \
      from 3.5*i+0,3 \
      to   3.5*i+3,3 \
      to   3.5*i+3,0+(i+1.)*scl \
      to   3.5*i+0,1+(i+1.)*scl \
      to   3.5*i+0,3 \
      fc rgb '#0000FF' fs solid
  }

  plot \
    NaN notitle
}

# sharp and smoothed surfaces
{
  reset

  xmin = 0.
  xmax = 6.5
  ymin = 0.
  ymax = 3.

  mrgn = 0.25

  lx = xmax - xmin
  ly = ymax - ymin

  set xrange [xmin-mrgn : xmax+mrgn]
  set yrange [ymin-mrgn : ymax+mrgn]

  set terminal epslatex standalone color size lx, ly font ',17.28'
  set output 'intercept1.tex'

  unset border
  set lmargin 0.
  set rmargin 0.
  set bmargin 0.
  set tmargin 0.

  unset xtics
  unset ytics

  set size ratio -1

  set style line 1 lc rgb '#000000' lw 3
  set style arrow 1 nohead ls 1

  ln = 3.
  do for [i = 0 : 1 : 1] {
    ox = 0.
    oy = 0.
    set arrow from 3.5*i+ox,    oy    to 3.5*i+ox+ln, oy    as 1
    set arrow from 3.5*i+ox,    oy    to 3.5*i+ox,    oy+ln as 1
    set arrow from 3.5*i+ox+ln, oy    to 3.5*i+ox+ln, oy+ln as 1
    set arrow from 3.5*i+ox,    oy+ln to 3.5*i+ox+ln, oy+ln as 1
  }

  # sharp
  scl = 0.5
  set object polygon \
    from 0,3 \
    to   3,3 \
    to   3,0+2.*scl \
    to   0,1+2.*scl \
    to   0,3 \
    fc rgb '#0000FF' fs solid

  # diffused
  nsub = 20
  do for [i = 0 : nsub : 1] {
    scl = 0.02
    c = 255 - 255 / nsub * i
    if (c < 0) {
      c = 0
    }
    if (c > 255) {
      c = 255
    }
    set object polygon \
      from 3.5+0,3 \
      to   3.5+3,3 \
      to   3.5+3,1+i*scl-0.5*scl*nsub \
      to   3.5+0,2+i*scl-0.5*scl*nsub \
      to   3.5+0,3 \
      fc rgb sprintf('#%02X%02XFF', c, c) fs solid
  }

  plot \
    NaN notitle
}
# beta effects

{
  reset

  xmin = 0.
  xmax = 5.
  ymin = 0.
  ymax = 3.5

  mrgn = 0.

  lx = xmax - xmin
  ly = ymax - ymin

  set terminal epslatex standalone color size lx, ly font ',17.28'
  set output 'intercept2.tex'

  set xlabel '$x$'
  set ylabel '$\hat{H}$'

  set xrange [-3.:+3.]
  set yrange [ 0.: 1.]

  set xtics 1.
  set ytics 0.5

  set style line 1 lc rgb '#FF0000' lw 3
  set style line 2 lc rgb '#0000FF' lw 3
  set style line 3 lc rgb '#33AA00' lw 3
  set style line 4 lc rgb '#000000' lw 3 dt 3

  set samples 1000

  f(x, b) = 0.5*(1.+tanh(b*x))

  set key right bottom

  plot \
    f(x, 1.) t '$1$'      ls 1 w l, \
    f(x, 2.) t '$2$'      ls 2 w l, \
    f(x, 4.) t '$4$'      ls 3 w l, \
    f(x, 99) t '$\infty$' ls 4 w l
}

{
  reset

  set terminal epslatex standalone color size 4.0, 4.0 font ',17.28'
  set output 'intercept3.tex'

  set xlabel '$x$'
  set ylabel '$y$'
  set zlabel '$\hat{H}$'

  set xrange [-0.5:+0.5]
  set yrange [-0.5:+0.5]
  set zrange [ 0.0:+1.0]

  unset xtics
  unset ytics
  set ztics 1.

  set format z '$% .0f$'

  beta = 4.

  nx = 1.
  ny = 2.

  nrm = (nx**2.+ny**2.)**0.5

  nx = nx / nrm
  ny = ny / nrm

  psi(x, y) = nx * x + ny * y

  f(x, y, d) = 0.5*(1.+tanh(beta*(psi(x, y)+d)))

  set xyplane 0

  set isosamples 40,40

  set style line 1 lc rgb '#0000FF' lw 2
  set style line 2 lc rgb '#0040DF' lw 2
  set style line 3 lc rgb '#0080BF' lw 2
  set style line 4 lc rgb '#00C09F' lw 2
  set style line 5 lc rgb '#00FF7F' lw 2

  set style line 1 lc rgb '#00007F' lw 2
  set style line 2 lc rgb '#0080FF' lw 2
  set style line 3 lc rgb '#7CFF79' lw 2
  set style line 4 lc rgb '#FF9400' lw 2
  set style line 5 lc rgb '#7F0000' lw 2

  set style line 1 lc rgb '#00FFFF' lw 2
  set style line 2 lc rgb '#40BFFF' lw 2
  set style line 3 lc rgb '#807FFF' lw 2
  set style line 4 lc rgb '#C03FFF' lw 2
  set style line 5 lc rgb '#FF00FF' lw 2

  splot \
    f(x, y, -0.6) notitle ls 1 w l, \
    f(x, y, -0.3) notitle ls 2 w l, \
    f(x, y,  0.0) notitle ls 3 w l, \
    f(x, y, +0.3) notitle ls 4 w l, \
    f(x, y, +0.6) notitle ls 5 w l
}

