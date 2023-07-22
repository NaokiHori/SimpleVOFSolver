
do for [case = 1 : 2 : 1] {

  reset

  array xf[3] = [0., 2.5, 6.]
  array yf[2] = [0., 3.]

  array xc[2]
  array yc[1]

  do for [i = 1 : 2 : 1] {
    xc[i] = 0.5*xf[i] + 0.5*xf[i+1]
  }

  do for [j = 1 : 1 : 1] {
    yc[j] = 0.5*yf[j] + 0.5*yf[j+1]
  }

  xmin = 0.
  xmax = 6.
  ymin = 0.
  ymax = 3.

  mrgn = 0.25

  lx = xmax - xmin
  ly = ymax - ymin

  set xrange [xmin-mrgn : xmax+mrgn]
  set yrange [ymin-mrgn : ymax+mrgn]

  set terminal epslatex standalone color size lx, ly font ',17.28'
  set output sprintf('flux%d.tex', case-1)

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

  ox = 0.
  oy = 0.
  do for [i = 1 : 3 : 1] {
    set arrow from ox+xf[i], oy+yf[1] to ox+xf[i], oy+yf[2] as 1
  }
  do for [j = 1 : 2 : 1] {
    set arrow from ox+xf[1], oy+yf[j] to ox+xf[3], oy+yf[j] as 1
  }

  nsub = 20
  scl = 0.02

  # left cell, diffused interface
  if(case == 1) {
    do for [i = 0 : nsub : 1] {
      c = 255 - 255 / nsub * i
      if (c < 0) {
        c = 0
      }
      if (c > 255) {
        c = 255
      }
      lft = 0.7
      rgt = 0.3
      set object polygon \
        from xf[1],yf[2] \
        to   xf[2],yf[2] \
        to   xf[2],(1.-rgt)*yf[1]+rgt*yf[2]+scl*i \
        to   xf[1],(1.-lft)*yf[1]+lft*yf[2]+scl*i \
        to   xf[1],yf[2] \
        fc rgb sprintf('#%02X%02XFF', c, c) fs solid back
    }
  }

  # right cell, diffused interface
  if(case == 2){
    do for [i = 0 : nsub : 1] {
      c = 255 - 255 / nsub * i
      if (c < 0) {
        c = 0
      }
      if (c > 255) {
        c = 255
      }
      lft = 0.35
      rgt = 0.50
      set object polygon \
        from xf[2],yf[2] \
        to   xf[3],yf[2] \
        to   xf[3],(1.-rgt)*yf[1]+rgt*yf[2]+scl*i \
        to   xf[2],(1.-lft)*yf[1]+lft*yf[2]+scl*i \
        to   xf[2],yf[2] \
        fc rgb sprintf('#%02X%02XFF', c, c) fs solid back
    }
  }

  # velocity
  arrlngt = 0.75
  set style line 2 lc rgb '#FF0000' lw 20
  set style arrow 2 head size 0.5, 15 nofilled front ls 2
  if(case == 1){
    # ->
    set arrow from xf[2]-arrlngt, yc[1] to xf[2]+arrlngt, yc[1] as 2
  }else{
    # <-
    set arrow from xf[2]+arrlngt, yc[1] to xf[2]-arrlngt, yc[1] as 2
  }
  set label '$\left. u_x \right|_{i-\frac{1}{2}, j}$' center at xf[2], yc[1]-arrlngt*0.5 front

  # phi
  set object circle center xc[case], yc[1] size 0.125*arrlngt fs solid fc rgb '#000000' lc rgb '#000000' front
  if(case == 1){
    i_idx = 'i-1'
  }else{
    i_idx = 'i'
  }
  set label sprintf('$\phi_{%s, j}$', i_idx) center at xc[case], yc[1]-arrlngt*0.5 front

  plot \
    NaN notitle
}

