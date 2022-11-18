do for [level = 0 : 2 : 1] {
  reset

  load 'config_normal.gp'

  set terminal epslatex standalone color size lx, ly font ',17.28'
  set output sprintf('normal%d.tex', level)

  set style line 1 lc rgb '#000000' lw 3
  set style line 2 lc rgb '#AA3333' lw 7
  set style line 3 lc rgb '#FF0000' lw 11
  set style arrow 1 nohead ls 1
  set style arrow 2 head size 0.1, 15. filled front ls 2
  set style arrow 3 head size 0.2, 15. filled front ls 3

  do for [i = 1 : 4 : 1] {
    set arrow from xf[i], yf[1] to xf[i], yf[4] as 1
  }

  do for [j = 1 : 4 : 1] {
    set arrow from xf[1], yf[j] to xf[4], yf[j] as 1
  }

  set label sprintf('$% .1f$', phi00) center at xc[1], yc[1]
  set label sprintf('$% .1f$', phi10) center at xc[2], yc[1]
  set label sprintf('$% .1f$', phi20) center at xc[3], yc[1]
  set label sprintf('$% .1f$', phi01) center at xc[1], yc[2]
  set label sprintf('$% .1f$', phi11) center at xc[2], yc[2]
  set label sprintf('$% .1f$', phi21) center at xc[3], yc[2]
  set label sprintf('$% .1f$', phi02) center at xc[1], yc[3]
  set label sprintf('$% .1f$', phi12) center at xc[2], yc[3]
  set label sprintf('$% .1f$', phi22) center at xc[3], yc[3]

  if (level > 0){
    #
    nxmm = 0.5*(+phi10-phi00+phi11-phi01)/(xc[2]-xc[1])
    nymm = 0.5*(-phi10-phi00+phi11+phi01)/(yc[2]-yc[1])
    nrm = (nxmm**2.+nymm**2.)**0.5
    nxmm = nxmm / nrm
    nymm = nymm / nrm
    set arrow from xf[2], yf[2] to xf[2]+nxmm, yf[2]+nymm as 2
    #
    nxpm = 0.5*(+phi20-phi10+phi21-phi11)/(xc[3]-xc[2])
    nypm = 0.5*(-phi20-phi10+phi21+phi11)/(yc[3]-yc[2])
    nrm = (nxpm**2.+nypm**2.)**0.5
    nxpm = nxpm / nrm
    nypm = nypm / nrm
    set arrow from xf[3], yf[2] to xf[3]+nxpm, yf[2]+nypm as 2
    #
    nxmp = 0.5*(+phi11-phi01+phi12-phi02)/(xc[2]-xc[1])
    nymp = 0.5*(-phi11-phi01+phi12+phi02)/(yc[2]-yc[1])
    nrm = (nxmp**2.+nymp**2.)**0.5
    nxmp = nxmp / nrm
    nymp = nymp / nrm
    set arrow from xf[2], yf[3] to xf[2]+nxmp, yf[3]+nymp as 2
    #
    nxpp = 0.5*(+phi21-phi11+phi22-phi12)/(xc[3]-xc[2])
    nypp = 0.5*(-phi21-phi11+phi22+phi12)/(yc[3]-yc[2])
    nrm = (nxpp**2.+nypp**2.)**0.5
    nxpp = nxpp / nrm
    nypp = nypp / nrm
    set arrow from xf[3], yf[3] to xf[3]+nxpp, yf[3]+nypp as 2
  }

  if(level > 1){
    nx = 0.25*(nxmm+nxpm+nxmp+nxpp)
    ny = 0.25*(nymm+nypm+nymp+nypp)
    nrm = (nx**2.+ny**2.)**0.5
    nx = nx / nrm
    ny = ny / nrm
    set arrow from xc[2], yc[2] to xc[2]+nx, yc[2]+ny as 3
  }

  plot \
    NaN notitle
}

{
  reset

  load 'config_normal.gp'

  set terminal epslatex standalone color size lx, ly font ',17.28'
  set output 'normal3.tex'

  set style line 1 lc rgb '#000000' lw 3
  set style line 2 lc rgb '#000000' lw 13
  set style line 3 lc rgb '#000000' lw 6.5
  set style line 4 lc rgb '#FF0000' lw 4
  set style line 5 lc rgb '#FF0000' lw 6.5
  set style arrow 1 nohead ls 1
  set style arrow 2 head size 0.30, 15. filled front ls 2
  set style arrow 3 head size 0.15, 15. filled front ls 3
  set style arrow 4 head size 0.10, 15. filled front ls 4
  set style arrow 5 head size 0.25, 15. filled front ls 5

  ## old

  set arrow from 0.25+xmin, ymin to 0.25+1.50, ymin as 1 dt 3
  set arrow from 0.25+xmin, ymin to 0.25+xmin, ymax as 1 dt 3
  set arrow from 0.25+1.50, ymin to 0.25+1.50, ymax as 1 dt 3
  set arrow from 0.25+xmin, ymax to 0.25+1.50, ymax as 1 dt 3

  set arrow from 0.25+xmin, ymin to 0.25+1.50, ymin as 2
  set arrow from 0.25+xmin, ymin to 0.25+xmin, 1.50 as 2
  set label '$\underline{e}_x$' center at 0.15+1.50, ymin+0.25
  set label '$\underline{e}_y$' center at 0.50+xmin, 1.40

  ox0 = 0.8
  oy0 = 1.7
  scl = 0.5
  ox1 = ox0+scl*(xf[3]-xf[2])
  oy1 = oy0+scl*(yf[3]-yf[2])
  set arrow from ox0, 0.5*(oy0+oy1) to ox1, 0.5*(oy0+oy1) as 3
  set arrow from 0.5*(ox0+ox1), oy0 to 0.5*(ox0+ox1), oy1 as 3
  set arrow from ox0, oy0 to ox1, oy0 as 1
  set arrow from ox0, oy0 to ox0, oy1 as 1
  set arrow from ox1, oy0 to ox1, oy1 as 1
  set arrow from ox0, oy1 to ox1, oy1 as 1
  set arrow \
    from 0.5*(ox0+ox1), 0.5*(oy0+oy1) \
    to   0.5*(ox0+ox1)+scl*nx, 0.5*(oy0+oy1)+scl*ny \
    as 4
  set label '$\underline{e}_X$' center at ox1, oy0-0.2
  set label '$\underline{e}_Y$' center at ox0-0.2, oy1

  ## new

  ox0 = 2.5
  oy0 = 0.5
  ox1 = 4.5
  oy1 = 2.5
  set arrow from ox0, 0.5*(oy0+oy1) to ox1, 0.5*(oy0+oy1) as 2
  set arrow from 0.5*(ox0+ox1), oy0 to 0.5*(ox0+ox1), oy1 as 2
  set arrow from ox0, oy0 to ox1, oy0 as 1
  set arrow from ox0, oy0 to ox0, oy1 as 1
  set arrow from ox1, oy0 to ox1, oy1 as 1
  set arrow from ox0, oy1 to ox1, oy1 as 1
  set label '$\underline{e}_X$' center at ox1-0.2, oy0+0.5*(oy1-oy0)-0.2
  set label '$\underline{e}_Y$' center at ox0+0.5*(ox1-ox0)-0.2, oy1-0.2

  nx = nx / (xf[3]-xf[2])
  ny = ny / (yf[3]-yf[2])
  nrm = (nx**2.+ny**2.)**0.5
  nx = nx / nrm
  ny = ny / nrm
  scl = 1.0
  set arrow \
    from 0.5*(ox0+ox1), 0.5*(oy0+oy1) \
    to   0.5*(ox0+ox1)+nx*scl, 0.5*(oy0+oy1)+ny*scl \
    as 5

  plot \
    NaN notitle
}
