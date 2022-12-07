;;
;;  Sub solar or Anti solar pointをプロットする
;;
function TS_plot_pat_point_ccd,Dis,Rv,l_v,V_c,Axes,lat,lon,ssc_lon,dtheta,anti=anti
  pi = !dpi
  Dis = Dis/Rv
  Rv = 1.
  rad2sec = 180./!dpi*3600.

  ;; Plot
  ;; Projection
  ;; For Visible check
  theta = asin(Rv/Dis)
  th_s = -sin(theta)
  
  if keyword_set(anti) eq 0 then begin
    tmplat = lat
    tmplon = lon
  endif else begin
    tmplat = -lat
    tmplon = lon+!dpi
  endelse
  tmplon -= ssc_lon
  
  Vx = cos(tmplat)*cos(tmplon)*Axes[0,0] + cos(tmplat)*sin(tmplon)*Axes[1,0] + sin(tmplat)*Axes[2,0] + V_c[0]
  Vy = cos(tmplat)*cos(tmplon)*Axes[0,1] + cos(tmplat)*sin(tmplon)*Axes[1,1] + sin(tmplat)*Axes[2,1] + V_c[1]
  Vz = cos(tmplat)*cos(tmplon)*Axes[0,2] + cos(tmplat)*sin(tmplon)*Axes[1,2] + sin(tmplat)*Axes[2,2] + V_c[2]
  
  pr = l_v[0]*(Vx-V_c[0])+l_v[1]*(Vy-V_c[1])+l_v[2]*(Vz-V_c[2])
  if pr gt th_s then begin
    return,0
  endif
  
  Vx /= -Vz
  Vy /= -Vz
  
  ;Creating a Filled Circle Symbol
  psy_a = FINDGEN(17) * (!PI*2./16.)
  plot_point = where(pr le th_s)
  
  if plot_point[0] ne -1 then begin
    siz = size(plot_point)
    if keyword_set(anti) eq 0 then begin
      ;; Pixel
      ;plots,Vx[plot_point[0]]/dtheta,Vy[plot_point[0]]/dtheta,psym=1,symsize=3,thick=3,color=0
      ;; Sec
      plots,Vx[plot_point[0]]*rad2sec,Vy[plot_point[0]]*rad2sec,psym=1,symsize=3,thick=3,color=0

    endif else begin
      ;; Pixel
      ;plots,Vx[plot_point[0]]/dtheta,Vy[plot_point[0]]/dtheta,psym=7,symsize=3,thick=3,color=0
      ;; sec
      plots,Vx[plot_point[0]]*rad2sec,Vy[plot_point[0]]*rad2sec,psym=7,symsize=3,thick=3,color=255
    endelse
  endif
  
  return,0
end
