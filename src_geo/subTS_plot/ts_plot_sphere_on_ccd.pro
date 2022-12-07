;;
;; CCD上に球面を描く
;;
function TS_plot_sphere_on_ccd,Dis,Rv,l_v,V_c,Axes,ssl_lat,ssl_lon,ssc_lon,dtheta $
  ,prange_x=prange_x, prange_y=prange_y $
  ,reverse_x=reverse_x, reverse_y=reverse_y $
  ,terminate = terminate,line_color=line_color
  pi = !dpi
  rad2sec = 180./!dpi*3600.
  ;Dis = Dis/Rv
  ;Rv = 1.
  
  ;; Resolution
  if n_elements(prange_x) eq 0 then begin
    ;; Nayoro
    prange_x = 0.39/3600.*!dpi/180./dtheta*200.
  endif
  if n_elements(prange_y) eq 0 then begin
    prange_y = pragne_x
  endif
  
  ;; ひっくり返しに対応　もともとはピンホール条件を基準にしている
  if keyword_set(reverse_x) then prange_x = -prange_x
  if keyword_set(reverse_y) then prange_y = -prange_y
  
  V_x0 = Rv*Axes[0,0]+V_c[0]
  V_y0 = Rv*Axes[0,1]+V_c[1]
  V_z0 = Rv*Axes[0,2]+V_c[2]
  
  V_x1 = Rv*Axes[1,0]+V_c[0]
  V_y1 = Rv*Axes[1,1]+V_c[1]
  V_z1 = Rv*Axes[1,2]+V_c[2]
  
  V_x2 = Rv*Axes[2,0]+V_c[0]
  V_y2 = Rv*Axes[2,1]+V_c[1]
  V_z2 = Rv*Axes[2,2]+V_c[2]
  
  x = 0.*findgen(10)
  y = 0.*findgen(10)
  z = 0.*findgen(10)
  
  ;; Plot flag
  ps = 0.
  if n_elements(line_color) eq 0 then line_color=0
  
  ;; For Visible check
  theta = asin(Rv/Dis)
  th_s = -sin(theta)
  eps = 1.0e-2
  
  
  ;; Plot Lon Lat lines
  ;; 経線のインターバル
  lon_int = 15 * 2d
  lon_r = 360./lon_int

  ;; 緯線のインターバル
  lat_int = 15 * 2d
  lat_r = 90./lat_int
  
  ;;;;;;;;;;;;;;;;;;;;
  ;; Longitude line ;;
  ;;;;;;;;;;;;;;;;;;;;
  pn = 360
  Vx = dblarr(pn)
  Vy = dblarr(pn)
  Vz = dblarr(pn)
  ;ssc_lon=7.5*!dpi/180.*0
  for j = 0,lon_r-1,1 do begin
    lon = double(j)*lon_int*pi/180. - ssc_lon ;; ssc_lonの符号に疑問あり
    for i=0,pn-1,1 do begin
      ;lat = (double(i)*180./(pn-1)-(90.-lat_int) < (90.-lat_int) )*pi/180.
      lat = (double(i)*180./(pn-1)-(90.) < (90.) )*pi/180.
      Vx[i] = cos(lat)*cos(lon)*Axes[0,0] + cos(lat)*sin(lon)*Axes[1,0] + sin(lat)*Axes[2,0] + V_c[0]
      Vy[i] = cos(lat)*cos(lon)*Axes[0,1] + cos(lat)*sin(lon)*Axes[1,1] + sin(lat)*Axes[2,1] + V_c[1]
      Vz[i] = cos(lat)*cos(lon)*Axes[0,2] + cos(lat)*sin(lon)*Axes[1,2] + sin(lat)*Axes[2,2] + V_c[2]
    endfor
    pr = l_v[0]*(Vx-V_c[0])+l_v[1]*(Vy-V_c[1])+l_v[2]*(Vz-V_c[2])
    Vx /= -Vz
    Vy /= -Vz

    ;; 線幅 ;;
    if double(j)*lon_int eq 90. then begin
      line_thick = 1
      ;stop
    endif else begin
      line_thick = 1
    endelse

    if min(pr) le th_s then begin
      if ps eq 0 then begin
        plot,Vx[where(pr le th_s)]*rad2sec,Vy[where(pr le th_s)]*rad2sec,thick=line_thick $
          ,xrange=[prange_x,-prange_x],yrange=[prange_y,-prange_y] $
          ,xstyle=5,ystyle=5,color=line_color,/noerase

        ps=1
      endif else if ps ne 0 then begin
        oplot,Vx[where(pr le th_s)]*rad2sec,Vy[where(pr le th_s)]*rad2sec,thick=line_thick,color=line_color
      endif
    endif

  endfor
  
  ;;;;;;;;;;;;;;;;;;;
  ;; Latitude line ;;
  ;;;;;;;;;;;;;;;;;;;
  pn = 360
  Vx = dblarr(pn)
  Vy = dblarr(pn)
  Vz = dblarr(pn)
  for j = -lat_r,lat_r,1 do begin
    lat = double(j)*lat_int*pi/180.
    for i=0,pn-1,1 do begin
      lon = (double(i-pn/2.)*360./(pn-1))*pi/180.
      Vx[i] = cos(lat)*cos(lon)*Axes[0,0] + cos(lat)*sin(lon)*Axes[1,0] + sin(lat)*Axes[2,0] + V_c[0]
      Vy[i] = cos(lat)*cos(lon)*Axes[0,1] + cos(lat)*sin(lon)*Axes[1,1] + sin(lat)*Axes[2,1] + V_c[1]
      Vz[i] = cos(lat)*cos(lon)*Axes[0,2] + cos(lat)*sin(lon)*Axes[1,2] + sin(lat)*Axes[2,2] + V_c[2]
    endfor
    pr = l_v[0]*(Vx-V_c[0])+l_v[1]*(Vy-V_c[1])+l_v[2]*(Vz-V_c[2])
    if min(pr) le th_s then begin
      Vx /= -Vz
      Vy /= -Vz
      if ps eq 0 then begin
        plot,Vx[where(pr le th_s)]*rad2sec,Vy[where(pr le th_s)]*rad2sec,thick=1 $
          ,xrange=[prange_x,-prange_x],yrange=[prange_y,-prange_y] $
          ,xstyle=1,ystyle=1,color=line_color,/noerase
        ps = 1.
      endif else begin
        oplot,Vx[where(pr le th_s)]*rad2sec,Vy[where(pr le th_s)]*rad2sec,thick=1,color=line_color
      endelse

      ;; Equator
      if j eq 0 and ps ne 0 then begin
        ;loadct,0,/silent
        oplot,Vx[where(pr le th_s)]*rad2sec,Vy[where(pr le th_s)]*rad2sec,thick=2,color=line_color;,color=0
        ;loadct,39,/silent
      endif
    endif
  endfor
  
  ;;;;;;;;;;;;;;;;;;;;;
  ;; Terminater Line ;; using quoternion
  ;;;;;;;;;;;;;;;;;;;;;
  ;; 太陽方向ベクトル=回転軸
  ssl_axis = dblarr(3)
  ssl_axis[0] = cos(ssl_lat)*cos(ssl_lon - ssc_lon)
  ssl_axis[1] = cos(ssl_lat)*sin(ssl_lon - ssc_lon)
  ssl_axis[2] = sin(ssl_lat)
  absv = sqrt(total(ssl_axis^2.))
  vec = ssl_axis/absv
  
  ; Original point :: y0 = 0 , i.e. 90 deg
  x0 = vec[2]/sqrt(vec[0]^2.+vec[2]^2.)
  z0 = -vec[0]/sqrt(vec[0]^2.+vec[2]^2.)
  ; Original Vevctor :: 回転の始点
  v0 = dblarr(3)
  v0[0]=x0
  v0[1]=0
  v0[2]=z0
  v0_vec = v0/sqrt(total(v0^2.))
  
  ;; Original Vevctor の Quaternion
  qt_v0 = qtcompose(v0_vec,!pi)
  
  ;; Rotation of Theta from Original Vector
  v1 = dblarr(361,3)
  theta = dindgen(361)*!dpi/180.
  for i=0,360,1 do begin
    qt = qtcompose(vec, theta[i])
    qt_in = qtcompose(vec, -theta[i])
    tmp_v1 = qtmult(qtmult(qt_in,qt_v0),qt)
    v1[i,0:2] = real_part(tmp_v1[0:2])
  endfor
  lat = asin(v1[*,2])
  lon = atan(v1[*,1],v1[*,0])
  
  Vx = cos(lat)*cos(lon)*Axes[0,0] + cos(lat)*sin(lon)*Axes[1,0] + sin(lat)*Axes[2,0] + V_c[0]
  Vy = cos(lat)*cos(lon)*Axes[0,1] + cos(lat)*sin(lon)*Axes[1,1] + sin(lat)*Axes[2,1] + V_c[1]
  Vz = cos(lat)*cos(lon)*Axes[0,2] + cos(lat)*sin(lon)*Axes[1,2] + sin(lat)*Axes[2,2] + V_c[2]
  pr = l_v[0]*(Vx-V_c[0])+l_v[1]*(Vy-V_c[1])+l_v[2]*(Vz-V_c[2])
  
  Vx /= -Vz
  Vy /= -Vz
  
  front_pos = where(pr lt 0)
  if front_pos[0] ne -1 then begin
    Vx = Vx[front_pos]
    Vy = Vy[front_pos]
    
    theta_s = shift(theta[front_pos],1)
    d_theta = abs(theta[front_pos]-theta_s)
    pr_siz = size(Vx)
    for i=1,pr_siz[1]-1,1 do begin
      if d_theta[i]*180./!dpi gt 2. then begin
        Vx = shift(Vx,-i)
        Vy = shift(Vy,-i)
        ;print,i
        break
      endif
    endfor
    
    ;; 昼夜境界線 ;;
    if keyword_set(terminate) then begin
      ;loadct,0,/silent
      oplot,Vx*rad2sec,Vy*rad2sec,thick=2,color=64,linestyle=2
      ;loadct,39,/silent
    endif
  endif
  
  return,0
end
