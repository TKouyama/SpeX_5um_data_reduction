;;
;; 与えられた幾何情報を基に画像に写る惑星を再現するプログラム
;; 幾何情報の処理の仕方、計算の仕方はKouyama PhD. thesisの2章を参照のこと
;; https://dl.dropboxusercontent.com/u/58145314/papers/D_thesis_Kouyama_pub_main_en.pdf
;;

;; quaternion関連のfunctionを先にコンパイルしておく必要あり
;; qtaxis.proとqtcompose.proとqtmult.pro
  
;; Sub-routine 群 ;;
; Forward_Function TS_plot_sphere_on_ccd
; Forward_Function TS_plot_pat_point_ccd
; Forward_Function TS_ccd_latlon

;; Forward_Function TS_plot_spheroid_on_ccd

;;
;; 金星のような球体とみなせる天体を想定したプログラム
;; 有限距離からの観測を考慮するとCCDには楕円体として本来は投影されているので
;; 視直径/2の半径を長軸として、どんな観測もリムの形状は楕円とみなして計算をしている。
;; とはいえ金星半径×100倍以上の距離から撮像した場合は円とみなしても問題ない。
;; 画像処理から画像中の楕円長軸長さを視直径/2の半径として持ってこられれば矛盾ない。
;;

;
; Main
;
pro subTS_plot_ex,ssl_lat_deg,ssl_lon_deg,ssc_lat_deg,ssc_lon_deg,NA_deg,ang_diam_deg,dtheta $
              ,ccd_geo $　;; <= ジオメトリ情報を収めたデータキューブのアウトプット
              ,pix_x=pix_x, pix_y=pix_y,cx=cx, cy=cy $
              ,sub_ssl_plot=sub_ssl_plot $
              ,reverse_x= reverse_x, reverse_y=reverse_y $
              ,line_only=line_only,noerase=noerase, noplot=noplot, line_color=line_color

  ;; 入力がなかった時の処理 (テスト用)
  if n_elements(ssl_lon_deg) eq 0 then ssl_lon_deg = 60.
  if n_elements(ssl_lat_deg) eq 0 then ssl_lat_deg = 0.
  if n_elements(ssc_lon_deg) eq 0 then ssc_lon_deg = 0.
  if n_elements(ssc_lat_deg) eq 0 then ssc_lat_deg = 30.
  if n_elements(NA_deg) eq 0 then NA_deg = 90.
  if n_elements(ang_diam_deg) eq 0 then ang_diam_deg = 150./3600. ;; 
  if n_elements(dtheta) eq 0 then dtheta = 0.39/3600*!dpi/180. ;; MSI
  if n_elements(cx) eq 0 then cx = 0
  if n_elements(cy) eq 0 then cy = 0

  pix2sec = 3600.*180./!dpi*dtheta
  rad2sec = 3600.*180./!dpi

  ;; 画像描画pixel数 的な ;;
  if n_elements(pix_x) eq 0 then begin
    pix_x = 512.
  endif
  if n_elements(pix_y) eq 0 then begin
    pix_y = 512.
  endif

  prange_x = pix_x/2.*pix2sec ;; Unit = arcsec 
  prange_y = pix_y/2.*pix2sec ;; Unit = arcsec 

  ; Sub solor latitude, longitude
  ssl_lat = ssl_lat_deg*!dpi/180.
  ssl_lon = ssl_lon_deg*!dpi/180.

  ; Sub spacecraft latitude
  ssc_lat = ssc_lat_deg*!dpi/180.
  ssc_lon = ssc_lon_deg*!dpi/180.
  ; North Azimuth defined by Kouyama Phd thesis
  Na =NA_deg*!dpi/180.

  ;; Slit中心と金星中心のRA, DEC差から求める. rad単位?
  ;; 楕円中心
  x0 = tan(-cx)
  y0 = tan(-cy)

  ;; 長軸の長さ :: 視半径から換算 :: ほぼ円なので
  ;semi_a = 100.*dtheta
  semi_a = tan(ang_diam_deg/2.*!dpi/180.)
  ;print,semi_a
  
  ;; CCD平面上での原点から楕円中心までの距離
  r = sqrt(x0^2.+y0^2.) ; = tan(pi/2.-ele)
  ;; 楕円の傾き
  lam = atan(y0,x0)
  
  ;; 角度に換算 focal_lenght = 1
  fl = 1.
  phi1 = atan((r-semi_a)/fl)
  phi2 = atan((r+semi_a)/fl)
  ;; 視線ベクトルと金星中心のなす角
  phi = (phi1+phi2)/2.
  theta = (phi2-phi1)/2.
  ;print,"Phi, Theta : ",phi*180./pi,theta*180./pi," ( Phi1, Phi2 : ",phi1*180./pi,phi2*180./pi,")"
  
  ;; 金星中心の投影点
  xc = tan(phi)*cos(lam)
  yc = tan(phi)*sin(lam)
  ;print,"Xc, Yc, Lamda : ",xc,yc,lam*180./pi
  
  semi_b = semi_a*sqrt(cos(phi)^2.-sin(phi)^2.*tan(theta)^2.)
  tmpRv = 6051. + 70. ;; Venus Radius with cloud
  tmpDis = tmpRv/sin(theta)
  
  Rv = tmpRv
  Crange = tmpDis
  ;print,"## Distance/Rvenus : ",Crange/Rv, Crange,Rv
  ;print,"## Eccentricity : ",sqrt(1-(semi_b/semi_a)^2.),semi_a,semi_b

  ;; lv vector
  l_v = dblarr(3)
  abs_lv = sqrt(xc^2.+yc^2.+1.^2.)
  l_v[0] = -xc/abs_lv
  l_v[1] = -yc/abs_lv
  l_v[2] = 1./abs_lv
  
  ;print,"Lv Vector : ",l_v[0],l_v[1],l_v[2] ;,sqrt(total(l_v^2.))

  ;; Normarization by Planet radius
  Dis = Crange/Rv
  Rv = 1.
  ;; Venus Center
  V_c = dblarr(3)
  V_c[0] = Dis*l_v[0]
  V_c[1] = Dis*l_v[1]
  V_c[2] = Dis*l_v[2]
  
  ;; Vector from Venus Center to North Pole
  alpha = atan((l_v[0]*cos(Na)+l_v[1]*sin(Na))/l_v[2])
  AA = l_v[2]
  BB = (l_v[0]*cos(Na)+l_v[1]*sin(Na))
  ;print,"Alpha, sqrt", alpha*180./pi,-sin(ssc_lat)/sqrt(AA^2.+BB^2.)
  ele = asin(-sin(ssc_lat)/sqrt(AA^2.+BB^2.))-alpha
  if abs(ele) gt !dpi/2. then begin
    ele = (!dpi-abs(ele))*ele/(abs(ele))
    Na = Na+!dpi
  endif

  ;tmp_al = -asin(sin(alpha)*sqrt(AA^2.+BB^2.))
  ;print,"North Vector Elevation in CCD Frame Coordinate",ele*180./pi,tmp_al*180./pi
  
  ;; Venus 3 axes
  Axes = dblarr(3,3)
  ;; Z :: From Center to North Pole
  Axes[2,0] = cos(ele)*cos(Na)
  Axes[2,1] = cos(ele)*sin(Na)
  Axes[2,2] = sin(ele)
  ;; Y = (-l_v) × Z
  Axes[1,0] = 1./cos(ssc_lat)*(-l_v[1]*Axes[2,2]+l_v[2]*Axes[2,1])
  Axes[1,1] = 1./cos(ssc_lat)*(-l_v[2]*Axes[2,0]+l_v[0]*Axes[2,2])
  Axes[1,2] = 1./cos(ssc_lat)*(-l_v[0]*Axes[2,1]+l_v[1]*Axes[2,0])
  ;; X = Z × Y
  Axes[0,0] = (Axes[2,1]*Axes[1,2]-Axes[2,2]*Axes[1,1])
  Axes[0,1] = (Axes[2,2]*Axes[1,0]-Axes[2,0]*Axes[1,2])
  Axes[0,2] = (Axes[2,0]*Axes[1,1]-Axes[2,1]*Axes[1,0])
  
  pr = l_v[0]*Axes[0,0]+l_v[1]*Axes[0,1]+l_v[2]*Axes[0,2]
  pr = pr > (-1) & pr = pr < 1.
  ;print,"E0 Vector : ",Axes[0,0],Axes[0,1],Axes[0,2],180.-acos(pr)*180./pi
  ;print,"E1 Vector : ",Axes[1,0],Axes[1,1],Axes[1,2]
  ;print,"E2 Vector : ",Axes[2,0],Axes[2,1],Axes[2,2]

  ;;;;;;;;;;;;;;;;;;;;;;;
  ;; 線だけ表示ならここで終了 ;;
  ;;;;;;;;;;;;;;;;;;;;;;;
  if keyword_set(line_only) then begin
    nn_x = pix_x
    nn_y = pix_y
    ;; 枠
    contour,indgen(nn_x,nn_y),(dindgen(nn_x)-(nn_x-1)/2)*pix2sec,(dindgen(nn_y)-(nn_y-1)/2)*pix2sec $
      ,xrange=[-prange_x,prange_x],yrange=[-prange_y,prange_y] $
      ,xstyle=1,ystyle=1, xtitle="arcsec",ytitle="arcsec" $
      ,/fill,c_color=[196,64],noerase=noerase,/nodata
      
    ;; Draw Lon, lat, terminater line
    ;; 緯度経度線重ねたくないときはここをコメントアウト
    error = TS_plot_sphere_on_ccd(Dis,Rv,l_v,V_c,Axes,ssl_lat,ssl_lon,ssc_lon,dtheta $
      ,prange_x=prange_x, prange_y=prange_y $
      ,reverse_x=reverse_x, reverse_y=reverse_y $
      ,line_color= line_color,/terminate)
    
    omega = findgen(360)*!dpi/180.
    x = cos(omega)
    y = sin(omega)
    tmp_x = x
    tmp_y = y
    x = semi_a*tmp_x*cos(lam)-semi_b*tmp_y*sin(lam)
    y = semi_a*tmp_x*sin(lam)+semi_b*tmp_y*cos(lam)
    x = (x + x0);/dtheta
    y = (y + y0);/dtheta

    oplot,x*rad2sec,y*rad2sec,thick=1,color=line_color
    
    ;; Plot sub solr point and anit-sub solr point
    if keyword_set(sub_ssl_plot)then begin
      error = TS_plot_pat_point_ccd(Dis,Rv,l_v,V_c,Axes,ssl_lat,ssl_lon,ssc_lon,dtheta)
      error = TS_plot_pat_point_ccd(Dis,Rv,l_v,V_c,Axes,ssl_lat,ssl_lon,ssc_lon,dtheta,/anti)
    endif
    return
  endif

  ;;;;;;;;;;;;;;;;;;;;
  ;; CCD simulation ;;
  ;;;;;;;;;;;;;;;;;;;;
  ;; Geometry calculation
  L_N = Axes[2,*]
  ccd_geo = TS_ccd_latlon(dtheta,pix_x,pix_y,Dis,Rv,L_v,L_N,ssl_lat,ssl_lon,ssc_lat,ssc_lon)
  ;ccd_geo = reverse(reverse(ccd_geo,1),2)
  inc = ccd_geo[*,*,2]
  emi = ccd_geo[*,*,3]
  bright = cos(inc)
  
  ;出射角依存性を見る
;  imu = cos(inc)
;  emu = cos(emi)
;  bright = 0.59/!dpi*(abs(emu*imu))^0.90/emu*(1-exp(-imu/0.0547))/(1-exp(-emu/0.0039))
;  bright /= cos(inc)
;  bright *= cos(inc) gt 0
;  bright /= max(bright)
;  tvscl,(bright> 0.5)
;  stop

  ;;;;;;;;;;;;
  ;; Output ;; 
  ;;;;;;;;;;;;
  if keyword_set(reverse_x) then begin
    bright = reverse(bright,1)
    ccd_geo = reverse(ccd_geo,1)
  endif
  if keyword_set(reverse_y) then begin
    bright = reverse(bright,2)
    ccd_geo = reverse(ccd_geo,2)
  endif

  ;wset,1

  ;;;;;;;;;;;;;;;;;;;;
  ;; Draw CCD Image ;;
  ;;;;;;;;;;;;;;;;;;;;
  omega = findgen(360)*!dpi/180.
  x = cos(omega)
  y = sin(omega)
  tmp_x = x
  tmp_y = y
  x = semi_a*tmp_x*cos(lam)-semi_b*tmp_y*sin(lam)
  y = semi_a*tmp_x*sin(lam)+semi_b*tmp_y*cos(lam)
  x = (x + x0);/dtheta
  y = (y + y0);/dtheta
  
  ;; Draw Dayside and Night side
  bright = 1.*(bright ge 0.) + 0.*(bright gt -1.)  -1*(bright le -1)
  nn_x = pix_x
  nn_y = pix_y
  
  if keyword_set(noplot) then begin
  endif else if keyword_set(noerase) then begin
    ;loadct,39,/silent
    ;; 枠
    contour,bright,(dindgen(nn_x)-(nn_x-1)/2)*pix2sec,(dindgen(nn_y)-(nn_y-1)/2)*pix2sec $
           ,xrange=[-prange_x,prange_x],yrange=[-prange_y,prange_y] $
           ,xstyle=1,ystyle=1, xtitle="arcsec",ytitle="arcsec" $
           ,/fill,c_color=[196,64],noerase=noerase,/nodata

    ;; Draw Lon, lat, terminater line
    ;; 緯度経度線重ねたくないときはここをコメントアウト
    error = TS_plot_sphere_on_ccd(Dis,Rv,l_v,V_c,Axes,ssl_lat,ssl_lon,ssc_lon,dtheta $
                                 ,prange_x=prange_x, prange_y=prange_y $
                                 ,reverse_x=reverse_x, reverse_y=reverse_y $
                                 ,line_color= line_color)

    oplot,x*rad2sec,y*rad2sec,thick=1,color=line_color

    ;; Plot sub solr point and anit-sub solr point
    if keyword_set(sub_ssl_plot)then begin
      error = TS_plot_pat_point_ccd(Dis,Rv,l_v,V_c,Axes,ssl_lat,ssl_lon,ssc_lon,dtheta)
      error = TS_plot_pat_point_ccd(Dis,Rv,l_v,V_c,Axes,ssl_lat,ssl_lon,ssc_lon,dtheta,/anti)
    endif

  endif else begin
    ;loadct,0,/silent
    ;; 日照面プロット
    contour,bright,(dindgen(nn_x)-(nn_x-1)/2)*pix2sec,(dindgen(nn_y)-(nn_y-1)/2)*pix2sec $
           ,xrange=[prange_x,-prange_x],yrange=[prange_y,-prange_y] $
           ,xstyle=5,ystyle=5 $
           ;,/fill,c_color=[196,64] ;,/noerase
           ,/fill,c_color=[254,128] 

    contour,bright,(dindgen(nn_x)-(nn_x-1)/2)*pix2sec,(dindgen(nn_y)-(nn_y-1)/2)*pix2sec $
           ,xrange=[prange_x,-prange_x],yrange=[prange_y,-prange_y] $
           ,xstyle=5,ystyle=5 $
           ;,/fill,c_color=[196,64] ;,/noerase
           ,/fill,c_color=[254,128] 
    ;; 枠
    contour,bright,(dindgen(nn_x)-(nn_x-1)/2)*pix2sec,(dindgen(nn_y)-(nn_y-1)/2)*pix2sec $
           ,xrange=[-prange_x,prange_x],yrange=[-prange_y,prange_y] $
           ,xstyle=1,ystyle=1, xtitle="arcsec",ytitle="arcsec" $
           ,/fill,c_color=[196,64] ,/noerase,/nodata
    ;loadct,39,/silent

    ;; Draw Lon, lat, terminater line
    error = TS_plot_sphere_on_ccd(Dis,Rv,l_v,V_c,Axes,ssl_lat,ssl_lon,ssc_lon,dtheta $
                                 ,prange_x=prange_x, prange_y=prange_y $
                                 ,reverse_x=reverse_x, reverse_y=reverse_y $
                                 ,/terminate)

    ;; For spheroid
    ;error = TS_plot_spheroid_on_ccd(Dis,Rv,Rv*0.5,l_v,V_c,Axes,ssl_lat,ssl_lon,ssc_lon,dtheta $
    ;                             ,prange_x=prange_x, prange_y=prange_y $
    ;                             ,reverse_x=reverse_x, reverse_y=reverse_y $
    ;                             ,/terminate)

    ;; Limb line
    oplot,x*rad2sec,y*rad2sec,thick=1,color=line_color

    ;; Plot sub solr point and anit-sub solr point
    if keyword_set(sub_ssl_plot)then begin
      error = TS_plot_pat_point_ccd(Dis,Rv,l_v,V_c,Axes,ssl_lat,ssl_lon,ssc_lon,dtheta)
      error = TS_plot_pat_point_ccd(Dis,Rv,l_v,V_c,Axes,ssl_lat,ssl_lon,ssc_lon,dtheta,/anti)
    endif
  endelse


  ;wset,0
  
  ;print,"---"
end
