;; CCD座標→Lon Lat
function TS_ccd_latlon,dtheta,pix_x,pix_y,crange,Rv,f_Lv,f_Nv,ssl_lat,ssl_lon,ssc_lat,ssc_lon

  ;; CCD 座標定義
  c_ccd_x = dblarr(pix_x,pix_y)
  for j=0,pix_y-1,1 do begin
    c_ccd_x[*,j] = dindgen(pix_x)
  endfor
  c_ccd_y = dblarr(pix_x,pix_y)
  for i=0,pix_x-1,1 do begin
    c_ccd_y[i,*] = dindgen(pix_y)
  endfor
  
  ;; CCD中心
  ccd_center_x = (pix_x-1.)/2.;255.5
  ccd_center_y = (pix_y-1.)/2.;255.5
  
  ;; 各CCD素子の視線ベクトル
  c_ccd_Lv = dblarr(pix_x,pix_y,3)
  c_ccd_Lv[*,*,0] = -(c_ccd_x-ccd_center_x)*dtheta
  c_ccd_Lv[*,*,1] = -(c_ccd_y-ccd_center_y)*dtheta
  c_ccd_Lv[*,*,2] = 1.
  
  ;; 単位ベクトル化 ;;
  tmp_c_ccd_Lv = c_ccd_Lv
  for i=0,2,1 do begin
    tmp_c_ccd_Lv[*,*,i] /= sqrt(c_ccd_Lv[*,*,0]^2.+c_ccd_Lv[*,*,1]^2.+c_ccd_Lv[*,*,2]^2.)
  endfor
  c_ccd_Lv = tmp_c_ccd_Lv
  
  ;; Lvベクトル
  if n_elements(f_Lv) eq 0 then f_Lv = [0,0,1]
  
  ;; 各CCDから出る視線ベクトルとLvベクトルのなす角
  c_ccd_Lv_ang = acos(c_ccd_Lv[*,*,0]*f_Lv[0]+c_ccd_Lv[*,*,1]*f_Lv[1]+c_ccd_Lv[*,*,2]*f_Lv[2])
  ;; リムとLvベクトルのなす角 (= 最大角度)
  limit_ang  = asin(Rv/Crange)
  
  ;; 金星表面の投影点までの距離
  sq_surface_dis = Rv^2.-Crange^2.*sin(c_ccd_Lv_ang)^2.
  surface_dis = (Crange*cos(c_ccd_Lv_ang)-sqrt(sq_surface_dis > 0))*(c_ccd_Lv_ang le Limit_ang)
  
  ;; 金星中心から見たときの座標
  V_center = f_Lv*Crange
  V_sx = (c_ccd_Lv[*,*,0]*surface_dis - V_center[0])*(c_ccd_Lv_ang le Limit_ang)
  V_sy = (c_ccd_Lv[*,*,1]*surface_dis - V_center[1])*(c_ccd_Lv_ang le Limit_ang)
  V_sz = (c_ccd_Lv[*,*,2]*surface_dis - V_center[2])*(c_ccd_Lv_ang le Limit_ang)
  abs_V_s = sqrt(V_sx^2.+V_sy^2.+V_sz^2.)
  
  pos = where(abs_V_s ne 0)
  ;print,mean(abs_V_s[pos])
  
  ;; 値のあるところだけ規格化
  V_sx[pos]/=abs_V_s[pos]
  V_sy[pos]/=abs_V_s[pos]
  V_sz[pos]/=abs_V_s[pos]
  
  ;; 基本3軸
  f_E1 = crossp(f_Lv, f_Nv)/cos(ssc_lat)
  f_E0 = crossp(f_E1, f_Nv)
  f_E2 = f_Nv
  
  ;; 行列化
  r_mat = dblarr(3,3)
  r_mat[0,*]=f_E0
  r_mat[1,*]=f_E1
  r_mat[2,*]=f_E2
  i_r_mat = invert(r_mat)
  
  ;; 座標変換→緯度経度
  tmp_vec = dblarr(3)
  c_ccd_lon = dblarr(pix_x,pix_y)
  c_ccd_lat = dblarr(pix_x,pix_y)
  c_ccd_inc = dblarr(pix_x,pix_y)
  c_ccd_emi = dblarr(pix_x,pix_y)
  c_ccd_azi = dblarr(pix_x,pix_y)
  c_ccd_phase = dblarr(pix_x,pix_y)
  c_ccd_bri = dblarr(pix_x,pix_y)
  
  ;; 太陽直下点
  l_ssl = dblarr(3)
  l_ssl[0] = cos(ssl_lat)*cos(ssl_lon)
  l_ssl[1] = cos(ssl_lat)*sin(ssl_lon)
  l_ssl[2] = sin(ssl_lat)
  ;; 衛星直下点
  l_ssc = dblarr(3)
  l_ssc[0] = cos(ssc_lat)*cos(ssc_lon)
  l_ssc[1] = cos(ssc_lat)*sin(ssc_lon)
  l_ssc[2] = sin(ssc_lat)
  
  ;; phase angle
  ;; ssl lon - ssc lon ベクトル
  tmp_l_ssl = dblarr(3)
  ;tmp_l_ssl[0] = cos(ssl_lat)*cos(ssl_lon-ssc_lon)
  ;tmp_l_ssl[1] = cos(ssl_lat)*sin(ssl_lon-ssc_lon)
  
  ;; 東経は左手系?だからこっち? (2013/02/04)
  tmp_l_ssl[0] = cos(ssl_lat)*cos(ssc_lon-ssl_lon)
  tmp_l_ssl[1] = cos(ssl_lat)*sin(ssc_lon-ssl_lon)
  tmp_l_ssl[2] = sin(ssl_lat)
  
  ;; Rotate
  ssl_vec = r_mat##tmp_l_ssl
  ssl_vec /= sqrt(total(ssl_vec^2.))
  
  ;; Phase angle
  c_ccd_phase = acos(-c_ccd_Lv[*,*,0]*ssl_vec[0] $
    -c_ccd_Lv[*,*,1]*ssl_vec[1] $
    -c_ccd_Lv[*,*,2]*ssl_vec[2])
    
  R2 = Rv^2.
  D2 = Crange^2.
  DR2 = 2.*Rv*Crange
  l_ccd = dblarr(3)
  for j=0, pix_y-1,1 do begin
    for i=0, pix_x-1,1 do begin
      if c_ccd_Lv_ang[i,j] lt Limit_ang then begin
        tmp_vec[0]=V_sx[i,j]
        tmp_vec[1]=V_sy[i,j]
        tmp_vec[2]=V_sz[i,j]
        latlon_vec = i_r_mat##tmp_vec
        c_ccd_lat[i,j]=asin(latlon_vec[2])
        ;; "-"をつける理由未解決(2012/02/25)
        ;; 東経は左手系?だから? (2013/02/04)
        c_ccd_lon[i,j]=-atan(latlon_vec[1],latlon_vec[0])+ssc_lon
        
        l_ccd[0] = cos(c_ccd_lat[i,j])*cos(c_ccd_lon[i,j])
        l_ccd[1] = cos(c_ccd_lat[i,j])*sin(c_ccd_lon[i,j])
        l_ccd[2] = sin(c_ccd_lat[i,j])
        ;; Incident Angle
        c_ccd_inc[i,j] = acos(total(l_ssl*l_ccd))
        ;; Emission Angle
        tmp_angle = acos(total(l_ssc*l_ccd))
        tmp_X = sqrt(R2+D2-DR2*cos(tmp_angle))
        tmp_angle2 = acos((R2+tmp_X^2.-D2)/(2.*Rv*tmp_X))
        c_ccd_emi[i,j] = !dpi-tmp_angle2
        ;; Bright
        imu = cos(c_ccd_inc[i,j])
        emu = cos(c_ccd_emi[i,j])
        if imu gt 0. then c_ccd_bri[i,j] = 0.59/!dpi*(abs(emu*imu))^0.90/emu*(1-exp(-imu/0.0547))/(1-exp(-emu/0.0039))
      endif
    endfor
  endfor
  c_ccd_lon = c_ccd_lon + 2*!dpi*(c_ccd_lon lt -!dpi) - 2*!dpi*(c_ccd_lon gt !dpi)
  ;c_ccd_lon = -c_ccd_lon
  
  ;; Azimuthal angle
  imu = cos(c_ccd_inc)
  emu = cos(c_ccd_emi)
  inu = sqrt(1.-imu^2.);sin(c_ccd_inc)
  enu = sqrt(1.-emu^2.);sin(c_ccd_emi)
  azi_pos = where(inu*enu ne 0)
  
  tmp_azi = dblarr(pix_x,pix_y)
  c_ccd_azi = dblarr(pix_x,pix_y)
  tmp_cos_phase = cos(c_ccd_phase)
  tmp_azi[azi_pos] = (tmp_cos_phase[azi_pos]-(imu[azi_pos]*emu[azi_pos])) $
    /(inu[azi_pos]*enu[azi_pos])
  c_ccd_azi[azi_pos] = acos(tmp_azi[azi_pos])
  
  c_ccd_lon = c_ccd_lon*(c_ccd_lon ne 0) - !dpi*(c_ccd_lon eq 0)
  c_ccd_lat = c_ccd_lat*(c_ccd_lat ne 0) - !dpi*(c_ccd_lat eq 0)
  c_ccd_inc = c_ccd_inc*(c_ccd_inc ne 0) - !dpi*(c_ccd_inc eq 0)
  c_ccd_emi = c_ccd_emi*(c_ccd_emi ne 0) - !dpi*(c_ccd_emi eq 0)
  
  c_ccd_geo = dblarr(pix_x,pix_y,6)
  c_ccd_geo[*,*,0]=c_ccd_lon
  c_ccd_geo[*,*,1]=c_ccd_lat
  c_ccd_geo[*,*,2]=c_ccd_inc
  c_ccd_geo[*,*,3]=c_ccd_emi
  c_ccd_geo[*,*,4]=c_ccd_azi
  c_ccd_geo[*,*,5]=c_ccd_phase
  
  error = 0.
  return,c_ccd_geo
  
end
