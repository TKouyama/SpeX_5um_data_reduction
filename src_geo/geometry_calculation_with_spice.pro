;;
;;
;;
function geometry_calculation_with_spice $
  ,fldname_kernel, metakr $
  ,utctim, obs_lat, obs_lon, obs_alt

  common obs_name,TARGET_BODY, OBS_BODY

  common obs_geo,Azi_Lv,Ele_Lv,ang_diam,lam $
    ,Na_deg_J2000, EV_distance, VS_distance $
    ,sub_solar_lon,sub_solar_lat,sub_earth_lon $
    ,sub_earth_lat,sub_venus_lat,angle_sun_venus_earth,angle_venus_earth_sun

  ;;;;;;;;;;;;;;;;;;;
  ;; Input kernels ;;
  ;;;;;;;;;;;;;;;;;;;
  cd, fldname_kernel, current = old ;; Kernel が置かれているディレクトリに移動
  CSPICE_FURNSH, METAKR ;; Kernel の読み込み
  cd, old ;; 元に戻る

  ;TARGET_BODY = 'Venus'
  ;OBS_BODY = 'Earth'

  case target_body of
    'EARTH': TARGET_BODY_F = 'ITRF93'
    'MOON': TARGET_BODY_F = 'MOON_ME'
    else:  TARGET_BODY_F = 'IAU_'+TARGET_BODY
  endcase

  case obs_body of
    'EARTH': OBS_BODY_F = 'ITRF93'
    'MOON': OBS_BODY_F = 'MOON_ME'
    else:  OBS_BODY_F = 'IAU_'+OBS_BODY
  endcase

  ;; UTC => Et
  CSPICE_STR2ET, utctim, et  ;; lsk\naif0010.tls

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; 気球周りのベクトルや気球を基準とするあれこれについて計算 ;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;; 気球のzベクトルを求める。これと気球から金星に向かうベクトルを基準に自転軸の傾きを求める
  ;; 地球偏平率を求める => pgrrecに必要
  cspice_bodvrd,OBS_BODY, 'RADII', 3, radii
  e_f = (radii[0]-radii[2])/radii[0]

  ;; 地球中心から気球までのベクトル
  cspice_pgrrec, OBS_BODY, obs_lon, obs_lat, obs_alt, radii[0], e_f, b_rectan
  ;; 地球中心から地表-10までのベクトル
  cspice_pgrrec, OBS_BODY, obs_lon, obs_lat, -10., radii[0], e_f, s_rectan

  ;; 気球直上ベクトル = zベクトルの計算  地表面から気球に向かうベクトルの向き
  z = (b_rectan-s_rectan)/sqrt(total((b_rectan-s_rectan)^2.))

  ;; 重力基準でzベクトルを考えるver ;;
  G = 6.673e-11
  r_e = 6.378e6
  r_p = 6.356e6
  r_lat = sqrt(total(b_rectan^2.))*1e3
  a_lat = sqrt(b_rectan[0]^2.+b_rectan[1]^2.)*1e3
  M_e = 5.974e24
  om_e = 7.269e-5
  g_lat = G*M_e/r_lat^2.
  c_lat = a_lat*om_e^2. ;; 遠心力
  alpha = 1.-c_lat/g_lat

  g_z = b_rectan/sqrt(total(b_rectan^2.))
  g_z[0] = alpha*g_z[0]
  g_z[1] = alpha*g_z[1]

  ;; 重力verに更新するなら以下を使う
  ;z = g_z/sqrt(total(g_z^2.))

  ;; 北極ベクトルとのなす角
  cspice_pgrrec, OBS_BODY, 0., !dpi/2., 0., radii[0], e_f, NP_rectan
  e_NP = NP_rectan/sqrt(total(NP_rectan^2.));[0,0,1.] ;; [0,0,1] = 地球中心から地球の北極へ向かうベクトル

  theta_NP_z = acos(total(z*e_NP))
  ;; 気球から見て東を向くベクトル Eastward vector
  if theta_NP_z ne 0. then begin
    EWv = crossp(e_NP,z)/abs(sin(theta_NP_z))
  endif else begin
    EWv = [1,0,0]
  endelse
  ;; 気球から見て北を向くベクトル Northward vector
  NWv = crossp(z,EWv)

  ;; 気球から金星へ向かうLvベクトルの計算
  ;; 地球から見た金星位置を求めておく + 光航時間 (地球基準座標)
  CSPICE_SPKPOS, TARGET_BODY,et,OBS_BODY_F,"LT+S", OBS_BODY,V_pos, ltime
  ;; 気球位置からのベクトルにし、規格化する
  Lv = (V_pos - b_rectan)
  Lv /= sqrt(total((V_pos - b_rectan)^2.))
  abs_Lv = abs(total(Lv^2.))

  ;; 金星からみた太陽位置 ;;
  CSPICE_SPKPOS, 'SUN',et-ltime,TARGET_BODY_F,"LT+S", OBS_BODY,SUN_pos_V, ltime_SUN

  ;; Target - Observer distance
  EV_distance = sqrt(total(V_pos^2.))
  ;; Target - Observer distance
  VS_distance = sqrt(total(SUN_pos_V^2.))

  ;; 気球から見た金星の方位角 北(0) - 東(90) - 南 (180) - 西(270) - 北(360)
  Azi_Lv = !dpi/2.-atan(total(Lv*NWv),total(Lv*EWv))
  if Azi_Lv lt 0 then  Azi_Lv = !dpi*2. + Azi_Lv

  ;; 気球から見た金星の仰角
  ;; zとLvのなす角度をもとめる。 この角度は 90°- Elevation と一致
  theta_z_Lv = acos(total(z*Lv))
  Ele_Lv = !dpi/2.-theta_z_Lv

  ;; Elevationが0°以下のとき観測不成立
  if Ele_Lv lt 0. then begin
    print, "!!! Not observable at given date, Venus is below horizon seen from observer. !!!"
    ;; 参考
    ;    Azi_Lv_deg = Azi_Lv*180./!dpi
    ;    Azi_Lv_min = (Azi_Lv_deg-long(Azi_Lv_deg))*60.
    ;    Azi_Lv_sec = (Azi_Lv_min-long(Azi_Lv_min))*60.
    ;
    ;    Ele_Lv_deg = Ele_Lv*180./!dpi
    ;    Ele_Lv_min = (Ele_Lv_deg-long(Ele_Lv_deg))*60.
    ;    Ele_Lv_sec = (Ele_Lv_min-long(Ele_Lv_min))*60.
    ;    print," Azimuth: ", long(Azi_Lv_deg),long(Azi_Lv_min), Azi_Lv_sec ;, " [deg]"
    ;    print," Elevation: ", long(Ele_Lv_deg),long(Ele_Lv_min), Ele_Lv_sec ; " [deg]"
    ;
    ;    return, -1
  endif

  ;; 外積計算を利用して投影画像面のx軸、y軸を求める
  ;; x軸はzとLvベクトル(=視線ベクトル)両者に垂直となるように定義
  ;; x = ±(z × Lv)/sin(theta)
  ;; このxは画像の水平方向と一致するはず
  ;; 一致しなくてもズレ角度は機器の調整からわかると思われる
  ;; (符号は機器のセットアップにより決定, ここでは＋を使ってみる)
  if theta_z_Lv ne 0. then begin
    x = crossp(z,Lv)/abs(sin(theta_z_Lv))
  endif else begin
    ;; z = Lvとなるときは
    ;; 角度が出意義出来ない特異点
    stop
  endelse
  ;; y = Lv × x
  ;; このyは画像の垂直方向と一致する
  y = crossp(Lv,x)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; 画像に写る金星自転軸に関わる計算;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;; 気球から金星に向かうLvベクトルがカメラの視線ベクトルと一致すると仮定する
  ;; 金星は地球から十分遠いので観測では金星に向かうベクトルとカメラの視線ベクトルは
  ;; だいたい一致することからこの仮定は成立する

  cspice_bodvrd,TARGET_BODY, 'RADII', 3, target_radii
  if TARGET_BODY eq "VENUS" or TARGET_BODY eq "Venus" then target_radii += 70.

  v_f = (target_radii[0]-target_radii[2])/target_radii[0]
  pix_scale = atan(0.175/3600.*!dpi/180.)

  ;; 時刻矛盾が無いように天体の北極ベクトルを算出する ;;
  ;; 金星表面を指す3つのベクトルをSPICEで求めて連立方程式を解いて求める ;;
  ;; 0: 直下点
  fov_0 = V_pos/sqrt(total(V_pos^2.))
  cspice_sincpt, "Ellipsoid", TARGET_BODY, et, TARGET_BODY_F, 'LT+S', $
    OBS_BODY, OBS_BODY_F, fov_0, V_spoint_0, trgepc, $
    srfvec_0, found_0
  cspice_recpgr, TARGET_BODY, V_spoint_0, target_radii[0], v_f, lon_0, lat_0, alt
  V_distance_0 = sqrt(total(srfvec_0^2.))
  V_pos_0 = fov_0*V_distance_0
  c_0 = [cos(lat_0)*cos(lon_0),cos(lat_0)*sin(lon_0),sin(lat_0)]

  V_pos = fov_0*(V_distance_0+target_radii[0])
  distance = sqrt(total(V_pos^2.))

  ;; 1: +x方向ずれ
  fov_1 = fov_0 +x*pix_scale*1.+y*pix_scale*0
  fov_1 /= sqrt(total(fov_1^2.))
  cspice_sincpt, "Ellipsoid", TARGET_BODY, et, TARGET_BODY_F, 'LT+S', $
    OBS_BODY, OBS_BODY_F, fov_1, V_spoint_1, trgepc, $
    srfvec_1, found_1
  cspice_recpgr, TARGET_BODY, V_spoint_1, target_radii[0], v_f, lon_1, lat_1, alt
  V_distance_1 = sqrt(total(srfvec_1^2.))
  V_pos_1 = fov_1*V_distance_1
  c_1 = [cos(lat_1)*cos(lon_1),cos(lat_1)*sin(lon_1),sin(lat_1)]

  ;; 2: +y方向ずれ
  fov_2 = fov_0 -x*pix_scale*0+y*pix_scale*1
  fov_2 /= sqrt(total(fov_2^2.))
  cspice_sincpt, "Ellipsoid", TARGET_BODY, et, TARGET_BODY_F, 'LT+S', $
    OBS_BODY, OBS_BODY_F, fov_2, V_spoint_2, trgepc, $
    srfvec_2, found_2
  cspice_recpgr, TARGET_BODY, V_spoint_2, target_radii[0], v_f, lon_2, lat_2, alt
  V_distance_2 = sqrt(total(srfvec_2^2.))
  V_pos_2 = fov_2*V_distance_2
  c_2 = [cos(lat_2)*cos(lon_2),cos(lat_2)*sin(lon_2),sin(lat_2)]

  b_vec_0 = V_pos_0 - V_pos
  b_vec_1 = V_pos_1 - V_pos
  b_vec_2 = V_pos_2 - V_pos

  a_arr = [[c_0[0],c_0[1],c_0[2]] $
    ,[c_1[0],c_1[1],c_1[2]] $
    ,[c_2[0],c_2[1],c_2[2]] ]

  x_b_vec = [b_vec_0[0],b_vec_1[0],b_vec_2[0]]/target_radii[0]
  y_b_vec = [b_vec_0[1],b_vec_1[1],b_vec_2[1]]/target_radii[0]
  z_b_vec = [b_vec_0[2],b_vec_1[2],b_vec_2[2]]/target_radii[0]

  i_a_arr = invert(a_arr)
  x_vec = i_a_arr ## x_b_vec
  y_vec = i_a_arr ## y_b_vec
  z_vec = i_a_arr ## z_b_vec
  np_vec = [x_vec[2],y_vec[2],z_vec[2]]
  np_vec /= sqrt(total(np_vec^2.))
  Ln = np_vec

  ;; 確認 ;;
  ;  np_pos = V_pos + np_vec*target_radii[0]*1.
  ;  cspice_sincpt, "Ellipsoid", TARGET_BODY, et, TARGET_BODY_F, 'LT+S', $
  ;                 OBS_BODY, OBS_BODY_F, np_pos, V_spoint_np, trgepc, $
  ;                 srfvec_np, found_np
  ;  cspice_recpgr, TARGET_BODY, V_spoint_np, target_radii[0], v_f, lon_np, lat_np, alt
  ;  print,lat_np*180./!dpi

  ;; 金星中心から金星北極へ向かうベクトル Ln をもとめる(まずは金星基準の座標系で)
  ;; 時刻矛盾が取れない? 3つのベクトルのうち1つでも表面を指さなかった場合に使う
  if found_1 eq 0 or found_2 eq 0 then begin
    v_np_lat = !dpi/2d
    v_np_lon = 0.
    v_np_alt = 0.
    cspice_pgrrec,TARGET_BODY, v_np_lon, v_np_lat, v_np_alt, target_radii[0], v_f, v_np_rectan

    ;; 地球基準の座標系での表現に変換する
    ;; 金星基準座標系から地球基準へ変換するRotate matrixを求める (ltime 巻き戻し)
    cspice_pxform, TARGET_BODY_F, OBS_BODY_F, et-ltime, rotate
    ;; 変換
    cspice_mxv,rotate,v_np_rectan,v_np_rectan_bs
    ;; 規格化
    Ln = v_np_rectan_bs/sqrt(total(v_np_rectan_bs^2.))
  endif

  ;; 観測者から見た北極の位置 ;;
  CSPICE_SPKPOS, TARGET_BODY,et,"J2000","LT+S", OBS_BODY,V_pos_J2000, ltime
  abs_Vpos_J2000=1.

  ;; 金星中心のra dec ;;
  cspice_recrad, V_pos_J2000, abd_Vpos_J2000, venus_ra, venus_dec

  ;; 金星北極のra dec ;;
  cspice_pxform, TARGET_BODY_F, "J2000", et-ltime, rotate_J2000
  tmp_Ln=[0,0,1]
  cspice_mxv,rotate_J2000,tmp_Ln,Ln_J2000

  ;; もとめたLnを使う ;;
  ;cspice_pxform, OBS_BODY_F, "J2000", et, rotate_J2000
  ;cspice_mxv,rotate_J2000,Ln,Ln_J2000

  Ln_real = V_pos_J2000 + Ln_J2000*target_radii[0]
  cspice_recrad, Ln_real, abs_Ln_real, venus_n_ra, venus_n_dec


  ;; 天の東からの角度(CCW)を求める
  Na_deg_J2000 = atan((venus_n_dec-venus_dec),(-venus_n_ra+venus_ra))*180./!dpi

  ;; RA, Dec ;;
  ra_h = long(venus_ra*180./!dpi/15.)
  ra_m = long((venus_ra*180./!dpi/15.-ra_h)*60.)
  ra_s = (((venus_ra*180./!dpi/15.-ra_h)*60.-ra_m)*60.)

  dec_d = long(venus_dec*180./!dpi)
  dec_m = long((venus_dec*180./!dpi-dec_d)*60.)
  dec_s = (((venus_dec*180./!dpi-dec_d)*60.-dec_m)*60.)

  print,"### Center: Ra[h, m, s]",ra_h, ra_m,ra_s
  print,"### Center: Dec[d, m, s]",dec_d, dec_m,dec_s
  print,"### Center: Ra[h], Dec[deg]",venus_ra*180./!dpi/15.,venus_dec*180./!dpi
  print,"### North pole: Ra[h], Dec[deg]",venus_n_ra*180./!dpi/15.,venus_n_dec*180./!dpi
  print,"### Ra, Dec difference (sec) between Center and North pole: ",(venus_n_ra-venus_ra)*180./!dpi*3600.,(venus_n_dec-venus_dec)*180./!dpi*3600.
  print,"### North azimuthal angle from celestial equator: ",Na_deg_J2000

  ;; 画像上の自転軸の傾き(x軸基準)
  lam = atan(total(Ln*y),total(Ln*x))

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; その他ジオメトリに関わる計算 ;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; 金星における太陽直下点を求める
  CSPICE_SUBSLR, "Intercept: ellipsoid",TARGET_BODY,et,TARGET_BODY_F, "LT+S", $
    OBS_BODY, subpoint_solar, subepc, subsrfvec_solar
  CSPICE_RECLAT, subpoint_solar, t_radius_solar, sub_solar_lon, sub_solar_lat

  ;; 地球直下点を求める
  CSPICE_SUBPNT, "Intercept: ellipsoid",TARGET_BODY,et,TARGET_BODY_F, "LT+S",OBS_BODY,$
    subpoint_earth, subepc, subsrfvec_earth
  CSPICE_RECLAT, subpoint_earth, t_radius_earth, sub_earth_lon, sub_earth_lat

  ;; 地球における金星直下点の緯度。
  CSPICE_SUBPNT, "Intercept: ellipsoid",OBS_BODY,et,OBS_BODY_F, "LT+S",TARGET_BODY,$
    subpoint_venus, subepc, subsrfvec_venus
  CSPICE_RECLAT, subpoint_venus, obs_radii, sub_venus_lon, sub_venus_lat

  ;; 太陽-金星-地球のなす角
  unit_solar_vector = subpoint_solar/t_radius_solar
  unit_earth_vector = subpoint_earth/t_radius_earth
  angle_sun_venus_earth = acos(total(unit_solar_vector*unit_earth_vector))

  ;; 観測者-Target-Sunのなす角
  CSPICE_SPKPOS, OBS_BODY,et,OBS_BODY_F,"LT+S", TARGET_BODY,V_pos_OBS, ltime
  CSPICE_SPKPOS, OBS_BODY,et,OBS_BODY_F,"LT+S", "SUN",S_pos_OBS, ltime
  unit_solar_vector2 = S_pos_OBS/sqrt(total(S_pos_OBS^2.))
  unit_venus_vector = V_pos_OBS/sqrt(total(V_pos_OBS^2.))
  angle_venus_earth_sun = acos(total(unit_solar_vector2*unit_venus_vector))
  print,"### Angle of Venus-Earth-Sun: ",angle_venus_earth_sun*180./!dpi

  ;; 気球から見た金星の視直径
  Distance = sqrt(total((V_pos - b_rectan)^2.)) ; [km]
  Rv = target_radii[0];6115. ; [km]
  ang_diam = asin(Rv/Distance)*2.

  ;; Free kernels
  CSPICE_KCLEAR
  error=0
  return, error
end

