;;
;; Program for simulating obsevation geometry at given observation date.
;;
;; usege:
;;
;; at first, please set IDL-path for the source directory
;; (or compiling all pro files in this example, including files in "lib" directory)
;; then, change directory where "location_files", "generic_kernels", and other files are in.
;;
;; IDL> cd, '/your/source/path/'
;;
;; then
;;
;; IDL> simulate_TO, /now, ccd_geo = ccd_geo
;;
;; or specify observation date like,
;;
;; IDL> simulate_TO, utctim='2022-02-04T12:00:00', ccd_geo = ccd_geo
;;
;; if you do not want to display calculation result,
;;
;; IDL> simulate_TO, utctim='2022-02-04T12:00:00', ccd_geo = ccd_geo, /noplot
;;
;; will work.
;;
;; "ccd_geo" is a 3D array (nx, ny, 6) containg geomery infomation
;; 0: Longitude (radian)
;; 1: Latitude (radian)
;; 2: Incident angle (radian)
;; 3: Emission angle (radian)
;; 4: Azimuthal angle (radian)
;; 5: Phase angle (radian)
;;

Forward_Function read_location_file
Forward_Function geometry_calculation_with_spice

pro simulate_TO,utctim = utctim $
               ,now = now $
               ,ifname_location = ifname_location $ ;; <= parameter file
               ,noplot = noplot $
               ,ofname_log = ofname_log $
               ,cx_p = cx_p $ ;; offset X pixels from center of image frame
               ,cy_p = cy_p $ ;; offset Y pixels from center of image frame
               ,reverse_x = reverse_x $ ;; flip the image in X direction
               ,reverse_y = reverse_y $ ;; flip the image in Y direction
               ,latlonline_only = latlonline_only $
               ;; Output
               ,ccd_geo = ccd_geo $
               ,obs_geo_struct = obs_geo_struct

  ;;
  ;; Obervation date ;;
  ;;
  if n_elements(utctim) eq 0 then begin
    utctim = '2020-10-13T00:00:00' ;; SpeX/NSFCam2 観測?
  endif

  ;; 現在時刻指定, 上記 utcimは破棄される ;;
  if keyword_set(now) then begin
    time = SYSTIME(/julian,/utc)
    CALDAT, time, M, D, Y, H, Min, S
    utctim = strcompress(string(Y)+'-'+string(M)+'-'+string(D)+'T'+ $
                       string(H)+':'+string(Min)+':'+string(long(S)),/remove)
  endif
  PRINT, utctim

  ;;;;;;;;;;;;;;;;;;;;;;;
  ;; Observer location ;;
  ;;;;;;;;;;;;;;;;;;;;;;;
  if n_elements(ifname_location) eq 0 then begin
    ;;
    ;; Spcify direcoty where location files are in
    ;;
    ifldname_location = './'

    ;; IRTF/SpeX ;;
    ;ifname_location = ifldname_location + 'locationfile_IRTF_SpeX_Image.txt'
    ifname_location = ifldname_location + 'locationfile_IRTF_SpeX_Image_NorthUp.txt'
  endif
  
  print, "Input location file: ",ifname_location
  ;print,"Change directory to: ",file_dirname(ifname_location)
  ;cd,file_dirname(ifname_location),current=old_dir

  error = read_location_file(ifname_location,obs_lat,obs_lon,obs_alt $
                           ,fldname_kernel,metakr $
                           ,Tel_res_sec,pix_x,pix_y $
                           ,Na_deg_mode, North_up $
                           ,image_rot_angle_deg=image_rot_angle_deg $
                           ,target_name = target_name $
                           ,geo_ifname = geo_ifname)

  print,fldname_kernel
  print,metakr

  obs_lat *= !dpi/180.
  obs_lon *= !dpi/180.

  ;;;;;;;;;;;;;;;;;;;;
  ;; 再現画像用パラメータ ;;
  ;;;;;;;;;;;;;;;;;;;;

  ;; 望遠鏡or検出器の視線ベクトル周りの回転角基準 ;; Location fileであたえる 
  if n_elements(Na_deg_mode) eq 0 then Na_deg_mode = 'Fix' ;; <= Y軸方向=天球北
  ; Na_deg_mode = 'Free' ;; <= 観測者基準。観測者から見て"真上方向"が画像の上になる

  ;; 画像上を北にするか否か ;; Location fileであたえる 
  ; North_up = 'Y'
  if n_elements(North_up) eq 0 then North_up = 'N'
  
  ;; Lat Lon線を引くか否か ;;
  latlonline_to_sim_image = 1

  ;; 任意の回転角度 ;; North_upがYのとき無視される。
  if n_elements(image_rot_angle_deg) eq 0 then image_rot_angle_deg=0.

  ;; 画像上における金星中心の位置 ;;
  ;; 観測画像を解析するなどして自分で用意する
  if n_elements(cx_p) eq 0 then begin
    cx_p = 0
  endif
  if n_elements(cy_p) eq 0 then begin
    cy_p = 0
  endif
  ;cx_p = 0 ;; [pixel] 0=中心
  ;cy_p = 0 ;; [pixel] 0=中心

  ;; シーイング ;; 0より大きい値を入れるとぼやけ具合を計算して画像をぼやけさせる
  FWHM = 1. ;; [pixel] ;; <= arcsecにしたほうが良い？

  ;;
  ;; Display topograhpy?
  ;;
  if n_elements(geo_ifname) eq 0 then begin
    topography_plot = 0l ;; No topography
  endif else begin
    ;topography_plot = 1l ;; Whole disk
    topography_plot = 2l ;; Night side + Dayside enhanced
    ;topography_plot = 3l ;; Dayside only    
  endelse


  common obs_name,TARGET_BODY, OBS_BODY
  if n_elements(target_name) eq 0 then begin
    TARGET_BODY = 'Venus'
  endif else begin
    TARGET_BODY = strtrim(Target_name,2)
  endelse
  OBS_BODY = 'Earth'
  
  ;; SPICEを使って各種ジオメトリ情報を得る。
  ;; 計算結果はCommonで定義した変数の中に入る
  common obs_geo,Azi_Lv,Ele_Lv,ang_diam,lam $
                ,Na_deg_J2000, EV_distance, VS_distance $
                ,sub_solar_lon,sub_solar_lat,sub_earth_lon,sub_earth_lat $
                ,sub_venus_lat,angle_sun_venus_earth,angle_venus_earth_sun
                
  error= geometry_calculation_with_spice($
         fldname_kernel, metakr $
         ,utctim, obs_lat, obs_lon, obs_alt)

  if error eq -1 then begin
    print,"Geometry Calculation error..."
    return
  endif

  ;; 金星のNorth pole angle ;;
  ;; 天球北からの角度 ;;Horizons表記に準拠
  Na_deg_J2000_from_north = Na_deg_J2000+270.
  if (Na_deg_J2000_from_north gt 360.) then Na_deg_J2000_from_north -= 360.
  if (Na_deg_J2000_from_north gt 180.) then Na_deg_J2000_from_north -= 360.

  ;; 観測者座標系の真上からの角度 ;;
  Na_deg_obs_from_north = lam*180./!dpi+270.
  if (Na_deg_obs_from_north gt 360.) then Na_deg_obs_from_north -= 360.
  if (Na_deg_obs_from_north gt 180.) then Na_deg_obs_from_north -= 360.

  ;;;;;;;;;;;;;;;
  ;; 構造体の作成 ;;
  ;;;;;;;;;;;;;;;
  tags = ['date', 'obs_latitude', 'obs_longitude', 'obs_altitude' $
    ,'trg_azimuth', 'trg_elevation', 'trg_airmass', 'trg_diameter' $
    ,'trg_distance', 'sun_distance' $
    ,'J2000_north_azimuth','view_north_azimuth' $
    ,'ssl_latitude', 'ssl_longitude', 'ssc_latitude', 'ssc_longitude' $
    ,'sst_latitude', 'phase_angle', 'dif_angle']
  obs_geo_struct = create_struct(tags $
    ,utctim,obs_lat*180./!dpi,obs_lon*180./!dpi, obs_alt $
    ,Azi_Lv*180./!dpi,Ele_Lv*180./!dpi,1./cos(!dpi/2.-Ele_Lv),ang_diam*180./!dpi*3600. $
    ,EV_distance, VS_distance $
    ,Na_deg_J2000_from_north,Na_deg_obs_from_north $
    ,sub_solar_lat*180./!dpi,sub_solar_lon*180./!dpi,sub_earth_lat*180./!dpi,sub_earth_lon*180./!dpi $
    ,sub_venus_lat*180./!dpi,angle_sun_venus_earth*180./!dpi,angle_venus_earth_Sun*180./!dpi $
    )

  ;;;;;;;;;;;;;;;;;;;;;
  ;; IDL コンソールに出力 ;;
  ;;;;;;;;;;;;;;;;;;;;;
  print, "-- Ground based observation --"
  print,"Date: ",utctim
  print," Observer Lat: ", obs_lat*180./!dpi , " [deg]"
  print," Observer Lon: ", obs_lon*180./!dpi , " [deg]"
  print," Observer Alt: ", obs_alt, " [km]"
  print, "From Observer --"
  print," Azimuth  :     ", Azi_Lv*180./!dpi, " [deg]"
  print," Elevation:     ", Ele_Lv*180./!dpi, " [deg]"
  print," Airmass  :     ", 1./cos(!dpi/2.-Ele_Lv)
  print," Angular Diam:  ",ang_diam*180./!dpi*3600., " [sec]"
  print,"In images --"
  print," North azimuth from Celestial North: ",Na_deg_J2000_from_north, " [deg]"
  print," North azimuth from Observer Viewing Perp: ",Na_deg_obs_from_north, " [deg]"
  print,"On ", TARGET_BODY, "--"
  print," Sub solar Lat: ",sub_solar_lat*180./!dpi, " [deg]"
  print," Sub solar Lon: ",sub_solar_lon*180./!dpi, " [deg]"
  print," Sub earth Lat: ",sub_earth_lat*180./!dpi, " [deg]"
  print," Sub earth Lon: ",sub_earth_lon*180./!dpi, " [deg]"
  print,"ETC --"
  print," Sub ",TARGET_BODY," Lat (on earth): ",sub_venus_lat*180./!dpi, " [deg]"
  print," Phase angle of Sun - ", TARGET_BODY," - Earth: ", angle_sun_venus_earth*180./!dpi, " [deg]"
  print," Angle between ", TARGET_BODY," - Earth - Sun : ", angle_venus_earth_Sun*180./!dpi, " [deg]"
  print," Distance from ", OBS_BODY, " To ", Target_Body,EV_distance, " [km]"

  ;;;;;;;;;;;;;;;;;;;
  ;; 指定ファイルに出力 ;;
  ;;;;;;;;;;;;;;;;;;;
  if n_elements(ofname_log) ne 0 then begin
    openw,2,ofname_log
    printf,2, "-- Ground based observation --"
    printf,2,"Date: ",utctim
    printf,2, "From Observer --"
    printf,2," Azimuth  :     ", Azi_Lv*180./!dpi, " [deg]"
    printf,2," Elevation:     ", Ele_Lv*180./!dpi, " [deg]"
    printf,2," Angular Diam:  ",ang_diam*180./!dpi*3600., " [sec]"
    printf,2,"In images --"
    printf,2," North azimuth from Observer Image Y axis: ",Na_deg_obs_from_north, " [deg]"
    printf,2," North azimuth from Celestial North: ",Na_deg_J2000_from_north, " [deg]"
    printf,2," Image ratation angle: ",image_rot_angle_deg, " [deg]"
    printf,2,"On ", TARGET_BODY, "--"
    printf,2," Sub solar Lat: ",sub_solar_lat*180./!dpi, " [deg]"
    printf,2," Sub solar Lon: ",sub_solar_lon*180./!dpi, " [deg]"
    printf,2," Sub earth Lat: ",sub_earth_lat*180./!dpi, " [deg]"
    printf,2," Sub earth Lon: ",sub_earth_lon*180./!dpi, " [deg]"
    printf,2,"ETC --"
    printf,2," Sub ",TARGET_BODY," Lat (on earth): ",sub_venus_lat*180./!dpi, " [deg]"
    printf,2," Phase angle of Sun - ", TARGET_BODY," - Earth: ", angle_sun_venus_earth*180./!dpi, " [deg]"
    printf,2," Angle between ", TARGET_BODY," - Sun: ", angle_venus_earth_Sun*180./!dpi, " [deg]"
    printf,2," Distance from ", OBS_BODY, " To ", Target_Body,EV_distance, " [km]"
    close,2
  endif
  
  ;;;;;;;;;;
  ;; おまけ ;;
  ;;;;;;;;;;
  ;; 計算結果を受けてプロットする ;; Structureに渡す？
  ssl_lon_deg = sub_solar_lon*180./!dpi
  ssl_lat_deg = sub_solar_lat*180./!dpi
  ssc_lon_deg = sub_earth_lon*180./!dpi
  ssc_lat_deg = sub_earth_lat*180./!dpi
  V_distance = EV_distance

  ang_diam_deg = ang_diam*180./!dpi

  ;;;;;;;;;;;;;;;;;;;;;;
  ;; 惑星の地軸向きの設定 ;;
  ;;;;;;;;;;;;;;;;;;;;;;
  if Na_deg_mode ne 'Free' and Na_deg_mode ne 'Fix' then begin
    print,'Invalid North deg mode: ', Na_deg_mode
    print,'North deg mode to be Fix'
    Na_deg_mode = 'Fix' 
  endif

  if Na_deg_mode eq 'Free' then begin
    if North_up eq 'Y' then begin
      image_rot_angle_deg = 90.-lam*180./!dpi + image_rot_angle_deg
      if (image_rot_angle_deg gt 360.) then image_rot_angle_deg -= 360.
      if (image_rot_angle_deg gt 180.) then image_rot_angle_deg -= 360.
    endif

    Na_deg = lam*180./!dpi + image_rot_angle_deg
  endif else if Na_deg_mode eq 'Fix' then begin
    if North_up eq 'Y' then begin
      image_rot_angle_deg = 360.-Na_deg_J2000_from_north + image_rot_angle_deg
      if (image_rot_angle_deg gt 360.) then image_rot_angle_deg -= 360.
      if (image_rot_angle_deg gt 180.) then image_rot_angle_deg -= 360.
    endif

    Na_deg = Na_deg_J2000 + image_rot_angle_deg
  endif

  ;;;;;;;;;;;;
  ;; Output ;;
  ;;;;;;;;;;;;
  ;; 解像度 ;;
  dtheta = Tel_res_sec/3600.*!dpi/180.

  if keyword_set(noplot) then begin
    ;; No plotの場合直接subTS_plot_exを呼び出しジオメトリだけ計算して終了 ;;
    cx = cx_p*dtheta
    cy = cy_p*dtheta
    subTS_plot_ex,ssl_lat_deg,ssl_lon_deg,ssc_lat_deg,ssc_lon_deg $
      ,NA_deg,ang_diam_deg,dtheta,ccd_geo $
      ,reverse_x= reverse_x, reverse_y=reverse_y $
      ,pix_x=pix_x,pix_y=pix_y,cx=cx,cy=cy,/noplot

  endif else if keyword_set(latlonline_only) then begin

    cx = cx_p*dtheta
    cy = cy_p*dtheta

    subTS_plot_ex,ssl_lat_deg,ssl_lon_deg,ssc_lat_deg,ssc_lon_deg $
      ,NA_deg,ang_diam_deg,dtheta,ccd_geo $
      ,reverse_x= reverse_x, reverse_y=reverse_y $
      ,pix_x=pix_x,pix_y=pix_y,cx=cx,cy=cy,/noplot

    subTS_plot_ex_eng,ssl_lat_deg,ssl_lon_deg,ssc_lat_deg,ssc_lon_deg $
      ,NA_deg,ang_diam_deg,dtheta,ccd_geo $
      ,pix_x=pix_x,pix_y=pix_y $
      ,cx=cx,cy=cy $
      ,reverse_x= reverse_x, reverse_y=reverse_y $
      ,/sub_ssl_plot,/line_only,/noerase,line_color=196

    plot,dindgen(10),/nodata,/noerase,xs=5,ys=5,title=utctim,charsize=1.

  endif else begin
    ;; Plotの場合、plot_simulated_imageから呼び出して出力 ;;
    tmp_position = !P.position
    !P.position = [0.15,0.15,0.95,0.95]
    
    ;ssc_lat_deg = 5d

    ;; 緯度経度線付をプロット ;;
    plot_simulated_image,ssl_lon_deg,ssl_lat_deg,ssc_lon_deg,ssc_lat_deg $
                          ,Na_deg,ang_diam_deg, Tel_res_sec, FWHM = FWHM $
                          ,ccd_geo=ccd_geo, sim_image = sim_image $
                          ,reverse_x= reverse_x, reverse_y=reverse_y $
                          ,pix_x = pix_x, pix_y = pix_y,cx_p = cx_p,cy_p = cy_p $
                          ,latlonline_to_sim_image=latlonline_to_sim_image,/sub_ssl_plot

    ;; 線なしをプロット ;;
    ;plot_simulated_image,ssl_lon_deg,ssl_lat_deg,ssc_lon_deg,ssc_lat_deg $
    ;  ,Na_deg,ang_diam_deg, Tel_res_sec, FWHM = FWHM $
    ;  ,ccd_geo=ccd_geo, sim_image = sim_image $
    ;  ,pix_x = pix_x, pix_y = pix_y,cx_p = cx_p,cy_p = cy_p ;$
    ;  ;,/reverse_y

    plot,dindgen(10),/nodata,/noerase,xs=5,ys=5,title=utctim,charsize=1.0

    !P.position = tmp_position    

  endelse

  ;;
  ;; Important: Flipping ccd_geo array in X and Y direction for display,
  ;; because subTS_plot_ex provides geometry with a (X, Y) reverse condition due to considering CCD projection.
  ;;
  ccd_geo = reverse(reverse(ccd_geo,1),2)

  ;;
  ;; 地形を表示するか否か
  ;;
  if topography_plot ne 0 then begin

    ;geo_ifldname = './'
    ;geo_ifname = geo_ifldname + 'Venus_Magellan_Topography_Global_4641m_gapfilled_v02.tif'

    ;; data読み込み (金星中心からの距離になる) ;;
    geo_data = read_tiff(geo_ifname) + 6051000.
    
    ;; 標高変換 ;;
    tmp_topo = (double(geo_data) - 6051800)/1000.

    topo_siz = size(tmp_topo)
    topo_res_x = topo_siz[1]
    topo_res_y = topo_siz[2]
    grid_per_deg = topo_res_x / 360.

    ;; 0.5degで平滑化 ;;
    ;tmp_topo = g_smooth(tmp_topo,grid_per_deg/2.,/cylindrical)

    ;; Y反転 ;;
    tmp_topo = reverse(tmp_topo,2)

    ;; 0degを左端に ;;
    tmp_topo = shift(tmp_topo, topo_res_x/2,0)

    tmp_topo = tmp_topo > (-3.)
    ;tvscl,congrid(tmp_topo,!d.x_size,!d.y_size) < 6.
    ;stop

    ;; 経度方向のループ対策 ;;
    topo = dblarr(topo_res_x+2,topo_res_y)
    topo[1:topo_res_x,*]=tmp_topo
    topo[0,*]=tmp_topo[topo_res_x-1,*]
    topo[topo_res_x+1,*]=tmp_topo[0,*]

    ;; 座標出し ;;
    ix = (ccd_geo[*,*,0]*180./!dpi)
    ;; 負の部分を正にしておく ;;
    ix_pos = where(ix lt 0 and ix gt -180.)
    ix[ix_pos] += 360.
    ix = ix * grid_per_deg > 0.

    iy = (ccd_geo[*,*,1]*180./!dpi)
    iy = (iy + 90.)*grid_per_deg > 0

    ;; 回転いる ? ;;
    ;ix = rot(ix,180)
    ;iy = rot(iy,180)
    ;;;;;;;;;;;;;;;;;;;
    ;; Interpolation ;;
    ;;;;;;;;;;;;;;;;;;;
    ;; Nearest neighbor ;;
    ;ccd_topo = topo[long(ix+1.+0.5-0.5),long(iy+0.5-0.5)]

    ;; Bilinear ::
    ccd_topo = bilinear(topo,ix+1.,iy) ;; <= loop分を意識, 0.5degずれは考慮する?

    max_topo = 6.
    min_topo = -3.

    ;ccd_topo *= -1.
    ;max_topo = 3.
    ;min_topo = -6.

    ccd_topo = ccd_topo*(ccd_geo[*,*,1] ge -!dpi/2.)
    inc = cos(ccd_geo[*,*,2])
    ;; 回転いる ?
    ;inc = rot(inc,180)

    topo_n_pos = where(ccd_geo[*,*,1] lt -!dpi/2.)
    ccd_topo[topo_n_pos] = min(ccd_topo)

    ;; Day side ;;
    ccd_topo_day = ccd_topo ;* (inc ge 0)
    ccd_topo_day[topo_n_pos] = min(ccd_topo)
    ;ccd_topo_day[where(inc lt 0)] -= 2; min_topo

    ;; Night side ;;
    ccd_topo_night = ccd_topo ;* (inc lt 0)
    ccd_topo_night[topo_n_pos] = min(ccd_topo)
    ccd_topo_night[where(inc ge 0)] += 2

    ;; 出力する場合
    if keyword_set(noplot) eq 0 then begin
      px = !x.window*!d.x_vsize
      py = !y.window*!d.y_vsize
      sx = px(1)-px(0)+1
      sy = py(1)-py(0)+1

      ;; 地形出力 ;; 2: Night(+Day enhanced), 3; Day, 1: Both
      if topography_plot eq 2 then begin
        output_geo = congrid(ccd_topo_night,sx,sy)

      endif else if topography_plot eq 3 then begin
        output_geo = congrid(ccd_topo_day,sx,sy)

      endif else begin
        output_geo = congrid(ccd_topo,sx,sy)
      endelse

      output_geo[0,0] = max_topo
      output_geo[0,1] = min_topo
      output_geo = (output_geo < output_geo[0,0]) > output_geo[0,1]
      tvscl,output_geo,px[0],py[0]

      ;; 枠を描く ;;
      tmp_position = !P.position
      !P.position = [0.15,0.15,0.95,0.95]
      ;    subTS_plot_ex,ssl_lat_deg,ssl_lon_deg,ssc_lat_deg,ssc_lon_deg $
      ;      ,NA_deg,ang_diam_deg,dtheta,ccd_geo $
      ;      ,pix_x=pix_x,pix_y=pix_y $
      ;      ,cx=cx,cy=cy $
      ;      ,reverse_x= reverse_x, reverse_y=reverse_y $
      ;      ,/sub_ssl_plot,/line_only,/noerase,line_color=196

      !P.position = tmp_position

    endif

  endif

  print,"Plot Keywords --"
  print," Target name: ", target_body
  print," North Azimuth Direction Mode: '", Na_deg_mode, "'"
  print," North UP: '",North_up,"'"
  print," Image rotation angle (CCW): ",image_rot_angle_deg, " [deg]"
  print," Angular Semi Diam:  ",ang_diam/2./dtheta, " [pix]"
  if ~keyword_set(noplot) and ~keyword_set(latlonline_only) then print," Bright Max: ",max(sim_image)

  print, "==="

end


