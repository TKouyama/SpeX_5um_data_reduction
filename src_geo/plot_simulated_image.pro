;;; qtaxis.proとqtcompose.proとqtmult.proとsubTS_plot.proを先にコンパイルしておく

;;
;;  SeeingによるPoint spread functionをガウシアンで表現
;;
function simulate_seeing, HWHM
  nr = (long(HWHM*6+0.5)+1)>7
  if nr mod 2 eq 0 then nr+=1
  
  psf = dblarr(nr,nr)
  psf_c = dblarr(nr,nr)
  sigma2 = ((HWHM/sqrt(2.*alog(2.))))^2.
  for j=0,nr-1, 1. do begin
    y = double(j)-double(nr-1)/2.
    for i=0, nr-1,1 do begin
      x =  double(i)-double(nr-1)/2.
      psf[i,j] = exp(-(x^2.+y^2.)/2./sigma2)
    endfor
  endfor
  return,psf
end

pro plot_simulated_image,ssl_lon_deg,ssl_lat_deg,ssc_lon_deg,ssc_lat_deg $
                   ,Na_deg,ang_diam_deg, Tel_res_sec,sim_image=sim_image $
                   ,distance = distance, FWHM = FWHM $
                   ,ccd_geo=ccd_geo $
                   ,pix_x = pix_x, pix_y = pix_y, cx_p = cx_p,cy_p = cy_p $
                   ,reverse_x = reverse_x, reverse_y = reverse_y $
                   ,latlonline_to_sim_image=latlonline_to_sim_image $
                   ,sub_ssl_plot=sub_ssl_plot, unset_window=unset_window

  ;; Test実行用フラグ。下をtest_run=1とすると単体で動く
  test_run = 0

  if test_run eq 1 then begin
  ;; 必要なパラメータを設定, SPICEからもらってきてもよい ;;
    ssl_lon_deg = -43.890 ; 太陽直下経度
    ssl_lat_deg = 2.601 ; 太陽直下緯度
    ssc_lon_deg = 59.393 ; 衛星直下経度
    ssc_lat_deg = 3.554 ; 衛星直下緯度
    Na_deg = 43.9 ; 画像のx軸方向から何度の方向に北極があるか。 90で真上方向

    ;; Akatsuki: 2015 12/07 ;;
    ssl_lon_deg =  141.476; 太陽直下経度
    ssl_lat_deg = -2.63397 ; 太陽直下緯度
    ssc_lon_deg = 94.8273 ; 衛星直下経度
    ssc_lat_deg = -5.52871 ; 衛星直下緯度
    Na_deg = 97.4763 ; 画像のx軸方向から何度の方向に北極があるか。 90で真上方向

    ;; Akatsuki: 2015 12/09 ;;
    ssl_lon_deg =  148.824; 太陽直下経度
    ssl_lat_deg = -2.63804 ; 太陽直下緯度
    ssc_lon_deg = 74.7316 ; 衛星直下経度
    ssc_lat_deg = -5.37198 ; 衛星直下緯度
    Na_deg = 99.4862 ; 画像のx軸方向から何度の方向に北極があるか。 90で真上方向

    ;; Akatsuki: 2016 05/06 ;;
    ;ssl_lon_deg =  148.824; 太陽直下経度
    ;ssl_lat_deg = -2.63804 ; 太陽直下緯度
    ;ssc_lon_deg = 74.7316 ; 衛星直下経度
    ;ssc_lat_deg = -5.37198 ; 衛星直下緯度
    ;Na_deg = 99.4862 ; 画像のx軸方向から何度の方向に北極があるか。 90で真上方向

    ;; 惑星の視直径 [deg], 直接指定する
    ;ang_diam_deg = 29.983/3600.
    ang_diam_deg =  (6371./342394.60)* 2. * 180./!dpi ;; <= Hayabusa2

    ; or 金星までの距離から計算する
    ;Distance = 1.630160e+06*2
    ;; 画像の解像度, 画像サイズ
    ;Tel_res_sec = 0.39; 0.39" <= MSIの解像度 [sec/pix]
    ;Tel_res_sec = 20.490; 20.490" ;; <= Hayabusa2/ONC [sec/pix]
    Tel_res_sec = 20.490; 20.490" ;; <= Hayabusa2/ONC [sec/pix]

    pix_x = 512; 659 ;; x画素数
    pix_y = 512; 494 ;; y画素数

    cx_p = 0 ;; 金星の中心がccdのどこに写っているか.(pixel) 0=中心
    cy_p = 0 ;; 金星の中心がccdのどこに写っているか.(pixel) 0=中心
  endif

  ;; 観測距離を指定した場合、そこから視直径を算出する
  if n_elements(distace) ne 0 then begin
    Rv = 6115.
    ang_diam_deg = asin(Rv/Distance)*2.*180./!dpi
  endif
  
  if n_elements(FWHM) eq 0 then begin
    FWHM = 0. ;; シーイング [pixel]
  endif

  if pix_x*pix_y gt 2100.^2. then begin
    ;; 計算破たん防止
    stop
    return
  endif  

  dtheta = Tel_res_sec/3600.*!dpi/180. ;; 画像の解像度 (radian)
  ;; Galileo/SSIの場合
  ;dtheta = 10.16e-6
  cx = cx_p*dtheta ;;(radian)
  cy = cy_p*dtheta ;;(radian)

  ;; パラメータ設定ここまで ;;

  ;; Plot windowの設定(適宜変更)
  ;xs = pix_x & ys = pix_y
  xs = 640 & ys = 640

  if keyword_set(unset_window) ne 1 then begin

    if FWHM gt 0. then window,2,xs=xs,ys=ys,xpos=700
    window,xs=xs,ys=ys,1 ;; 緯度経度線
    window,xs=xs,ys=ys,0 ;; 模擬画像
    device, decomposed=0
    !p.background = 255
    !p.color = 0
    loadct,0

  endif

  ;; subTS_plotを呼び出す, ccd_geoは計算結果が入っている(画像上での緯度経度入射角出射角)
  if FWHM eq 0. then begin
    wset,0
    subTS_plot_ex,ssl_lat_deg,ssl_lon_deg,ssc_lat_deg,ssc_lon_deg $
               ,NA_deg,ang_diam_deg,dtheta $
               ,ccd_geo,pix_x=pix_x,pix_y=pix_y $
               ,cx=cx,cy=cy,reverse_x=reverse_x,reverse_y=reverse_y

    ;; 惑星模擬画像のTest 出力 ;;
    inc = ccd_geo[*,*,2]
    bright = cos(inc) > 0

    output = dblarr(pix_x,pix_y,7)
    output[*,*,0:5] = ccd_geo
    output[*,*,6] = bright
    ;output = reverse(reverse(output,1),2) ;; ←x, yひっくり返すときに使う

    ;; Fits file 出力
    ;ofname = "C:\work\ENVI_IDL\sample\subaru_reduction\tmp.fits"
    ;writefits,ofname,float(output)

    bright = reverse(reverse(bright,1),2) ;; ccd_geoはx, yひっくり返って入っている
    s_bright = bright

  endif else if FWHM gt 0. then begin

    ;; 自力で Convolution
    tmp_pix_x = pix_x + FWHM*2.
    tmp_pix_y = pix_y + FWHM*2.

    subTS_plot_ex,ssl_lat_deg,ssl_lon_deg,ssc_lat_deg,ssc_lon_deg $
      ,NA_deg,ang_diam_deg,dtheta $
      ,ccd_geo,pix_x=tmp_pix_x,pix_y=tmp_pix_y $
      ,cx=cx,cy=cy,reverse_x=reverse_x,reverse_y=reverse_y,sub_ssl_plot=sub_ssl_plot
    
    ;; 惑星模擬画像のTest 出力 ;;
    output = dblarr(pix_x,pix_y,7)
    output[*,*,0:5] = ccd_geo[FWHM:FWHM+pix_x-1,FWHM:FWHM+pix_y-1,0:5]

    inc = ccd_geo[*,*,2]
    bright = cos(inc) > 0

    output[*,*,6] = bright[FWHM:FWHM+pix_x-1,FWHM:FWHM+pix_y-1]

    ;; 畳み込みを利用したpsfによる模擬画像の平滑化。
    ;; 数学的な内容はwikiなどを参照
    psf = simulate_seeing(FWHM)
    
    ;; 組み込みconvol関数で畳み込みするとき
    ;s_bright = convol(bright,psf)/total(psf)

    ;; PSF の配列サイズを金星画像と同じにする。拡張した部分の値には0を入れておく
    r_psf = dblarr(tmp_pix_x,tmp_pix_y) ;; <= Brightはちょっと大きい 
    siz_psf = size(psf)
    r_psf[0:siz_psf[1]-1,0:siz_psf[2]-1] = psf
    ;; 座標[0,0]にpsfのピークを持ってくる
    r_psf = shift(r_psf,-siz_psf[1]/2, -siz_psf[2]/2)
    ;; Kernel の規格化
    r_psf /= total(r_psf)

    ;; 本来ここで配列のフリッピングが必要だが
    ;; 今回はKernelがフリッピングに対して対象なので必要ない。

    ;; FFTを利用したConvolution
    fft_psf = fft(r_psf)
    fft_bright = fft(bright)
    fft_s_bright = fft_bright*fft_psf
    tmp_s_bright = real_part(fft(fft_s_bright,/inverse))*tmp_pix_x*tmp_pix_y ;*pix_x*pix_y
    s_bright = tmp_s_bright[FWHM:FWHM+pix_x-1,FWHM:FWHM+pix_y-1]

    s_bright = reverse(reverse(s_bright,1),2) ;; ccd_geoはx, yひっくり返って入っている

    ;; 自力 Deconvolution !! Noiseのない画像にだけ使える !!
    ;s_bright += randomn(seed, pix_x, pix_y)*1e-10
    ;fft_psf = fft(r_psf)
    ;fft_s_bright = fft(s_bright)
    ;fft_d_bright = fft_s_bright/fft_psf
    ;d_bright = real_part(fft(fft_d_bright,/inverse))/(pix_x*pix_y)
  endif
  
  ;; subTS_plot_ex.pro でreverse処理済
  ;if keyword_set(reverse_x) then s_bright = reverse(s_bright,1)
  ;if keyword_set(reverse_y) then s_bright = reverse(s_bright,2)

  ;; 外に計算結果を引き継ぐ ;;
  sim_image = s_bright

  ;;;;;;;;;;;;;;;;;
  ;; 模擬画像出力 ;;
  ;;;;;;;;;;;;;;;;;
  wset,1
  ;; 画像表示枠の大きさを得る
  contour,dindgen(100,100),/nodata,/noerase,xstyle=5,ystyle=5
  px = !x.window*!d.x_vsize
  py = !y.window*!d.y_vsize
  sx = px(1)-px(0)+1
  sy = py(1)-py(0)+1

  loadct,0,/silent
  erase

  tvscl,congrid(s_bright,sx,sy,/interp),px(0),py(0)  
  ;; 緯度経度線
  if keyword_set(latlonline_to_sim_image) then begin
    subTS_plot_ex,ssl_lat_deg,ssl_lon_deg,ssc_lat_deg,ssc_lon_deg $
               ,NA_deg,ang_diam_deg,dtheta $
               ,ccd_geo,pix_x=pix_x,pix_y=pix_y $
               ,cx=cx,cy=cy,reverse_x=reverse_x,reverse_y=reverse_y $
               ,/line_only,/noerase,line_color=196,sub_ssl_plot=sub_ssl_plot
  endif

  wset,0

  if FWHM gt 0. then begin
    ;window,xs=250,ys=250,2,xpos=700 ;; Point spread function
    wset,2
    erase
    ;shade_surf,psf,thick=1,SHADES=BYTSCL(psf, TOP = !D.TABLE_SIZE) 
    surface,psf,thick=1,/noerase
    wset,0

    ;; Deconvolution test
    ;window,xs=xs,ys=ys,3,xpos=700 ;; Deconvolution image
    ;wset,3
    ;loadct,0,/silent
    ;tvscl, d_bright
    ;loadct,39,/silent
    ;wset,0
  endif
 
  return
end