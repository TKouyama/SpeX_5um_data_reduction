;;
;; gauss_smooth  を高速処理で行うプログラム
;; アルゴリズムはx方向ガウシアンスムース実行後にy方向のガウシアンスムースを行う
;; 

;;
;; limb: ignore value (or 0)をハイパスに使わない
;; cylindrical: x方向にデータがつながっていることを仮定(緯度経度展開図など)
;;

;;
;; Author: Toru Kouyama (t.kouyama@aist.go.jp)
;;

function g_smooth,data,sigma,limb=limb,cylindrical=cylindrical $
                 ,ignore_value=ignore_value,mergin_w=mergin_w
  if sigma le 0. then begin
    return, data
  endif

  im_siz = size(data)

  if n_elements(mergin_w) eq 0 then begin
    mergin_w = 6d ;; <= 6sigma まで取るということ    
  endif

  ;; 配列サイズの制約 ;;
  weight_w = long(mergin_w*sigma+0.5) > 8
  weight_w = weight_w < (min([im_siz[1],im_siz[2]]) - 3.)

  weight_y = dblarr(1,weight_w+1)

  distance_y = (dindgen(weight_w+1) - weight_w/2)

  weight_y[0,*]= exp(-distance_y^2./(2.*sigma^2.))
  ;weight_y[0,*]= 1./exp(-distance_y^2./(2.*sigma^2.))
  ;weight_y[0,*]= exp(-distance_y^2./(sigma^2.))

  ;plot,weight_y
  ;stop

  ;; work領域を使う ;;
  if keyword_set(cylindrical) then begin
    tmp_i_siz = size(data)
    tmp_im_x = tmp_i_siz[1]
    tmp_im_y = tmp_i_siz[2]

    tmp_data = dblarr(tmp_im_x+2*weight_w, tmp_im_y)
    tmp_data[0:weight_w-1,*] = data[tmp_im_x-1-(weight_w-1):tmp_im_x-1,*]
    tmp_data[weight_w:weight_w+tmp_im_x-1,*] = data 
    tmp_data[tmp_im_x+weight_w:tmp_im_x-1+2*weight_w,*] = data[0:weight_w-1,*]
  endif else begin
    tmp_data = data
  endelse
  
  i_siz = size(tmp_data)
  im_x = i_siz[1]
  im_y = i_siz[2]

  sumdata = dblarr(im_x,im_y)
  sumcount = dblarr(im_x,im_y)
  tmp_g = dblarr(1,im_y)
  tmp_gV = dblarr(1,im_y)
  
  ;;
  ;; 境界を強調する
  ;;
  if keyword_set(limb) then begin
    if n_elements(ignore_value) eq 0 then ignore_value=0.

    ;d_check = 1.*(tmp_data gt ignore_value)
    d_check = 1.*(tmp_data ne ignore_value)

    for gi=0,im_x-1,1 do begin
      tmp_g[0,*] = tmp_data[gi,*] * (d_check[gi,*] ne 0)
      tmp_gV[0,*] = d_check[gi,*]
      sumdata[gi,*]  = convol(tmp_g,weight_y,/EDGE_TRUNCATE)*(tmp_gV[0,*] gt 0.)
      sumcount[gi,*] = convol(1.*(tmp_gV[0,*] gt 0.),weight_y,/EDGE_TRUNCATE)*(tmp_gV[0,*] gt 0.)
    endfor

    tmp_avedata = sumdata
    tmp_count = sumcount
  
    tmp_g = dblarr(1,im_x)
    tmp_gV = dblarr(1,im_x)
    ;tmp_g = dblarr(1,im_y)
    ;tmp_gV = dblarr(1,im_y)
    for gj=0,im_y-1,1 do begin
      tmp_g[0,*] = tmp_avedata[*,gj]
      tmp_gV[0,*] = tmp_count[*,gj]
      sumdata[*,gj]  = convol(tmp_g,weight_y,/EDGE_TRUNCATE)*(tmp_gV[0,*] gt 0.)
      sumcount[*,gj] = convol(tmp_gV,weight_y,/EDGE_TRUNCATE)*(tmp_gV[0,*] gt 0.)
    endfor

  endif else begin

    for gi=0,im_x-1,1 do begin
      tmp_g[0,*] = tmp_data[gi,*]
      tmp_gV[0,*] = replicate(1d,im_y)
      sumdata[gi,*]  = convol(tmp_g,weight_y,/EDGE_TRUNCATE)
      sumcount[gi,*] = convol(tmp_gV,weight_y,/EDGE_TRUNCATE)
    endfor
    
    tmp_avedata = sumdata
    tmp_count = sumcount
    
    tmp_g = dblarr(1,im_x)
    tmp_gV = dblarr(1,im_x)
    for gj=0,im_y-1,1 do begin
      tmp_g[0,*] = tmp_avedata[*,gj]
      tmp_gV[0,*] = tmp_count[*,gj]
      sumdata[*,gj]  = convol(tmp_g,weight_y,/EDGE_TRUNCATE)
      sumcount[*,gj] = convol(tmp_gV,weight_y,/EDGE_TRUNCATE)
    endfor
    
  endelse

  ;; 規格化 ;;
  smooth_pos = where(sumcount gt 0)
  smooth_data = dblarr(im_x,im_y)
  smooth_data[smooth_pos] = sumdata[smooth_pos]/sumcount[smooth_pos]
  ;smooth_data = (sumdata/(sumcount > 1))
  
  if n_elements(ignore_value) ne 0 then begin
    n_smooth_pos = where(sumcount eq 0)
    smooth_data[n_smooth_pos] = ignore_value
  endif

  if keyword_set(cylindrical) then begin
    smooth_data = smooth_data[weight_w:weight_w+tmp_im_x-1,*]
  endif

  return, smooth_data
end