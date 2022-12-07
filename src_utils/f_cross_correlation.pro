;;
;; Fast correlation ;;
;;
;; fftを使った相関計算
function f_cross_correlation,array,array2,csx,csy,cyclic=cyclic
  siz = size(array)
  tsx = siz[1] ;; template size x
  tsy = siz[2] ;; template size y

  siz2 = size(array2)
  ssx = siz2[1] ;; search area size x
  ssy = siz2[2] ;; search area size y

  ;; Array1 の平均値と二乗和
  m_array = mean(array,/double)
  s_array = total(array^2.d)

  ;; Array2 の平均値と二乗和
  m_array2 = f_summation(array2,tsx,tsy)
  m_array2 /= tsx*tsy
  s_array2 = f_summation((array2)^2.d,tsx,tsy)

  ;; Covarriance の計算
  tmp_array = dblarr(siz2[1],siz2[2])
  tmp_array[0:siz[1]-1,0:siz[2]-1]=array ;; テンプレートサイズ以上のところは0

  flip_array2 = flip_array(array2)

  tmp_cov_fft = fft(tmp_array,/double)*fft(flip_array2,/double)
  tmp_cov_r = real_part(fft(tmp_cov_fft,/double,/inverse))*ssx*ssy
  cor_v = flip_array(tmp_cov_r) - tsx*tsy*m_array*m_array2

  ;; 分母に来る係数 ;;
  cor_a = sqrt(s_array-tsx*tsy*m_array^2.d)
  cor_b = sqrt(s_array2-tsx*tsy*m_array2^2.d)

  ;; 相関の計算 ;;
  if keyword_set(cyclic) then begin    
    csx = ssx
    csy = ssy
  endif else begin
    csx = ssx-tsx
    csy = ssy-tsy
  endelse

  result=dblarr(csx,csy)
  tmp_cor_v = cor_v[0:csx-1,0:csy-1]
  tmp_cor_b = cor_b[0:csx-1,0:csy-1]

  pos = where(tmp_cor_b ne 0) ;; 0で割ることを回避するため、分母が0以外のところを探す
  if pos[0] ne -1 then begin
    result[pos] = tmp_cor_v[pos]/cor_a/tmp_cor_b[pos]
  endif

  ;; Debug
  ;  tmp_array = array[0:tsx-1,0:tsy-1]
  ;  tmp_array2 = dblarr(ssx,ssy)
  ;  result2 = dblarr(csx,csy)
  ;  for j_shift=0,csy-1,1 do begin
  ;    for i_shift=0,csx-1,1 do begin
  ;      tmp_array2 = array2[i_shift:i_shift+tsx-1,j_shift:j_shift+tsy-1]
  ;      result2[i_shift,j_shift] = correlate(tmp_array,tmp_array2)
  ;      ;print,i_shift,j_shift,result[i_shift,j_shift]
  ;      ;stop
  ;    endfor
  ;    print,j_shift
  ;  endfor
  ;  print,max(result2-result),min(result2-result);max(tmp_cov_r),max(cor_a*cor_b)

  return,result
end



;;
;;
;;
function f_cross_correlation_bug,array1, array2, tx,ty

  s_array1 = f_summation(double(array1),tx,ty)/double(tx*ty)
  s2_array1 = f_summation(double(array1)^2d,tx,ty)

  ;s_array1 = total(double(array1))/double(tx*ty)
  ;s2_array1 = total(double(array1)^2d)

  s_array2 = f_summation(double(array2),tx,ty)/double(tx*ty)
  s2_array2 = f_summation(double(array2)^2d,tx,ty)

  ;; Covariance ;;
  cov_image = f_summation(double(array1)*double(array2),tx,ty)

  C1 = sqrt(s2_array1 - s_array1^2d * double(tx*ty))
  C2 = sqrt(s2_array2 - s_array2^2d * double(tx*ty))

  C0 = cov_image - s_array2*s_array1*tx*ty
  cor_image = C0*0d
  pos = [where(C0 ne 0)]
  cor_image[pos] = C0[pos]/(C1[pos]*C2[pos])
  cor_image = shift(cor_image,tx/2,ty/2)

  return,cor_image  

end
