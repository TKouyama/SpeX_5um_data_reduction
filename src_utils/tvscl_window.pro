;;
;;
;;
pro tvscl_window,im_array,min_value, max_value, center=center,interp=interp, true=true, pos_n=pos_n
  if keyword_set(true) eq 0 then begin

    if n_elements(min_value) gt 0 or n_elements(min_value) gt 0 then begin
      tmp_out = congrid(im_array,!d.x_size,!d.y_size,center=center,interp=interp)
      if n_elements(min_value) gt 0 then begin
        tmp_out[0,0] = min_value
      endif

      if n_elements(min_value) gt 0 then begin
        tmp_out[0,1] = max_value
      endif

      tvscl,tmp_out > min_value < max_value        

    endif else begin      

      tvscl,congrid(im_array,!d.x_size,!d.y_size,center=center,interp=interp)
      
    endelse


  endif else begin
    ;; pngなどカラー表示 ;;
    tvscl,congrid(im_array,3,!d.x_size,!d.y_size,center=center,interp=interp),true=true      
  endelse
  
  return
end