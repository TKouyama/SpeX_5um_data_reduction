;;
;; Linear interpolation
;;

;; History
;; 2017.09.27: add a keyword: remove_cycle

function d_shift,data,dx,dy,remove_cycle=remove_cycle, remove_value=remove_value

  siz = size(data)
  nx = siz[1]
  ny = siz[2]
  cx = dindgen(nx)
  cy = dindgen(ny)
  
  if dx ge 0 then begin
    i_dx = long(dx)
  endif else begin
    i_dx = long(dx)-1
  endelse
  
  if dy ge 0 then begin
    i_dy = long(dy)
  endif else begin
    i_dy = long(dy)-1
  endelse
  
  tmp_data = shift(data,i_dx,i_dy)

  ;; Cyclicになった部分はremove valudeで埋める ;;
  if n_elements(remove_value) eq 0 then begin
    remove_value = 0.
  endif
  
  if keyword_set(remove_cycle) then begin
    if i_dx gt 0 then begin
      tmp_data[0:i_dx-1,*] = remove_value[0]
    endif else begin
      tmp_data[nx+i_dx-1:nx-1,*] = remove_value[0]
    endelse

    if i_dy gt 0 then begin
      tmp_data[*,0:i_dy-1] = remove_value[0]
    endif else begin
      tmp_data[*,ny+i_dy-1:ny-1] = remove_value[0]
    endelse
  endif

  d_dx = dx-i_dx
  d_dy = dy-i_dy

  ax = d_dx
  ay = d_dy
  bx = 1.-d_dx
  by = 1.-d_dy
  tmp_data1 = shift(tmp_data,1,0)
  tmp_data2 = shift(tmp_data,0,1)
  tmp_data3 = shift(tmp_data,1,1)
  d_shift_data = bx*by*tmp_data+ax*by*tmp_data1+bx*ay*tmp_data2+ax*ay*tmp_data3

;  d_shift_data = tmp_data*0
;  for j=1,ny-2,1 do begin
;    for i=1,nx-2,1 do begin
;      d_shift_data[i,j] = tmp_data[i,j]*bx*by+tmp_data[i-1,j]*ax*by $
;                      +tmp_data[i,j-1]*bx*ay+tmp_data[i-1,j-1]*ax*ay
;    endfor
;  endfor
  
  return,d_shift_data
end
