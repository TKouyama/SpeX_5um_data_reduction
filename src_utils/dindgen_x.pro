;;
;; 一様にふえる座標を作る
;;
function dindgen_x,nx,ny,y=y
  array = dblarr(nx,ny)

  if keyword_set(y) then begin

    for j=0, ny-1, 1 do begin
      array[*,j] = double(j)
    endfor

  endif else begin    
    for i=0, nx-1, 1 do begin
      array[i,*] = double(i)
    endfor
  endelse

  
  return,array
end
