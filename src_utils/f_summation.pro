;;
;;
;;
function f_summation,array,tx,ty
  ; array: 入力配列, tx, ty: 平滑化に用いるボックスの大きさ
  siz_array = size(array)
  sx = siz_array[1]
  sy = siz_array[2]
  s_array = dblarr(sx,sy)
  s_array[0,0]=total(array[0:tx-1,0:ty-1])
  
  for i=1,sx-tx-1,1 do begin
    s_array[i,0]=s_array[i-1,0]+total(array[i+tx-1,0:ty-1])-total(array[i-1,0:ty-1])
  endfor
  for j=1,sy-ty-1,1 do begin
    s_array[0,j]=s_array[0,j-1]+total(array[0:tx-1,j+ty-1])-total(array[0:tx-1,j-1])
  endfor
  tmp_array = dblarr(sx)
  for i=0,sx-1,1 do begin
    tmp_array[i]=total(array[i,0:ty-1])
  endfor
  
  for j=1,sy-ty-1,1 do begin
    for i=0,sx-1,1 do begin
      tmp_array[i]=tmp_array[i]+array[i,j+ty-1]-array[i,j-1]
    endfor
    for i=1,sx-tx-1,1 do begin
      s_array[i,j]=s_array[i-1,j]+(tmp_array[i+tx-1]-tmp_array[i-1])
    endfor
  endfor
  
  ;;;Debug
  ;     for j=0,sy-ty-1,1 do begin
  ;             for i=0,sx-tx-1,1 do begin
  ;                     ;tmp=total(array[i:i+tx-1,j:j+ty-1])-s_array[i,j]
  ;                     ;print,i,j,tmp_x,s_array[i,j]
  ;                     tmp=total(array[i:i+tx-1,j:j+ty-1])
  ;                     s_array[i,j]=tmp
  ;             endfor
  ;     endfor
  
  return,s_array
end
