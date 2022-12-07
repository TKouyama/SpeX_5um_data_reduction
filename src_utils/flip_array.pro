;;
;; 2次元行列を必要な形にひっくり返す関数 ;;
;;
function flip_array,array
  siz_array = size(array)
  sx = siz_array[1]
  sy = siz_array[2]

  tmp_array = array[1:sx-1,*]
  rev_tmp_array = reverse(tmp_array)

  f_array = array
  f_array[1:sx-1,0:sy-1]=rev_tmp_array

  tmp_array = f_array[*,1:sy-1]
  rev_tmp_array = reverse(tmp_array,2)

  f_array[0:sx-1,1:sy-1]=rev_tmp_array

  return,f_array
end
