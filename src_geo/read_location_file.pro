;;
;; Inputファイルから計算設定などよみこみ
;;
function read_location_file,ifname,obs_lat,obs_lon,obs_alt $
  ,fldname_kernel,metakr $
  ,Tel_res_sec,pix_x,pix_y $
  ,North_deg_mode, North_up $
  ,image_rot_angle_deg=image_rot_angle_deg $
  ,Target_name = Target_name $
  ,geo_ifname = geo_ifname

  openr,1,ifname, ERROR = err
  if err ne 0 then begin
    close,1
    return,err
  endif

  tmp_line=''
  while(not eof(1)) do begin
    readf,1,tmp_line
    ;print,tmp_line
    if strmid(tmp_line,0,1) eq '#' then begin
    endif else begin
      tmp_line = strcompress(tmp_line,/remove)
      result = strpos(tmp_line,'=')
      if result ne -1 then begin
        tmp_str = strmid(tmp_line,0,result)
        tmp_len = strlen(tmp_line)
        tmp_value = strmid(tmp_line,result+1,tmp_len-result)
        case tmp_str of
          'Kernel_dir': fldname_kernel = tmp_value
          'Kernel_list': metakr = tmp_value
          'Topo_data_file_name': geo_ifname = tmp_value
          'Observer_latitude': obs_lat = double(tmp_value)
          'Observer_longitude': obs_lon = double(tmp_value)
          'Observer_altitude': obs_alt = double(tmp_value)
          'Resolution': Tel_res_sec = double(tmp_value)
          'Pixel_width': pix_x = double(tmp_value)
          'Pixel_height': pix_y = double(tmp_value)
          'North_direction_mode': North_deg_mode = tmp_value
          'North_up': North_up = tmp_value
          'Image_rotation_angle' : image_rot_angle_deg = double(tmp_value)
          'Target_name': Target_name = tmp_value
          else: print,"!!",tmp_value
        endcase
      endif
    endelse

  endwhile
  close,1

  return, err
end
