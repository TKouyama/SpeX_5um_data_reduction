;;
;; IDL code for data reduction for 5um Venus images obtained by IRTF/SpeX
;;
;; 1. Generating Sky image by averaging sky frames, which is used for background subtraction
;; 2. Generating Sky flat from the Sky image with reducing line noises
;; 3. Preparing a reference Venus image for template matching
;; 4. Stacking Venus images by shift-and-add procedure based on template matching
;; (option) 5. adding geometry information

;;
;; Toru Kouyama :: t.kouyama@aist.go.jp
;;
pro IRTF_SpeX_5um_processing
  ;;
  src_directory = '/path/to/copy_directory/' ;; <= Modify this path to your copy directory
  cd,src_directory

  ;; use three IDL windows
  window_Set,3,ct=0
  
  n_obs = 1000l
  obs_date    = strarr(n_obs)
  obs_seq_number = strarr(n_obs) ;; Observation sequence number in each day
  st_n        = lonarr(n_obs)
  en_n        = lonarr(n_obs)
  ref_frame_i = lonarr(n_obs)

  ;;
  ;; you can add setting of any observation sequence here.
  ;; Each observation sequence is expected to be composed of
  ;;  "Venus imaging" and "Sky imaging".
  ;;
  obs_n = 0l

;  obs_date[obs_n]       = '20200912' ;; observation date and directory name
;  obs_seq_number[obs_n] = "1"
;  st_n[obs_n]           = 000
;  en_n[obs_n]           = 581
;  ref_frame_i[obs_n]    = 281
;  obs_n++
;
;  obs_date[obs_n]       = '20220220'
;  obs_seq_number[obs_n] = "1" 
;  st_n[obs_n]           = 000
;  en_n[obs_n]           = 300
;  ref_frame_i[obs_n]    = 138
;  obs_n++ 

  ;
  obs_date[obs_n]       = '20220220'
  obs_seq_number[obs_n] = "1"
  st_n[obs_n]           = 000
  en_n[obs_n]           = 300
  ref_frame_i[obs_n]    = 138
  obs_n++ ;;

  obs_date[obs_n]       = '20220220'
  obs_seq_number[obs_n] = "2"
  st_n[obs_n]           = 301
  en_n[obs_n]           = 600
  ref_frame_i[obs_n]    = 353
  obs_n++ ;;

  ;;
  ;; Test
  ;;
  ;test_n = 0
  test_n = 1

  ;; Badpixel map
  ifldname_badpix = './data/20210326/'
  ifname_badpix = ifldname_badpix + 'badPixMask.fits'
  badpix_image = readfits(ifname_badpix)
  badpix_image[264:275,*] = 0d
  ;tvscl,badpix_image

  ;; Directory path for 5.1 um images
  ;;  Folder for each obsearvation day is expected to locate under {ifldname_data},
  ;;  and simple name indicating observatoin date, such as, "20220220"

  ifldname_data = './data/' ;; <= Replace your data direcotry path
  ifldname = ifldname_data + obs_date[test_n]+'/'

  ;; Output directory
  ofldname = './results/'
  result = file_test(ofldname,/directory)
  if result eq 0 then begin
    file_mkdir,ofldname
    print,"MKDIR: ",ofldname
  endif

  exec_IRTF_SpeX_5um_processing,obs_date[test_n] $
                               ,obs_seq_number[test_n] $
                               ,st_n[test_n] $
                               ,en_n[test_n] $
                               ,ref_frame_i[test_n] $
                               ,badpix_image $
                               ,ifldname $
                               ,ofldname

  return
end

;;
;; main
;;
pro exec_IRTF_SpeX_5um_processing,obs_date $
                               ,obs_seq_number $
                               ,st_n $
                               ,en_n $
                               ,ref_frame_i $
                               ,badpix_image $
                               ,ifldname $
                               ,ofldname

  
  ;;
  ;; Data load
  ;;
  ifnames = file_search(ifldname,'*.fits',count=f_count)
  
  ;;
  ;; Data screening thresholds
  ;;  
  thres_venus_count = 100d
  ;thres_cor = 0.975d
  thres_cor = 0.995d
  plot_r = 30d

  ;;
  ;; SpeX specifications
  ;; 
  ;; SpeX plate scale
  tmp_dtheta_sec = 0.1175 ;; arcsec

  ;; image width and hight
  nx = 512l
  ny = 512l
  upper_line_y = 465l

  ;; X and Y 2d coordinates
  im_x = dindgen_x(nx,ny)
  im_y = dindgen_x(nx,ny,/y)

  ;;;;;;;;;;;;;;;;;;
  ;; 1. Sky image ;;
  ;;;;;;;;;;;;;;;;;;
  print,"Sky..."
  sky_im = dblarr(nx,ny)
  valid_im = dblarr(nx,ny)
  sky_im_c = 0l
  badpix_thres = 47500d

  for i=st_n, en_n-1, 1 do begin
    tmp_head = headfits(ifnames[i])
    tmp_flt = fits_header_keyword_value(tmp_head,'GFLT')
    tmp_com_S = fits_header_keyword_value(tmp_head,"COMMENT = 'S")

    if strmid(strtrim(tmp_com_S,2),1,3) eq "Sky" then begin
      raw_image = double(readfits(ifnames[i],/silent))
      
      ;; revresing count order
      tmp_sky_im = 65535d - raw_image
      sky_im += tmp_sky_im * (1d - badpix_image)
      
      ;; for eliminating dead pixels, raw_image is used
      valid_im += raw_image * (1d - badpix_image)
      
      ;; counting sky image number
      sky_im_c++
    endif
  endfor
  
  ;; Check Valid position ;;
  mean_valid_im = valid_im/sky_im_c
  valid_pos_im = double(mean_valid_im gt badpix_thres) * (sky_im ne 0 )

  ;; Sky image
  mean_sky_im = sky_im/sky_im_c
  s_sky_im = g_smooth(mean_sky_im, 20, ignore_value=0)
  mean_sky_im[where(valid_pos_im eq 0)] = s_sky_im[where(valid_pos_im eq 0)]

  ;; 0 filling for out of FOV
  mean_sky_im[0:24,*] = 0
  mean_sky_im[*,upper_line_y:511] = 0

  wset,0
  tvscl_window,mean_sky_im > 5000d
  wset,1
  tvscl_window,valid_pos_im
  wset,0
  print,mean(mean_sky_im),max(mean_sky_im),min(mean_sky_im)
  stop

  ;;;;;;;;;;;;;;;;;;
  ;; 2. Sky flat  ;; with reducing line noises (readout noise)
  ;;;;;;;;;;;;;;;;;;
  ;; extract valid region
  tmp_sky_flat = mean_sky_im[25:511,0:upper_line_y -1] / mean(mean_sky_im[25:511,0:upper_line_y -1])
  f_sky_flat = fft(tmp_sky_flat)

  ;; highest frequency components in x and y directions should be due to readout noise
  ;; thus filling 0

  ;; Reducing y-line noise
  f_sky_flat[1:486,0] = complex(0.,0)

  ;; Reducing x-line noise (only high frequency components)
  ;f_sky_flat[0,1:upper_line_y -1] = complex(0.,0)
  f_sky_flat[0,upper_line_y/2-10:upper_line_y/2+10] = complex(0.,0)

  i_f_sky_im_flat = real_part(fft(f_sky_flat,/inverse))
  tmp_sky_flat = i_f_sky_im_flat / mean(i_f_sky_im_flat[242-64:242+64,237-64:237+64])

  sky_flat = dblarr(nx,ny)
  sky_flat[*,*] = 1d
  sky_flat[25:511,0:upper_line_y -1] = tmp_sky_flat

  wset,0
  tmp_out = sky_flat
  tmp_out[where(sky_flat eq 1d)] = 0.7
  tvscl_window,tmp_out > 0.7 < 1.3
  wset,0
  stop

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; 3. Make a reference frame ;; which will be used for stacking
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  template_im_raw = 65535d - double(readfits(ifnames[ref_frame_i],/silent))
  template_im = template_im_raw - mean_sky_im
  ;; flat correction ;;
  template_im = template_im / sky_flat * valid_pos_im

  ;; Replacing dead pixel value to mean value of surrounding pixels
  s_template_im = g_smooth(template_im,5, ignore_Value=0)
  template_im[where(valid_pos_im eq 0)] = s_template_im[where(valid_pos_im eq 0)]

  ;; Simple background level estimation, need to be modified
  for j=0, upper_line_y-1, 1 do begin
    tmp_line = template_im[*,j]
    ;bg_level = mean([tmp_line[30:40],tmp_line[490:500]])
    MEANCLIP,[tmp_line[30:40],tmp_line[490:500]], bg_level, Sigma
    template_im[*,j] -= bg_level
  endfor
  template_im = template_im * valid_pos_im

  wset,0
  tvscl_window,template_im < 1500 > (-500)
  wset,0

  ;;
  ;; Reducing readout noise in y direction
  ;;
  sigma_y = 25d
  ;; component of higher frequency in x direction
  tmp_g_smooth1 = g_smooth_diff_sig(template_im[25:511,0:upper_line_y -1], 0.01, sigma_y,/limb)
  ;; component of mid-higher frequency  in x direction
  tmp_g_smooth2 = g_smooth_diff_sig(template_im[25:511,0:upper_line_y -1], 5, sigma_y,/limb)
  
  ;; By subtracting mid-higher from higher, then only target higher components will remain
  line_noise = tmp_g_smooth1 - tmp_g_smooth2

  tmp_reduce_line = dblarr(nx,ny)
  tmp_reduce_line[25:511,0:upper_line_y-1] = template_im[25:511,0:upper_line_y-1] - line_noise
  template_im = tmp_reduce_line

  template_im[0:30,*] = 0.
  template_im[*,upper_line_y:511] = 0.

  wset,1
  tvscl_window,template_im < 1500 > (-500)
  wset,0

  ;;
  ;; Centering Venus position in a frame (just using a bright center, not so accurate)
  ;;
  
  ;; Finding Venus disk region => "1"
  template_im_c = (template_im gt thres_venus_count) * valid_pos_im

  ;; Estimating brighness center
  template_im_c_pos = where(template_im_c eq 1,v_count)
  template_br_cx = mean(im_x[template_im_c_pos] * template_im_c[template_im_c_pos])
  template_br_cy = mean(im_y[template_im_c_pos] * template_im_c[template_im_c_pos])

  ;; subpixel shifting ;;
  ;template_im = shift(tmp_Venus_im,255.5-br_cx,255.5-br_cy)
  template_im_center = d_shift(template_im,255.5-template_br_cx,255.5-template_br_cy)

  ;; Extracting only valid Venus disk region ;;
  template_im_orig = template_im_center
  
  ;; clipping only venus disk region
  template_im_center = template_im_orig[50:511-50, 50 : ny-1-50] > 0
  offset_x = 50
  offset_y = 50

  wset,2
  tvscl_window,template_im_center < 1500 > (-500)
  wset,0
  stop

  ;;;;;;;;;;;;;;;;;;;;;;
  ;; 4. Shift and Add ;; based on template matching
  ;;;;;;;;;;;;;;;;;;;;;;
  Venus_im   = dblarr(nx, ny)
  Venus_im_c = dblarr(nx, ny)
  thres_venus = 50d
  used_image_n = 0l
  ;loop_end = 25l ;; for debug

  for i=st_n, en_n-1, 1 do begin
  ;for i=st_n, st_n+loop_end, 1 do begin

    tmp_head = headfits(ifnames[i])
    tmp_com_V = fits_header_keyword_value(tmp_head,"COMMENT = 'V")
    tmp_mode = long(fits_header_keyword_value(tmp_head,'RO_MODE'))
    tmp_flt = fits_header_keyword_value(tmp_head,'GFLT')

    if strmid(strtrim(tmp_com_V,2),1,3) eq "Ven" and tmp_mode eq 0 then begin

      ;; Header から必要な情報を取り出す ;;
      Date_obs=fits_header_keyword_value(tmp_head,'DATE_OBS')
      Time_obs=fits_header_keyword_value(tmp_head,'TIME_OBS')
      Date_obs = strmid(Date_obs,1,10)
      Time_obs = strmid(Time_obs,1,15)
      obs_date_utc = Date_obs+'T'+Time_obs
      obs_date_utc = reform(obs_date_utc)

      ;;
      ;; Load image ;;
      ;;
      tmp_Venus_im_raw = 65535d - double(readfits(ifnames[i],/silent))
      tmp_Venus_im = tmp_Venus_im_raw - mean_sky_im
      tmp_Venus_im = tmp_Venus_im / sky_flat * valid_pos_im

      ;; replacing values at invalid pixels with mean value of surrouding pixels ;;
      tmp_s_Venus_im = g_smooth(tmp_Venus_im,5)
      tmp_Venus_im[where(valid_pos_im eq 0)] = tmp_s_Venus_im[where(valid_pos_im eq 0)]

      ;; Filling 0 in out of FOV region
      tmp_Venus_im[0:24,*] = 0.
      tmp_Venus_im[*,upper_line_y:511] = 0.

      ;; Simple background level estimation, need to be modified
      for j=0, upper_line_y-1, 1 do begin
        tmp_line = tmp_Venus_im[*,j]
        ;bg_level = mean([tmp_line[30:40],tmp_line[501:511]])
        ;MEANCLIP,[tmp_line[30:40],tmp_line[490:500]], bg_level, Sigma
        MEANCLIP,[tmp_line[490:500]], bg_level, Sigma
        tmp_Venus_im[*,j] -= bg_level
      endfor
      tmp_Venus_im = tmp_Venus_im * valid_pos_im

      ;;
      ;; Reducing readout noise in y direction with asymmetric Gausssian filter
      ;;
      sigma_y = 25d
      ;; component of higher frequency in x direction
      tmp_g_smooth1 = g_smooth_diff_sig(tmp_Venus_im[25:511,0:upper_line_y -1], 0.01, sigma_y,/limb)
      ;; component of mid-higher frequency  in x direction
      tmp_g_smooth2 = g_smooth_diff_sig(tmp_Venus_im[25:511,0:upper_line_y -1], 5, sigma_y,/limb)

      ;; By subtracting mid-higher from higher, then only target higher components will remain
      line_noise = tmp_g_smooth1 - tmp_g_smooth2

      tmp_reduce_line = dblarr(nx,ny)
      tmp_reduce_line[25:511,0:upper_line_y-1] = tmp_Venus_im[25:511,0:upper_line_y-1] - line_noise
      tmp_Venus_im = tmp_reduce_line

      tmp_Venus_im[0:30,*] = 0.
      tmp_Venus_im[*,upper_line_y:511] = 0.
      ;; for debugging
      ;wset,1
      ;tvscl_window,line_noise
      ;wset,0
      
      ;;
      ;; Reducing readout noise in x direction
      ;; comment: Reducing line noise in x direction seems to be also needed.
      ;;
      sigma_x = 25d
      tmp_g_smooth1 = g_smooth_diff_sig(tmp_Venus_im[50:511,0:upper_line_y-1], sigma_x, 0.01,/limb)
      tmp_g_smooth2 = g_smooth_diff_sig(tmp_Venus_im[50:511,0:upper_line_y-1], sigma_x, 5,/limb)
      line_noise = tmp_g_smooth1 - tmp_g_smooth2
      tmp_n_line = dblarr(nx,ny)
      tmp_n_line[50:511,0:upper_line_y-1] = tmp_Venus_im[50:511,0:upper_line_y-1] - line_noise
      tmp_Venus_im = tmp_n_line

      ;;
      ;; Performing shift and add based on template matching
      ;;
      tmp_venus_c = (tmp_Venus_im gt thres_venus) * valid_pos_im

      ;; Estimating bright center for determining Venus position rougly ;;
      tmp_venus_c_pos = where(tmp_venus_c eq 1,v_count)
      br_cx = mean(im_x[tmp_venus_c_pos] * tmp_venus_c[tmp_venus_c_pos])
      br_cy = mean(im_y[tmp_venus_c_pos] * tmp_venus_c[tmp_venus_c_pos])

      thres_v_c = 50d*50d*!dpi ;; Reference Venus disk size for eliminating error ;;
      if sqrt((br_cx-255.5)^2.+(br_cy-255.5)^2.) lt 250d $
         and br_cy lt 400 $
         and v_count gt thres_v_c then begin

        ;; keeping original frame
        tmp_Venus_im_ori = tmp_Venus_im

        ;; Shifting the Venus position to reduce computation cost for estimating correlation surface ;;
        tmp_Venus_im = shift(tmp_Venus_im,255.5-br_cx,255.5-br_cy)
        tmp_Venus_c = shift(tmp_Venus_c,255.5-br_cx,255.5-br_cy)

        ;; Estimating correlation surface
        result_cor = f_cross_correlation(template_im_center,tmp_Venus_im > 0,csx,csy)
        ;stop

        ;; Finding maximum position of the correlation surface (integer)
        max_cor=max(result_cor,tmp_sub);,min(result)
        tmp_center_pos = [tmp_sub mod csx, tmp_sub/csx]

        ;; Fiding maximum position of the coorelation surface with sub-pixel level
        sub_center_pos = subpix_vec_1t_ellipse_surf(result_cor,tmp_center_pos)
        s_delay_x = sub_center_pos[0] - csx/2 ;;- (offset_x)
        s_delay_Y = sub_center_pos[1] - csy/2 ;;- (offset_y)

        tmp_Venus_im = d_shift(tmp_Venus_im,-s_delay_x,-s_delay_y)
        tmp_Venus_c = d_shift(tmp_Venus_c,-s_delay_x,-s_delay_y)
        ;stop

        ;; shift後の画像でcを作る
        ;; validな値が入っているところをチェック ;;
        tmp_venus_c_pos = where(tmp_venus_c gt 0.5)

        if max_cor gt thres_cor then begin
          ;; valid 判定のところだけ像に活かす ;;
          Venus_im[tmp_venus_c_pos] += tmp_Venus_im[tmp_venus_c_pos]
          venus_im_c[tmp_venus_c_pos] += 1.
        endif

        tmp_distance = sqrt((br_cx-template_br_cx)^2.+(br_cy-template_br_cy)^2.)
        measured_distance = sqrt((-(sub_center_pos[0] - csx/2) + br_cx-template_br_cx)^2. $
                               +(-(sub_center_pos[1] - csy/2) + br_cy-template_br_cy)^2.)

        used_image_n++

        print,i ,used_image_n, strmid(obs_date_utc,0,19), float(max_cor), tmp_distance, measured_distance $
          ,format='(I5, " ", I5, " ", a19, " ", f11.8, " ", 2f11.3)'

        wset,0
        tmp_out = tmp_Venus_im
        tmp_out[0,0] = -100
        tmp_out[0,1] = 1300
        tvscl_window,tmp_out > tmp_out[0,0] < tmp_out[0,1]
        out_string = strmid(obs_date_utc,0,19)+" "+string(float(max_cor))
        xyouts,50,600,out_string, charsize=1.5,color=255, /device
        wset,0

        wset,1
        tmp_out = tmp_Venus_im_ori
        tmp_out[0,0] = -100
        tmp_out[0,1] = 1300
        tvscl_window,tmp_out > tmp_out[0,0] < tmp_out[0,1]
        xyouts,50,600,out_string, charsize=1.5,color=255, /device
        wset,0

        ;; ql output for debug
        ;ofname_png = ofldname + file_basename(ifnames[i],'.fits')+'.png'
        ;tmp_out_png = byte(double(tmp_out-tmp_out[0,0])/(tmp_out[0,1]-tmp_out[0,0])*255d > 0 < 255)
        ;write_png,ofname_png,tmp_out_png
        wait,0.1

        ;stop
      endif
    endif

    if i mod 100 eq 0 then print,i
    wait,0.01

  endfor

  ;; averaging
  venus_im /= venus_im_c > 1

  ;;
  ;; QL & output ;;
  ;;
  wset,0
  tvscl_window,Venus_im > 0 < 1500

  ofname_png =  ofldname + obs_date + "_"+ obs_seq_number + '.png'
  tmp_out = venus_im
  tmp_out[0,0] = 0 & tmp_out[0,1] = max(Venus_im)
  tmp_out_png = byte(double(tmp_out-tmp_out[0,0])/(tmp_out[0,1]-tmp_out[0,0])*255d > 0 < 255)
  write_png,ofname_png,tmp_out_png

  wset,0
  
  ;;
  ;; bzero check and output
  ;;
  SXADDPAR, tmp_head, 'BZERO', 0
  ofname_fits = ofldname + obs_date + "_"+ obs_seq_number + '.fits'
  writefits,ofname_fits,Venus_im,tmp_head

  stop
  ;; if you want to add geometry information, comment-out below line.
  ;return

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; 5. Adding geometry infromation (not completed) ;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; need appropriate setting for running simulate_TO
  ;; need upside down

  ;;
  ;; Below is Kouyama's setting ;;
  ;;
  loc_file = './locationfile_IRTF_SpeX_NorthUp.txt'
  cx_p = 0 & cy_p = 0
  simulate_TO,utctim=obs_date_utc[0]  $
    ,ifname_location=loc_file      $
    ,cx_p=cx_p,cy_p=cy_p            $
    ,obs_geo_struct=obs_geo_struct  $
    ,ccd_geo=ccd_geo,/noplot

  ;; north azimuth of Venus seeing from IRTF position ;;
  Na_Deg = obs_geo_struct.J2000_NORTH_AZIMUTH

  ;; Frame rotation angle from header
  POSANGLE = fits_header_keyword_value(tmp_head,'POSANGLE')

  wset,0
  tvscl_window,Venus_im > 0 < 1500

  wset,1
  Venus_im_rot = rot(Venus_im,POSANGLE - Na_Deg,/interp)
  tvscl_window, Venus_im_rot > (-100) < 1300

  ofname_png =  ofldname + obs_date + "_"+ obs_seq_number + '_rot.png'
  tmp_out = Venus_im_rot
  tmp_out[0,0] = 0 & tmp_out[0,1] = max(Venus_im)
  tmp_out_png = byte(double(tmp_out-tmp_out[0,0])/(tmp_out[0,1]-tmp_out[0,0])*255d > 0 < 255)
  write_png,ofname_png,tmp_out_png


  wset,2
  ;; high passed
  g_smooth_r = 10d ;; for 2020 06, 07
  Venus_im_rot_high_passed = rot((Venus_im-g_smooth(Venus_im, g_smooth_r,/limb)),POSANGLE - Na_Deg,/interp)
  tvscl_window, Venus_im_rot_high_passed > (-plot_r) < plot_r
  wset,0

  ;;
  ;; bzero check
  ;;
  SXADDPAR, tmp_head, 'BZERO', 0
  ;; original
  out_image = dblarr(nx,ny,2)
  out_image[*,*,0] = Venus_im_rot
  out_image[*,*,1] = Venus_im_rot_high_passed

  ofname_fits = ofldname + obs_date + "_"+ obs_seq_number + '_rot.fits'
  writefits,ofname_fits,out_image,tmp_head

  return
end
