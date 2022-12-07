;;
;; Finding sub-pixel peak position of correlation surface
;;
function subpix_vec_1t_ellipse_surf,cor,distance
  distance_sub = double(distance)

  td0 = distance[0]
  td1 = distance[1]
  cor_siz = size(cor)
  ndata = 9
  xdat = dblarr(ndata)
  ydat = dblarr(ndata)
  fdat = dblarr(ndata)

  ;; Peakが端の場合は終了
  if td0 lt 2  and td0 ge cor_siz[1]-2 and td0 lt 2  and td0 ge cor_siz[1]-2 then  return, distance_sub

  xdat[0]=td0;
  xdat[1]=td0-1;
  xdat[2]=td0+1;
  xdat[3]=td0-1;
  xdat[4]=td0+1;
  xdat[5]=td0-1;
  xdat[6]=td0+1;
  xdat[7]=td0;
  xdat[8]=td0;

  ydat[0]=td1;
  ydat[1]=td1-1;
  ydat[2]=td1-1;
  ydat[3]=td1+1;
  ydat[4]=td1+1;
  ydat[5]=td1;
  ydat[6]=td1;
  ydat[7]=td1-1;
  ydat[8]=td1+1;

  fdat=cor[xdat,ydat]

  tmp_xdat = dblarr(6)
  tmp_ydat = dblarr(6)
  tmp_fdat = dblarr(6)

  tmp_xm=dblarr(4)
  tmp_ym=dblarr(4)
  tmp_fm=dblarr(4)

  xm=0. & ym = 0. & fm = 0. & mcount=0.

  ;; Pattern 1 x type
  tmp_xdat[0:4] = xdat[0:4]
  tmp_ydat[0:4] = ydat[0:4]
  tmp_fdat[0:4] = fdat[0:4]

  for i=0, 4-1,1 do begin
    tmp_xdat[5] = xdat[5+i]
    tmp_ydat[5] = ydat[5+i]
    tmp_fdat[5] = fdat[5+i]

    fm0 = qubicmax_invert(tmp_xdat, tmp_ydat,tmp_fdat,xm0,ym0)
    ;  print,"##",distance[0],distance[1]
    ;  print,"##",fdat
    ;  print,"##",xm0,ym0,ym0
    ;  stop

    if abs(xm0 - td0) le 1 and abs(ym0-td1) le 1 then begin
      xm+=xm0
      ym+=ym0
      fm+=xm0
      mcount++
    endif
  endfor

  ;; Pattern 2 + type
  tmp_xdat[1:4] = xdat[5:8]
  tmp_ydat[1:4] = ydat[5:8]
  tmp_fdat[1:4] = fdat[5:8]

  for i=0, 4-1,1 do begin
    tmp_xdat[5] = xdat[i+1]
    tmp_ydat[5] = ydat[i+1]
    tmp_fdat[5] = fdat[i+1]

    fm0 = qubicmax_invert(tmp_xdat, tmp_ydat,tmp_fdat,xm0,ym0)
    if abs(xm0 - td0) le 1 and abs(ym0-td1) le 1 then begin
      xm+=xm0
      ym+=ym0
      fm+=xm0
      mcount++
    endif
  endfor

  ;; パターン3    □_ type
  tmp_xdat[1] = xdat[1] & tmp_ydat[1] = ydat[1] & tmp_fdat[1] = fdat[1];
  tmp_xdat[2] = xdat[4] & tmp_ydat[2] = ydat[4] & tmp_fdat[2] = fdat[4];
  tmp_xdat[3] = xdat[5] & tmp_ydat[3] = ydat[5] & tmp_fdat[3] = fdat[5];
  tmp_xdat[4] = xdat[7] & tmp_ydat[4] = ydat[7] & tmp_fdat[4] = fdat[7];

  tmp_xdat[5] = xdat[2] & tmp_ydat[5] = ydat[2] & tmp_fdat[5] = fdat[2];

  fm0 = qubicmax_invert(tmp_xdat,tmp_ydat,tmp_fdat,xm0,ym0);
  if abs(xm0 - td0) le 1 and abs(ym0-td1) le 1 then begin
    xm+=xm0 & ym+=ym0 &  fm+=xm0 & mcount++
  endif

  tmp_xdat[5] = xdat[3] & tmp_ydat[5] = ydat[3] & tmp_fdat[5] = fdat[3]
  fm0 = qubicmax_invert(tmp_xdat,tmp_ydat,tmp_fdat,xm0,ym0);
  if abs(xm0 - td0) le 1 and abs(ym0-td1) le 1 then begin
    xm+=xm0 & ym+=ym0 &  fm+=xm0 & mcount++
  endif

  tmp_xdat[5] = xdat[6] & tmp_ydat[5] = ydat[6] & tmp_fdat[5] = fdat[6] &
  fm0 = qubicmax_invert(tmp_xdat,tmp_ydat,tmp_fdat,xm0,ym0) &
  if abs(xm0 - td0) le 1 and abs(ym0-td1) le 1 then begin
    xm+=xm0 & ym+=ym0 &  fm+=xm0 & mcount++
  endif

  tmp_xdat[5] = xdat[8] & tmp_ydat[5] = ydat[8] & tmp_fdat[5] = fdat[8] &
  fm0 = qubicmax_invert(tmp_xdat,tmp_ydat,tmp_fdat,xm0,ym0);
  if abs(xm0 - td0) le 1 and abs(ym0-td1) le 1 then begin
    xm+=xm0 & ym+=ym0 &  fm+=xm0 & mcount++
  endif


  ; パターン4    _□ type
  tmp_xdat[1] = xdat[2] & tmp_ydat[1] = ydat[2] & tmp_fdat[1] = fdat[2];
  tmp_xdat[2] = xdat[3] & tmp_ydat[2] = ydat[3] & tmp_fdat[2] = fdat[3];
  tmp_xdat[3] = xdat[6] & tmp_ydat[3] = ydat[6] & tmp_fdat[3] = fdat[6];
  tmp_xdat[4] = xdat[7] & tmp_ydat[4] = ydat[7] & tmp_fdat[4] = fdat[7];

  tmp_xdat[5] = xdat[1] & tmp_ydat[5] = ydat[1] & tmp_fdat[5] = fdat[1];
  fm0 = qubicmax_invert(tmp_xdat,tmp_ydat,tmp_fdat,xm0,ym0);
  if abs(xm0 - td0) le 1 and abs(ym0-td1) le 1 then begin
    xm+=xm0 & ym+=ym0 &  fm+=xm0 & mcount++
  endif

  tmp_xdat[5] = xdat[4] & tmp_ydat[5] = ydat[4] & tmp_fdat[5] = fdat[4]
  fm0 = qubicmax_invert(tmp_xdat,tmp_ydat,tmp_fdat,xm0,ym0)
  if abs(xm0 - td0) le 1 and abs(ym0-td1) le 1 then begin
    xm+=xm0 & ym+=ym0 &  fm+=xm0 & mcount++
  endif

  tmp_xdat[5] = xdat[5] & tmp_ydat[5] = ydat[5] & tmp_fdat[5] = fdat[5]
  fm0 = qubicmax_invert(tmp_xdat,tmp_ydat,tmp_fdat,xm0,ym0)
  if abs(xm0 - td0) le 1 and abs(ym0-td1) le 1 then begin
    xm+=xm0 & ym+=ym0 &  fm+=xm0 & mcount++
  endif
  tmp_xdat[5] = xdat[8] & tmp_ydat[5] = ydat[8] & tmp_fdat[5] = fdat[8]
  fm0 = qubicmax_invert(tmp_xdat,tmp_ydat,tmp_fdat,xm0,ym0)
  if abs(xm0 - td0) le 1 and abs(ym0-td1) le 1 then begin
    xm+=xm0 & ym+=ym0 &  fm+=xm0 & mcount++
  endif

  ; パターン5    □~ type
  tmp_xdat[1] = xdat[2] & tmp_ydat[1] = ydat[2] & tmp_fdat[1] = fdat[2]
  tmp_xdat[2] = xdat[3] & tmp_ydat[2] = ydat[3] & tmp_fdat[2] = fdat[3];
  tmp_xdat[3] = xdat[5] & tmp_ydat[3] = ydat[5] & tmp_fdat[3] = fdat[5];
  tmp_xdat[4] = xdat[8] & tmp_ydat[4] = ydat[8] & tmp_fdat[4] = fdat[8];

  tmp_xdat[5] = xdat[1] & tmp_ydat[5] = ydat[1] & tmp_fdat[5] = fdat[1];
  fm0 = qubicmax_invert(tmp_xdat,tmp_ydat,tmp_fdat,xm0,ym0);
  if abs(xm0 - td0) le 1 and abs(ym0-td1) le 1 then begin
    xm+=xm0 & ym+=ym0 &  fm+=xm0 & mcount++
  endif
  tmp_xdat[5] = xdat[4] & tmp_ydat[5] = ydat[4] & tmp_fdat[5] = fdat[4];
  fm0 = qubicmax_invert(tmp_xdat,tmp_ydat,tmp_fdat,xm0,ym0);
  if abs(xm0 - td0) le 1 and abs(ym0-td1) le 1 then begin
    xm+=xm0 & ym+=ym0 &  fm+=xm0 & mcount++
  endif
  tmp_xdat[5] = xdat[7] & tmp_ydat[5] = ydat[7] & tmp_fdat[5] = fdat[7];
  fm0 = qubicmax_invert(tmp_xdat,tmp_ydat,tmp_fdat,xm0,ym0);
  if abs(xm0 - td0) le 1 and abs(ym0-td1) le 1 then begin
    xm+=xm0 & ym+=ym0 &  fm+=xm0 & mcount++
  endif
  tmp_xdat[5] = xdat[6] & tmp_ydat[5] = ydat[6] & tmp_fdat[5] = fdat[6];
  fm0 = qubicmax_invert(tmp_xdat,tmp_ydat,tmp_fdat,xm0,ym0);
  if abs(xm0 - td0) le 1 and abs(ym0-td1) le 1 then begin
    xm+=xm0 & ym+=ym0 &  fm+=xm0 & mcount++
  endif

  ; パターン6    ~□ type
  tmp_xdat[1] = xdat[1] & tmp_ydat[1] = ydat[1] & tmp_fdat[1] = fdat[1];
  tmp_xdat[2] = xdat[4] & tmp_ydat[2] = ydat[4] & tmp_fdat[2] = fdat[4];
  tmp_xdat[3] = xdat[6] & tmp_ydat[3] = ydat[6] & tmp_fdat[3] = fdat[6];
  tmp_xdat[4] = xdat[8] & tmp_ydat[4] = ydat[8] & tmp_fdat[4] = fdat[8];

  tmp_xdat[5] = xdat[2] & tmp_ydat[5] = ydat[2] & tmp_fdat[5] = fdat[2];
  fm0 = qubicmax_invert(tmp_xdat,tmp_ydat,tmp_fdat,xm0,ym0);
  if abs(xm0 - td0) le 1 and abs(ym0-td1) le 1 then begin
    xm+=xm0 & ym+=ym0 &  fm+=xm0 & mcount++
  endif
  tmp_xdat[5] = xdat[3] & tmp_ydat[5] = ydat[3] & tmp_fdat[5] = fdat[3];
  fm0 = qubicmax_invert(tmp_xdat,tmp_ydat,tmp_fdat,xm0,ym0);
  if abs(xm0 - td0) le 1 and abs(ym0-td1) le 1 then begin
    xm+=xm0 & ym+=ym0 &  fm+=xm0 & mcount++
  endif
  tmp_xdat[5] = xdat[5] & tmp_ydat[5] = ydat[5] & tmp_fdat[5] = fdat[5];
  fm0 = qubicmax_invert(tmp_xdat,tmp_ydat,tmp_fdat,xm0,ym0);
  if abs(xm0 - td0) le 1 and abs(ym0-td1) le 1 then begin
    xm+=xm0 & ym+=ym0 &  fm+=xm0 & mcount++
  endif
  tmp_xdat[5] = xdat[7] & tmp_ydat[5] = ydat[7] & tmp_fdat[5] = fdat[7];
  fm0 = qubicmax_invert(tmp_xdat,tmp_ydat,tmp_fdat,xm0,ym0);
  if abs(xm0 - td0) le 1 and abs(ym0-td1) le 1 then begin
    xm+=xm0 & ym+=ym0 &  fm+=xm0 & mcount++
  endif

  if mcount ge 1 then begin
    xm/=mcount
    ym/=mcount
    fm/=mcount
  endif
  ;  print,distance[0],distance[1]
  ;  print,xm,ym,mcount
  ;  stop

  distance_sub[0]=xm
  distance_sub[1]=ym

  return, distance_sub
end
