;;
;; Fitting a curve surface to given 2d data
;;
function qubicmax_invert, xdat, ydat,fdat,xm,ym
  error = 0.
  Adata = dblarr(6,6)
  Adata[0,*] = [1.,1.,1.,1.,1.,1.]
  Adata[1,*] = xdat
  Adata[2,*] = ydat
  Adata[3,*] = xdat^2.
  Adata[4,*] = xdat*ydat
  Adata[5,*] = ydat^2.

  iAdata = invert(Adata,/double)
  tmp_coef = iAdata##transpose(fdat)

  det = tmp_coef[4]*tmp_coef[4]-4.*tmp_coef[3]*tmp_coef[5]
  xm = (2*tmp_coef[1]*tmp_coef[5]-tmp_coef[2]*tmp_coef[4])/det;
  ym = (2*tmp_coef[2]*tmp_coef[3]-tmp_coef[1]*tmp_coef[4])/det;
  fm = tmp_coef[0]+tmp_coef[1]*xm+tmp_coef[2]*ym+tmp_coef[3]*xm*xm+tmp_coef[4]*xm*ym+tmp_coef[5]*ym*ym;

  return, fm
end
