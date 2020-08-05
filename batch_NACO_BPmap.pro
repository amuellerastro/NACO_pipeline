function one_gauss, x, p

  z1 = (x-p[1])/p[2]
  g1 = p[0]*exp(-z1^2./2.d0)
  fit = g1

  return, fit

end


function two_gauss, x, p

  z1 = (x-p[1])/p[2]
  z2 = (x-p[4])/p[5]
  g1 = p[0]*exp(-z1^2./2.d0)
  g2 = p[3]*exp(-z2^2./2.d0)
  fit = g1+g2

  return, fit

end


pro batch_NACO_BPmap, path, display

fileo = file_search(path+'NACO*fits', count=nfiles)

out = strarr(nfiles)
object = strarr(nfiles)
type = strarr(nfiles)
catg = strarr(nfiles)
opti1 = strarr(nfiles)
opti3 = strarr(nfiles)
opti6 = strarr(nfiles)
opti7 = strarr(nfiles)
dit = strarr(nfiles)
naxis1 = strarr(nfiles)
naxis2 = strarr(nfiles)
airmass = strarr(nfiles)

st = strarr(nfiles)

for i=0,nfiles-1 do begin

  pos1 = strpos(fileo[i], '/', /reverse_search)
  pos2 = strpos(fileo[i], '.fits', /reverse_search)
  out[i] = strmid(fileo[i], pos1+1, pos2-pos1-1)

  hdr = headfits(fileo[i], exten=0, /silent)

  object[i] = strcompress(get_eso_keyword(hdr,'OBJECT'),/rem)
  naxis1[i] = strcompress(get_eso_keyword(hdr,'NAXIS1'),/rem)
  naxis2[i] = strcompress(get_eso_keyword(hdr,'NAXIS2'),/rem)
  type[i] = strcompress(get_eso_keyword(hdr,'HIERARCH ESO DPR TYPE'),/rem)
  catg[i] = strcompress(get_eso_keyword(hdr,'HIERARCH ESO DPR CATG'),/rem)
  dit[i] = strcompress(get_eso_keyword(hdr,'HIERARCH ESO DET DIT'),/rem)
  dit[i] = sigfig(dit[i],2)
  opti1[i] = strcompress(get_eso_keyword(hdr,'HIERARCH ESO INS OPTI1 NAME'),/rem)	;mask
  opti3[i] = strcompress(get_eso_keyword(hdr,'HIERARCH ESO INS OPTI3 NAME'),/rem)	;pupil stop
  opti6[i] = strcompress(get_eso_keyword(hdr,'HIERARCH ESO INS OPTI6 NAME'),/rem)	;filter
  opti7[i] = strcompress(get_eso_keyword(hdr,'HIERARCH ESO INS OPTI7 NAME'),/rem)	;Objective
  airmass[i] = strcompress(get_eso_keyword(hdr,'AIRMASS'),/rem)
  airmass[i] = sigfig(airmass[i],5)

  st[i] = dit[i]+'s_'+naxis1[i]+'px'

endfor

idx = where(catg eq 'CALIB' and type eq 'FLAT,SKY')
fff = fileo[idx]
object = object[idx]
type = type[idx]
catg = catg[idx]
opti1 = opti1[idx]
opti3 = opti3[idx]
opti6 = opti6[idx]
opti7 = opti7[idx]
dit = dit[idx]
naxis1 = naxis1[idx]
naxis2 = naxis2[idx]
airmass = airmass[idx]
st = st[idx]

;select low and high airmass
dum = min(airmass, idxmin)
dum = max(airmass, idxmax)

fff = [fff[idxmin],fff[idxmax]]

; idx = uniq(st, sort(st))
; ust = st[idx]
; udit = dit[idx]
; unaxis1 = naxis1[idx]
; nu = n_elements(ust)
; 
; if (nu ne 1) then stop	;there should be only one parameter set of sky,flats

bgf = file_search(path+'../Reduced/'+'master_dark*fits', count=nd)


;===============================================================================================

for xx=0,n_elements(bgf)-1 do begin

  ;read in of data

  bg = (mrdfits(bgf[xx], 0, hdrbg, /silent))[*,*,0]
  exptime = strcompress(get_eso_keyword(hdrbg, 'EXPTIME'),/rem)

  ff1 = mrdfits(fff[0], 0, hdrff, /silent)
  ff2 = mrdfits(fff[1], 0, /silent)
  
  fsx = float(get_eso_keyword(hdrff, 'HIERARCH ESO DET WIN STARTX'))
  fsy = float(get_eso_keyword(hdrff, 'HIERARCH ESO DET WIN STARTY'))
  fnx = float(get_eso_keyword(hdrff, 'HIERARCH ESO DET WIN NX'))
  fny = float(get_eso_keyword(hdrff, 'HIERARCH ESO DET WIN NY'))

  ;cut out FF
  naxis1 = strcompress(get_eso_keyword(hdrbg,'NAXIS1'),/rem)
  nx = float(get_eso_keyword(hdrbg, 'HIERARCH ESO DET WIN NX'))
  ny = float(get_eso_keyword(hdrbg, 'HIERARCH ESO DET WIN NY'))
  ny = nx
  sx = float(get_eso_keyword(hdrbg, 'HIERARCH ESO DET WIN STARTX'))
  sy = float(get_eso_keyword(hdrbg, 'HIERARCH ESO DET WIN STARTY'))

  ff1 = ff1[(sx-fsx):(sx-fsx)+nx-1,(sy-fsy):(sy-fsy)+ny-1]
  ff2 = ff2[(sx-fsx):(sx-fsx)+nx-1,(sy-fsy):(sy-fsy)+ny-1]

  ; ff1 = ff1+32768.
  ; ff2 = ff2+32768.

  ;===============================================================================================
; 
;   ;create a bad pixel mask form the dark by fitting 2 Gaussians to the histogram of counts
; 
;   bpbg = intarr(nx,ny)
; 
;   idxs = sort(bg)
;   bgs = bg[idxs]
;   bgx = dindgen(n_elements(bgs))
; 
;   hd = histogram(bgs, binsize=100., locations=hx)
; 
;   ; window, 2
;   nbins = floor((max(bgs)-min(bgs))/100.+1.)
;   ; hx = dindgen(nbins)*100.d0+0.05d0
; 
;   if (n_elements(hx) ne nbins) then begin
; 
;     print, ''
;     print, 'Number of bins not computed correclty. Stop.'
;     stop
; 
;   endif
; 
;   dumerr = n_elements(hx)
;   dumerr[*] = 1.
;   peak = max(hd, idxmax)
; 
; ;=========================================================================
; idx = where(hx gt -1.d4 and hx lt 1d4)
; t1 = hd[idx]
; t2 = hx[idx]
; peak = max(t1, idxmax)
; start_val = [peak, t2[idxmax], 1.d3]
; hd = t1
; hx = t2
; ;=========================================================================
;     start_val = [peak, hx[idxmax], 1.d3]
; 
;     result = mpfitfun('one_gauss', hx, hd, dumerr, start_val, weights=sqrt(hd), yfit=yfit, /quiet);, parinfo=pi)
; 
;     fwhm1 = 2.d0*sqrt(2.d0*alog(2.d0))*result[2]
; 
;   ;   x1 = abs(result[1]-fwhm1)
;     x1 = (result[1]-1.5*fwhm1)	;did this for J-band observations. Gaussian around 0
;     x2 = (result[1]+1.5*fwhm1)
;     xr = [x1,x2]
;     x0 = min(xr)
;     x1 = max(xr)
; 
;   ; endelse
; 
;   if keyword_set(display) then begin
; 
;     window, 0, xs=1000, ys=600
;     !p.multi = [0,1,2]
;     plot, bgx, bgs, psym=3, xst=1, charsize=2, yst=1, thick=3;, yr=[-5.*median(bgs),5.*median(bgs)]
;     plot, hx, hd, xst=1, charsize=2, thick=3
;       oplot, hx, yfit, color=cgcolor('red')
;       g1 = one_gauss(hx, result[0:2])
;       ;if (ditbg gt 32.d0) then g2 = one_gauss(hx, result[3:5])
;       oplot, hx, g1, color=cgcolor('yellow')
;       ;if (ditbg gt 32.d0) then oplot, hx, g2, color=cgcolor('yellow')
;       plots, x0*[1,1], !y.crange, color=cgcolor('green')
;       plots, x1*[1,1], !y.crange, color=cgcolor('green')
; 
;     !p.multi = [0,1,0]
; 
;   endif
; 
;   print, ''
;   print, 'Counts considered: ', x0, x1
;   print, ''
; 
;   ; quest = ''
;   ; read, 'Limits OK? y/n: ', quest
;   ; ; quest = 'n'
;   ; if (quest ne 'y') then begin
;   ; 
;   ;   read, 'Lower Limit: ', x0
;   ;   read, 'Upper Limit: ', x1
;   ; ; x0 = -12.d3
;   ; ; x1 = -3.d3
;   ; 
;   ;   plots, x0*[1,1], !y.crange, color=cgcolor('red')
;   ;   plots, x1*[1,1], !y.crange, color=cgcolor('red')
;   ; 
;   ; endif
; 
;   idx = where(bg lt x0 or bg gt x1)
; 
;   idxbad = array_indices(bpbg, idx)
;   for i=0L,n_elements(idxbad[0,*])-1 do bpbg[idxbad[0,i], idxbad[1,i]] = 1

  ;===============================================================================================

  bp1 = intarr(nx,ny)
  bp1[*,*] = 1
  bp2 = bp1

  ;compute ratio of both FFs and reject outliers
  ; ratio = ff2/ff1
  ratio = ff2/ff1

  rm = robust_mean(ratio, 4., sigma, numrej, goodind=good)
  ;bp1[*,*] = 0
  idxgood = array_indices(bp1, good)
  for i=0L,n_elements(idxgood[0,*])-1 do bp1[idxgood[0,i], idxgood[1,i]] = 0

  ;make a selection of bad pixels (in this case pixels which do not get illuminated) based on flux
  rm = robust_mean(ff2, 4., sigma, numrej, goodind=good)
  idxgood = array_indices(bp2, good)
  for i=0L,n_elements(idxgood[0,*])-1 do bp2[idxgood[0,i], idxgood[1,i]] = 0

  ;combine all BP maps

  bp = bp1+bp2;+bpbg
  idx0 = where(bp gt 0.)
  idx = array_indices(bp, idx0)
  for i=0L,n_elements(idx[0,*])-1 do bp[idx[0,i], idx[1,i]] = 1

  ;===============================================================================================

  if (n_elements(bp[*,0]) lt n_elements(bp[1,*])) then bp = bp[*,0:n_elements(bp[*,0])-1]

  writefits, path+'../Reduced/'+'static_badpixels_'+sigfig(exptime,2)+'s_'+naxis1+'px.fits', bp

  ;===============================================================================================

  if keyword_set(display) then begin

    window, 2, xs=nx, ys=ny
    cgimage, bp, stretch=2

  endif

  print, ''
  print, 'Percantage of bad Pixel: ', 100.*double(n_elements(where(bp eq 1.)))/double(n_elements(bp))
  print, ''

endfor


end
