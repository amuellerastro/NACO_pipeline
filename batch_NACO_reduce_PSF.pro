pro batch_NACO_reduce_PSF, star, path, fwhm, display, filter

restore, path+'../Reduced/PSF_sky_ident.sav', /v

files = file_search(path+'NACO*fits', count=nfiles)

; out = strarr(nfiles)
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
ndit = strarr(nfiles)
nexp = uintarr(nfiles)
st = strarr(nfiles)

for i=0,nfiles-1 do begin

  pos1 = strpos(files[i], '/', /reverse_search)
  pos2 = strpos(files[i], '.fits', /reverse_search)
;   out[i] = strmid(fileo[i], pos1+1, pos2-pos1-1)

  hdr = headfits(files[i], exten=0, /silent)

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
  ndit[i] = strcompress(get_eso_keyword(hdr, 'HIERARCH ESO DET NDIT'),/rem)
  nexp[i] = uint(get_eso_keyword(hdr, 'HIERARCH ESO SEQ NEXPO'))
  st[i] = dit[i]+'s_'+naxis1[i]+'px'

endfor

if (filter eq 'l') then tmpsz = 256.
if (filter eq 'm') then tmpsz = 128.

idx = where(catg eq 'SCIENCE' and naxis1 le tmpsz)
rawf = files[idx]
scidit = dit[idx[0]]
scinx = naxis1[idx[0]]

;==================================================================================================

;calibration data

;bad pixel
bpf = path+'../Reduced/'+'static_badpixels_'+scidit+'s_'+strcompress(scinx,/rem)+'px.fits'
;dark
dkf = path+'../Reduced/'+'master_dark_'+scidit+'s_'+strcompress(scinx,/rem)+'px.fits'
;flat field
fff = path+'../Reduced/'+'instrument_flat.fits'

;sky background
; skf = file_search(path+'../Reduced/'+'PSF_sky_background_*fits', count=nsk)

;==================================================================================================

nx = uint(scinx)
ny = nx
dim = nx

;read in files

ff = mrdfits(fff, 0, hdrff, /silent)
bp = mrdfits(bpf, 0, hdrbg, /silent)
dk = mrdfits(dkf, 0, hdrdk, /silent)

naxis1 = strcompress(get_eso_keyword(hdrbg,'NAXIS1'),/rem)

nx = float(naxis1)
fnx = strcompress(get_eso_keyword(hdrff,'NAXIS1'),/rem)
ny = nx

ff = ff[fnx/2.-nx/2:fnx/2.+nx/2-1,fnx/2.-nx/2:fnx/2.+nx/2-1]

;==================================================================================================

quad = uniq(flag, sort(flag))
step = n_elements(quad)

for qu=0,step-1 do begin

  idx = where(flag eq qu+1)
  rawf = path+rfname[idx]
  nfiles = n_elements(rawf)

  ;sky background
  skf = file_search(path+'../Reduced/'+'PSF_sky_'+strcompress(qu+1,/rem)+'_background_*fits', count=nsk)

  ;extract date from sky files to find the closest one to science later
  jdsk = dblarr(n_elements(skf))
  for i=0,n_elements(skf)-1 do begin

    head = headfits(skf[i], exten=0, /silent)
    jdsk[i] = double(get_eso_keyword(head, 'JD'))

  endfor

  ;read in all sky frames for later interpolation
  hdrtmp = headfits(skf[0], exten=0, /silent)
  sknx = float(get_eso_keyword(hdrtmp, 'NAXIS1'))
  skny = float(get_eso_keyword(hdrtmp, 'NAXIS2'))
  sktmp = dblarr(sknx, skny, n_elements(skf))
  for j=0,n_elements(skf)-1 do sktmp[*,*,j] = mrdfits(skf[j], 0, /silent)

  ;==================================================================================================

  skyflag = 'y'

  ;read in cubes, do cosmetics, average

  nfiles = n_elements(rawf)
  jd = dblarr(nfiles)
  outname = strarr(nfiles)
  used_frames = intarr(2,nfiles)
  date = strarr(nfiles)

  for i=0,nfiles-1 do begin
  ;for i=1,1 do begin

    print, ''
    print, 'Reducing Cube '+strcompress(i+1,/rem)+' / '+strcompress(nfiles,/rem)
    print, 'Reading File'
    cube = mrdfits(rawf[i], 0, hdrraw, /silent)
    if (i eq 0.) then hdrsciref = hdrraw

    szcube = size(cube)
    if (szcube[1] lt szcube[2]) then cube = cube[*,0:nx-1,*]
    if (n_elements(szcube) gt 5) then cube = cube[*,*,30:n_elements(cube[0,0,*])-2]	;last frame is already averaged, remove it
    nframes = n_elements(cube[0,0,*])
    ;cube = cube

    ;----------------------------------------------------------------------------------

    ;output name
    pos1 = strpos(rawf[i], '/', /reverse_search)
    pos2 = strpos(rawf[i], '.fits', /reverse_search)
    outname[i] = strmid(rawf[i], pos1+1, pos2-pos1-1)

    exptime = double(get_eso_keyword(hdrraw, 'EXPTIME'))
    ndit = double(get_eso_keyword(hdrraw, 'HIERARCH ESO DET NDIT'))
    dit = double(get_eso_keyword(hdrraw,'HIERARCH ESO DET DIT'))
    jd[i] = double(get_eso_keyword(hdrraw, 'MJD-OBS'))+2400000.5d0
    date[i] = get_eso_keyword(hdrraw,'DATE')

    ;----------------------------------------------------------------------------------

    ;find corresponding sky frame
    if (n_elements(jdsk) gt 1) then begin

      idx = closest(jdsk, jd[i])

      if (jd[i] lt jdsk[0]) then sk = mrdfits(skf[idx],0,/silent)	;first obs with sky afterwards, no interpolation needed

      if (jd[i] gt jdsk[n_elements(jdsk)-1]) then sk = mrdfits(skf[idx],0,/silent)	;last observation is target and no sky, again no interpolation

      if (jd[i] gt jdsk[0] and jd[i] lt jdsk[n_elements(jdsk)-1]) then begin

	for xx=0,sknx-1 do begin
	  for yy=0,skny-1 do begin
	    sk[xx,yy] = interpol(sktmp[xx,yy,*], jdsk, jd[i])
	  endfor
	endfor

      endif

      sk = sk[*,0:nx-1]

    endif else begin

      sk = mrdfits(skf[0],0,/silent)
      sk = sk[*,0:nx-1]

    endelse

    ;----------------------------------------------------------------------------------

    im = fltarr(nx,ny,nframes)
    corim = fltarr(nx,ny,nframes)	;cosmetically corrected images

    for j=0,nframes-1 do begin

      ;cosmetics

      if (skyflag eq 'y') then im[*,*,j] = (cube[*,*,j]-sk)/ff $
	  else im[*,*,j] = (cube[*,*,j]-bg)/ff
  ;       if (skyflag eq 'y') then im[*,*,j] = (cube[*,*,j]-(sk*(median(cube[*,*,j])/median(sk))))/ff $
  ; 	else im[*,*,j] = (cube[*,*,j]-bg)/ff

      tmp = im[*,*,j]
      fixpix_mod, tmp, bp, outim, npix=24, /weight, /silent
      corim[*,*,j] = outim

  ;       tmp = corim[*,*,j]
  ;       la_cosmic_MOD, tmp, outim, path, skyval=median(sk+32768.,/even), gain=100.;, gain=9.8, readn=4.4
  ;       corim[*,*,j] = tmp

      ;----------------------------------------------------------------------------------

      ;remove additional horizontal additive pattern by subtracting the median of each row (NACO manuel, p62)
      ;we have 4 quadrants, hence define windows for each quadrant

      if (nx gt 500) then ws = 40 else ws = 20
      ll = median(corim[0:ws-1,*,j], dimension=1)	;left side, left window
      lr = median(corim[(nx/2)-ws:(nx/2)-1,*,j], dimension=1)	;left side, right window
      rl = median(corim[nx/2:nx/2+ws-1,*,j], dimension=1)	;right side, left window
      rr = median(corim[nx-1-ws:*,*,j], dimension=1)	;right side, right window

      for k=0,ny-1 do begin

	corim[0:nx/2-1,k,j] = corim[0:nx/2-1,k,j]-mean([ll[k],lr[k]])
	corim[nx/2:*,k,j] = corim[nx/2:*,k,j]-mean([rl[k],rr[k]])

      endfor

      ;----------------------------------------------------------------------------------


      proceeding_text,loop=nframes, i=j, prompt='> Frame   '+string(j+1,form='(I4)')

    endfor

    cube = 0
    imdk = 0
    im  = 0

    ;----------------------------------------------------------------------------------

    ;estimate star position using median combined image and do quality control

    flux = dblarr(nframes)
    mflux = dblarr(nframes)
    sdev = dblarr(nframes)

    if (n_elements(size(corim)) gt 5) then medtmp = median(corim,dim=3,/even) $
      else medtmp = corim
    dum = max(smooth(medtmp,3, /edge_truncate), idx)	;smooth to be more resistant against left over BP
    idxmax = array_indices(medtmp, idx)
    txc = idxmax[0]
    tyc = idxmax[1]
    medtmp = 0
    dum = 0

  ;   mask_t = shift(dist(nx), txc, tyc)
  ;   mask = mask_t ge 1.*fwhm and mask_t le 4.*fwhm
  ;   mask = mask[*]	;convert to 1D vector
  ;   idxmask = where(mask eq 1)
  ; 
  ;   for j=0,nframes-1 do begin
  ; 
  ;     ;measure stddev between 1-4 lambda/D, see Absil 2013
  ;     tmp = corim[*,*,j]
  ;     tmp = tmp[*]
  ;     sdev[j] = stddev(tmp[idxmask])
  ; 
  ;     ;measure flux, additional criteria, is more robust, e.g. betaPic/RAW/NACO.2013-02-01T03:46:19.205.fits
  ;     tmp = corim[*,*,j]
  ;     circint_MOD, tmp, txc, tyc, 2.*fwhm, tot, mtot, meantot, maxpx, sdtot, npix, totdif, npixdif, t8
  ;     flux[j] = tot
  ;     mflux[j] = mtot
  ; 
  ;     proceeding_text,loop=nframes, i=j, prompt='> Quality   '+string(j+1,form='(I4)')
  ; 
  ;   endfor

    sim = dblarr(nframes)
    medim = median(corim, dim=3)
    for j=0,nframes-1 do sim[j] = stddev(corim[*,*,j]-medim)
    resistant_mean, sim, 3, t1, t2, nbad, /double, goodvec=idxgood, badvec=idxbad

    if keyword_set(display) then begin

      window, 0, title=star
      xa = dindgen(nframes)
      plot, xa, sim, /yn, title=star
      oplot, xa[idxgood], sim[idxgood], psym=2, color=cgcolor('green')
      if (idxbad[0] ne -1) then oplot, xa[idxbad], sim[idxbad], psym=2, color=cgcolor('red')

    endif

    ;----------------------------------------------------------------------------------

    if (n_elements(size(corim)) gt 5) then begin

  ;     ;reject bad frames
  ;     idxgood1 = where(sdev lt median(sdev)+stddev(sdev))
  ;     idxgood2 = where(flux gt median(flux)-stddev(flux))
  ; 
  ;     flag = intarr(nframes)
  ; 
  ;     for j=0,nframes-1 do begin
  ; 
  ;       idx1 = where(j eq idxgood1)
  ;       if (idx1[0] eq -1) then flag[j] = 1
  ;       idx2 = where(j eq idxgood2)
  ;       if (idx2[0] eq -1) then flag[j] = 1
  ; 
  ;     endfor
  ; 
  ;     idxgood = where(flag eq 0)

      used_frames[*,i] = [n_elements(idxgood), nframes+1]	;+1 because we remove the first frame

      print, strcompress(uint(nframes-n_elements(idxgood)),/rem)+' / '+strcompress(uint(nframes),/rem)+' frames rejected'

    ;   window, 0
    ;   !p.multi=[0,1,2]
    ;   plot, sdev, /yn, psym=2, xtitle='Frame', ytitle='SDEV'
    ;     oplot, !x.crange, median(sdev)*[1,1]
    ;     oplot, !x.crange, (median(sdev)+stddev(sdev))*[1,1], color=cgcolor('red')
    ;   plot, mflux, /yn, psym=2, xtitle='Frame', ytitle='Median Flux'
    ;     oplot, !x.crange, median(mflux)*[1,1]
    ;     oplot, !x.crange, (median(flux)-stddev(flux))*[1,1], color=cgcolor('red')
    ;   !p.multi=[0,1,0]

      nframes = n_elements(idxgood)
      corim = corim[*,*,idxgood]

    endif

    ;----------------------------------------------------------------------------------

    ;shift corrected images to geometric center, find center by fitting a Gaussian

    scorim = fltarr(nx,ny,nframes)
    xpos = dblarr(nframes)
    ypos = dblarr(nframes)

    if (dim gt 128.) then radius = 25.
    if (dim eq 64.) then radius = 7.
    if (dim eq 128.) then radius = 15.

    if (star eq 'HD19668' and strmatch(date, '*2015-12-19*') eq 1) then radius = 15.
    if (star eq 'HD131835' and strmatch(date, '*2016-05-04*') eq 1) then radius = 15.
    if (star eq 'HD22049' and strmatch(date, '*2015-12-17*') eq 1) then radius = 5.
    if (star eq 'HD121617') then radius = 13.
    if (star eq 'HD145263' and strmatch(date, '*2018-06-19*') eq 1) then radius = 15.
    if (star eq 'CPD-722713' and strmatch(date, '*2018-06-19*') eq 1) then radius = 15.
    if (star eq 'HD121191') then radius = 15.
    ;if (star eq 'HD103234' and strmatch(date, '*2018-06-04*') eq 1) then radius = 9.

    
    if (txc-radius le 0. or txc+radius ge dim or tyc-radius le 0 or tyc+radius ge dim) then begin
    
        txc = dim/2.
        tyc = dim/2.    
    
    endif

    for j=0,nframes-1 do begin

      cutim = corim[txc-radius:txc+radius,tyc-radius:tyc+radius,j]
      xa = dindgen(2*radius+1)+1.d0 & ya = xa

      estimates = [median(cutim), max(cutim), fwhm, fwhm, radius, radius, 0., 1.]

      ;dist_circle, weight, 2.*radius+1
      ;weight = 1./(weight+1.)

;       weights = cutim
;       idx1 = where(cutim eq 0.)	;e.g. beta Pic close to the center which is masked out
;       if (idx1[0] ne -1) then begin
; 	idx2 = array_indices(weights, idx1)
; 	for k=0,n_elements(idx2[0,*])-1 do weights[idx2[0,k],idx2[1,k]] = 1.d5
;       endif
; 
;   ;       ;work around for having not signum becauseof IDL < v8.3
;       signdum = fltarr(n_elements(cutim[*,0]),n_elements(cutim[0,*]))
;       for isign=0,n_elements(cutim[*,0])-1 do begin
; 	for jsign=0,n_elements(cutim[0,*])-1 do begin
; 
; 	  signdum[isign,jsign] = sign((reform(cutim[isign,jsign]))[0])
; 
; 	endfor
;       endfor
;       sign = signdum
;   ;      sign = signum(cutim)
; 
;       idx1 = where(sign eq -1.)
;       if (idx1[0] ne -1) then begin
; 	idx2 = array_indices(weights, idx1)
; 	for k=0,n_elements(idx2[0,*])-1 do weights[idx2[0,k],idx2[1,k]] = 1.d5
;       endif
;       weights = 1./sqrt(weights)

      yfit = mpfit2dpeak(cutim, A, xa, ya, /moffat, estimates=estimates, dof=dof, chisq=chisq, perror=perror, sigma=sigma, /tilt, /quiet);, weights=weights)

      xpos[j] = A[4]-radius+txc
      ypos[j] = A[5]-radius+tyc

      ;scorim[*,*,j] = shift_sub(corim[*,*,j], nx/2.-xpos[j]+1, ny/2.-ypos[j]+1)
    ;to shift the image using FFT set a frame around the image because of ringing
      wframe = 150.	;width of frame
      tmp = dblarr(nx+2.*wframe, ny+2.*wframe)
      tmp[wframe:wframe+nx-1, wframe:wframe+nx-1] = corim[*,*,j]
      stmp = fftshift(tmp, nx/2.-xpos[j]+1, ny/2.-ypos[j]+1)
      scorim[*,*,j] = stmp[(nx+2.*wframe)/2.-nx/2:(nx+2.*wframe)/2.+nx/2.-1,(nx+2.*wframe)/2.-nx/2:(nx+2.*wframe)/2.+nx/2.-1]


      proceeding_text,loop=nframes, i=j, prompt='> Fitting Frame   '+string(j+1,form='(I4)')

    endfor

    if (n_elements(size(scorim)) gt 5) then begin	;it can happen that all frames are rejected, except for one...skip entire cube

      smedcorim = median(scorim, dim=3, /even)

      scorim = 0
      corim = 0

    ;----------------------------------------------------------------------------------

      ;cut out frame, write out result

      csmedcorim = smedcorim[nx/2.-radius:nx/2.+radius, ny/2.-radius:ny/2.+radius]

  ;     mwrfits, csmedcorim, path+'PSF_'+date[i]+'.fits', /silent
      writefits, path+'../Reduced/'+'PSF_'+date[i]+'.fits', csmedcorim

    endif else begin

      smedcorim = scorim
      csmedcorim = smedcorim[nx/2.-radius:nx/2.+radius, ny/2.-radius:ny/2.+radius]
      writefits, path+'../Reduced/'+'PSF_'+date[i]+'.fits', csmedcorim

    endelse

  endfor

endfor


; endif else begin	;no AGPM
; 
;   file = file_search(path+'img_Lp_dc.fits')
; 
;   cube = mrdfits(file, 0, /silent)
;   psf = median(cube, dim=3, /even)
; 
;   szpsf = size(psf)
; 
;   radius = 25.
;   cutpsf = psf[szpsf[1]/2.-radius:szpsf[1]/2.+radius,szpsf[2]/2.-radius:szpsf[2]/2.+radius]
; 
;   mwrfits, cutpsf, path+'PSF_Lp.fits', /silent
; 
; endelse


;=============================================================================================

; print, ''
; print, '*******************************************'
; print, 'DELETE CORRUPTED PSF FILES AND PRESS ENTER!'
; print, '*******************************************'
; spawn, 'ds9 '+path+'PSF_2*.fits'
; hak

file = file_search(path+'../Reduced/'+'PSF_2*.fits', count=npsf)

psf = fltarr(2.*radius+1,2.*radius+1,npsf)

for i=0,npsf-1 do psf[*,*,i] = mrdfits(file[i],0,/silent)

sim = dblarr(npsf)
medim = median(psf, dim=3)
for i=0,npsf-1 do sim[i] = stddev(psf[*,*,i]-medim)
resistant_mean, sim, 3, t1, t2, nbad, /double, goodvec=idxgood, badvec=idxbad

print, 'Removed '+strcompress(round(n_elements(idxbad)))+' / '+strcompress(round(npsf))+' PSF frames'

psf_final = median(psf[*,*,idxgood], dim=3, /even)

; mwrfits, psf_final, path+'PSF_Lp.fits', /silent
if (filter eq 'l') then writefits, path+'../Reduced/'+'PSF_Lp.fits', psf_final, hdrsciref
if (filter eq 'm') then writefits, path+'../Reduced/'+'PSF_Mp.fits', psf_final, hdrsciref

end
