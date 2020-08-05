function get_parangle, jd, epoch, ra, dec, pma, pmd, radvel, plx, lat, lon, altitude, temp, pres

  caldat, jd, mo, day, yr, hh, mm, ss

  ;get date of the form yyyy-mm-ddThh:mm:ss
  date = dblarr(6)
  date[0] = yr
  date[1] = mo
  date[2] = day
  date[3] = hh
  date[4] = mm
  date[5] = ss

  ;get year of observations for new RA, DEC, i.e. epoch
  hour = hh+mm/60.d0+ss/3600.d0
  if ((double(yr) mod 4.) eq 0.) then date2 = yr+mo/12.d0+day/366.d0+hour/8766.0d0 $
    else date2 = yr+mo/12.d0+day/365.d0+hour/8766.0d0	;takes leap year into account

  ;compute current cordinates and parallactic angle at time of observation
  dt = date2-epoch
  eq2hor_MOD, ra, dec, pma, pmd, radvel, plx, dt, jd, alt, az, hatmp, lat=lat, lon=lon, altitude=altitude, pres=pres, temp=temp, outra=outra, outdec=outdec;, /verbose
  ha = hatmp
  ha = ha/15.

  parang = parangle(ha, outdec, lat)

  if (outdec gt lat and parang lt 0.) then parang = parang+360.d0

  return, {parang:parang, ha:ha}

end


pro batch_NACO_reduce_Science, star, path, agpmflag, filestar, onlymed, fwhm, qintsky, display

restore, filestar, /verbose
ra_st = ra
dec_st = dec

;==================================================================================================
;science files

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
mjd = dblarr(nfiles)

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
  mjd[i] = double(get_eso_keyword(hdr,'MJD-OBS'))

endfor

if keyword_set(agpmflag) then qagpm = 'y' else qagpm = 'n'

if keyword_set(agpmflag) then begin

  ;this will include the raw sky frames as science frames too
  idx = where(catg eq 'SCIENCE' and naxis1 ge 512 and nexp ge median(nexp))

  ;after this MJD we used the 'new' AGPM OBs where we sky is recorded from the OB with proper keyword instead doing it manual
  if (min(mjd[idx]) lt 57693.d0) then begin

    rawf = files[idx]
    scidit = dit[idx[0]]
    scinx = naxis1[idx[0]]

  endif else begin

    idx = where(catg eq 'SCIENCE' and naxis1 ge 512 and nexp ge median(nexp) and type eq 'OBJECT')

    rawf = files[idx]
    scidit = dit[idx[0]]
    scinx = naxis1[idx[0]]

  endelse


endif else begin

  idx = where(catg eq 'SCIENCE' and naxis1 ge 512)
  rawf = files[idx]

  scidit = dit[idx[0]]
  scinx = naxis1[idx[0]]

endelse

;==================================================================================================

;calibration data

;bad pixel
bpf = path+'../Reduced/'+'static_badpixels_'+scidit+'s_'+strcompress(scinx,/rem)+'px.fits'
;dark
dkf = path+'../Reduced/'+'master_dark_'+scidit+'s_'+strcompress(scinx,/rem)+'px.fits'
;flat field
fff = path+'../Reduced/'+'instrument_flat.fits'

;sky background
skf = file_search(path+'../Reduced/'+'sky_background_*fits', count=nsk)
; cnx = uintarr(nsk)
; cdit = fltarr(nsk)
; 
; for i=0,nsk-1 do begin
; 
;   hdr = headfits(cfile[i], exten=0)
;   cnx[i] = get_eso_keyword(hdr, 'NAXIS1')
;   cdit[i] = get_eso_keyword(hdr, 'DETDIT')
; 
; endfor
;obsolete because PSF sky is named differently
; idx = where(cnx eq scinx and cdit eq scidit)
; if (idx[0] eq -1) then stop
; skf = cfile[idx]
; skf = cfile

;==================================================================================================

nx = uint(scinx)
ny = nx
dim = float(nx)

;read in files

ff = mrdfits(fff, 0, hdrff, /silent)
bp = mrdfits(bpf, 0, hdrbg, /silent)
dk = mrdfits(dkf, 0, hdrdk, /silent)

sx = get_eso_keyword(hdrdk, 'HIERARCH ESO DET WIN STARTX')
sy = get_eso_keyword(hdrdk, 'HIERARCH ESO DET WIN STARTY')

ff = ff[sx-1:sx+nx-2,sy-1:sy+ny-2,*]

;==================================================================================================

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
nsky = n_elements(skf)
for j=0,nsky-1 do sktmp[*,*,j] = mrdfits(skf[j], 0, /silent)

;==================================================================================================

;read in cubes, do cosmetics, average

nfiles = n_elements(rawf)
jd = dblarr(nfiles)
used_frames = intarr(2,nfiles)
outname = strarr(nfiles)

; if (nx gt 380) then csmedcorim = fltarr(380,380,nfiles) $	;cut out version
;   else csmedcorim = fltarr(nx,ny,nfiles)

flagcut = intarr(nfiles)

;==================================================================================================

;load all sky frames and compute PCs and Eigenvalues for later

if (qintsky eq 'p') then begin

  truncate_pca = nsky
  obj = sktmp
  nobj = (size(obj))[3]

  for i=0,nobj-1 do begin

    dum = 0
    dum = reform(obj[0:dim-1,0:dim-1,i])
    dum = dum-mean(dum)
    obj[0:dim-1,0:dim-1,i] = dum

  endfor

  ;PCA
  data = transpose(reform(obj,dim*dim,nobj))	;[0:dim-1,0:dim-1,*]
  covMatrix = matrix_multiply(data, data, /btranspose)

  eigenval = la_eigenql(covMatrix, EIGENVECTORS=eigenvect, range=[nobj-truncate_pca,nobj-1], /DOUBLE)
  eigenval = reverse(eigenval)

  eigenvect = reverse(eigenvect,2)
  pc_orig = matrix_multiply(eigenvect,data,/atranspose)
  pc = pc_orig
  for k=0,nsky-1 do pc[k,*] = pc_orig[k,*]/(eigenval[k])

endif

;==================================================================================================


for i=0,nfiles-1 do begin
;for i=55,55 do begin

  print, ''
  print, 'Reducing Cube '+strcompress(i+1,/rem)+' / '+strcompress(nfiles,/rem)
  print, 'Reading File'

  cube = mrdfits(rawf[i], 0, hdrraw, /silent)
  if (i eq 0) then hdrsciref = hdrraw

  if keyword_set(onlymed) then cube = cube[*,*,n_elements(cube[0,0,*])-1]  

  szcube = size(cube)
  if (szcube[1] lt szcube[2]) then cube = cube[*,0:szcube[1]-1,*]	;assumption: X additional rows in y-dir

  if (n_elements(szcube) gt 5) then cube = cube[*,*,30:n_elements(cube[0,0,*])-2]	;last frame is already averaged, remove it, 1st frame looks bad as well

  nframes = n_elements(cube[0,0,*])
;   cube = cube;+32768.

  ;----------------------------------------------------------------------------------

  tel = get_eso_keyword(hdrraw, 'TELESCOP')
  ;Observatory parameters
  ;from http://www.eso.org/sci/facilities/paranal/astroclimate/site.html

  if (strmatch(tel, '*U4*') eq 1) then begin

    lat = ten(-24.d0, 37.d0, 31.00d0)	;for UT4
    lon = ten(-70.d0, 24.d0, 8.00d0)	;for UT4

  endif

  if (strmatch(tel, '*U1*') eq 1) then begin

    lat = ten(-24.d0, 37.d0, 33.117d0)	;for UT1
    lon = ten(-70.d0, 24.d0, 11.642d0)	;for UT1

  endif

  altitude = 2635.43d0

  ;----------------------------------------------------------------------------------

  ;output name
  pos1 = strpos(rawf[i], '/', /reverse_search)
  pos2 = strpos(rawf[i], '.fits', /reverse_search)
  outname[i] = strmid(rawf[i], pos1+1, pos2-pos1-1)

  exptime = double(get_eso_keyword(hdrraw, 'EXPTIME'))
  ndit = double(get_eso_keyword(hdrraw, 'HIERARCH ESO DET NDIT'))
  dit = double(get_eso_keyword(hdrraw,'HIERARCH ESO DET DIT'))
  jd[i] = double(get_eso_keyword(hdrraw, 'MJD-OBS'))+2400000.5d0

  ;jitter
  offx = double(get_eso_keyword(hdrraw, 'HIERARCH ESO SEQ CUMOFFSETX'))	;in px
  offy = double(get_eso_keyword(hdrraw, 'HIERARCH ESO SEQ CUMOFFSETY'))	;in px

  pres = double(get_eso_keyword(hdrraw, 'HIERARCH ESO TEL AMBI PRES START'))
  temp = double(get_eso_keyword(hdrraw, 'HIERARCH ESO TEL AMBI TEMP'))+273.15d0
  epoch = double(get_eso_keyword(hdrraw, 'HIERARCH ESO TEL TARG EPOCH'))
  pupilpos = double(get_eso_keyword(hdrraw, 'HIERARCH ESO ADA PUPILPOS'))

  ;----------------------------------------------------------------------------------

  ;find/create corresponding sky frame

  if (qintsky eq 'n') then begin

    idx = closest(jdsk, jd[i])
    sk = mrdfits(skf[idx],0,/silent)

  endif

  if (qintsky eq 'y') then begin

    if (n_elements(jdsk) gt 1) then begin

      idx = closest(jdsk, jd[i])

      if (jd[i] lt jdsk[0]) then sk = mrdfits(skf[idx],0,/silent)	;first obs with sky afterwards, no interpolation needed

      if (jd[i] gt jdsk[n_elements(jdsk)-1]) then sk = mrdfits(skf[idx],0,/silent)	;last observation is target and no sky, again no interpolation
	
      if (jd[i] gt jdsk[0] and jd[i] lt jdsk[n_elements(jdsk)-1]) then begin

	sk = dblarr(nx,ny)

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

  endif

  ;----------------------------------------------------------------------------------

  im = fltarr(nx,ny,nframes)
  corim = fltarr(nx,ny,nframes)	;cosmetically corrected images
  imdk = fltarr(nx,ny,nframes)
  corimdk = fltarr(nx,ny,nframes)
  sdev = dblarr(nframes)
  flux = dblarr(nframes)
  mflux = dblarr(nframes)

  ;----------------------------------------------------------------------------------

  if (qagpm eq 'y') then begin

    if (i eq 0) then begin

;      device, cursor_standard=2
;      window, 0, xs=1000, ys=1000

;       if (n_elements(size(cube)) gt 5) then tmp = median(cube[*,*,1:n_elements(cube[0,0,*])-2], dimension=3)-sktmp[*,*,0] else tmp = cube-sktmp[*,*,0]
      if (n_elements(size(cube)) gt 5) then tmp = median(cube[*,*,30:n_elements(cube[0,0,*])-2], dimension=3)-dk else tmp = cube-dk

      ;-----------------------------------------------------------------------------

      u = dim/2.+40.
      v = dim/2.

      trad = 50.
      tim = tmp[u-trad:u+trad-1,v-trad:v+trad-1]
      tbp = bp[u-trad:u+trad-1,v-trad:v+trad-1]

      fixpix_mod, tim, tbp, outim, npix=24, /weight, /silent
      tim = outim

      xa = dindgen(2.*trad);+1.d0 
      ya = xa
      dum = max(tim, idx0)
      idx1 = array_indices(tim, idx0)
      estimates = [median(tim), max(tim), 5., 5., idx1[0], idx1[1], 0., 1.]
      pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]},8)
      pi[1].fixed[0] = 1
      pi[2].fixed[0] = 1
      pi[3].fixed[0] = 1
      yfit = mpfit2dpeak(tim, A, xa, ya, /moffat, estimates=estimates, dof=dof, chisq=chisq, perror=perror, sigma=sigma, parinfo=pi, /quiet)

      txc = u+a[4]-trad
      tyc = v+a[5]-trad

      ;-----------------------------------------------------------------------------
      ;old manual way

;       cgimage, scale_image_am(tmp), /axis;, stretch=8, /axis;, minvalue=min(tmp), maxvalue=max(tmp)
;       print, ''
;       print, 'Click at the center of the AGPM'
;       cursor, xclick, yclick, 3, /data; /device
; 
;       txc = xclick
;       tyc = yclick

      ;-----------------------------------------------------------------------------

;       mask_t = shift(dist(nx), txc, tyc)
;       mask = mask_t ge 1.*fwhm and mask_t le 4.*fwhm
;       mask = mask[*]	;convert to 1D vector
;       idxmask = where(mask eq 1)

    endif

  endif

  ;----------------------------------------------------------------------------------

;   ;mask for quality control using AGPM
;   if (qagpm eq 'y') then begin
; 
;     if (jd[i] gt julday(6., 1., 2015, 0., 0., 0.)) then begin
; 
;       txc = nx/2.+40.	;estimated positions
;       tyc = ny/2.-8.
; 
;     endif else begin
; 
;       txc = nx/2.-44.	;estimated positions
;       tyc = ny/2.-9.
; 
;     endelse
; 
;     mask_t = shift(dist(nx), txc, tyc)
;     mask = mask_t ge 1.*fwhm and mask_t le 4.*fwhm
;     mask = mask[*]	;convert to 1D vector
;     idxmask = where(mask eq 1)
; 
;   endif

  for j=0,nframes-1 do begin

    ;cosmetics

    imdk[*,*,j] = (cube[*,*,j]-dk)/ff	;subtract dark
    tmp = imdk[*,*,j]
    fixpix_mod, tmp, bp, outimdk, npix=24, /weight, /silent
    corimdk[*,*,j] = outimdk

    if (qintsky eq 'y' or qintsky eq 'n') then begin

      im[*,*,j] = (cube[*,*,j]-sk)/ff
  ;     im[*,*,j] = (cube[*,*,j]-(sk*(median(cube[*,*,j])/median(sk))))/ff	;subtract sky
      tmp = im[*,*,j]
      fixpix_mod, tmp, bp, outim, npix=24, /weight, /silent
      corim[*,*,j] = outim

    endif

    if (qintsky eq 'p') then begin

      tim = cube[*,*,j]
      fixpix_mod, tim, bp, outim, npix=24, /weight, /silent
      outim = outim-mean(outim)

      data2 = transpose(reform(outim,dim*dim))
      s1 = matrix_multiply(pc_orig, data2, /btranspose)
      sk = reform(matrix_multiply(s1, pc), dim, dim)

      corim[*,*,j] = (outim-sk)/ff

    endif

    ;----------------------------------------------------------------------------------

    ;remove additional horizontal additive pattern by subtracting the median of each row (NACO manual, p62 or Sec.8.6: Cube mode)
    ;we have 4 quadrants, hence define windows for each quadrant

    if (qagpm eq 'n') then begin

      if (nx gt 500) then ws = 40 else ws = 20
      ll = median(corim[0:ws-1,*,j], dimension=1)	;left det., left window
      lr = median(corim[(nx/2)-ws:(nx/2)-1,*,j], dimension=1)	;left det., right window
      rl = median(corim[nx/2:nx/2+ws-1,*,j], dimension=1)	;right det., left window
      rr = median(corim[nx-1-ws:*,*,j], dimension=1)	;right det., right window

      for k=0,ny-1 do begin

	corim[0:nx/2-1,k,j] = corim[0:nx/2-1,k,j]-mean([ll[k],lr[k]])
	corim[nx/2:*,k,j] = corim[nx/2:*,k,j]-mean([rl[k],rr[k]])

      endfor

    endif else begin	;if AGPM is used the removal of horizontal lines on the upper left detector introduces black sripes because of the star altering the median level of the row

      if (nx gt 500) then ws = 40 else ws = 20
      ll = median(corim[0:ws-1,*,j], dimension=1)	;left det., left window
      ;lr = median(corim[(nx/2)-ws:(nx/2)-1,*,j], dimension=1)	;left det., right window
      ;rl = median(corim[nx/2:nx/2+ws-1,*,j], dimension=1)	;right det., left window
      rr = median(corim[nx-1-ws:*,*,j], dimension=1)	;right det., right window

      for k=0,ny-1 do begin

	corim[0:nx/2-1,k,j] = corim[0:nx/2-1,k,j]-mean(ll[k])
	corim[nx/2:*,k,j] = corim[nx/2:*,k,j]-mean(rr[k])

      endfor

    endelse

    ;----------------------------------------------------------------------------------


;     ;mask for quality control if jitter is used
;     if (qagpm eq 'n') then begin
; 
;       txc = nx/2.+offx	;estimated positions
;       tyc = ny/2.+offy
; 
;       mask_t = shift(dist(nx), txc, tyc)
;       mask = mask_t ge 1.*fwhm and mask_t le 4.*fwhm
;       mask = mask[*]	;convert to 1D vector
;       idxmask = where(mask eq 1)
; 
;     endif


;     if (qagpm eq 'y') then begin
; 
;       ;measure stddev between 1-4 lambda/D, see Absil 2013
;       tmp = corim[*,*,j]
;       tmp = tmp[*]
; 
;       sdev[j] = stddev(tmp[idxmask])
; 
;       ;measure flux, additional criteria, is more robust, e.g. betaPic/RAW/NACO.2013-02-01T03:46:19.205.fits
;       tmp = corim[*,*,j]
;       circint_MOD, tmp, txc, tyc, 4.*fwhm, tot, mtot, meantot, maxpx, sdtot, npix, totdif, npixdif, t8
;       flux[j] = tot
;       mflux[j] = mtot
; 
;     endif

    proceeding_text,loop=nframes, i=j, prompt='> Reducing Frame   '+string(j+1,form='(I4)')

  endfor

  ;I am doing this because its more robust computing a mean image to find the star position than relying on offset positions in the header
  if (qagpm eq 'n') then begin

    if (n_elements(size(corim)) gt 5) then begin

      ;cutting the image because of bad bad pixel correction at the edges
      if (nx eq 1024) then begin

	mtmp = median(corim[512-384:512+383,512-384:512+383,*], dim=3)
	flag1024 = 1

      endif else begin

	mtmp = median(corim, dim=3)
	flag1024 = 0

      endelse

      dum = max(mtmp, idxmax)
      idx0 = array_indices(mtmp, idxmax)
      ;estimated OBJECT positions
      txc = idx0[0]
      tyc = idx0[1]

      if (flag1024 eq 1) then begin

	txc = txc+128
	tyc = tyc+128

      endif

;       mask_t = shift(dist(nx), txc, tyc)
;       mask = mask_t ge 1.*fwhm and mask_t le 4.*fwhm
;       mask = mask[*]	;convert to 1D vector
;       idxmask = where(mask eq 1)
; 
;       for j=0,nframes-1 do begin
; 
; 	;measure stddev between 1-4 lambda/D, see Absil 2013
; 	tmp = corim[*,*,j]
; 	tmp = tmp[*]
; 
; 	sdev[j] = stddev(tmp[idxmask])
; 
; 	;measure flux, additional criteria, is more robust, e.g. betaPic/RAW/NACO.2013-02-01T03:46:19.205.fits
; 	tmp = corim[*,*,j]
; 	circint_MOD, tmp, txc, tyc, 4.*fwhm, tot, mtot, meantot, maxpx, sdtot, npix, totdif, npixdif, t8
; 	flux[j] = tot
; 	mflux[j] = mtot
; 
;       endfor

    endif else begin

      dum = max(corim, idxmax)
      idx0 = array_indices(corim, idxmax)
      ;estimated OBJECT positions
      txc = idx0[0]
      tyc = idx0[1]

    endelse

  endif

  cube = 0
  imdk = 0
  im  = 0

  ;----------------------------------------------------------------------------------

  sim = dblarr(nframes)
  medim = median(corim, dim=3)
;   for j=0,nframes-1 do sim[j] = stddev(corim[*,*,j]-medim)
  for j=0,nframes-1 do sim[j] = stddev(corim[nx/2.-174:nx/2.+175,nx/2.-174:nx/2.+175,j]-medim)
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
;     idxgood2 = where(mflux lt median(mflux)+stddev(mflux))
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

    used_frames[*,i] = [n_elements(idxgood), nframes+1]	;+1 because we remove the first frame, we remove the last frame as well, but this is the mean combined one

    print, strcompress(uint(nframes-n_elements(idxgood)),/rem)+' / '+strcompress(uint(nframes),/rem)+' frames rejected'

  ;   window, 0
  ;   !p.multi=[0,1,2]
  ;   plot, sdev, /yn, psym=2, xtitle='Frame', ytitle='SDEV'
  ;     oplot, !x.crange, median(sdev)*[1,1]
  ;     oplot, !x.crange, (median(sdev)+stddev(sdev))*[1,1], color=cgcolor('red')
  ;   plot, mflux, /yn, psym=2, xtitle='Frame', ytitle='Median Flux'
  ;     oplot, !x.crange, median(mflux)*[1,1]
  ;     oplot, !x.crange, (median(mflux)+stddev(mflux))*[1,1], color=cgcolor('red')
  ;   !p.multi=[0,1,0]

    nframes = n_elements(idxgood)
    corim = corim[*,*,idxgood]
    corimdk = corimdk[*,*,idxgood]

  ;   medcorim[*,*,i] = median(corim, dim=3, /even)
    medcorimdk = median(corimdk, dim=3, /even)

  endif

  ;----------------------------------------------------------------------------------

  ;compute parallactic angle for exposure start and end of a cube

  ra = [double(strmid(ra_st,0,2)),double(strmid(ra_st,2,2)),double(strmid(ra_st,4,6))]
  sign1 = strmid(dec_st,0,1)
  if (sign1 eq '-') then sign1 = -1.d0 else sign1 = 1.d0
  dec = [sign1*double(strmid(dec_st,1,2)), double(strmid(dec_st,3,2)), double(strmid(dec_st,5,6))]
  ra = ten(ra)*15.
  dec = ten(dec)

  if (finite(radvel[0] ne 1)) then radvel = 0.
  if (finite(plx[0] ne 1)) then plx = 1./100.

  if (n_elements(size(corim)) gt 5) then begin

    ;+1 because I remove the first frame anyway
    st = get_parangle(jd[i]+(dit*(idxgood[0]+30.))/86400.d0, epoch, ra, dec, pma, pmd, radvel, plx, lat, lon, altitude, temp, pres)
    parang_start = st.parang

    st = get_parangle(jd[i]+(dit*(idxgood[n_elements(idxgood)-1]+30.))/86400.d0, epoch, ra, dec, pma, pmd, radvel, plx, lat, lon, altitude, temp, pres)
    parang_end = st.parang

    parang = (parang_start+parang_end)/2.

  endif else begin

    st = get_parangle(jd[i], epoch, ra, dec, pma, pmd, radvel, plx, lat, lon, altitude, temp, pres)
    parang = st.parang

  endelse

  ha = st.ha
  if (ha gt 12.) then ha = ha-24.
  if (ha lt -12.) then ha = ha+24.


  rotptoff = 90.d0+(89.44-pupilpos)	;user manual page 71

  PA_onsky = parang-rotptoff

  ;----------------------------------------------------------------------------------

  ;shift corrected images to geometric center, find center by fitting a Gaussian to only dark subtracted images instead of sky subtracted images (to have a gaussian signal instead of donut)
  scorim = fltarr(nx, ny, nframes)	;shifted images
  xpos = dblarr(nframes)
  ypos = dblarr(nframes)

  for j=0,nframes-1 do begin

    radius = 35.
    if (txc-radius lt 0. or txc+radius gt dim-1.) then txc = dim/2.
    if (tyc-radius lt 0. or tyc+radius gt dim-1.) then tyc = dim/2.

    cutim = corimdk[txc-radius:txc+radius,tyc-radius:tyc+radius,j]
    xa = dindgen(2*radius+1)+1.d0 & ya = xa

    estimates = [median(cutim), max(cutim), fwhm, fwhm, radius, radius, 0., 1.]

    ;dist_circle, weight, 2.*radius+1
    ;weight = 1./(weight+1.)

;     weights = cutim
;     idx1 = where(cutim eq 0.)	;e.g. beta Pic close to the center which is masked out
;     if (idx1[0] ne -1) then begin
; 
;       idx2 = array_indices(weights, idx1)
;       for k=0,n_elements(idx2[0,*])-1 do weights[idx2[0,k],idx2[1,k]] = 1.d5
; 
;     endif
; 
;     ;work around for having not signum because of IDL < v8.3
;     signdum = fltarr(n_elements(cutim[*,0]),n_elements(cutim[0,*]))
;     for isign=0,n_elements(cutim[*,0])-1 do begin
; 
;       for jsign=0,n_elements(cutim[0,*])-1 do begin
; 
; 	signdum[isign,jsign] = sign((reform(cutim[isign,jsign]))[0])
; 
;       endfor
; 
;     endfor
;     sign = signdum
;     ;sign = signum(cutim)
; 
;     idx1 = where(sign eq -1.)
;     if (idx1[0] ne -1) then begin
; 
;       idx2 = array_indices(weights, idx1)
;       for k=0,n_elements(idx2[0,*])-1 do weights[idx2[0,k],idx2[1,k]] = 1.d5
; 
;     endif
;     weights = 1./sqrt(weights)

    yfit = mpfit2dpeak(cutim, A, xa, ya, /moffat, estimates=estimates, dof=dof, chisq=chisq, perror=perror, sigma=sigma, /tilt, /quiet);, weights=weights)

    xpos[j] = A[4]-radius+txc
    ypos[j] = A[5]-radius+tyc

    ;scorim[*,*,j] = shift_sub(corim[*,*,j], txc-xpos[j], tyc-ypos[j])

    ;to shift the image using FFT set a frame around the image because of ringing
    wframe = 150.	;width of frame
    tmp = dblarr(nx+2.*wframe, ny+2.*wframe)
    tmp[wframe:wframe+nx-1, wframe:wframe+nx-1] = corim[*,*,j]
    stmp = fftshift(tmp, txc-xpos[j], tyc-ypos[j])
    scorim[*,*,j] = stmp[(nx+2.*wframe)/2.-nx/2:(nx+2.*wframe)/2.+nx/2.-1,(nx+2.*wframe)/2.-nx/2:(nx+2.*wframe)/2.+nx/2.-1]

    proceeding_text,loop=nframes, i=j, prompt='> Fitting Frame   '+string(j+1,form='(I4)')

  endfor

  if (n_elements(size(corim)) gt 5) then smedcorim = median(scorim, dim=3, /even) else smedcorim = scorim

;   scorim = 0
;   corim = 0

; fit averaged corimdk as well for comparison
; shift corim frames and cut (380x380) and average


;   for j=0,nframes-1 do begin
; 
;     ;weights as a function of distance
;     dist_circle, weight, nx/2
;     tmp = shift(weight, txc, tyc)
;     weight = 1./(tmp+1.)
; 
;     estimates = [2000., 400., fwhm, fwhm, txc, tyc, 0.]
; 
;     yfit = mpfit2dpeak(corim[*,*,j], A, xa, ya, /gaussian, estimates=estimates, dof=dof, chisq=chisq, perror=perror, sigma=sigma, /tilt, weights=weight);, /quiet
; 
; 
; stop
;   endfor

  ;----------------------------------------------------------------------------------

  ;cut out frames

  if (qagpm eq 'y') then begin

    if (nx eq 512) then cutrad = 190
    if (nx eq 768) then cutrad = 250
    if (nx eq 1024) then cutrad = 400

;       x0 =txc-190-1
;       x1 = txc+188
;       y0 = tyc-190-1
;       y1 = tyc+188
    x0 = txc-cutrad-1
    x1 = txc+cutrad-2
    y0 = tyc-cutrad-1
    y1 = tyc+cutrad-2

    if (x0 lt szcube[1] and x1 lt szcube[1] and y0 lt szcube[1] and y1 lt szcube[1]) then begin

      csmedcorim = smedcorim[x0:x1,y0:y1]
      flagcut[i] = 1

    endif else begin

      flagcut[i] = 0

    endelse

  endif else begin

    if (nx eq 512) then cutrad = 100
    if (nx eq 768) then cutrad = 150
    if (nx eq 1024) then cutrad = 200

    x0 = txc-cutrad-1
    x1 = txc+cutrad-2
    y0 = tyc-cutrad-1
    y1 = tyc+cutrad-2

    if (x0 lt szcube[1] and x1 lt szcube[1] and y0 lt szcube[1] and y1 lt szcube[1] and x0 ge 0 and x1 ge 0 and y0 ge 0 and y0 ge 0) then begin

      csmedcorim = smedcorim[x0:x1,y0:y1]
      flagcut[i] = 1

    endif else begin

      flagcut[i] = 0

    endelse

  endelse


  if (flagcut[i] eq 1) then begin


    ;I leave this out for the moment. We average over many frames -> cosmics shouldn't be a problem

;     tmp = csmedcorim
;     la_cosmic_MOD, tmp, outim, path, skyval=median(sk+32768.,/even), gain=100.;, gain=9.8, readn=4.4	;NACO manual, table 5-2, p.29
;     ;the gain=100 value was manually set
; 
;     ;sometimes frames having NANs from la_cosmic, the following removes them
;     idx = where(finite(outim) ne 1)
;     if (idx[0] ne -1) then begin
; 
;       idx0 = array_indices(outim, idx)
;       tmp_bp = tmp
;       tmp_bp[*,*] = 0.
;       for xx=0,n_elements(idx)-1 do tmp_bp[idx0[0,xx], idx0[1,xx]] = 1.
;       fixpix_mod, outim, tmp_bp, tmp_outim, npix=24, /weight, /silent
;       outim = tmp_outim
; 
;     endif
; 
;     csmedcorim = outim

    ;----------------------------------------------------------------------------------

    ;output
    writefits, path+'../Reduced/'+'cube_'+outname[i]+'_reduced.fits', csmedcorim, hdrraw
    writefits, path+'../Reduced/'+'cube_'+outname[i]+'_paral.fits', [PA_onsky]
    writefits, path+'../Reduced/'+'cube_'+outname[i]+'_paral.fits', [ha], /append

  endif

endfor

;=============================================================================================

;output

;find reduced frames and merge them into one file

file = file_search(path+'../Reduced/'+'cube*reduced.fits', count=nred)
if (n_elements(size(corim)) gt 5) then filenosky = file_search(path+'../Reduced/'+'cube*NoSky.fits', count=nnosky)
filepa = file_search(path+'../Reduced/'+'cube*paral.fits', count=npa)

hdr = headfits(file[0], exten=0, /silent)
rnx = strcompress(get_eso_keyword(hdr,'NAXIS1'),/rem)
rny = strcompress(get_eso_keyword(hdr,'NAXIS2'),/rem)
redim = fltarr(rnx,rny,nred)
pa = dblarr(npa)
ha = dblarr(npa)

for i=0,nred-1 do begin

  im = mrdfits(file[i], 0, /silent)
  redim[*,*,i] = im

  ang = mrdfits(filepa[i], 0, /silent)
  pa[i] = ang

  tmp = mrdfits(filepa[i], 1, /silent)
  ha[i] = tmp


endfor

writefits, path+'../Reduced/'+'img_Lp_dc.fits', redim, hdrsciref
writefits, path+'../Reduced/'+'vec_Lp_paral.fits', pa
writefits, path+'../Reduced/'+'vec_Lp_paral.fits', ha, /append

;----------------------------------------------------------------------------------

if (n_elements(size(corim)) gt 5) then begin

  ;write out the number of used frames

  openw, lun, path+'../Reduced/'+'used_frames.txt', width=2000, /get_lun

    ratio = 100.*(double(used_frames[0,*])/double(used_frames[1,*]))

    idx = where(flagcut eq 1)

    printf, lun, '                          File Used Total Ratio%'
    for i=0,nfiles-1 do begin

	if (flagcut[i] eq 1) then printf, lun, outname[i], strcompress(used_frames[0,i],/rem), strcompress(used_frames[1,i],/rem), sigfig(ratio[i],4), format='(a30,2a5,a8)' $
	else printf, lun, outname[i], 0, strcompress(used_frames[1,i],/rem), 0, format='(a30,2a5,a8)'

	printf, lun, ''
	printf, lun, 'Percentage of used frames: '+sigfig(100.*(double(total(used_frames[0,idx]))/double(total(used_frames[1,*]))),4)

    endfor

  close, lun
  free_lun,lun

endif


end
