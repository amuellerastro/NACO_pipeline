pro batch_NACO_SkyBG, path, agpmflag, filter

fileo = file_search(path+'NACO*fits', count=nfiles)

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

st = strarr(nfiles)

for i=0,nfiles-1 do begin

  pos1 = strpos(fileo[i], '/', /reverse_search)
  pos2 = strpos(fileo[i], '.fits', /reverse_search)
;   out[i] = strmid(fileo[i], pos1+1, pos2-pos1-1)

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
  ndit[i] = strcompress(get_eso_keyword(hdr, 'HIERARCH ESO DET NDIT'),/rem)

  st[i] = dit[i]+'s_'+naxis1[i]+'px'

endfor

;===============================================================================================

;NO AGPM used, sky from science frames itself

if (agpmflag eq 0) then begin

  idx = where(catg eq 'SCIENCE')
  file = fileo[idx]
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
  st = st[idx]

  idx = uniq(st, sort(st))
  ust = st[idx]
  nu = n_elements(idx)
  if (nu gt 2) then begin

    print, ''
    print, 'More than 2 sets of different PSF files. There should be only 1 set of PSF and 1 set of science! Stop.'
    stop

  endif

  nx_all = naxis1[uniq(naxis1, sort(naxis1))]

  for xx=0,nu-1 do begin

    nx = nx_all[xx]

    ;PSF frames always have <=256px, science >256px
    idx = where(naxis1 eq nx)
    rawf = file[idx]
    nfiles = n_elements(rawf)
    ny = nx

    dateim = strarr(nfiles)
    jd = dblarr(nfiles)
    reloffx = dblarr(nfiles)
    reloffy = dblarr(nfiles)

    for i=0,nfiles-1 do begin

      hdr = headfits(rawf[i], exten=0, /silent)
      reloffx[i] = get_eso_keyword(hdr, 'HIERARCH ESO SEQ RELOFFSETX')
      reloffy[i] = get_eso_keyword(hdr, 'HIERARCH ESO SEQ RELOFFSETY')
      dateim[i] = get_eso_keyword(hdr, 'DATE')
      jd[i] = double(get_eso_keyword(hdr, 'MJD-OBS'))+2400000.5d0
      if (i eq 0) then ditsk = sigfig(get_eso_keyword(hdr,'HIERARCH ESO DET DIT'),2)

    endfor

    reloffs = strcompress(reloffx,/rem)+strcompress(reloffy,/rem)
    ureloffs = reloffs[uniq(reloffs, sort(reloffs))]

    bpf = path+'../Reduced/'+'static_badpixels_'+ditsk+'s_'+strcompress(nx,/rem)+'px.fits'
    bp = mrdfits(bpf, 0, /silent)

    pos = strpos(rawf, '/', /reverse_search)
    rfname = strmid(rawf, pos[0]+1, strlen(rawf[0])-pos[0]-1)
    flag = uintarr(n_elements(rawf))

    ;in case of PSF frames we cannot do PCA subtraction because of lack of frames. Therefore we create an average sky frame for each PSF frame depending on detector position.

    if (filter eq 'l') then thsz = 256.
    if (filter eq 'm') then thsz = 128.
    if (nx gt thsz) then begin

      for io=0,n_elements(ureloffs)-1 do begin

        idx = where(reloffs ne ureloffs[io])

        idxq = where(reloffs eq ureloffs[io])
        flag[idxq] = io+1

        for i=0,n_elements(idx)-1 do begin

        im = mrdfits(rawf[idx[i]],0,head,/silent)

        im = im[*,0:n_elements(im[*,0,0])-1,*]

        if (n_elements(size(im)) gt 5) then medim = median(im[*,*,30:n_elements(im[0,0,*])-2], dim=3, /even) else medim = im
        tmpdate = dateim[idx[i]]

        tmp = medim
        fixpix_mod, tmp, bp, outim, npix=24, /weight, /silent

        sxaddpar, head, 'JD', jd[idx[i]], format='(f18.10)'

        skyname = path+'../Reduced/'+'sky_'+strcompress(io+1,/rem)+'_background_'+tmpdate+'.fits'
        writefits, skyname, outim, head

        endfor

      endfor

      fn = path+'../Reduced/'+'sky_ident.sav'
      save, rfname, flag, filename=fn

    endif else begin

      medim = dblarr(nx,ny,nfiles)

      for i=0,nfiles-1 do begin

        im = mrdfits(rawf[i],0,/silent)
        im = im[*,0:n_elements(im[*,0,0])-1,30:n_elements(im[0,0,*])-2]
        medim[*,*,i] = median(im, dim=3, /even)

        tmp = reform(medim[*,*,i])
        fixpix_mod, tmp, bp, outim, npix=24, /weight, /silent

        medim[*,*,i] = outim

      endfor

      if (nfiles/2. eq 2.) then begin

        flag = [1,2,3,4]

        for i=0,nfiles-1 do begin

        head = headfits(rawf[i], exten=0, /silent)
        sxaddpar, head, 'JD', jd[i], format='(f18.10)'
        skyname = path+'../Reduced/'+'PSF_sky_'+strcompress(flag[i],/rem)+'_background_'+dateim[i]+'.fits'

        if (i eq 0) then writefits, skyname, medim[*,*,1], head
        if (i eq 1) then writefits, skyname, medim[*,*,0], head
        if (i eq 2) then writefits, skyname, medim[*,*,3], head
        if (i eq 3) then writefits, skyname, medim[*,*,2], head

        endfor

      endif

      if (nfiles/2. eq 3.) then begin

        flag = [1,2,3,4,5,6]

        for i=0,nfiles-1 do begin

        head = headfits(rawf[i], exten=0, /silent)
        sxaddpar, head, 'JD', jd[i], format='(f18.10)'
        skyname = path+'../Reduced/'+'PSF_sky_'+strcompress(flag[i],/rem)+'_background_'+dateim[i]+'.fits'

        if (i eq 0) then writefits, skyname, median(medim[*,*,1:2],dim=3, /even), head
        if (i eq 1) then writefits, skyname, median(medim[*,*,[0,2]],dim=3, /even), head
        if (i eq 2) then writefits, skyname, median(medim[*,*,0:1],dim=3, /even), head
        if (i eq 3) then writefits, skyname, median(medim[*,*,4:5],dim=3, /even), head
        if (i eq 4) then writefits, skyname, median(medim[*,*,[3,5]],dim=3, /even), head
        if (i eq 5) then writefits, skyname, median(medim[*,*,3:4],dim=3, /even), head

        endfor

      endif

      if (nfiles eq 3) then begin

        flag = [1,2,3]

        for i=0,nfiles-1 do begin

        head = headfits(rawf[i], exten=0, /silent)
        sxaddpar, head, 'JD', jd[i], format='(f18.10)'
        skyname = path+'../Reduced/'+'PSF_sky_'+strcompress(flag[i],/rem)+'_background_'+dateim[i]+'.fits'

        if (i eq 0) then writefits, skyname, median(medim[*,*,1:2],dim=3, /even), head
        if (i eq 1) then writefits, skyname, median(medim[*,*,[0,2]],dim=3, /even), head
        if (i eq 2) then writefits, skyname, median(medim[*,*,0:1],dim=3, /even), head

        endfor

      endif

      if (nfiles/2. eq 4.) then begin

	flag = [1,2,3,4,5,6,7,8]

	for i=0,nfiles-1 do begin

	  head = headfits(rawf[i], exten=0, /silent)
	  sxaddpar, head, 'JD', jd[i], format='(f18.10)'
	  skyname = path+'../Reduced/'+'PSF_sky_'+strcompress(flag[i],/rem)+'_background_'+dateim[i]+'.fits'

	  if (i eq 0) then writefits, skyname, median(medim[*,*,1:3],dim=3, /even), head
	  if (i eq 1) then writefits, skyname, median(medim[*,*,[0,[2:3]]],dim=3, /even), head
	  if (i eq 2) then writefits, skyname, median(medim[*,*,[[0:1],3]],dim=3, /even), head
	  if (i eq 3) then writefits, skyname, median(medim[*,*,0:2],dim=3, /even), head
	  if (i eq 4) then writefits, skyname, median(medim[*,*,5:7],dim=3, /even), head
	  if (i eq 5) then writefits, skyname, median(medim[*,*,[4,[6:7]]],dim=3, /even), head
	  if (i eq 6) then writefits, skyname, median(medim[*,*,[[4:5],7]],dim=3, /even), head
	  if (i eq 7) then writefits, skyname, median(medim[*,*,4:6],dim=3, /even), head

	endfor

      endif

      fn = path+'../Reduced/'+'PSF_sky_ident.sav'
      save, rfname, flag, filename=fn

    endelse

  endfor

endif

;===============================================================================================
;===============================================================================================

;WITH AGPM
;AGPM data have manual offset applied (hence not present in fits header)

if (agpmflag eq 1) then begin

  idx = where(catg eq 'SCIENCE')
  file = fileo[idx]
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
  st = st[idx]

  idx = uniq(st, sort(st))
  ust = st[idx]
  nu = n_elements(idx)
  if (nu gt 2) then begin

    print, ''
    print, 'More than 2 sets of different PSF files. There should be only 1 set of PSF and 1 set of science! Stop.'
    stop

  endif

  ;-----------------------------------------------------------------------------------------------
  ;sky PSF frames
  ;the problem here is that the offsets are applied manually and are not written in the header
  ;So I am assuming that each PSF cube has the star at a different quadrant
  ;ASSUMPTION: naxis1 <= 256px

  idx = where(naxis1 le 256)

  if (idx[0] ne -1) then begin

    rawf = file[idx]
    nx = naxis1[idx[0]]
    nfiles = n_elements(rawf)
  ;   if ((nfiles mod 2) ne 0) then begin
  ; 
  ;     print, ''
  ;     print, 'The number of PSF frames is odd! This should not happen. Stop.'
  ;     stop
  ; 
  ;   endif

    ny = nx

    dateim = strarr(nfiles)
    jd = dblarr(nfiles)
  ;   reloffx = dblarr(nfiles)
  ;   reloffy = dblarr(nfiles)

    for i=0,nfiles-1 do begin

      hdr = headfits(rawf[i], exten=0, /silent)
  ;     reloffx[i] = get_eso_keyword(hdr, 'HIERARCH ESO SEQ RELOFFSETX')
  ;     reloffy[i] = get_eso_keyword(hdr, 'HIERARCH ESO SEQ RELOFFSETY')
      dateim[i] = get_eso_keyword(hdr, 'DATE')
      jd[i] = double(get_eso_keyword(hdr, 'MJD-OBS'))+2400000.5d0
      if (i eq 0) then ditsk = sigfig(get_eso_keyword(hdr,'HIERARCH ESO DET DIT'),2)

    endfor

  ;   reloffs = strcompress(reloffx,/rem)+strcompress(reloffy,/rem)
  ;   ureloffs = reloffs[uniq(reloffs, sort(reloffs))]


    bpf = path+'../Reduced/'+'static_badpixels_'+ditsk+'s_'+strcompress(nx,/rem)+'px.fits'
    bp = mrdfits(bpf, 0, /silent)

    pos = strpos(rawf, '/', /reverse_search)
    rfname = strmid(rawf, pos[0]+1, strlen(rawf[0])-pos[0]-1)
    flag = uintarr(n_elements(rawf))

    medim = dblarr(nx,ny,nfiles)

    for i=0,nfiles-1 do begin

      im = mrdfits(rawf[i],0,/silent)
      im = im[*,0:n_elements(im[*,0,0])-1,30:n_elements(im[0,0,*])-2]
      medim[*,*,i] = median(im, dim=3, /even)

      tmp = reform(medim[*,*,i])
      fixpix_mod, tmp, bp, outim, npix=24, /weight, /silent

      medim[*,*,i] = outim

    endfor

    if (nfiles/2. eq 2.) then begin

      flag = [1,2,3,4]

      for i=0,nfiles-1 do begin

	head = headfits(rawf[i], exten=0, /silent)
	sxaddpar, head, 'JD', jd[i], format='(f18.10)'
	skyname = path+'../Reduced/'+'PSF_sky_'+strcompress(flag[i],/rem)+'_background_'+dateim[i]+'.fits'

	if (i eq 0) then writefits, skyname, medim[*,*,1], head
	if (i eq 1) then writefits, skyname, medim[*,*,0], head
	if (i eq 2) then writefits, skyname, medim[*,*,3], head
	if (i eq 3) then writefits, skyname, medim[*,*,2], head

      endfor

    endif

    if (nfiles/2. eq 3.) then begin

      flag = [1,2,3,4,5,6]

      for i=0,nfiles-1 do begin

	head = headfits(rawf[i], exten=0, /silent)
	sxaddpar, head, 'JD', jd[i], format='(f18.10)'
	skyname = path+'../Reduced/'+'PSF_sky_'+strcompress(flag[i],/rem)+'_background_'+dateim[i]+'.fits'

	if (i eq 0) then writefits, skyname, median(medim[*,*,1:2],dim=3, /even), head
	if (i eq 1) then writefits, skyname, median(medim[*,*,[0,2]],dim=3, /even), head
	if (i eq 2) then writefits, skyname, median(medim[*,*,0:1],dim=3, /even), head
	if (i eq 3) then writefits, skyname, median(medim[*,*,4:5],dim=3, /even), head
	if (i eq 4) then writefits, skyname, median(medim[*,*,[3,5]],dim=3, /even), head
	if (i eq 5) then writefits, skyname, median(medim[*,*,3:4],dim=3, /even), head

      endfor

    endif

    if (nfiles eq 3) then begin

      flag = [1,2,3]

      for i=0,nfiles-1 do begin

	head = headfits(rawf[i], exten=0, /silent)
	sxaddpar, head, 'JD', jd[i], format='(f18.10)'
	skyname = path+'../Reduced/'+'PSF_sky_'+strcompress(flag[i],/rem)+'_background_'+dateim[i]+'.fits'

	if (i eq 0) then writefits, skyname, median(medim[*,*,1:2],dim=3, /even), head
	if (i eq 1) then writefits, skyname, median(medim[*,*,[0,2]],dim=3, /even), head
	if (i eq 2) then writefits, skyname, median(medim[*,*,0:1],dim=3, /even), head

      endfor

    endif

    if (nfiles/2. eq 4.) then begin

      flag = [1,2,3,4,5,6,7,8]

      for i=0,nfiles-1 do begin

	head = headfits(rawf[i], exten=0, /silent)
	sxaddpar, head, 'JD', jd[i], format='(f18.10)'
	skyname = path+'../Reduced/'+'PSF_sky_'+strcompress(flag[i],/rem)+'_background_'+dateim[i]+'.fits'

	if (i eq 0) then writefits, skyname, median(medim[*,*,1:3],dim=3, /even), head
	if (i eq 1) then writefits, skyname, median(medim[*,*,[0,[2:3]]],dim=3, /even), head
	if (i eq 2) then writefits, skyname, median(medim[*,*,[[0:1],3]],dim=3, /even), head
	if (i eq 3) then writefits, skyname, median(medim[*,*,0:2],dim=3, /even), head
	if (i eq 4) then writefits, skyname, median(medim[*,*,5:7],dim=3, /even), head
	if (i eq 5) then writefits, skyname, median(medim[*,*,[4,[6:7]]],dim=3, /even), head
	if (i eq 6) then writefits, skyname, median(medim[*,*,[[4:5],7]],dim=3, /even), head
	if (i eq 7) then writefits, skyname, median(medim[*,*,4:6],dim=3, /even), head

      endfor

    endif

    fn = path+'../Reduced/'+'PSF_sky_ident.sav'
    save, rfname, flag, filename=fn

  endif

  ;-----------------------------------------------------------------------------------------------
  ;Sky for AGPM frames, can be identified using nexp, offsets not written in header because manual applied
  ;since P98 SKY keyword can be present if correct OB used

  idx = where(catg eq 'SCIENCE' and naxis1 ge 512 and type eq 'SKY')
  rawf = file[idx]
  naxis1 = naxis1[idx]


  nexp = uintarr(n_elements(rawf))
  type = strarr(n_elements(rawf))
  for i=0,n_elements(rawf)-1 do begin

    hdr = headfits(rawf[i], exten=0, /silent)
    nexp[i] = uint(get_eso_keyword(hdr, 'HIERARCH ESO SEQ NEXPO'))
    type[i] = strcompress(get_eso_keyword(hdr,'HIERARCH ESO DPR TYPE'),/rem)

  endfor

  idx = where(nexp lt median(nexp) or type eq 'SKY')
  
  nsk = n_elements(idx)
  skf = rawf[idx]
  nx = naxis1[idx[0]]
  ny = nx
;   nsk = n_elements(idx)
;   skf = rawf
;   nx = 768.;naxis1[idx[0]]
;   ny = nx

  sky = fltarr(nx,ny,nsk)
  date = strarr(nsk)
  jd = dblarr(nsk)
  ditsk = dblarr(nsk)
  nsets = nsk

  for i=0,nsk-1 do begin

    tmp = mrdfits(skf[i], 0, hdr, /silent)
    tmp = tmp[*,0:nx-1,*]
    jd[i] = double(get_eso_keyword(hdr, 'MJD-OBS'))+2400000.5d0
    date[i] = get_eso_keyword(hdr, 'DATE')
    ditsk[i] = double(get_eso_keyword(hdr,'HIERARCH ESO DET DIT'))
    if (n_elements(tmp[0,0,*]) gt 80) then tmp = tmp[*,*,30:n_elements(tmp[0,0,*])-2] else tmp = tmp[*,*,1:n_elements(tmp[0,0,*])-2]	;remove 1st and last frame

    bpf = path+'../Reduced/'+'static_badpixels_'+sigfig(ditsk[i],2)+'s_'+strcompress(nx,/rem)+'px.fits'
    bp = mrdfits(bpf, 0, /silent)

    sky[*,*,i] = median(tmp, dimension=3, /even)	;+32768.

    tmp = sky[*,*,i]
    fixpix_mod, tmp, bp, outim, npix=24, /weight, /silent
    sky[*,*,i] = outim

    mkhdr, head, sky[*,*,i]
    sxaddpar, head, 'JD', jd[i], format='(f18.10)'
    sxaddpar, head, 'DETDIT', ditsk[i]
    writefits, path+'../Reduced/'+'sky_background_'+date[i]+'.fits', sky[*,*,i], head

    proceeding_text,loop=nsk, i=i, prompt='> Sky Frames Science  '+string(i+1,form='(I4)')

  endfor

  ;-----------------------------------------------------------------------------------------------


endif

;===============================================================================================


end
