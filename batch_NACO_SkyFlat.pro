pro batch_NACO_SkyFlat, path, display

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

idx = uniq(st, sort(st))
if (n_elements(idx) gt 1) then stop
ust = st[idx]

;===============================================================================================

;get BP map

bpf = file_search(path+'../Reduced/'+'static_badpixels*fits', count=nbp)
match = strmatch(bpf, '*'+ust+'*')
idx = where(match eq 1)
if (n_elements(idx) gt 1 or idx[0] eq -1) then stop
bpf = bpf[idx]

;===============================================================================================

bp = mrdfits(bpf, 0, /silent)

nff = n_elements(fff)

airmass = dblarr(nff)
alt = dblarr(nff)
dit = dblarr(nff)
jd = dblarr(nff)

for i=0,nff-1 do begin

  hdr = headfits(fff[i], exten=0)
  airmass[i] = double(get_eso_keyword(hdr, 'AIRMASS'))
  alt[i] = double(get_eso_keyword(hdr, 'HIERARCH ESO TEL ALT'))
  dit[i] = strcompress(get_eso_keyword(hdr,'HIERARCH ESO DET DIT'),/rem)
  jd[i] = double(get_eso_keyword(hdr, 'MJD-OBS'))

  if (i eq 0) then begin

    ffnx = float(get_eso_keyword(hdr, 'NAXIS1'))
    ffny = float(get_eso_keyword(hdr, 'NAXIS2'))

  endif

endfor

;sort frames w.r.t. ascending airmass
ff = dblarr(ffnx,ffny,nff)

idxsort = sort(airmass)
fff = fff[idxsort]
airmass = airmass[idxsort]
alt = alt[idxsort]
jd = jd[idxsort]

if (total(airmass eq 0.)) then begin

  flagair = 0
  idxsort = sort(jd)
  jd = jd[idxsort]
  fff = fff[idxsort]

endif else begin

  flagair = 1

endelse

for i=0,nff-1 do ff[*,*,i] = mrdfits(fff[i], 0, hdrff, /silent)
nx = float(get_eso_keyword(hdrff, 'HIERARCH ESO DET WIN NX'))
; ny = float(get_eso_keyword(hdrff, 'HIERARCH ESO DET WIN NY'))
ny = nx
; ff = ff[sx-1:sx+nx-2,sy-1:sy+ny-2,*]
ff = ff[0:nx-1,0:nx-1,*]


;===============================================================================================

xa = dindgen(nff)
slope = dblarr(nx,ny)

dumerr = dblarr(nff)
dumerr[*] = 1.

; window, 0
for xx=0,nx-1 do begin

  for yy=0,ny-1 do begin

    if (bp[xx,yy] ne 1.) then begin

      if (flagair eq 1) then begin

	sixlin, airmass, ff[xx,yy,*], a, siga, b, sigb
	slope[xx,yy] = b[0]

      endif

      if (flagair eq 0) then begin

	sixlin, jd, ff[xx,yy,*], a, siga, b, sigb
	slope[xx,yy] = b[0]

      endif

    endif else begin

      slope[xx,yy] = -99.

    endelse

  endfor

  proceeding_text,loop=(nx), i=xx, prompt='> x direction   '+string(xx+1,form='(I4)')

endfor

;===============================================================================================

;normalisation

;median slope
idx = where(bp eq 0.d0)
msl = median(slope[idx], /even)
; msl = median(slope, /even)

ff_final = slope/msl

;===============================================================================================

tmp = ff_final
fixpix_mod, tmp, bp, outim, npix=24, /weight, /silent
ff_final = outim

;===============================================================================================

; mwrfits, ff_final, path+'instrument_flat.fits', /silent
writefits, path+'../Reduced/'+'instrument_flat.fits', ff_final

;===============================================================================================

if keyword_set(display) then begin

  window, 0
  !p.multi=[0,1,0]
  ; idx = where(bp eq 0.d0)
  ; plot, slope[idx], psym=3, /yn, xst=1
  plot, slope, psym=3, /yn, xst=1
  !p.multi=[0,1,0]

  window, 2, xs=nx, ys=ny
  cgimage, ff_final, stretch=2

  window, 1, xs=500, ys=500
  cghistoplot, ff_final, xr=[0.8,1.2], bin=0.003, yr=[0,5.d4]

endif


end
