@batch_NACO_Dark.pro
@batch_NACO_BPmap.pro
@batch_NACO_SkyFlat.pro
@batch_NACO_SkyBG.pro
; @batch_NACO_SkyBG_v2.pro
@batch_NACO_reduce_Science.pro
@batch_NACO_reduce_Science_noAGPM.pro
@batch_NACO_reduce_PSF.pro
@batch_NACO_QC.pro
@circint_MOD.pro
@eq2hor_MOD.pro
@sign.pro
@readcol.pro
@remchar.pro
@gettok.pro
@headfits.pro
@fxposit.pro
@mrd_hread.pro
@get_eso_keyword.pro
@sigfig.pro
@mrdfits.pro
@fxpar.pro
@valid_num.pro
@mrd_skip.pro
@arrdelete.pro
@writefits.pro
@check_fits.pro
@fxaddpar.pro
@sxdelpar.pro
@sxpar.pro
@fits_add_checksum.pro
@checksum32.pro
@n_bytes.pro
@is_ieee_big.pro
@host_to_ieee.pro
@get_date.pro
@daycnv.pro
@sxaddpar.pro
@fits_ascii_encode.pro
@mpfitfun.pro
@mpfit.pro
@cgcolor.pro
@cggetcolorstate.pro
@cgsnapshot.pro
@cgcolor24.pro
@array_indices.pro
@robust_mean.pro
@avg.pro
@mkhdr.pro
@cgimage.pro
@image_dimensions.pro
@cgdefcharsize.pro
@setdefaultvalue.pro
@cgdefaultcolor.pro
@cgerase.pro
@cgsetcolorstate.pro
@clipscl.pro
@cgresizeimage.pro
@sixlin.pro
@proceeding_text.pro
@cgscalevector.pro
@fpufix.pro
@cghistoplot.pro                                                                                                                                                                     
@convert_to_type.pro                                                                                                                                                                 
@cgcheckforsymbols.pro     
@dist.pro 
@ten.pro
@closest.pro
@gaussscl.pro
@cgplot.pro
@cgbitget.pro
@colorsareidentical.pro
@fixpix_mod.pro
@dist_circle.pro
@caldat.pro
@precess.pro
@premat.pro
@co_nutate.pro
@nutate.pro
@poly.pro
@cirrange.pro
@isarray.pro
@co_aberration.pro
@sunpos.pro
@addpm.pro
@ct2lst.pro
@hadec2altaz.pro
@co_refract.pro
@parangle.pro
@mpfit2dpeak.pro
@mpfit2dfun.pro
@fftshift.pro
@interpol.pro
@resistant_mean.pro
@scale_image_am.pro
@sky.pro
@mmm.pro
@asinh.pro
@remove.pro
@ts_diff.pro
@batch_NACO_reduce_Science_noAGPM_Mband.pro

pro batch_NACO_reduce, display=display

;===============================================================================================

readcol, '/home/amueller/work/IDLlibs/AO/NACO/datapaths.txt', tmp, format='a', /silent
print, ''
for i=0,n_elements(tmp)-1 do print, strcompress(i+1, /rem), ' ', tmp[i]

print, ''
read, 'Select Path: ', selp
; selp = 1
path = tmp[selp-1]+'RAW/'
path2 = tmp[selp-1]

; pathr = path+'../Reduced/'
spawn, 'mkdir -p '+path+'../Reduced/'

;===============================================================================================

;select star and get its needed properties

filestar = file_search('/home/amueller/work/IDLlibs/AO/TargetProperties/Targets/*.sav', count=nstars)
stars = strarr(nstars)
for i=0,nstars-1 do begin

  stars[i] = strmid(filestar[i], 56, strlen(filestar[i])-56)
  stars[i] = strmid(stars[i], 0, strlen(stars[i])-4)

endfor

nstars = n_elements(stars)

flag = uintarr(nstars)
for i=0,nstars-1 do begin

  match = strmatch(path, '*'+stars[i]+'*')
  flag[i] = total(match)

endfor

idx = where(flag eq 1)
if (idx[0] ne -1) then begin

  star = stars[idx]
  filestar = filestar[idx]

endif else begin

  print, ''
  print, 'No corresponding stellar parameter set found! Stop.'
  stop

endelse

; print, ''
; for i=0,nstars-1 do print, strcompress(i+1, /rem), ' ', stars[i]
; read, 'Select Star?: ', selstar
; ; selstar = 216
; star = stars[selstar-1]
; filestar = filestar[selstar-1]

if (strmatch(path, '*'+star+'*') ne 1) then begin

  print, ''
  print, 'Wrong star selected! Stop.'
  stop

endif

; qintsky = ''
; read, 'SkyPCA (p) or Interpolate Sky frames (y/n): ', qintsky
qintsky = 'p'
;For PSF extraction no PCA is used for sky subtraction because there is only one sky for each set

;===============================================================================================

filter = ''
read, 'Band l/m: ', filter

; plsc = 13.19d-3
plsc = 27.19d-3
diam = 8.2

if (filter eq 'l') then lambda = 3.8d-6	;L'
if (filter eq 'm') then lambda = 4.78d-6	;L'

fwhm = lambda/diam*206265.d0/plsc

;===============================================================================================

;DARK
batch_NACO_Dark, path

;===============================================================================================

;Bad Pixel map
batch_NACO_BPmap, path, display

;===============================================================================================

;Sky Flat
batch_NACO_SkyFlat, path, display

;===============================================================================================

;Sky extraction
agpmflag = strmatch(path, '*noAGPM*')
agpmflag = abs(agpmflag-1)	;no agpm = 0, w/ agpm = 1
batch_NACO_SkyBG, path, agpmflag, filter

;===============================================================================================

;reduce cubes
onlymed = 0	;if set to 1 only median combined images are reduced instead of all single frames in cube
agpmflag = strmatch(path, '*noAGPM*')
agpmflag = abs(agpmflag-1)	;no agpm = 0, w/ agpm = 1
;I leave batch_NACO_reduce_Science untouched in case I need to switch back to the original noAGPM sky subtraction for whatever reason
if (filter eq 'l') then begin
    if (agpmflag eq 1) then batch_NACO_reduce_Science, star, path, agpmflag, filestar, onlymed, fwhm, qintsky, display
    if (agpmflag eq 0) then batch_NACO_reduce_Science_noAGPM, star, path, filestar, onlymed, fwhm, qintsky, display
endif

if (filter eq 'm') then begin

    if (agpmflag eq 0) then batch_NACO_reduce_Science_noAGPM_Mband, star, path, filestar, onlymed, fwhm, qintsky, display
    
endif

;===============================================================================================

;reduce PSF frames
batch_NACO_reduce_PSF, star, path, fwhm, display, filter

;===============================================================================================

;QC
batch_NACO_QC, star, path, fwhm, display, filter

;===============================================================================================

if (filter eq 'l') then begin

    spawn, 'cp '+path+'../Reduced/'+'img_Lp_dc.fits '+path2
    spawn, 'cp '+path+'../Reduced/'+'PSF_Lp.fits '+path2
    spawn, 'cp '+path+'../Reduced/'+'vec_Lp_paral.fits '+path2

endif

if (filter eq 'm') then begin

    spawn, 'cp '+path+'../Reduced/'+'img_Mp_dc.fits '+path2
    spawn, 'cp '+path+'../Reduced/'+'PSF_Mp.fits '+path2
    spawn, 'cp '+path+'../Reduced/'+'vec_Mp_paral.fits '+path2

endif
    
;===============================================================================================


print, ''
print, '*****************************************************'
print, 'CHECK REDUCED IMAGE CUBE, REMOVE BAD FRAMES IF NEEDED'
print, '*****************************************************'
print, ''
print, 'Reduction of '+star+' finished!'
print, ''
stop
end
