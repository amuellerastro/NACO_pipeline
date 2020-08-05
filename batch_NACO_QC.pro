pro batch_NACO_QC, star, path, fwhm, display, filter

if (filter eq 'l') then begin

    im_orig = mrdfits(path+'../Reduced/'+'img_Lp_dc.fits', 0, hdrref, /silent)
    paral = mrdfits(path+'../Reduced/'+'vec_Lp_paral.fits', 0, /silent)
    ha  = mrdfits(path+'../Reduced/'+'vec_Lp_paral.fits', 1, /silent)

endif

if (filter eq 'm') then begin

    im_orig = mrdfits(path+'../Reduced/'+'img_Mp_dc.fits', 0, hdrref, /silent)
    paral = mrdfits(path+'../Reduced/'+'vec_Mp_paral.fits', 0, /silent)
    ha  = mrdfits(path+'../Reduced/'+'vec_Mp_paral.fits', 1, /silent)

endif

;=====================================================================================

im = im_orig

nframes = n_elements(im[0,0,*])
medim = median(im, dim=3)

; val = dblarr(nframes)
; for i=0,nframes-1 do val[i] = ms_ssim(medim, im[*,*,i])
; resistant_mean, val, 3, t1, t2, nbad, /double, goodvec=idxgood, badvec=idxbad
; 

sim = dblarr(nframes)
for i=0,nframes-1 do sim[i] = stddev(im[*,*,i]-medim)

resistant_mean, sim, 3, t1, t2, nbad, /double, goodvec=idxgood, badvec=idxbad

idxbad1 = idxbad

index = intarr(nframes)
index[idxbad1] = 1

;-------------------------------------------------------------------------------------

mval = dblarr(n_elements(im[0,0,*]))
dim = n_elements(im[*,0,0])
;for i=0,nframes-1 do mval[i] = max(im[*,*,i])
for i=0,nframes-1 do mval[i] = max(im[dim/2,dim/2,i])
resistant_mean, mval, 3, t1, t2, nbad, /double, goodvec=idxgood, badvec=idxbad

idxbad2 = idxbad
index[idxbad2] = 1


idxgood = where(index eq 0)
idxbad = where(index eq 1)

;-------------------------------------------------------------------------------------

if keyword_set(display) then begin

  window, 0, title=star, xs=700, ys=600
  !p.multi=[0,1,2]

  x = dindgen(nframes)
  plot, x, sim, /yn, title=star, xst=1
  oplot, x[idxgood], sim[idxgood], psym=2, color=cgcolor('green')
  if (idxbad[0] ne -1) then oplot, x[idxbad], sim[idxbad], psym=2, color=cgcolor('red')

  x = dindgen(nframes)
  plot, x, mval, /yn, title=star, xst=1
  oplot, x[idxgood], mval[idxgood], psym=2, color=cgcolor('green')
  if (idxbad[0] ne -1) then oplot, x[idxbad], mval[idxbad], psym=2, color=cgcolor('red')

  !p.multi=[0,1,2]

endif

if (filter eq 'l') then begin

    if (idxbad[0] ne -1) then begin

    writefits, path+'../Reduced/'+'img_BAD.fits', im[*,*,idxbad]

    spawn, 'cp '+path+'../Reduced/'+'img_Lp_dc.fits '+path+'../Reduced/'+'img_Lp_dc.fits.orig'
    spawn, 'cp '+path+'../Reduced/'+'vec_Lp_paral.fits '+path+'../Reduced/'+'vec_Lp_paral.fits.orig'

    writefits, path+'../Reduced/'+'img_Lp_dc.fits', im[*,*,idxgood], hdrref
    writefits, path+'../Reduced/'+'vec_Lp_paral.fits', paral[idxgood]
    writefits, path+'../Reduced/'+'vec_Lp_paral.fits', ha[idxgood], /append

    endif

endif

if (filter eq 'm') then begin

    if (idxbad[0] ne -1) then begin

    writefits, path+'../Reduced/'+'img_BAD.fits', im[*,*,idxbad]

    spawn, 'cp '+path+'../Reduced/'+'img_Mp_dc.fits '+path+'../Reduced/'+'img_Mp_dc.fits.orig'
    spawn, 'cp '+path+'../Reduced/'+'vec_Mp_paral.fits '+path+'../Reduced/'+'vec_Mp_paral.fits.orig'

    writefits, path+'../Reduced/'+'img_Mp_dc.fits', im[*,*,idxgood], hdrref
    writefits, path+'../Reduced/'+'vec_Mp_paral.fits', paral[idxgood]
    writefits, path+'../Reduced/'+'vec_Mp_paral.fits', ha[idxgood], /append

    endif

endif


; ;========================================================================================
; 
; nframes = n_elements(im[0,0,*])
; norig = nframes
; nx = n_elements(im[*,0,0])
; ny = nx
; 
; flag = uintarr(nframes)
; flag[*] = 1
; 
; ;check for NaNs
; for i=0,nframes-1 do begin
; 
;   idx = where(finite(im[*,*,i]) ne 1)
;   if (idx[0] ne -1) then begin
; 
;     flag[i] = 0
;     print, 'NANs detected in frame '+strcompress(i,/rem)+' (IDL index)'
; 
;   endif
; 
; endfor
; 
; idx0 = where(flag eq 0)
; if (idx0[0] ne -1) then imbad = im[*,*,idx0]
; idx1 = where(flag eq 1)
; im = im[*,*,idx1]
; paral = paral[idx1]
; 
; nframes = n_elements(idx1)
; flag = uintarr(nframes)
; 
; 
; ;=====================================================================================
; 
; sumim = dblarr(nframes)	;simple summation of all pixels in 
; medim = dblarr(nframes)
; sdevim = dblarr(nframes)
; sumc = dblarr(nframes)
; medc = dblarr(nframes)
; maxc = dblarr(nframes)
; sdc = dblarr(nframes)
; 
; for i=0,nframes-1 do begin
; 
;   sumim[i] = total(im[*,*,i])/nx^2.
;   medim[i] = median(im[*,*,i])
;   sdevim[i] = stddev(im[*,*,i])
; 
;   tmp = im[*,*,i]
;   circint_MOD, tmp, nx/2., ny/2., 2.*fwhm, tot, mtot, meantot, maxpx, sdtot, npix, totdif, npixdif, t8
;   sumc[i] = tot 
;   medc[i] = mtot
;   maxc[i] = maxpx
;   sdc[i] = sdtot
; 
; endfor
; 
; 
; if (median(sumim) lt 0.) then begin
; 
;   print, ''
;   print, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
;   print, 'Average of sum of intensity is NEGATIVE!'
;   print, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
;   print, ''
; 
; endif
; 
; resistant_mean, medc, 1., t1, t2, nbad, /double, GoodVec=good
; sumc = sumc[good]
; medc = medc[good]
; maxc = maxc[good]
; sdc = sdc[good]
; sumim = sumim[good]
; medim = medim[good]
; sdevim = sdevim[good]
; flag[good] = 1
; 
; nframes = n_elements(good)
; 
; 
; idx0 = where(flag eq 0)
; if (idx0[0] ne -1) then begin
; 
; ;   t1 = imbad
; ;   t2 = im[*,*,idx0]
; ; 
; ;   imbad = dblarr(nx,ny,n_elements(t1[0,0,*])+n_elements(t2[0,0,*]))
; ;   imbad[*,*,0:]
; ; 
; ; stop
;   writefits, path+'img_BAD.fits', im[*,*,idx0]
; 
;   spawn, 'cp '+path+'img_Lp_dc.fits '+path+'img_Lp_dc.fits.orig'
;   spawn, 'cp '+path+'vec_Lp_paral.fits '+path+'vec_Lp_paral.fits.orig'
; 
;   idx1 = where(flag eq 1)
;   writefits, path+'img_Lp_dc.fits', im[*,*,idx1]
;   writefits, path+'vec_Lp_paral.fits', paral[idx1]
; 
; endif
; 
; print, ''
; print, 'Removed '+strcompress(norig-nframes,/rem)+' / '+strcompress(norig, /rem)+' frames'
; 
; ;========================================================================================

print, ''
print, 'Removed '+strcompress(n_elements(idxbad),/rem)+' / '+strcompress(nframes, /rem)+' frames'


end
