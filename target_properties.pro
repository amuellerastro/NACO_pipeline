@QuerySimbad_mod.pro
@queryvizier.pro
@readcol.pro
@remchar.pro
@gettok.pro
@repstr.pro
@webget.pro
@adstring.pro
@radec.pro
@repchr.pro
@interpol.pro
@strnumber.pro
@proceeding_text.pro


pro target_properties

dir = '/home/amueller/work/IDLlibs/AO/TargetProperties/'
resdir = dir+'Targets/'
namefile = 'target_ids.txt'
readcol, namefile, name, format='a';, numline=37;, skipline=2;, /silent

;------------------------------------------------------------------------------

;check for double entries
idxs = sort(name)
name = name[idxs]

idxu = uniq(name)
name = name[idxu]


;------------------------------------------------------------------------------

nid = n_elements(name)

; id = strarr(nid)
; spt = strarr(nid)
; ra = strarr(nid)
; dec = strarr(nid)
; pma = dblarr(nid)
; pmd = dblarr(nid)
; jmag = dblarr(nid)
; hmag = dblarr(nid)
; kmag = dblarr(nid)
; lmag = dblarr(nid)
; ejmag = dblarr(nid)
; ehmag = dblarr(nid)
; ekmag = dblarr(nid)
; elmag = dblarr(nid)
; plx = dblarr(nid)
; eplx = dblarr(nid)
; age = dblarr(nid)
; eage = dblarr(nid)
; radvel = dblarr(nid)

;=====================================================================================================

;check if we have already the target in saved

file = file_search(resdir+'*.sav', count=nfiles)

if (nfiles gt 0) then begin

  id_present = strarr(nfiles)

  pos1 = strpos(file, '/', /reverse_search)
  pos2 = strpos(file, '.sav', /reverse_search)

  for i=0,nfiles-1 do id_present[i] = strmid(file[i], pos1[i]+1, pos2[i]-pos1[i]-1)


  flag = uintarr(nid)

  for i=0,nid-1 do begin

    idx = where(name[i] eq id_present)
    if (idx[0] ne -1) then flag[i] = 1

  endfor

endif else begin

  flag = uintarr(nid)

endelse


;=====================================================================================================

;Allen's Astrophysical Quantities, Arthur N. Cox, 4th edition
;assuming we have only class V objects

spt_ref = ['B0','B1','B2','B3','B4','B5','B6','B7','B8','B9','A0','A2','A5','A7','F0','F2','F5','F7','G0','G2','G4','G6','K0','K2','K4','K5','K7','M0','M1','M2','M3','M4','M5','M6']
spt_refidx = [1,2,3,4,5,6,7,8,9,10,11,13,16,18,21,23,26,28,31,33,35,37,41,43,45,46,48,51,52,53,54,55,56,57]
delta_ref = [0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.01,0.02,0.03,0.03,0.03,0.04,0.04,0.05,0.05,0.05,0.05,0.06,0.07,0.10,0.11,0.13,0.17,0.21,0.23,0.32,0.37,0.42,0.48]

spt_all = ['B0','B1','B2','B3','B4','B5','B6','B7','B8','B9','A0','A1','A2','A3','A4','A5','A6','A7','A8','A9','F0','F1','F2','F3','F4','F5','F6','F7','F8','F9','G0','G1','G2','G3','G4','G5','G6','G7','G8','G9','K0','K1','K2','K3','K4','K5','K6','K7','K8','K9','M0','M1','M2','M3','M4','M5','M6','M7','M8','M9']

spt_idxall = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60]

;=====================================================================================================

for i=0,nid-1 do begin
;for i=107,107 do begin

  if (flag[i] eq 0) then begin

    print, 'Retrieving information for '+name[i]

    tmpname = strcompress(name[i],/rem)
	if (tmpname eq '49Cet') then tmpname = 'V* 49 Cet'
	if (tmpname eq '51Oph') then tmpname = 'V* 51 Oph'
	if (tmpname eq 'ABAur') then tmpname = 'V* AB Aur'
	if (tmpname eq 'AKSco') then tmpname = 'V* AK Sco'
	if (tmpname eq 'ATPyx') then tmpname = 'V* AT Pyx'
	if (tmpname eq 'betaPic') then tmpname = 'beta Pic'
	if (tmpname eq 'BFOri') then tmpname = 'V* BF Ori'
	if (tmpname eq 'DXCha') then tmpname = 'V* DX Cha'
	if (tmpname eq 'EM*SR21A') then tmpname = 'EM* SR 21A'
	if (tmpname eq 'EM*SR24S') then tmpname = 'EM* SR 24S'
	if (tmpname eq 'GGTau') then tmpname = 'V* GG Tau'
	if (tmpname eq 'GMAur') then tmpname = 'V* GM Aur'
	if (tmpname eq 'HLTau') then tmpname = 'V* HL Tau'
	if (tmpname eq 'IMLup') then tmpname = 'V* IM Lup'
	if (tmpname eq 'KKOph') then tmpname = 'V* KK Oph'
	if (tmpname eq 'MYLup') then tmpname = 'V* MY Lup'
	if (tmpname eq 'NXPup') then tmpname = 'V* NX Pup'
	if (tmpname eq 'RCrA') then tmpname = 'V* R CrA'
	if (tmpname eq 'RYLup') then tmpname = 'V* RY Lup'
	if (tmpname eq 'RYTau') then tmpname = 'V* RY Tau'
	if (tmpname eq 'TCha') then tmpname = 'V* T Cha'
	if (tmpname eq 'TCrA') then tmpname = 'V* T CrA'
	if (tmpname eq 'TWHya') then tmpname = 'V* TW Hya'
	if (tmpname eq 'UXOri') then tmpname = 'V* UX Ori'
	if (tmpname eq 'UXTauA') then tmpname = 'HBC 43'
    
    

    ;----------------------------------------------------------------------------------

    tjmag = 0.
    thmag = 0.
    tkmag = 0.
    tspt = ''
    trv = 0.
    jmag = 0.
    hmag = 0.
    kmag = 0.
    spt = ''
    rv = 0.
    tra = 0.
    tdec = 0.
    tprimid = ''
    tepmra = 0.
    tepmde = 0.

    QuerySimbad_mod, tmpname, tra, tdec, tprimid, found=tfound, Jmag=tJmag, Hmag=tHmag, Kmag=tKmag, pmra=tpmra, pmde=tpmdec, epmra=tepmra, epmde=tepmde, parallax=tparallax, eparallax=teparallax, spt=tspt, rv=trv, result=result, status=status;, /verbose

    if (status eq -1) then begin

      read, 'Provide different target name: ', tmpname
      QuerySimbad_mod, tmpname, tra, tdec, tprimid, found=tfound, Jmag=tJmag, Hmag=tHmag, Kmag=tKmag, pmra=tpmra, pmde=tpmdec, epmra=tepmra, epmde=tepmde, parallax=tparallax, eparallax=teparallax, spt=tspt, rv=trv, result=result, status=status, /verbose

      if (status eq -1) then begin

	print, 'Something is fucked up...'
	stop

      endif

    endif

    id = name[i]
    spt = tspt[0]
    radvel = trv

    pma = tpmra/1.d3
    pmd = tpmdec/1.d3

    epma = tepmra/1.d3
    epmd = tepmde/1.d3


    coord = adstring([tra, tdec], precision=3)
    ra = strcompress(strmid(coord, 0, 15),/rem)
    dec = strcompress(strmid(coord, 15, 15),/rem)

    plx = tparallax/1.d3
    eplx = teparallax/1.d3

    ;----------------------------------------------------------------------------------

  ;   dis = 0.15
  ;   if (name[i] eq 'Sirius') then dis = 0.5
  ; 
  ;   ;HIP2
  ;   st = Queryvizier('I/311/HIP2', tmpname, dis);, /verbose)
  ; 
  ; ;   if (size(st,/dimensions) gt 0) then begin
  ; 
  ;     if (size(st, /type) eq 8) then begin
  ; 
  ;       if (n_elements(st.plx) gt 1) then begin
  ; 	print, 'More than one result. Stop.'
  ; 	stop
  ;       endif
  ; 
  ; ;       pma = st.pmra/1.d3
  ; ;       pmd = st.pmde/1.d3
  ;       plx = st.plx/1.d3
  ;       eplx = st.e_plx/1.d3
  ; 
  ;     endif else begin
  ; 
  ;       print, ''
  ;       print, tmpname
  ;       print, 'Provide following values manually!'
  ;       read, 'PM RA [mas]: ', pma
  ;       read, 'PM DEC [mas]: ', pmd
  ;       read, 'Parallax [mas]: ', plx
  ;       read, 'Error Parallax [mas]: ', eplx
  ; 
  ;       pma = pma/1.d3
  ;       pmd = pmd/1.d3
  ;       plx = plx/1.d3
  ;       eplx = eplx/1.d3
  ; 
  ;     endelse

  ;   endif

;     if (name[i] eq 'SAO150676') then plx = 1./78.	;vlaue from ISPY table
;     if (name[i] eq 'TWHya') then plx = 17.8571/1d3	;56pc
;     if (name[i] eq 'TCha') then begin
;       plx = 1./108.	;Torres2008
;       eplx = 7.8d-4	;+/-9pc
;     endif

    ;----------------------------------------------------------------------------------


    ;2MASS
    dis = 0.15
    st = Queryvizier('II/246', tmpname, dis);, /verbose)

    if (size(st, /type) eq 8) then begin

      if (n_elements(st.kmag) gt 1) then begin

	print, ''
	print, tmpname
	print, 'More than one result for K band magnitude. Select correct one:'
	for xx=0,n_elements(st.jmag)-1 do print, strcompress(xx+1,/rem), st[xx].kmag, format='(a3, a10)'
	read, 'Select correct value: ', cv

	jmag = st[cv-1].jmag
	hmag = st[cv-1].hmag
	kmag = st[cv-1].kmag
	ejmag = st[cv-1].e_jmag
	ehmag = st[cv-1].e_hmag
	ekmag = st[cv-1].e_kmag

  ;     print, 'More than one result. Stop.'
  ;     stop
      endif else begin

      jmag = st.jmag
      hmag = st.hmag
      kmag = st.kmag
      ejmag = st.e_jmag
      ehmag = st.e_hmag
      ekmag = st.e_kmag

      endelse

    endif else begin

      print, ''
      print, tmpname
      print, 'Provide following values manually!'
      read, 'J mag: ', jmag
      read, 'J mag error: ', ejmag
      read, 'H mag: ', hmag
      read, 'H mag error: ', ehmag
      read, 'K mag: ', kmag
      read, 'K mag error: ', ekmag

    endelse


    ;----------------------------------------------------------------------------------

    ;get L' mag
    tmpspt = strmid(spt, 0, 2)
    idx = where(tmpspt eq spt_all)
    if (idx[0] eq -1) then begin

      print, ''
      print, tmpname
      print, 'Check SpT. Using e.g. Simbad Measurements or Wiki or Astrid.'
      print, spt
      spt = ''
      read, 'Enter new SpT: ', spt

      idx = where(spt eq spt_all)

      spt_idx = spt_idxall[idx]
      delta = interpol(delta_ref, spt_refidx, spt_idx)
      lmag = kmag-delta
      elmag = ekmag

  ;     print, 'Check SpT. Stop.'
  ;     stop

    endif else begin

      spt_idx = spt_idxall[idx]
      delta = interpol(delta_ref, spt_refidx, spt_idx)
      lmag = kmag-delta
      elmag = ekmag

      if (name[i] eq 'betaPic') then begin
	lmag = 3.454	;Absil2013
	elmag = 0.003
      endif

    endelse

    ;----------------------------------------------------------------------------------

    age = 0.
    eage = 0.

    if (name[i] eq 'EM*SR24S') then name[i] = 'EMSR24S'
    if (name[i] eq 'EM*SR21A') then name[i] = 'EMSR21A'

    fn = resdir+name[i]+'.sav'
    save, id, ra, dec, spt, pma, pmd, epma, epmd, plx, eplx, jmag, ejmag, hmag, ehmag, kmag, ekmag, lmag, elmag, radvel, age, eage, filename=fn
;     print,  id, ra, dec, spt, pma, pmd, epma, epmd, plx, eplx, jmag, ejmag, hmag, ehmag, kmag, ekmag, lmag, elmag, radvel, age, eage
;     stop
    proceeding_text,loop=nid, i=i, prompt='> Target   '+string(i+1,form='(I4)')

  endif

endfor


stop
end
