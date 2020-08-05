pro batch_NACO_Dark, path

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

  st[i] = dit[i]+'_'+naxis1[i]

endfor

idx = where(catg eq 'CALIB' and strmatch(type, '*DARK*') eq 1)
fileo = fileo[idx]
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
udit = dit[idx]
unaxis1 = naxis1[idx]
nu = n_elements(ust)

for xx=0,nu-1 do begin

  idx = where(udit[xx] eq dit and unaxis1[xx] eq naxis1)

  file = fileo[idx]
  tmp = mrdfits(file[0], 0, hdr, /silent)
  sz = size(tmp)
  if (n_elements(sz) eq 5) then sz3 = 1 else sz3 = sz[3]

  tmpdark = dblarr(sz[1], sz[2], n_elements(file))


  for i=0,n_elements(file)-1 do begin

    tmp = mrdfits(file[i], 0, /silent)
    if (sz3 gt 1) then tmpdark[*,*,i] = median(tmp, dimension=3, /even) else tmpdark[*,*,i] = tmp

    ;tmpdark[*,*,i] = tmpdark[*,*,i]+32768.

  endfor

  if (n_elements(file) gt 1) then dark = median(tmpdark, dimension=3, /even) else dark = tmpdark

  if (n_elements(dark[*,0]) lt n_elements(dark[1,*])) then dark = dark[*,0:n_elements(dark[*,0])-1]

  exptime = strcompress(get_eso_keyword(hdr, 'EXPTIME'),/rem)
  exptime = sigfig(exptime,2)

  if (sz3 eq 1) then begin

    writefits, path+'../Reduced/'+'master_dark_'+exptime+'s_'+strcompress(unaxis1[xx])+'px.fits', dark, hdr

  endif else begin

    for j=0,n_elements(hdr)-1 do begin

      if (strmatch(hdr[j], '*NAXIS3*') eq 1) then begin

	newhdr = arrdelete(hdr, at=j, length=1)
	flag = 1

      endif

    endfor

    if (flag eq 1) then hdr = newhdr

    writefits, path+'../Reduced/'+'master_dark_'+exptime+'s_'+strcompress(unaxis1[xx])+'px.fits', dark, hdr

  endelse


endfor


end