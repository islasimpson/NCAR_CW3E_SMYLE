pro u120e_predictability_1980_2020_between_low_high


fid1=NCDF_OPEN('G:\My Drive\Isla\Isla_delta_predict_wind_section_120E_result.nc')
varid=ncdf_varid(fid1,'high_cor')
ncdf_varget,fid1,varid,high_cor
varid=ncdf_varid(fid1,'low_cor')
ncdf_varget,fid1,varid,low_cor
varid=ncdf_varid(fid1,'dif_cor')
ncdf_varget,fid1,varid,dif_cor
varid=ncdf_varid(fid1,'lat')
ncdf_varget,fid1,varid,lat
varid=ncdf_varid(fid1,'level')
ncdf_varget,fid1,varid,level
NCDF_CLOSE,fid1

high_cor1=fltarr(91,34)
for i =0,90 do begin
  for j=0,33 do begin
    high_cor1[i,j]=mean(high_cor[*,i,j])
  endfor
endfor

low_cor1=fltarr(91,34)
for i =0,90 do begin
  for j=0,33 do begin
    low_cor1[i,j]=mean(low_cor[*,i,j])
  endfor
endfor

dif_cor1=fltarr(91,34)
for i =0,90 do begin
  for j=0,33 do begin
    dif_cor1[i,j]=mean(dif_cor[*,i,j])
  endfor
endfor




fid1=NCDF_OPEN('G:\My Drive\Isla\Isla_delta_predict_wind_section_120E_quantile95.nc')
varid=ncdf_varid(fid1,'min95_low_cor')
ncdf_varget,fid1,varid,min95_low_cor
varid=ncdf_varid(fid1,'max95_low_cor')
ncdf_varget,fid1,varid,max95_low_cor
varid=ncdf_varid(fid1,'min95_high_cor')
ncdf_varget,fid1,varid,min95_high_cor
varid=ncdf_varid(fid1,'max95_high_cor')
ncdf_varget,fid1,varid,max95_high_cor
varid=ncdf_varid(fid1,'min95_dif_cor')
ncdf_varget,fid1,varid,min95_dif_cor
varid=ncdf_varid(fid1,'max95_dif_cor')
ncdf_varget,fid1,varid,max95_dif_cor
varid=ncdf_varid(fid1,'lat')
ncdf_varget,fid1,varid,lat
varid=ncdf_varid(fid1,'level')
ncdf_varget,fid1,varid,level
NCDF_CLOSE,fid1



low_cor2=fltarr(91,34)
for i =0,90 do begin
  for j=0,33 do begin
    if min95_low_cor[i,j] gt 0 or max95_low_cor[i,j] lt 0 then low_cor2[i,j]=1
  endfor
endfor

high_cor2=fltarr(91,34)
for i =0,90 do begin
  for j=0,33 do begin
    if min95_high_cor[i,j] gt 0 or max95_high_cor[i,j] lt 0 then high_cor2[i,j]=1
  endfor
endfor



dif_cor2=fltarr(91,34)
for i =0,90 do begin
  for j=0,33 do begin
    if min95_dif_cor[i,j] gt 0 or max95_dif_cor[i,j] lt 0 then dif_cor2[i,j]=1
  endfor
endfor


help,level
print,level[18],level[33]

PSOPEN,/portrait,  YOFFSET=2000, SPACING=1800, SPACE3=400,CHARSIZE=100,TCHARSIZE=100 ,/letter,FILE='G:\My Drive\Isla\u120e_predictability_1980_2020_between_low_high.eps'
yvals=[1000,  500, 300,200,  150, 100,70, 50, 30, 20,10,5]
xvals=[-30,-15,0,15,30]
POS, XOFFSET=2000,YOFFSET=6000,XSIZE=4000,YSIZE=9000
CS, SCALE=1,NCOLS=22
GSET, XMIN=-30, XMAX=30, YMIN=1000, YMAX=5,/YLOG
LEVS, MIN=-1, MAX=1, STEP=0.1,NDECS=1

CONC,  X=lat, Y=level, FIELD=high_cor1,/NOLINES,/NOCOLBAR,/NOaxes ,TITLE='high'
AXES, XSTEP=10, YVALS=yvals, YLABELS=yvals,XTITLE='Latitude',YTITLE='Pressure (hPa)'

xpts=[-10, -10, 10,10,-10]
ypts=[150, 200, 200, 150, 150]
GPLOT, X=xpts, Y=ypts,  THICK=100


ypts=[100, 100, 100.2, 100.2,100]
xpts=[-30, 30, 30, -30, -30]
GPLOT, X=xpts, Y=ypts,  THICK=300

for i=0,15 do begin
  for j=0,20 do begin
    if high_cor2[2*i+30,j] ne 1  then GPLOT, X=lat[2*i+30], Y=level[j], SYM=3, SIZE=20, COL=1;, /NOLINES, /DEVICE
  endfor
endfor

for i=0,15 do begin
  for j=21,33 do begin
    if high_cor2[2*i+30,j] ne 1 and (j mod 3) eq 0 then GPLOT, X=lat[2*i+30], Y=level[j], SYM=3, SIZE=20, COL=1;, /NOLINES, /DEVICE
  endfor
endfor


POS, XOFFSET=8000,YOFFSET=6000,XSIZE=4000,YSIZE=9000
CS, SCALE=1,NCOLS=22
GSET, XMIN=-30, XMAX=30, YMIN=1000, YMAX=5,/YLOG
LEVS, MIN=-1, MAX=1, STEP=0.1,NDECS=1

CONC,  X=lat, Y=level, FIELD=low_cor1,/NOLINES,/NOCOLBAR,/NOaxes ,TITLE='low'
AXES, XSTEP=10, YVALS=yvals, YLABELS=yvals,XTITLE='Latitude';,YTITLE='Pressure (hPa)'

xpts=[-10, -10, 10,10,-10]
ypts=[150, 200, 200, 150, 150]
GPLOT, X=xpts, Y=ypts,  THICK=100


ypts=[100, 100, 100.2, 100.2,100]
xpts=[-30, 30, 30, -30, -30]
GPLOT, X=xpts, Y=ypts,  THICK=300


for i=0,15 do begin
  for j=0,20 do begin
    if low_cor2[2*i+30,j] ne 1  then GPLOT, X=lat[2*i+30], Y=level[j], SYM=3, SIZE=20, COL=1;, /NOLINES, /DEVICE
  endfor
endfor

for i=0,15 do begin
  for j=21,33 do begin
    if low_cor2[2*i+30,j] ne 1 and (j mod 3) eq 0 then GPLOT, X=lat[2*i+30], Y=level[j], SYM=3, SIZE=20, COL=1;, /NOLINES, /DEVICE
  endfor
endfor



POS, XOFFSET=14000,YOFFSET=6000,XSIZE=4000,YSIZE=9000
CS, SCALE=1,NCOLS=22
GSET, XMIN=-30, XMAX=30, YMIN=1000, YMAX=5,/YLOG
LEVS, MIN=-1, MAX=1, STEP=0.1,NDECS=1

CONC,  X=lat, Y=level, FIELD=dif_cor1,/NOLINES,/NOCOLBAR,/NOaxes ,TITLE='high-low'
AXES, XSTEP=10, YVALS=yvals, YLABELS=yvals,XTITLE='Latitude';,YTITLE='Pressure (hPa)'

xpts=[-10, -10, 10,10,-10]
ypts=[150, 200, 200, 150, 150]
GPLOT, X=xpts, Y=ypts,  THICK=100


ypts=[100, 100, 100.2, 100.2,100]
xpts=[-30, 30, 30, -30, -30]
GPLOT, X=xpts, Y=ypts,  THICK=300


for i=0,15 do begin
  for j=0,20 do begin
    if dif_cor2[2*i+30,j] ne 1  then GPLOT, X=lat[2*i+30], Y=level[j], SYM=3, SIZE=20, COL=1;, /NOLINES, /DEVICE
  endfor
endfor

for i=0,15 do begin
  for j=21,33 do begin
    if dif_cor2[2*i+30,j] ne 1 and (j mod 3) eq 0 then GPLOT, X=lat[2*i+30], Y=level[j], SYM=3, SIZE=20, COL=1;, /NOLINES, /DEVICE
  endfor
endfor

COLBAR, COORDS=[20000, 6000, 20200, 15000], /NOLINES, LABELS=['-1.0','-0.5','0','0.5', '1'];,/TEXTPOS,TITLE='(m/s)'



PSCLOSE,/noview 

end