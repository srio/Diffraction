a = rascii('test_aperture.dat')
;xplot,a
z_min = min(a[0,*],max=z_max)
z_min = -0.4
z_max = 0.4

Iray = reform(a[1,*])
wIray = wave_setscale(Iray,a[0,0],a[0,1]-a[0,0],/P)

      xplot,make_set(wIray.x,wIray.w),$
        xtitle='-1', ytitle='-1', $
        title='Intensity distribution at screen ', $
        coltitles=['x ','Intensity']


npoint=100000L
print,"limits: ",z_min,z_max

plane = Dcomplexarr(npoint)
plane[*] = Dcomplex(1.0D,0.0D)
wplane = wave_setscale(plane,z_min,z_max,/I)

help,/str,wplane

intscale = wave_findvalue(wIray,wplane.x)
wplane.w[*] *= sqrt(intscale[*])

      tmp1 = (abs(wplane.w))^2
      tmp2 = atan(imaginary(wplane.w)/real_part(wplane.w))*180d0/!dpi
      tmpplot = make_set(wplane.x, tmp1,tmp2)
      coltitles=[' coord','Intensity(screen)','Phase [deg](screen)']
xplot,tmpplot,coltitles=coltitles

  ;apply ideal lens
 focallength = 100.0
  knum = 2.0D*!dpi/waLe 
  ;wplane.w[*] *= exp(Dcomplex(0.0D,-knum*wplane.x[*]*wplane.x[*]/focallength/2.0D))

; 1) compute FFT of the field
fft_delta = 1.0D/double(wplane.xsize)/wplane.xdelta
If (npoint mod 2) eq 1 then begin
    fft_offset = -fft_delta*double(npoint-1)/2.0D
endif else begin
    fft_offset = -fft_delta*double(npoint)/2.0D
endelse
help,fft_delta,fft_offset
if float(!version.release) GT 8.1 then begin 
      wfou_fft = wave_setscale(fft(wplane.w,/double,/center),fft_offset,fft_delta,/P)
endif else begin  ; CENTER kw not accepted
      wfou_fft = wave_setscale(shift(fft(wplane.w,/double),npoint/2),fft_offset,fft_delta,/P)
endelse

;for divergence, propagate to focallength 
; (results are wplane and wInt_angle)
wale = 5000.0e-8 ; wevelength in user units
wfou_fft.w[*] *= exp(Dcomplex(0.0D,-!dpi*waLe*focallength*wfou_fft.x[*]*wfou_fft.x[*])) 
wplane.w = fft(wfou_fft.w,/double,/inverse)
;angle intensity profile as a function of tan(dx) 
;wInt_angle = wave_setscale(abs(wplane.w)*abs(wplane.w),wplane.xoffset/distance,wplane.xdelta/distance,/P)
wInt_angle = wave_setscale(abs(wplane.w)*abs(wplane.w),wplane.xoffset/focallength,wplane.xdelta/focallength,/P)

;help,/str,wInt_angle,wplane
xplot,make_set(wInt_angle.x,wplane.x,wInt_angle.w),$
        ;xrange=[-2,2]*1d-6, xtitle='-1', ytitle='-1', $
        xtitle='-1', ytitle='-1', $
        title='Wavefront propagation: Intensity at focallength ('+$
        strcompress(focallength,/remove_all)+' cm)', $
        coltitles=['Angle [rads]','position [cm]','Intensity']
xplot,make_set(wplane.x,abs(wplane.w)^2),$
        xtitle='-1', ytitle='-1', $
        title='Wavefront propagation: Intensity at IMAGE ('+$
        strcompress(focallength,/remove_all)+' cm)', $
        coltitles=['position [cm]','Intensity']
end

