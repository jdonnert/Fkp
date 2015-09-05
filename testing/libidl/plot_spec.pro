pro plot_spec, snap

    @set_cgs
    setcolor, 1

    ngrid = 128

    fname = 'spec_128_'+strn(snap,len=3,padc='0')

    readspec, fname, spec, p=p

    plot, p, 10d^spec,/ylog,/xlog, psym=10, yrange=[1d6,1d13] $
        ,xtitle='p/m!De!N /c', ytitle='N(p)', ystyle=1

    fname = './spec_8096_'+strn(snap,len=3,padc='0')

   ; readspec, fname, np, p=p, head=head

    ;oplot, p, 10d^np, psym=10, col=color(0)
    
    return
end 

pro plot_txt_spec

 	@set_cgs
	common globals, gadget, cosmo

    setcolor, 1
	
	ngrid = 128

	for snap=1, 99 do begin

		fname = 'cassano/spec_'+strn(snap,len=3,padc='0')+'_128.txt'

		readcol, fname, p, spec

		plot, p/me/c, 10d^spec,/ylog,/xlog, psym=10, yrange=[1d6,1d13] $
        	,xtitle='p/m!De!N /c', ytitle='N(p)', ystyle=1


		fname = 'spec_'+strn(snap,len=3,padc='0')+'_128.txt'

		readcol, fname, p, spec
	
		oplot, p/me/c, 10D^spec, psym=10, col=color(0)

		wait, 0.1
	end

	return
end
