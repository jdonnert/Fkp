pro plot_expansion_spectra

    @set_cgs

    alpha = -1
    gamma = 2
    N0 = 1 
    L = 1d2
    H = 1d3
    t0 = 0

    dt = 0.1

    path = 'expansion/spec_'

    plot, [1.,1.], /nodata, xrange=[1,1e3], yrange=[1d-8, 1d-3], /xlog, $
        /ylog, xtitle='p', ytitle='n(p)', xstyle=1, ystyle=1

    for i=1, 39, 10 do begin
        
        fname = path + strn(i, len=3, padc='0')+'_128.txt'

        readcol, fname, p, np

        np = 10D^np

        oplot, p/me/c, np, col=color( floor(i/10.))

        pp = p/me/c ; analytic solution

        t = i * dt

        np_ana = N0 * pp^(-gamma)*exp(alpha*(t-t0)) $
            * (1 - pp/H*exp(-alpha*(t-t0))) * (1 - L/pp*exp(alpha*(t-t0)))

        oplot, pp, np_ana, linestyle=2
    end

    t = 0
    np_ana = N0 * pp^(-gamma)*exp(alpha*(t-t0)) $
            * (1 - pp/H*exp(-alpha*(t-t0))) * (1 - L/pp*exp(alpha*(t-t0)))

        oplot, pp, np_ana
    
    return
end
