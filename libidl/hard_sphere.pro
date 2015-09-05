pro plot_hard_sphere_steady_state, mode, snap

    path = './hard_sphere_'+strn(mode, len=1)

    ; analytic solution
    fname = path+'/solution.txt'
    readcol, fname, p, np

    plot, p, np, /xlog, /ylog, xrange=[1d-4, 1d4], $
        yrange=[1d-9, 1d1], ystyle=1, xtitle='p', ytitle='n(p)'

    ; numeric solution
    fname = path+'/spec_'+strn(snap, len=3, padc='0')+'_128.txt'
    readcol, fname, p_num, np_num

    oplot, p_num, 10.^np_num, psym=10, color=color(4)


    return
end

pro hard_sphere_evolution, mode

    setcolor

    path = './hard_sphere_'+strn(mode, len=1)

    fname = path+'/solution.txt'
    readcol, fname, p, np

    plot, p, np, /xlog, /ylog, xrange=[1d-4, 1d4],  $
        yrange=[1d-9, 1d1], ystyle=1, xtitle='p', ytitle='n(p)'

    delta = alog10(100./1)/10

    for i = 1, 10 do begin
        
        snap = round(1 * 10D^(delta * i) )
        
        print, i , snap

        fname = path+'/spec_'+strn(snap, len=3, padc='0')+'_128.txt'
        readcol, fname, p_num, np_num

        oplot, p_num, 10.^np_num, psym=10, color=color(i)

    end

    oplot, p, np, thick=2

    return
end

pro hard_sphere_evolution_5

    setcolor

    path = './hard_sphere_5'

    fname = path+'/solution_0.3.txt'
    readcol, fname, p, np

    plot, p, np, /xlog, /ylog, xrange=[1d-3, 1d3], xstyle=1,  $
        yrange=[1d-8, 1d2], ystyle=1, xtitle='p', ytitle='n(p)', $
        xtickformat='exponent'

    fname = path+'/solution_3.0.txt'
    readcol, fname, p, np

    oplot, p, np

    fname = path+'/solution_30.txt'
    readcol, fname, p, np

    oplot, p, np*1e11

    snaps=[3,30,300]

    for i = 0, 2 do begin
        
        snap = snaps[i]
        
        print, i , snap

        fname = path+'/spec_'+strn(snap, len=3, padc='0')+'_128.txt'
        readcol, fname, p_num, np_num

        oplot, p_num, 10D^np_num, psym=10, color=color(i)

        if i eq 2 then $
            oplot, p_num, 1d11*10D^np_num, psym=10, color=color(i)

    end

    xyouts, 40, 0.01, 't = 0.3', col=color(0), charsize=0.7
    xyouts, 40, 0.1, 't = 3', col=color(1), charsize=0.7
    xyouts, 40, 1, 't = 30', col=color(2), charsize=0.7

    return
end


pro make_pdfs
    outpath = './'

    ; 0 and 2 are broken, P+P96 are faulty

    tops, 'hard_sphere_evolution, 1', outpath+'hs_test_1', /pdf
    tops, 'hard_sphere_evolution, 3', outpath+'hs_test_3', /pdf
    tops, 'hard_sphere_evolution, 4', outpath+'hs_test_4', /pdf
    tops, 'hard_sphere_evolution_5', outpath+'hs_test_5', /pdf

    return
end
