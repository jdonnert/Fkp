; Contains compression testing routines, 
; this acts on Rossellas data
; cubic hermite splines

pro test, snap

    @set_cgs
    ngrid = 128

    fname = 'cassano/spec_'+strn(snap,len=3,padc='0')+'_'+strn(ngrid)+'.txt'

    readcol, fname, p, spec
    
    p /= me*c
    p = alog10(p)
    np = alog10(spec)

    np[1] = 0
    
    dp = make_array(ngrid,/double)
    for i=1, ngrid -2 do $
        dp[i] =  p[i] - p[i-1]

    dL = make_array(ngrid, /double)
    ddL = make_array(ngrid, /double)
    dR = make_array(ngrid, /double)
    ddR = make_array(ngrid, /double)
    thres = make_array(ngrid, /double)

    point = {   idx         : 0L, $
                dx_left     : 0D, $
                dx_right    : 0D, $
                d2x_left    : 0D, $
                d2x_right   : 0D  $
            }

    idx = make_array(19, /double)
    dydx = make_array(2, 19, /double)
    dyddx = make_array(2, 19, /double)

    ; find maximum
    maxnp = 0
    for i=1, ngrid-2 do $
        if maxnp lt np[i] then begin
            maxnp = np[i]
            maxi = i
        end
    print, maxi, maxnp, p[maxi]

    np /= maxnp
    prec = 0.2
    
    plot, p, np,psym=10 ,yrange=[-2,3], xrange=[0,6],xtitle='p/m!De!N /c', $
        ytitle='N(p)', ystyle=1

    ; find first & second derivative
    for i=1, ngrid-2 do begin
        dL[i+1] =  ( np[i+1] - np[i] ) / (0.5*(dp[i]+dp[i+1]))
        ddL[i+1] = (np[i+1] - 2*np[i] + np[i-1]) /  (dp[i]+dp[i+1])
        dR[i-1] = dL[i]
        ddR[i-1] = ddL[i]
    end

    bad = where(finite(dL) ne 1)
    dL[bad] = 0
    bad = where(finite(dR) ne 1)
    dR[bad] = 0

    bad = where(finite(ddL) ne 1)
    ddL[bad] = 0
    bad = where(finite(ddR) ne 1)
    ddR[bad] = 0

    j = 0
    ; find minima & maxima
	print, "Find MIN and MAX"
    for i=1, ngrid -2 do begin

        if (np[i-1] lt prec) and (np[i] ge prec) then begin
            idx[j] = i
            dydx[1,j] = mean(dR[i])
            dyddx[0,j] = 0
            print, 'start @', j,i-1, np[idx[j]], np[idx[j]+1], np[idx[j]]
            j++
        end

        if (np[i+1] lt prec) and (np[i] ge prec) then begin
            idx[j] = i
            dydx[0,j] = mean(dL[i])
            dyddx[0,j] = 0
            jmax = j
            print, 'stop @', j,i, np[i], np[i+1], np[i]
            break
        end

        if j gt 0 and dL[i] * dR[i] le 0 then begin

            print, 'max/min @', j, i, np[i]
            idx[j] = i
            dydx[0,j] = 0
            dydx[1,j] = 0
            dyddx[0,j] = mean(ddL[i])
            dyddx[1,j] = mean(ddR[i])
            j++
        end
        
    end

    ; draw curve
    Curve = make_array(2,ngrid,/double)
    Curve = make_array(2,ngrid,/double)
    for j=0, jmax do begin

        P0 = [p[idx[j]], np[idx[j]]]
        P1 = [p[idx[j+1]], np[idx[j+1]]]

        M0 = make_array(2,/double)
        M1 = make_array(2,/double)
        
        M0[0] = -1/6.*dyddx[0,j+1] - 1/3.*dyddx[1,j] - P0[0] + P1[0]
        M0[1] = M0[0] * dydx[1,j]

        M1[0] = 1/6*dyddx[1,j] + 1./3. * dyddx[0,j+1] -P0[0]+P1[0]
        M1[1] = M1[0] * dydx[0,j+1]
        
        A = 2*P0 - 2*P1 + M0 + M1
        B = -3.0*P0 + 3.0*P1 - 2.0*M0 - M1
        C = M0
        D = P0

        oplot, [1.,1]*(P0[0] + M0[0]), [1.,1]*(P0[1] + M0[1]) ,psym=4,$
			col=color(0)
        oplot, [1.,1]*(P1[0] + M1[0]), [1.,1]*(P1[1] + M1[1]) ,psym=4,$
			col=color(0)

        pmin = (where(p eq P0[0]))[0]
        pmax = (where(p eq P1[0]))[0]

        for i=pmin,pmax-1 do begin
            t = (p[i]-P0[0]) / (P1[0] - P0[0])
            Curve[0,i] = t^3*A[0] + t^2*B[0] + t*C[0] + D[0]
            Curve[1,i] = t^3*A[1] + t^2*B[1] + t*C[1] + D[1]
        end

        print, 'M0', j, idx[j], M0[0], M0[1], dyddx[1,j],  dydx[1,j]
        print, 'M1', j+1, idx[j+1], M1[0], M1[1], dyddx[0,j+1],  dydx[0,j+1]

    end

    error = (Curve[1,*]-np)/np * p
    error[idx[jmax]:*] = 0

    k=1
    
    ; add points from error measure
    derror = make_array(ngrid,/double)
	print, "ERROR points"
    for i=idx[0], ngrid-2 do begin
        derror[i] = ( error[i+1] - error[i] ) / (0.5*(dp[i]+dp[i+1]))
        
        if derror[i] * derror[i-1] lt 0 and abs(error[i]) gt 0.01 then begin

            print, 'err', i, jmax+k, derror[i]
            idx[jmax+k] = i
            dydx[0,jmax+k] = mean(dL[i])
            dydx[1,jmax+k] = mean(dR[i])
            dyddx[0,jmax+k] = mean(ddL[i])
            dyddx[1,jmax+k] = mean(ddR[i])
            k++
        end
        
    end

    jmax += k
    idx = idx[0:jmax-1]
    dydx = dydx[*,0:jmax-1]
    dyddx = dyddx[*,0:jmax-1]

    srt_i = sort(idx)
    idx = idx[srt_i]
    dydx[0,*] = dydx[0,srt_i]
    dydx[1,*] = dydx[1,srt_i]
    dyddx[0,*] = dyddx[0,srt_i]
    dyddx[1,*] = dyddx[1,srt_i]

    x = make_array(ngrid,/float)

    ; redraw curve
    for j=0, n_elements(srt_i)-2 do begin

        P0 = [p[idx[j]], np[idx[j]]]
        P1 = [p[idx[j+1]], np[idx[j+1]]]

        M0 = make_array(2,/double)
        M1 = make_array(2,/double)

        M0[0] = -1/6.*dyddx[0,j+1] - 1/3.*dyddx[1,j] - P0[0] + P1[0]
        M0[1] = M0[0] * dydx[1,j]
        M1[0] = 1/6*dyddx[1,j] + 1./3. * dyddx[0,j+1] -P0[0]+P1[0]
        M1[1] = M1[0] * dydx[0,j+1]
    
        print, j, idx[j], M0[0], M0[1], dyddx[0,j], dR[idx[j]]
        print, j+1, idx[j+1], M1[0], M1[1], ddR[idx[j+1]], dR[idx[j+1]]


        A = 2*P0 - 2*P1 + M0 + M1
        B = -3.0*P0 + 3.0*P1 - 2.0*M0 - M1
        C = M0
        D = P0

        oplot, [1.,1]*(P0[0] + M0[0]), [1.,1]*(P0[1] + M0[1]) ,psym=4, $
			col=color(0)
        oplot, [1.,1]*(P1[0] + M1[0]), [1.,1]*(P1[1] + M1[1]) ,psym=4, $
			col=color(0)

        for i=idx[j],idx[j+1] do begin
            t = (p[i]-P0[0]) / (P1[0] - P0[0])
            Curve[0,i] = t^3*A[0] + t^2*B[0] + t*C[0] + D[0]
            Curve[1,i] = t^3*A[1] + t^2*B[1] + t*C[1] + D[1]
            ;if i eq 22 then $
           ; print, i,t, Curve[0,i],Curve[1,i], A[1], B[1], C[1], D[1]
            ;print, i,t, Curve[1,i], P0[1], P1[1], M0[1], M1[1] 
            ;print, i, t, M0[1],  M0[0], dydx[1,j], 
        end
    end

    oplot, p[idx] ,np[idx] ,psym=2, symsize=3
    ;oplot, Curve[0,idx[0]:idx[jmax-1]], Curve[1,idx[0]:idx[jmax-1]], col=color(1)

 ;   oplot, p[idx[0]:idx[jmax-1]], Curve[1,idx[0]:idx[jmax-1]], col=color(1)
    oplot, x[idx[0]:idx[jmax-1]], Curve[1,idx[0]:idx[jmax-1]], col=color(2)
    
    readcol, 'spec_uncompressed_'+strn(snap, len=3, padc='0'), i, lp, c, sp
    oplot, lp, c, col=color(4), thick=1
    
    readcol, 'spec_compressed_'+strn(snap, len=3, padc='0'), i, lp, c
    oplot, lp, c, col=color(2), thick=2

	return
end


pro polyn, snap

    @set_cgs
    ngrid = 100

    fname = 'spec_'+strn(snap,len=3,padc='0')+'_'+strn(ngrid)+'.txt'

    readcol, fname, p, spec
    
    p /= me*c
    p = alog10(p)
    np = spec

    np[1] = 0

    dp = make_array(ngrid,/double)
    for i=1, ngrid -2 do $
        dp[i] =  p[i] - p[i-1]

    dL = make_array(ngrid, /double)
    ddL = make_array(ngrid, /double)
    dR = make_array(ngrid, /double)
    ddR = make_array(ngrid, /double)
    thres = make_array(ngrid, /double)

    idx = make_array(19, /double)
    dydx = make_array(2, 19, /double)
    dyddx = make_array(2, 19, /double)

    ; find maximum
    maxnp = 0
    for i=1, ngrid-2 do $
        if maxnp lt np[i] then begin
            maxnp = np[i]
            maxi = i
        end
    print, maxi, maxnp, p[maxi]

    prec = 4/maxnp
    np /= maxnp
    
    plot, p, np,psym=10 ,yrange=[6,16]/maxnp $
        ,xtitle='p/m!De!N /c', ytitle='N(p)'

    ; find first & second derivative
    for i=1, ngrid-2 do begin
        dL[i+1] =  ( np[i+1] - np[i] ) / (0.5*(dp[i]+dp[i+1]))
        ddL[i+1] = (np[i+1] - 2*np[i] + np[i-1]) /  (dp[i]+dp[i+1])
        dR[i-1] = dL[i]
        ddR[i-1] = ddL[i]
        thres[i] = 1 - np[i]
    end


    bad = where(finite(dL) ne 1)
    dL[bad] = 0
    bad = where(finite(dR) ne 1)
    dR[bad] = 0

    bad = where(finite(ddL) ne 1)
    ddL[bad] = 0
    bad = where(finite(ddR) ne 1)
    ddR[bad] = 0

    j = 0
    ; find minima & maxima
    for i=1, ngrid -2 do begin

        if (thres[i-1] gt prec) and (thres[i] lt prec) then begin
            idx[j] = i-1
            dydx[1,j] = mean(dR[i:i+1])
            dyddx[0,j] = 0
            print, 'start', j,i-1, thres[idx[j]], thres[idx[j]+1], np[idx[j]]
            j++
        end

        if (thres[i+1] gt prec) and (thres[i] lt prec) then begin
            idx[j] = i
            dydx[0,j] = mean(dL[i-1:i])
            dyddx[0,j] = 0
            jmax = j
            print, 'stop', j,i, thres[i], thres[i+1], np[i]
            break
        end

        if j gt 0 and dL[i] * dR[i] le 0 then begin

            if  idx[j-1] eq i-1 or $
                idx[j-1] eq i-2 or $
                idx[j-1] eq i-3 then continue
            print, 'max/min', j, np[i]
            idx[j] = i
            dydx[0,j] = 0
            dydx[1,j] = 0
            dyddx[0,j] = mean(ddL[i-1:i])
            dyddx[1,j] = mean(ddR[i:i+1])
            j++
        end
        
    end

    oplot, p[idx[0:jmax]], np[idx[0:jmax]], psym=2, symsize=4

;draw curve
    curve = make_array(ngrid, /double)
    for k=0, jmax-1 do begin
        
        fac = p[idx[k+1]] - p[idx[k]]

        f = np[idx[k]]
        e = dydx[1,k] *fac
        d = dyddx[1,k]/2*fac^2


        alpha =  np[idx[k+1]] - (d+e+f)
        beta = dydx[0,k+1]*fac - (2*d+e)
        gamma = dyddx[0,k+1] *fac^2 * 0.5 - d 

        c = gamma- 4*beta + 10*alpha
        b = 7*beta - 15*alpha -2*gamma
        a = 6*alpha - 3*beta + gamma

        for i=idx[k],idx[k+1] do begin
            t = (p[i]-p[idx[k]])/ fac

            curve[i] = a*t^5 + b*t^4 + c*t^3 + d*t^2 + e*t + f
        end
    end

    oplot, p, curve, col=color(0) 

    j=jmax+1
stop
    while j lt 7 do begin 

        delta = (curve-np)/np

        bad = where(finite(delta) ne 1)
        delta[bad] = 0

        delta[0:idx[0]] = 0
        delta[idx[j-1]:*] = 0

        i = where(abs(delta) eq max(abs(delta)))
    
        print, 'error', j,i, np[i]
        
        idx[j] = i
        dydx[0,j] = mean(dL[i-1:i+1])
        dydx[1,j] = mean(dL[i-1:i+1])
        dyddx[0,j] = mean(ddL[i-1:i+1])
        dyddx[1,j] = mean(ddL[i-1:i+1])

        srt = sort(idx[0:j])
        idx[0:j] = (idx[0:j])[srt]
        dydx[0,0:j] = (dydx[0,0:j])[srt]
        dydx[1,0:j] = (dydx[1,0:j])[srt]
        dyddx[0,0:j] = (dyddx[0,0:j])[srt]
        dyddx[1,0:j] = (dyddx[1,0:j])[srt]

        j++

        curve = make_array(ngrid, /double)
        for k=0, j-1 do begin
        
            fac = p[idx[k+1]] - p[idx[k]]

          ;  f = np[idx[k]]
          ;  e = dydx[1,k] *fac
          ;  d = dyddx[1,k]/2*fac^2

          ;  alpha =  np[idx[k+1]] - (d+e+f)
          ;  beta = dydx[0,k+1]*fac - (2*d+e)
          ;  gamma = dyddx[0,k+1] *fac^2 * 0.5 - d 

          ;  c = gamma- 4*beta + 10*alpha
          ;  b = 7*beta - 15*alpha -2*gamma
          ;  a = 6*alpha - 3*beta + gamma

            d =  np[idx[k]]
            c = dydx[1,k] *fac

            alpha =   np[idx[k+1]] - d -c
            beta =  dydx[0,k+1]*fac - c

            b = 3*alpha - beta
            a = beta - 2*alpha
stop
            for i=idx[k],idx[k+1] do begin
                t = (p[i]-p[idx[k]])/ fac
                curve[i] = a*t^3 + b*t^2 + c*t + d
          ;      curve[i] = a*t^5 + b*t^4 + c*t^3 + d*t^2 + e*t + f
            end

            oplot, p[idx[0:j]], np[idx[0:j]], psym=2, symsize=4
        good = where(curve ne 0)
        oplot, p[good], curve[good], col=color(j-1)
            stop
        end
         

    end

stop
    return
end

pro plot_compress_debug_output, snap

	@set_cgs

	fname = 'spec_uncompressed_'+strn(snap, len=3, padc='0')

	readcol, fname, idx, logp, spec, np

	p = 10d^logp

	plot, p, np, xrange=[1,1d6], yrange=[1d0,1d13], ystyle=1, /ylog, psym=10 $
		, thick=1, /xlog

	fname = 'spec_knots_'+strn(snap, len=3, padc='0')

	readcol,fname, num,  idx, P1, Mleftx, Mlefty, Mrightx, Mrighty, normal

	nKnots = n_elements(idx)

	P1 = 10D^(P1*normal)
	Mrightx = 10D^Mrightx
	Mleftx = 10D^Mleftx
	Mrighty = 10D^MrightY
	MleftY = 10D^MleftY

	oplot, p[idx], P1, psym=1,symsize=2, thick=2

	oplot, p[idx]*mleftx, P1*mlefty, psym=4, col=color(0) 
	oplot, p[idx]*mrightx, P1*mrighty, psym=4, col=color(1)

	fname = 'spec_compressed_'+strn(snap, len=3, padc='0')

	readcol, fname, idx, logp, reconstructed_curve

	reconstructed_curve = 10D^(reconstructed_curve*normal[0])

	oplot, 10^logp, reconstructed_curve, col=color(2), thick=2

	oplot, p, 10d^spec, col=color(3), thick=1
	
	return
end

