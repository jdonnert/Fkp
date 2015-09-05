function read_header, fp , info=info

    head = { snapnum    : ulong64(0),  	$
             nbins      : ulong64(0),  	$
             nall       : ulong64(0),  	$
             npart      : ulong64(0),  	$
             nfiles     : ulong64(0),  	$
             startID    : ulong64(0),	$
             pmin       : double(0), 	$
             pmax       : double(0), 	$	
             plow       : double(0), 	$
             phigh      : double(0),  	$
			 compressed : long64(0),	 	$
			 specsize	: long64(0), 	 	$
			 fill		: bytarr(88)	$
            }

    readu, fp, head

    return, head
end

function find_files, fbase

    maxfiles = 1e4

    if file_test(fbase) eq 1 then $
        return, 1

    for i=0, maxfiles-1 do begin
    
        fname = fbase+'.'+strn(i)

        if file_test(fname) eq 0 then $
            break
    end

    if i eq maxfiles-1 then $
        print, "Found "+strn(maxfiles)+" files ???"

    return, i 
end


pro readspec, fin, spec, head=head, p=p, info=info, old=old, idlist=idlist

    if not n_params() ge 2 then begin
        print, "readspec, fin, spec, head=head, p=p, info=info, "
        print, "                     old=old, idlist=idlist"
        return
    end

    if keyword_set(old) then begin
        print, 'Reading an old spectrum'
        readspec_old, fin, spec, head=head, p=p, info=info
        return
    end

    if keyword_set(spec) then $
        undefine, spec

    if keyword_set(idlist) then begin

        if keyword_set(info) then $
            print, 'Sorting IDlist'

        idlist = idlist[sort(idlist)]
    end
    
    ; handle multiple files
    nfiles = find_files(fin)

    if (nfiles eq 0) then begin
        print, 'File not found '+fin
        return
    end else if keyword_set(info) then begin
        print, 'Found '+strn(nfiles)+' file(s)'
    end

    fp = 7

    for file = 0,nfiles-1 do begin
    
        fname = fin
        if nfiles gt 1 then $
            fname += '.'+strn(file) 

		if keyword_set(info) then $
			print, "Reading ", fname

        openr, fp, fname

        head = read_header(fp, info=info)
        
		if keyword_set(info) && file eq 0 then begin
			print, "Header :"
        	print, "   Snapshot "+strn(head.snapnum)
         	print, "   Nbins "+strn( head.nbins)
	        print, "   nAll "+strn( head.nall)
    	    print, "   nPartFile "+strn( head.npart)
         	print, "   nFiles "+strn( head.nfiles)
         	print, "   First ID "+strn( head.startID)
         	print, "   Pmin "+strn( head.pmin )
         	print, "   Pmax "+strn( head.pmax)
         	print, "   Plow "+strn( head.plow)
         	print, "   Phigh "+strn(head.phigh)
         	print, "   compressed "+strn(head.compressed)
         	print, "   SpecSize "+strn(head.specsize)
         	print, "   unit of P: log10(cgs/me/c)"
    	end

    	pstep = alog10(Head.Pmax/Head.Pmin)/(head.nbins-1)
    	p = Head.Pmin * 10D^(lindgen(head.nbins) * pstep)

		if (head.compressed eq 0) then $
        	tmp = make_array(head.Npart * head.nbins,/float) $
        else $
			tmp = make_array(head.Npart * (head.SpecSize + 4 ), /byte)

        readu, fp, tmp
        
        close, fp
    
		if (head.compressed eq 0) then $
        	tmp = transpose(reform(tmp,  head.Nbins,head.npart,/overw)) $
        else $
            tmp = transpose(reform(tmp, head.SpecSize+4, head.Npart, /overw))

        if keyword_set(IDlist) then begin ; take only the subset

            IDrange = head.startID+[0,head.npart]

            good = where(IDlist ge IDrange[0] $
                     AND IDlist lt IDrange[1], nIDs)

            if nIDs gt 0 then begin

                idx  = IDlist[good] - head.startID

                tmp = tmp[idx,*]            
 
                if head.compressed eq 1 then $
			        tmp = uncompress(tmp, head, p=p, info=info)
            
                if  not keyword_set(spec)  then begin

                    spec = tmp 
                    
                end else $
                    spec = [spec, tmp]  
                
                heap_gc

            end 

        end else begin

		    if head.compressed eq 1 then $
			    tmp = uncompress(tmp, head, p=p, info=info)

            if file eq 0 then $
               spec = tmp $
            else $
                spec = [spec, tmp]  
            end

        tmp = 0
    end

    ; construct global header
    if keyword_set(IDlist) then begin
        head.npart = n_elements(IDlist) 
        head.startID = IDlist[0]
    end else begin
        head.npart = head.nall
        head.startID = 1
    end
    



;    plot, p, 10D^spec[0,*],/ylog,/xlog $
;        ,xrange=[head.pmin,head.pmax], yrange=minmax(10D^spec[0,*])
    
    return

end



; old version of the script for backward compatiblility
pro readspec_old, fname, ncre, head=head, p=p, info=info

    on_error, 2

    if not n_params() ge 2 then begin
        print, "readspec_old, fname, ncre, head=head, p=p, info=info"
        return
    end

    head = { snapnum    : long(0),  $
             nbins      : long(0),  $
             npart      : long(0),  $
             pmin       : float(0), $
             pmax       : float(0), $
             plow       : float(0), $
             phigh      : float(0)  $
            }

    ltmp = long64(0)
    ftmp = float(0)

    fp = 7

    openr, fp, fname

    readu, fp, ltmp
    head.snapnum = ltmp

    readu, fp, ltmp
    head.nbins = ltmp

    readu, fp, ltmp
    head.npart = ltmp
    
    readu, fp, ftmp
    head.pmin = ftmp

    readu, fp, ftmp
    head.pmax = ftmp

    readu, fp, ftmp 
    head.plow = ftmp

    readu, fp, ftmp
    head.phigh = ftmp
    
    ncre = make_array(head.npart * head.nbins,/float)

    readu, fp, ncre 

    ncre = transpose(reform(ncre,  head.nbins,head.npart,/overwrite))

    close, fp

    pstep = alog10(Head.Pmax/Head.Pmin)/(head.nbins-1)
    
    p = Head.Pmin * 10D^(lindgen(head.nbins) * pstep)

    if keyword_set(info) then begin
        print, "Reading ", fname
        print, "Snapshot", head.snapnum
        print, "# of spectral bins", head.nbins
        print, "# of particles", head.npart
        print, "Pmin", head.pmin, ""
        print, "Pmax", head.pmax
        print, "Plow", head.plow
        print, "Phigh", head.phigh
        print, "unit of P: log10(cgs/me/c)"
        plot, p, 10D^ncre[0,*],/ylog,/xlog $
            ,xrange=[head.pmin,head.pmax], yrange=minmax(10D^ncre[0,*])
    end

    return

end

; this is rather arcane, the c code is more clear ...
function uncompress, work, head, p=p, info=info

    if keyword_set(info) then $
		print, "Found "+strn(head.npart)+"  Compressed Spectra"

	npart = ( size(work) )[1] 

	specsize = head.specsize ; per spectrum

	log_p = alog10(p)

	spec = make_array( npart, head.nbins)

	for i = 0, npart-1 do begin

		norm = ( float(work[i, 0:3], 0, 1) )[0]
		
        debug = 0

        if i eq info  then begin
            debug = 1 
            print, "Knots of particle "+strn(i)
        end

		Knots = uncompress_knots_binary(work[i, 4:4+specsize-1], $
            nKnots=nKnots, debug=debug)
		
		reconstruct_control_points, Knots, p, nKnots
		
		spec[i,*] = norm * draw_curve(Knots, log_p, head.nbins)
	end

	; oplot, p, 10d^spec, psym=10, col=color(0)

	return, spec
end

function draw_curve, K, log_p, nbins

	KNOT_STOP = 4

	curve = make_array(nbins, /double, val=-1000)

	nKnots = (where(K.type eq KNOT_STOP))[0]

	for i=0, nKnots-1 do begin

		c0 = 2*K[i].P - 2*K[i+1].P + K[i].Mright + K[i+1].Mleft

        c1 = -3*K[i].P+3*K[i+1].P-2*K[i].Mright-K[i+1].Mleft

        c2 = K[i].Mright

        c3 = K[i].P

		;print, i, K[i].idx, K[i].P[0], K[i].P[1], c0[1],c1[1], c2[1],c3[1]

		; this is wrong
		t = (log_p[K[i].idx:K[i+1].idx]-K[i].P[0])/(K[i+1].P[0]-K[i].P[0])

		curve[K[i].idx:K[i+1].idx] = c0[1]*t^3+c1[1]*t^2+c2[1]*t+c3[1]

	end

	return, curve
end

function uncompress_knots_binary, bytearr, nKnots=nKnots, debug=debug

	KNOT_EMPTY = 0
	KNOT_START = 1
	KNOT_MINMAX = 2
	KNOT_FULL = 3
	KNOT_STOP = 4

	Knot = { 	type 	: fix(KNOT_EMPTY), 	$
				idx		: 0L,		$
				P		: [0.,0.],	$
				Mleft	: [0.,0.],	$
				Mright	: [0.,0.]}

	K = replicate(Knot, 10)
		
	sizeof_hp = 5
	sizeof_fp = 9
	
	b = 0L ; byte counter
	i = 0L ; knot counter

	K[i].type = KNOT_START
	K[i].idx = uint(bytearr[b++])
	K[i].P[1] = uncompressFloat16bit(bytearr[b++:b++])
	K[i].Mright[0] = uncompressFloat8bit(bytearr[b++])
	K[i++].Mright[1] = uncompressFloat8bit(bytearr[b++])

	if keyword_set(debug) then $
		print, i-1, " START", K[i-1].idx, K[i-1].P[1], $
							  K[i-1].Mright[0], K[i-1].Mright[1]
	
	if K[i-1].idx eq 0 then $ ; no spectrum here
		return, K 

	nMinMaxKnots = 0L

	while 1 do begin 
		
		if uint(bytearr[b]) ge 128 then begin ; marker bit set

			K[i].type = KNOT_FULL

			K[i].idx = uint(bytearr[b++]) AND 127 ; remove marker bit
			K[i].P[1] = uncompressFloat16bit(bytearr[b++:b++])

			K[i].Mleft[0] = uncompressFloat8bit(bytearr[b++]) 
			K[i].Mright[0] = uncompressFloat8bit(bytearr[b++]) 

			K[i].Mleft[1] = uncompressFloat16bit(bytearr[b++:b++])
			K[i].Mright[1] = uncompressFloat16bit(bytearr[b++:b++])

			if keyword_set(debug) then $
				print, i, " FULL", K[i].idx, K[i].P[1], $
					K[i].Mleft[0], K[i].Mleft[1], $
                    K[i].Mright[0], K[i].Mright[1]

		end else begin ; marker bit not set 
			
			K[i].type = KNOT_MINMAX
			
			K[i].idx = uint(bytearr[b++])
			K[i].P[1] = uncompressFloat16bit(bytearr[b++:b++])

			K[i].Mleft[0] = uncompressFloat8bit(bytearr[b++]) 
			K[i].Mleft[1] = 0;K[i].P[1] 

			K[i].Mright[0] = uncompressFloat8bit(bytearr[b++]) 
			K[i].Mright[1] = 0;K[i].P[1] 

			nMinMaxKnots++;

			if keyword_set(debug) then $
				print, i, " MINMAX", K[i].idx, K[i].P[1], $
					K[i].Mleft[0], K[i].Mright[0]
		end

		i++

		if  (nMinMaxKnots mod 2) and $
		  (b gt n_elements(bytearr) - sizeof_hp - sizeof_fp)  then $
			break;
	end
	
	K[i].type = KNOT_STOP
	K[i].idx = uint(bytearr[b++])
	K[i].P[1] = uncompressFloat16bit(bytearr[b++:b++])
	K[i].Mleft[0] = uncompressFloat8bit(bytearr[b++])
	K[i].Mleft[1] = uncompressFloat8bit(bytearr[b++])

	if (K[i].idx eq 0) then begin ; we overshot

		i--
		K[i].type = KNOT_STOP
		K[i].Mleft[1] = K[i].Mright[0]
	end
		

	if keyword_set(debug) then $
		print, i, " STOP", K[i].idx, K[i].P[1], K[i].Mleft[0], K[i].Mleft[1]
	
    nKnots = i+1

	return, K 
end


pro reconstruct_control_points, K, p, nKnots
	
	for i = 0, nKnots-1 do begin

		idx = K[i].idx

		K[i].P[0] = alog10(p[idx])

	end

	;i = 0 ; KNOT_START
	;K[i].Mright[0] =  K[i+1].P[0] - K[i].P[0]
	;K[i].Mright[1] = K[i].Mright[0] * K[i].dR

	;print, i, " KNOT_START ", K[i].idx, K[i].P[0], K[i].P[1], $
	;	K[i].Mleft[0], K[i].Mleft[1], K[i].Mright[0], K[i].Mright[1]

	;i = nKnots-1 ; KNOT_STOP
	;K[i].Mleft[0] = K[i].P[0] - K[i-1].P[0] 
	;K[i].Mleft[1] = K[i].Mleft[0] * K[i].dL

	;print, i, " KNOT_STOP", K[i].idx, K[i].P[0], K[i].P[1], $
;		K[i].Mleft[0], K[i].Mleft[1], K[i].Mright[0], K[i].Mright[1]

	return
end	

function uncompressFloat16bit, input

	number = long(input[1]) * 256L + long(input[0]) 

	return, float( number) * 0.000244140625 - 8.0
end

function uncompressFloat8bit, input

	return, uint( input ) * 0.015625 - 1.0
end
