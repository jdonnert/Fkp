pro makesnap_cassano

    @set_cgs
    gadget = obj_new('GADGETCODEOBJECT')

    readcol, 'cassano05.dat', z, dpp, nth, B

    rho = nth/gadget.density(1, /el)

    bfld  = B * [1.,1.,1.] / sqrt(3) 

    head = gadget.makehead()

    np = 1000

    head.npart[0] = np
    head.parttotal = head.npart
    head.num_files = 1
    head.boxsize = 20000
    
    for i=0,n_elements(z)-1 do  begin
        head.redshift = z[i]
        head.time = 1.0/(1+z[i])

        fname = 'snap_'+strn(i,length=3, padchar='0')

        gadget.writehead, fname, head

        tmp = np-findgen(np)
        gadget.addblock, fname, ULONG(tmp), 'ID'

        tmp[*] = gadget.t2u(1e8)
        gadget.addblock, fname,float(tmp), 'U'
        
		tmp[*] = 10
        gadget.addblock, fname,float(tmp), 'HSML'

        tmp[*] = dpp[i]
        gadget.addblock, fname, float(tmp), 'DPP'
      
        tmp[*] = rho[i]
        gadget.addblock, fname, float(tmp), 'RHO'

        tmp = make_array(3, np, /float)
        tmp[0,*] = bfld[0]
        tmp[1,*] = bfld[1]
        tmp[2,*] = bfld[2]
        gadget.addblock, fname, float(tmp), 'BFLD'

    end

    return 
end
