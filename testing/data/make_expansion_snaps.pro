pro makesnap_expansion
    @set_cgs

    gadget = obj_new('GADGETCODEOBJECT')
    head = gadget.makehead()

    np = 1

    head.npart[0] = np
    head.parttotal = head.npart
    head.num_files = 1
    head.boxsize = 20000
    
    nsnap = 100

    rho0 = 5.72977e-07

    for i=0, nsnap do begin

        a = i/float(nsnap)
        z = 1/a - 1

        fname = 'snap_'+strn(i,length=3, padchar='0')

        head.time = a
        head.redshift = z
        gadget.writehead, fname, head

        tmp = np-findgen(np)
        gadget.addblock, fname, ulong(tmp), 'ID'

        tmp[*] = gadget.T2U(1e8)
        gadget.addblock, fname,float(tmp), 'U'

        tmp[*] = 1e-18
        gadget.addblock, fname, float(tmp), 'DPP'

        tmp[*] = rho0 * exp(30 * a) ; 1e-4 nth
        gadget.addblock, fname, float(tmp), 'RHO'

    end


    return
end
