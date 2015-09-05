pro makesnap_primary_injection

	@set_cgs

    gadget = obj_new('GADGETCODEOBJECT')
    head = gadget.makehead()

    npart = 1

    head.npartart[0] = npart
    head.parttotal = head.npart
    head.num_files = 1
    head.boxsize = 20000
    
    nsnap = 100

    rho0 = 2d-28/gadget.density(1, h=1)
	bfld0 = 6e-6
	mach = 4.4
	hsml = 50D
	shock_speed = cs * Mach
	shock_rho = rho0
	shock_pressure = T * 

    for i=0, nsnap do begin

        a = (i+0.1)/float(nsnap)
        z = 1/a - 1

        fname = 'snap_'+strn(i,length=3, padchar='0')

        head.time = a
        head.redshift = z
        gadget.writehead, fname, head

        tmp = findgen(npart) + 1
        gadget.addblock, fname, ulong(tmp), 'ID'

        tmp[*] = gadget.T2U(1e8)
        gadget.addblock, fname, float(tmp), 'U'

        tmp[*] = 1e-30
        gadget.addblock, fname, float(tmp), 'DPP'

        tmp[*] = rho0 
        gadget.addblock, fname, float(tmp), 'RHO'

		tmp3d = make_array(3, npart, /float, val=0)
		tmp3d[0,*] = bfld0

        gadget.addblock, fname, float(tmp3d), 'BFLD'

        tmp[*] = mach
        gadget.addblock, fname, float(tmp), 'MACH'
        
        tmp[*] = hsml
		gadget.addblock, fname, float(tmp), 'HSML'

        tmp[*] = shock_speed
		gadget.addblock, fname, float(tmp), 'SHSP' 

        tmp[*] = shock_rho
		gadget.addblock, fname, float(tmp), 'SHRH' 

        tmp[*] = shock_pressure
		gadget.addblock, fname, float(tmp), 'SHPR' 
    end


    return
end
