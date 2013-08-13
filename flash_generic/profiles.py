def density_profile_1D_evolution(files,outdir):
    import matplotlib.pyplot as plt
    import yt.mods as ytm
    from numpy import linspace,max
    from toolbox import select_scale

    print 'Producing a time-evolution plot of the radial density profile.'

    pc   = 3.08568025e18
    AU   = 1.49598e13
    Rsun = 6.955e10
    scale = select_scale(6.17e18)
    
    #ts = ytm.TimeSeriesData.from_filenames(files)

    fig = plt.figure()
    ax = plt.subplot(1,1,1)
    ax.set_xlabel(r'Radius [{0}]'.format(r'$R_{\odot}$' if scale == 'Rsun' else scale))
    ax.set_ylabel(r'Density [g/cm$^3$]')
    ax.grid(True)
   
    numfiles = len(files)
    greys = linspace(0.8,0,numfiles)
    colors = [[greys[i],greys[i],greys[i]] for i in range(numfiles)]
    
    for i, file in enumerate(files):
        print 'Processing file {0} of {1}: {2}'.format(i+1,len(files),file)
        pf = ytm.load(file)
        a = ytm.PlotCollection(pf,center=[0.5,0.5,0.5]/pf["unitary"]).add_profile_sphere(0.5, "unitary",["Radius","Density"], weight="CellMassMsun")
        radii, densities = a.data['Radius']*pf[scale], a.data['Density']
        if i == numfiles - 1:
            ax.semilogy(radii,densities,color=colors[i],label="Radial average density profile")
            ax.legend()
        else:
            ax.semilogy(radii,densities,color=colors[i])
        if i == 1:
            ax.set_xlim((0,0.5*pf[scale]/pf["unitary"]))
        plt.savefig(outdir+'/'+'temp_radial_density_profile.png')

#    for sto,pf in ts.piter(storage=storage):
#        print 'Processing file {0} of {1}: {2}'.format(i+1,len(files),pf.basename)
#        a = ytm.PlotCollection(pf,center=[0.5,0.5,0.5]/pf["unitary"]).add_profile_sphere(0.5, "unitary",["Radius","Density"], weight="CellMassMsun")
#        sto.result = a.data
#        i+=1
#
#    for i in storage:
#        scale = select_scale(6.17e18)
#        radii = storage[i]['Radius']*pf[scale]
#        densities = storage[i]['Density']
#        if i == len(storage) - 1:
#            ax.semilogy(radii,densities,color=colors[i],label="Radial density profile")
#        else:
#            ax.semilogy(radii,densities,color=colors[i])

    formats = ['png','eps','pdf']
    for type in formats:
        plt.savefig(outdir+'/'+'radial_density_profile.'+type,format=type)

