def density_profile_1D_evolution(files,outdir):
    import matplotlib.pyplot as plt
    import yt.mods as ytm
    from numpy import linspace,max
    from toolbox import select_scale

    print 'Producing a time-evolution plot of the radial density profile.'

    pc   = 3.08568025e18
    AU   = 1.49598e13
    Rsun = 6.955e10
    
    storage = {}

    ts = ytm.TimeSeriesData.from_filenames(files)

    for sto,pf in ts.piter(storage=storage):
        a = ytm.PlotCollection(pf,center=[0.5,0.5,0.5]/pf["unitary"]).add_profile_sphere(0.5, "unitary",["Radius","Density"], weight="CellMassMsun")
        sto.result = a.data
    
    numfiles = len(files)
    greys = linspace(0.8,0,numfiles)
    colors = [[greys[i],greys[i],greys[i]] for i in range(numfiles)]

    fig = plt.figure()
    ax = plt.subplot(1,1,1)
    for i in storage:
        scale = select_scale(6.17e18)
        radii = storage[i]['Radius']*pf[scale]
        densities = storage[i]['Density']
        if i == len(storage) - 1:
            ax.semilogy(radii,densities,color=colors[i],label="Radial density profile")
        else:
            ax.semilogy(radii,densities,color=colors[i])

    if scale == 'Rsun': scale = r'$R_{\odot}$'
    ax.set_xlim((0,0.5*pf[scale]/pf["unitary"]))
    ax.set_xlabel(r'Radius [{0}]'.format(scale))
    ax.set_ylabel(r'Density [g/cm$^3$]')
    ax.grid(True)
    ax.legend()
   
    formats = ['png','eps','pdf']
    for type in formats:
        plt.savefig(outdir+'/'+'radial_density_profile.'+type,format=type)

