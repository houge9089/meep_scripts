import meep as mp
import math
import argparse

def main(args):


    resolution = 30

    nSiO2 = 1.4
    SiO2 = mp.Medium(index=nSiO2)

    # conversion factor for eV to 1/um
    eV_um_scale = 1/1.23984193

    # W, from Rakic et al., Applied Optics, vol. 32, p. 5274 (1998)
    W_eps_inf = 1
    W_plasma_frq = 13.22*eV_um_scale
    W_f0 = 0.206
    W_frq0 = 1e-10
    W_gam0 = 0.064*eV_um_scale
    W_sig0 = W_f0*W_plasma_frq**2/W_frq0**2
    W_f1 = 0.054
    W_frq1 = 1.004*eV_um_scale # 1.235 um
    W_gam1 = 0.530*eV_um_scale
    W_sig1 = W_f1*W_plasma_frq**2/W_frq1**2
    W_f2 = 0.166
    W_frq2 = 1.917*eV_um_scale # 0.647
    W_gam2 = 1.281*eV_um_scale
    W_sig2 = W_f2*W_plasma_frq**2/W_frq2**2

    W_susc = [ mp.DrudeSusceptibility(frequency=W_frq0, gamma=W_gam0, sigma=W_sig0),
               mp.LorentzianSusceptibility(frequency=W_frq1, gamma=W_gam1, sigma=W_sig1),
               mp.LorentzianSusceptibility(frequency=W_frq2, gamma=W_gam2, sigma=W_sig2) ]

    W = mp.Medium(epsilon=W_eps_inf, E_susceptibilities=W_susc)

    # crystalline Si, from M.A. Green, Solar Energy Materials and Solar Cells, vol. 92, pp. 1305-1310 (2008)
    # fitted Lorentzian parameters, only for 600nm-1100nm
    Si_eps_inf = 9.14
    Si_eps_imag = -0.0334
    Si_eps_imag_frq = 1/1.55
    Si_frq = 2.2384
    Si_gam = 4.3645e-02
    Si_sig = 14.797/Si_frq**2

    Si = mp.Medium(epsilon=Si_eps_inf,
                   D_conductivity=2*math.pi*Si_eps_imag_frq*Si_eps_imag/Si_eps_inf,
                   E_susceptibilities=[ mp.LorentzianSusceptibility(frequency=Si_frq, gamma=Si_gam, sigma=Si_sig) ])

    a = args.a                         # lattice periodicity
    cone_r = args.cone_r               # cone radius
    cone_h = args.cone_h               # cone height
    wire_w = args.wire_w               # metal-grid wire width
    wire_h = args.wire_h               # metal-grid wire height
    trench_w = args.trench_w           # trench width
    trench_h = args.trench_h           # trench height
    Np = args.Np                       # number of periods in supercell

    dair = 1.0                         # air gap thickness
    dmcl = 1.7                         # micro lens thickness
    dsub = 3000.0                         # substrate thickness
    dpml = 1.0                         # PML thickness

    sxy = Np*a
    sz = dpml+dair+dmcl+dsub+dpml
    cell_size = mp.Vector3(sxy, sxy, sz)

    boundary_layers = [ mp.PML(dpml, direction=mp.Z, side=mp.High),
                        mp.Absorber(dpml, direction=mp.Z, side=mp.Low) ]

    geometry = []

    if args.substrate:
        geometry = [ mp.Sphere(material=SiO2, radius=dmcl, center=mp.Vector3(0,0,0.5*sz-dpml-dair-dmcl)),
                     mp.Block(material=Si, size=mp.Vector3(mp.inf,mp.inf,dsub+dpml), center=mp.Vector3(0,0,-0.5*sz+0.5*(dsub+dpml))),
                     mp.Block(material=W, size=mp.Vector3(mp.inf, wire_w, wire_h), center=mp.Vector3(0,-0.5*sxy+0.5*wire_w,-0.5*sz+dpml+dsub+0.5*wire_h)),
                     mp.Block(material=W, size=mp.Vector3(mp.inf, wire_w, wire_h), center=mp.Vector3(0,+0.5*sxy-0.5*wire_w,-0.5*sz+dpml+dsub+0.5*wire_h)),
                     mp.Block(material=W, size=mp.Vector3(wire_w, mp.inf, wire_h), center=mp.Vector3(-0.5*sxy+0.5*wire_w,0,-0.5*sz+dpml+dsub+0.5*wire_h)),
                     mp.Block(material=W, size=mp.Vector3(wire_w, mp.inf, wire_h), center=mp.Vector3(+0.5*sxy-0.5*wire_w,0,-0.5*sz+dpml+dsub+0.5*wire_h)) ]

    if args.substrate and args.texture:
        for nx in range(Np):
            for ny in range(Np):
                cx = -0.5*sxy+(nx+0.5)*a
                cy = -0.5*sxy+(ny+0.5)*a
                geometry.append(mp.Cone(material=SiO2, radius=0, radius2=cone_r, height=cone_h, center=mp.Vector3(cx,cy,0.5*sz-dpml-dair-dmcl-0.5*cone_h)))

    if args.substrate:
        geometry.append(mp.Block(material=SiO2, size=mp.Vector3(mp.inf, trench_w, trench_h), center=mp.Vector3(0, -0.5*sxy+0.5*trench_w, 0.5*sz-dpml-dair-dmcl-0.5*trench_h)))
        geometry.append(mp.Block(material=SiO2, size=mp.Vector3(mp.inf, trench_w, trench_h), center=mp.Vector3(0, +0.5*sxy-0.5*trench_w, 0.5*sz-dpml-dair-dmcl-0.5*trench_h)))
        geometry.append(mp.Block(material=SiO2, size=mp.Vector3(trench_w, mp.inf, trench_h), center=mp.Vector3(-0.5*sxy+0.5*trench_w, 0, 0.5*sz-dpml-dair-dmcl-0.5*trench_h)))
        geometry.append(mp.Block(material=SiO2, size=mp.Vector3(trench_w, mp.inf, trench_h), center=mp.Vector3(+0.5*sxy-0.5*trench_w, 0, 0.5*sz-dpml-dair-dmcl-0.5*trench_h)))

    k_point = mp.Vector3(0,0,0)

    lambda_min = 0.7        # minimum source wavelength
    lambda_max = 1.0        # maximum source wavelength
    fmin = 1/lambda_max
    fmax = 1/lambda_min
    fcen = 0.5*(fmin+fmax)
    df = fmax-fmin

    sources = [ mp.Source(mp.GaussianSource(fcen, fwidth=df), component=mp.Ex, center=mp.Vector3(0, 0, 0.5*sz-dpml-0.5*dair), size=mp.Vector3(sxy, sxy, 0)) ]

    sim = mp.Simulation(resolution=resolution,
                        cell_size=cell_size,
                        boundary_layers=boundary_layers,
                        geometry=geometry,
                        dimensions=3,
                        k_point=k_point,
                        sources=sources)

    nfreq = 50
    refl = sim.add_flux(fcen, df, nfreq, mp.FluxRegion(center=mp.Vector3(0,0,0.5*sz-dpml),size=mp.Vector3(sxy,sxy,0)))
    trans_grid = sim.add_flux(fcen, df, nfreq, mp.FluxRegion(center=mp.Vector3(0,0,-0.5*sz+dpml+dsub+wire_h),size=mp.Vector3(sxy,sxy,0)))
    trans_sub_top = sim.add_flux(fcen, df, nfreq, mp.FluxRegion(center=mp.Vector3(0,0,-0.5*sz+dpml+dsub),size=mp.Vector3(sxy,sxy,0)))
    trans_sub_bot = sim.add_flux(fcen, df, nfreq, mp.FluxRegion(center=mp.Vector3(0,0,-0.5*sz+dpml),size=mp.Vector3(sxy,sxy,0)))

    sim.run(mp.at_beginning(mp.output_epsilon),until=0)

    if args.substrate:
        sim.load_minus_flux('refl-flux', refl)

    sim.run(until_after_sources=mp.stop_when_fields_decayed(50, mp.Ex, mp.Vector3(0,0,-0.5*sz+dpml+0.5*dsub), 1e-9))

    if not args.substrate:
        sim.save_flux('refl-flux', refl)

    sim.display_fluxes(refl, trans_grid, trans_sub_top, trans_sub_bot)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-substrate', action='store_true', default=False, help="add the substrate? (default: False)")
    parser.add_argument('-texture', action='store_true', default=False, help="add the texture? (default: False)")
    parser.add_argument('-a', type=float, default=0.4, help='lattice periodicity (default: 0.4 um)')
    parser.add_argument('-cone_r', type=float, default=0.2, help='cone radius (default: 0.2 um)')
    parser.add_argument('-cone_h', type=float, default=0.28247, help='cone height (default: 0.28247 um)')
    parser.add_argument('-wire_w', type=float, default=0.1, help='metal-grid wire width (default: 0.1 um)')
    parser.add_argument('-wire_h', type=float, default=0.2, help='metal-grid wire height (default: 0.2 um)')
    parser.add_argument('-trench_w', type=float, default=0.1, help='trench width (default: 0.1 um)')
    parser.add_argument('-trench_h', type=float, default=2.0, help='trench height (default: 2.0 um)')
    parser.add_argument('-Np', type=int, default=3, help='number of periods in supercell (default: 3)')
    args = parser.parse_args()
    main(args)
