import meep as mp
import math
import numpy as np
import h5py

cc=300
ff=300
frequency=ff/cc
f150=150*(1e-6)/cc
fmin=frequency-f150
fmax=frequency+f150

def fun(t):
    return math.cos(2*math.pi*frequency*t)

cell=mp.Vector3(12,12,12)

pml_layers = [mp.PML(1.0)]

sources=[mp.Source(mp.CustomSource(src_func=fun,start_time=0),component=mp.Ez,center=mp.Vector3(-4,0,0))]

sim=mp.Simulation(cell_size=cell,
                  boundary_layers=pml_layers,
                  sources=sources,
                  resolution=resolution)

freqz=sim.add_dft_fields([mp.Ez], fmin, fmax, 11,center=mp.Vector3(1,2,2),size=mp.Vector3(0,2,2))
freqx=sim.add_dft_fields([mp.Ex], fmin, fmax, 11,center=mp.Vector3(1,2,2),size=mp.Vector3(0,2,2))
freqy=sim.add_dft_fields([mp.Ey], fmin, fmax, 11,center=mp.Vector3(1,2,2),size=mp.Vector3(0,2,2))


sim.run(until=20)

sim.output_dft(freqz,'ddmz')
sim.output_dft(freqx,'ddmx')
sim.output_dft(freqy,'ddmy')
