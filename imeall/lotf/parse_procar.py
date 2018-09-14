import re
import numpy as np

from argparse import ArgumentParser
from quippy import set_fortran_indexing

#For the PROCAR file we split by k-points, then bands, then ions.
#For this iron system there are 1480 bands, 251 ions 27 kpoints
#we'll store a tuple with the correct bits and pieces.
#(energy, s, py, pz, px, dxy, dyz, dz2, dxz, dx2)

set_fortran_indexing(False)

def w0gauss(x,n=0):
    """
    ! Copyright (C) 2001 PWSCF group
    We only implement the simplest gaussian mixing
    n=0. If you wish to extend check q-e/Modules/dos.f90
    """
    sqrtpm1 = 1.0/np.pi
    arg = min (200.0, x**2) 
    return np.exp(-arg)*sqrtpm1

def gen_dos(nions=207, nks=2, n_bands=1260, weight=0.5, Emin = -15.00, Emax = 4.0):
    deltaE = 0.0329
    degauss = 0.1
    ndos = int((Emax - Emin)/deltaE + 0.500000)
    print 'Computing grid of ', ndos, ' points.'
    kpoint_regex = re.compile("k-point\s+[0-9]+\s+:\s+([+-]?([0-9]*[.])?[0-9]+\s){3}\s+weight\s+=\s+([0-9\.]+)")
    bands_regex = re.compile("band\s+([0-9]+|\*\*\*)\s#\s+energy\s+([+-]?[0-9\.]+)\s#\socc.")
    weight_regex = re.compile("weight = ([0-9\.]+)")
    
    with open('PROCAR','r') as f:
        procar_str = f.read()
    
    #build a matrix integrated over bands and kpoints of the projected density 
    #of states.
    energies = []
    N_E = np.zeros([ndos, nions, 10])
   
    weights = weight_regex.findall(procar_str)
    print "k-point weights:", weights
    kpoints = re.split(kpoint_regex, procar_str)
    for n in range(1, ndos):
        print 'step:', n
        E = Emin + (n - 1)*deltaE
        for k in range(1,nks):
            bands = re.split(bands_regex, kpoints[k*4])[1:]
            weight = float(weights[k-1])
            print 'weight', weight
            for n_b in range(n_bands):
                e_nk = float(bands[n_b*3+1])
                energies.append(e_nk)
                ions = bands[n_b*3+2].split('\n')[3:-3]
                #for last band in k-point block there is an extra new line:
                if len(ions)==(nions+1): ions=ions[:-1]
                for n_ion, line in enumerate(ions):
                    line = np.array(map(float, line.split()[1:]))
                    N_E[n, n_ion,:] += weight*w0gauss((E-e_nk)/degauss)*line
    
    with open('pdos.npy', 'w') as f:
        np.save(f, N_E)

#n_at in the index of the.
def print_dos(Emin = -15.00, Emax = 4.0, n_at=207):
    with open('pdos.npy','r') as f:
        N_E = np.load(f)
    deltaE = 0.0329
    ndos = int((Emax - Emin)/deltaE + 0.500000)
    with open('pdos_{}.dat'.format(n_at), 'w') as g:
        for n in range(ndos):
            E = Emin + (n - 1)*deltaE
            print >> g, E, ' '.join(map(str, N_E[n, n_at-1,:]))

if __name__=='__main__':
    parser = ArgumentParser()
    parser.add_argument("-g", "--gen_dos", action="store_true")
    parser.add_argument("-p", "--print_dos", action="store_true")
    parser.add_argument("-a", "--ats", type=int, nargs='+', help="atom position number in POSCAR file for projected dos")
    parser.add_argument("--Emin", default = -10.0, type=float)
    parser.add_argument("--Emax", default = 4.0, type=float)
    #Example nions=207, nks=2, n_bands=1260
    parser.add_argument("--nions", type=int, help="number of ions in system.", default=207)
    parser.add_argument("--nks", type=int, help="number of k points in system.", default=3)
    parser.add_argument("--n_bands", type=int, help="number of n_bands in system.", default=1260)

    args = parser.parse_args()

    if args.gen_dos:
        gen_dos(nions=args.nions, nks=args.nks, n_bands=args.n_bands) 

    if args.print_dos:
        for at in args.ats:
            print_dos(n_at=at, Emin=args.Emin, Emax=args.Emax)


