# -*- coding: utf-8 -*-
"""
Materials, simplified form for XAFSmass
---------------------------------------

Defines f1 and f2 scattering factors.
"""
__author__ = "Konstantin Klementiev, Roman Chernikov"
__date__ = "1 Mar 2016"
import os
import math
import struct
import numpy as np

ch = 12398.4186  # {5}   {c*h[eV*A]}
twoPi = math.pi * 2.
chbar = ch / twoPi  # {c*hbar[eV*A]}
r0 = 2.817940285e-5  # A
avogadro = 6.02214199e23  # atoms/mol

elementsList = (
    'none', 'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
    'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V',
    'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br',
    'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag',
    'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr',
    'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu',
    'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi',
    'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U')

elementsMass = {
    'H': 1.00797, 'He': 4.0026, 'Li': 6.939, 'Be': 9.0122, 'B': 10.811,
    'C': 12.01115, 'N': 14.0067, 'O': 15.9994, 'F': 18.9984, 'Ne': 20.179,
    'Na': 22.9898, 'Mg': 24.305, 'Al': 26.9815, 'Si': 28.086, 'P': 30.9738,
    'S': 32.064, 'Cl': 35.453, 'Ar': 39.948, 'K': 39.102, 'Ca': 40.08,
    'Sc': 44.956, 'Ti': 47.90, 'V': 50.942, 'Cr': 51.996, 'Mn': 54.9380,
    'Fe': 55.847, 'Co': 58.9332, 'Ni': 58.71, 'Cu': 63.546, 'Zn': 65.37,
    'Ga': 69.72, 'Ge': 72.59, 'As': 74.9216, 'Se': 78.96, 'Br': 79.904,
    'Kr': 83.80, 'Rb': 85.4678, 'Sr': 87.62, 'Y': 88.905, 'Zr': 91.22,
    'Nb': 92.906, 'Mo': 95.94, 'Tc': 99, 'Ru': 101.07, 'Rh': 102.905,
    'Pd': 106.42, 'Ag': 107.868, 'Cd': 112.40, 'In': 114.82, 'Sn': 118.69,
    'Sb': 121.75, 'Te': 127.60, 'I': 126.9044, 'Xe': 131.30, 'Cs': 132.905,
    'Ba': 137.34, 'La': 138.906, 'Ce': 140.12, 'Pr': 140.907, 'Nd': 144.24,
    'Pm': 147, 'Sm': 150.35, 'Eu': 151.95, 'Gd': 157.25, 'Tb': 158.924,
    'Dy': 162.50, 'Ho': 164.930, 'Er': 167.26, 'Tm': 168.934, 'Yb': 173.04,
    'Lu': 174.967, 'Hf': 178.49, 'Ta': 180.948, 'W': 183.85, 'Re': 186.207,
    'Os': 190.2, 'Ir': 192.22, 'Pt': 195.08, 'Au': 196.967, 'Hg': 200.59,
    'Tl': 204.383, 'Pb': 207.19, 'Bi': 208.980, 'Po': 209, 'At': 210,
    'Rn': 222, 'Fr': 223, 'Ra': 226.025, 'Ac': 227.028, 'Th': 232.038,
    'Pa': 231.036, 'U': 238.039}


def read_atomic_data(elem):
    """
    Reads atomic data from ``AtomicData.dat`` file adopted from XOP [XOP]_.
    It has the following data:
    0  AtomicRadius[A]  CovalentRadius[A]  AtomicMass  BoilingPoint[K]
    MeltingPoint[K]  Density[g/ccm]  AtomicVolume
    CoherentScatteringLength[1E-12cm]  IncoherentX-section[barn]
    Absorption@1.8A[barn]  DebyeTemperature[K]  ThermalConductivity[W/cmK]

    In :meth:`read_atomic_data` only the mass is inquired. The user may
    extend the method to get the other values by simply adding the
    corresponding array elements to the returned value."""
    if isinstance(elem, str):
        Z = elementsList.index(elem)
    elif isinstance(elem, int):
        Z = elem
    else:
        raise NameError('Wrong element')
    dataDir = os.path.dirname(__file__)
    with open(os.path.join(dataDir, 'data', 'AtomicData.dat')) as f:
        for li in f:
            fields = li.split()
            if int(fields[0]) == Z:
                atomicData = [float(x) for x in fields]
                break
    return atomicData[3]


class Element(object):
    """This class serves for accessing the scattering factors f0, f1 and f2 of
    a chemical element. It can also report other atomic data listed in
    ``AtomicData.dat`` file adopted from XOP [XOP]_.
    """
    def __init__(self, elem, table='Chantler'):
        """
        The element can be specified by its name (case sensitive) or its
        ordinal number. At the time of instantiation the tabulated scattering
        factors are read which are then interpolated at the requested **q**
        value and energy. *table* can be 'Henke' (10 eV < *E* < 30 keV)
        [Henke]_, 'Chantler' (11 eV < *E* < 405 keV) [Chantler]_ or 'BrCo'
        (30 eV < *E* < 509 keV) [BrCo]_.

        The tables of f2 factors consider only photoelectric cross-sections.
        The tabulation by Chantler can optionally have *total* absorption
        cross-sections. This option is enabled by *table*='Chantler total'.

        .. [Henke] http://henke.lbl.gov/optical_constants/asf.html
           B.L. Henke, E.M. Gullikson, and J.C. Davis, *X-ray interactions:
           photoabsorption, scattering, transmission, and reflection at
           E=50-30000 eV, Z=1-92*, Atomic Data and Nuclear Data Tables
           **54** (no.2) (1993) 181-342.

        .. [Chantler] http://physics.nist.gov/PhysRefData/FFast/Text/cover.html
           http://physics.nist.gov/PhysRefData/FFast/html/form.html
           C. T. Chantler, *Theoretical Form Factor, Attenuation, and
           Scattering Tabulation for Z = 1 - 92 from E = 1 - 10 eV to E = 0.4 -
           1.0 MeV*, J. Phys. Chem. Ref. Data **24** (1995) 71-643.

        .. [BrCo] http://www.bmsc.washington.edu/scatter/periodic-table.html
           ftp://ftpa.aps.anl.gov/pub/cross-section_codes/
           S. Brennan and P.L. Cowan, *A suite of programs for calculating
           x-ray absorption, reflection and diffraction performance for a
           variety of materials at arbitrary wavelengths*, Rev. Sci. Instrum.
           **63** (1992) 850-853.
        """
        if isinstance(elem, str):
            self.name = elem
            self.Z = elementsList.index(elem)
        elif isinstance(elem, int):
            self.name = elementsList[elem]
            self.Z = elem
        else:
            raise NameError('Wrong element')
        self.E, self.f1, self.f2 = self.read_f1f2_vs_E(table=table)
#        self.mass = read_atomic_data(self.Z)
        self.mass = elementsMass[self.name]

    def read_f1f2_vs_E(self, table):
        """Reads f1 and f2 scattering factors from the given *table* at the
        instantiation time."""
        dataDir = os.path.dirname(__file__)
        table_fn = table.split()[0]
        pname = os.path.join(dataDir, 'data', table_fn+'.npz')
        f2key = '_f2tot' if 'total' in table else '_f2'
        with open(pname, 'rb') as f:
            res = np.load(f)
            ef1f2 = (np.array(res[self.name+'_E']),
                     np.array(res[self.name+'_f1']),
                     np.array(res[self.name+f2key]))
        return ef1f2

    def get_f1f2(self, E):
        """Calculates (interpolates) f1 and f2 for the given array *E*."""
        if (E < self.E[0]) or (E > self.E[-1]):
            return (r'E={0} is out of the data table range [{1}, {2}]! ' +
                    r'Use another table.').format(E, self.E[0], self.E[-1])
        f1 = np.interp(E, self.E, self.f1)
        f2 = np.interp(E, self.E, self.f2)
        return f1 + 1j*f2
