# -*- coding: utf-8 -*-
"""
Materials, simplified form for XAFSmass
---------------------------------------

Defines f1 and f2 scattering factors.
"""
__author__ = "Konstantin Klementiev, Roman Chernikov"
__date__ = "30 Nov 2023"
import os
import math
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

elemental = dict(
    H=dict(name='Hydrogen', elements='H', rho=8.375e-05),
    He=dict(name='Helium', elements='He', rho=1.663e-04),
    Li=dict(name='Lithium', elements='Li', rho=5.340e-01),
    Be=dict(name='Beryllium', elements='Be', rho=1.848e+00),
    B=dict(name='Boron', elements='B', rho=2.370e+00),
    C=dict(name='Carbon', elements='C', rho=1.700e+00),
    N=dict(name='Nitrogen', elements='N', rho=1.165e-03),
    O=dict(name='Oxygen', elements='O', rho=1.332e-03),
    F=dict(name='Fluorine', elements='F', rho=1.580e-03),
    Ne=dict(name='Neon', elements='Ne', rho=8.385e-04),
    Na=dict(name='Sodium', elements='Na', rho=9.710e-01),
    Mg=dict(name='Magnesium', elements='Mg', rho=1.740e+00),
    Al=dict(name='Aluminum', elements='Al', rho=2.699e+00),
    Si=dict(name='Silicon', elements='Si', rho=2.330e+00),
    P=dict(name='Phosphorus', elements='P', rho=2.200e+00),
    S=dict(name='Sulfur', elements='S', rho=2.000e+00),
    Cl=dict(name='Chlorine', elements='Cl', rho=2.995e-03),
    Ar=dict(name='Argon', elements='Ar', rho=1.662e-03),
    K=dict(name='Potassium', elements='K', rho=8.620e-01),
    Ca=dict(name='Calcium', elements='Ca', rho=1.550e+00),
    Sc=dict(name='Scandium', elements='Sc', rho=2.989e+00),
    Ti=dict(name='Titanium', elements='Ti', rho=4.540e+00),
    V=dict(name='Vanadium', elements='V', rho=6.110e+00),
    Cr=dict(name='Chromium', elements='Cr', rho=7.180e+00),
    Mn=dict(name='Manganese', elements='Mn', rho=7.440e+00),
    Fe=dict(name='Iron', elements='Fe', rho=7.874e+00),
    Co=dict(name='Cobalt', elements='Co', rho=8.900e+00),
    Ni=dict(name='Nickel', elements='Ni', rho=8.902e+00),
    Cu=dict(name='Copper', elements='Cu', rho=8.960e+00),
    Zn=dict(name='Zinc', elements='Zn', rho=7.133e+00),
    Ga=dict(name='Gallium', elements='Ga', rho=5.904e+00),
    Ge=dict(name='Germanium', elements='Ge', rho=5.323e+00),
    As=dict(name='Arsenic', elements='As', rho=5.730e+00),
    Se=dict(name='Selenium', elements='Se', rho=4.500e+00),
    Br=dict(name='Bromine', elements='Br', rho=7.072e-03),
    Kr=dict(name='Krypton', elements='Kr', rho=3.478e-03),
    Rb=dict(name='Rubidium', elements='Rb', rho=1.532e+00),
    Sr=dict(name='Strontium', elements='Sr', rho=2.540e+00),
    Y=dict(name='Yttrium', elements='Y', rho=4.469e+00),
    Zr=dict(name='Zirconium', elements='Zr', rho=6.506e+00),
    Nb=dict(name='Niobium', elements='Nb', rho=8.570e+00),
    Mo=dict(name='Molybdenum', elements='Mo', rho=1.022e+01),
    Tc=dict(name='Technetium', elements='Tc', rho=1.150e+01),
    Ru=dict(name='Ruthenium', elements='Ru', rho=1.241e+01),
    Rh=dict(name='Rhodium', elements='Rh', rho=1.241e+01),
    Pd=dict(name='Palladium', elements='Pd', rho=1.202e+01),
    Ag=dict(name='Silver', elements='Ag', rho=1.050e+01),
    Cd=dict(name='Cadmium', elements='Cd', rho=8.650e+00),
    In=dict(name='Indium', elements='In', rho=7.310e+00),
    Sn=dict(name='Tin', elements='Sn', rho=7.310e+00),
    Sb=dict(name='Antimony', elements='Sb', rho=6.691e+00),
    Te=dict(name='Tellurium', elements='Te', rho=6.240e+00),
    I=dict(name='Iodine', elements='I', rho=4.930e+00),
    Xe=dict(name='Xenon', elements='Xe', rho=5.485e-03),
    Cs=dict(name='Cesium', elements='Cs', rho=1.873e+00),
    Ba=dict(name='Barium', elements='Ba', rho=3.500e+00),
    La=dict(name='Lanthanum', elements='La', rho=6.154e+00),
    Ce=dict(name='Cerium', elements='Ce', rho=6.657e+00),
    Pr=dict(name='Praseodymium', elements='Pr', rho=6.710e+00),
    Nd=dict(name='Neodymium', elements='Nd', rho=6.900e+00),
    Pm=dict(name='Promethium', elements='Pm', rho=7.220e+00),
    Sm=dict(name='Samarium', elements='Sm', rho=7.460e+00),
    Eu=dict(name='Europium', elements='Eu', rho=5.243e+00),
    Gd=dict(name='Gadolinium', elements='Gd', rho=7.900e+00),
    Tb=dict(name='Terbium', elements='Tb', rho=8.229e+00),
    Dy=dict(name='Dysprosium', elements='Dy', rho=8.550e+00),
    Ho=dict(name='Holmium', elements='Ho', rho=8.795e+00),
    Er=dict(name='Erbium', elements='Er', rho=9.066e+00),
    Tm=dict(name='Thulium', elements='Tm', rho=9.321e+00),
    Yb=dict(name='Ytterbium', elements='Yb', rho=6.730e+00),
    Lu=dict(name='Lutetium', elements='Lu', rho=9.840e+00),
    Hf=dict(name='Hafnium', elements='Hf', rho=1.331e+01),
    Ta=dict(name='Tantalum', elements='Ta', rho=1.665e+01),
    W=dict(name='Tungsten', elements='W', rho=1.930e+01),
    Re=dict(name='Rhenium', elements='Re', rho=2.102e+01),
    Os=dict(name='Osmium', elements='Os', rho=2.257e+01),
    Ir=dict(name='Iridium', elements='Ir', rho=2.242e+01),
    Pt=dict(name='Platinum', elements='Pt', rho=2.145e+01),
    Au=dict(name='Gold', elements='Au', rho=1.932e+01),
    Hg=dict(name='Mercury', elements='Hg', rho=1.355e+01),
    Tl=dict(name='Thallium', elements='Tl', rho=1.172e+01),
    Pb=dict(name='Lead', elements='Pb', rho=1.135e+01),
    Bi=dict(name='Bismuth', elements='Bi', rho=9.747e+00),
    Po=dict(name='Polonium', elements='Po', rho=9.320e+00),
    At=dict(name='Astatine', elements='At', rho=8.91e+00),
    Rn=dict(name='Radon', elements='Rn', rho=9.066e-03),
    Fr=dict(name='Francium', elements='Fr', rho=2.48e+00),
    Ra=dict(name='Radium', elements='Ra', rho=5.000e+00),
    Ac=dict(name='Actinium', elements='Ac', rho=1.007e+01),
    Th=dict(name='Thorium', elements='Th', rho=1.172e+01),
    Pa=dict(name='Protactinium', elements='Pa', rho=1.537e+01),
    U=dict(name='Uranium', elements='U', rho=1.895e+01),
    )

compounds = dict(
    SilverBromide=dict(formula='AgBr', elements=['Ag', 'Br'],
                       quantities=[1, 1], rho=6.473),
    AluminumArsenide=dict(formula='AlAs', elements=['Al', 'As'],
                          quantities=[1, 1], rho=3.81),
    Sapphire=dict(formula='Al2O3', elements=['Al', 'O'],
                  quantities=[2.0, 3.0], rho=3.97),
    AluminumPhosphide=dict(formula='AlP', elements=['Al', 'P'],
                           quantities=[1, 1], rho=2.42),
    BoronOxide=dict(formula='B2O3', elements=['B', 'O'],
                    quantities=[2.0, 3.0], rho=3.11),
    BoronCarbide=dict(formula='B4C', elements=['B', 'C'],
                      quantities=[4.0, 1], rho=2.52),
    BerylliumOxide=dict(formula='BeO', elements=['Be', 'O'],
                        quantities=[1, 1], rho=3.01),
    BoronNitride=dict(formula='BN', elements=['B', 'N'],
                      quantities=[1, 1], rho=2.25),
    Polyimide=dict(formula='C22H10N2O5', elements=['C', 'H', 'N', 'O'],
                   quantities=[22.0, 10.0, 2.0, 5.0], rho=1.43),
    Polypropylene=dict(formula='C3H6', elements=['C', 'H'],
                       quantities=[3.0, 6.0], rho=0.9),
    PMMA=dict(formula='C5H8O2', elements=['C', 'H', 'O'],
              quantities=[5.0, 8.0, 2.0], rho=1.19),
    Polycarbonate=dict(formula='C16H14O3', elements=['C', 'H', 'O'],
                       quantities=[16.0, 14.0, 3.0], rho=1.2),
    Kimfol=dict(formula='C16H14O3', elements=['C', 'H', 'O'],
                quantities=[16.0, 14.0, 3.0], rho=1.2),
    Mylar=dict(formula='C10H8O4', elements=['C', 'H', 'O'],
               quantities=[10.0, 8.0, 4.0], rho=1.4),
    Teflon=dict(formula='C2F4', elements=['C', 'F'],
                quantities=[2.0, 4.0], rho=2.2),
    ParyleneC=dict(formula='C8H7Cl', elements=['C', 'H', 'Cl'],
                   quantities=[8.0, 7.0, 1], rho=1.29),
    ParyleneN=dict(formula='C8H8', elements=['C', 'H'],
                   quantities=[8.0, 8.0], rho=1.11),
    Fluorite=dict(formula='CaF2', elements=['Ca', 'F'],
                  quantities=[1, 2.0], rho=3.18),
    CadmiumTungstate=dict(formula='CdWO4', elements=['Cd', 'W', 'O'],
                          quantities=[1, 1, 4.0], rho=7.9),
    CadmiumSulfide=dict(formula='CdS', elements=['Cd', 'S'],
                        quantities=[1, 1], rho=4.826),
    CadmiumTelluride=dict(formula='CdTe', elements=['Cd', 'Te'],
                          quantities=[1, 1], rho=5.85),
    CobaltSilicide=dict(formula='CoSi2', elements=['Co', 'Si'],
                        quantities=[1, 2.0], rho=5.3),
    Cromium3Oxide=dict(formula='Cr2O3', elements=['Cr', 'O'],
                       quantities=[2.0, 3.0], rho=5.21),
    CesiumIodide=dict(formula='CsI', elements=['Cs', 'I'],
                      quantities=[1, 1], rho=4.51),
    CopperIodide=dict(formula='CuI', elements=['Cu', 'I'],
                      quantities=[1, 1], rho=5.63),
    IndiumNitride=dict(formula='InN', elements=['In', 'N'],
                       quantities=[1, 1], rho=6.88),
    Indium3Oxide=dict(formula='In2O3', elements=['In', 'O'],
                      quantities=[2.0, 3.0], rho=7.179),
    IndiumAntimonide=dict(formula='InSb', elements=['In', 'Sb'],
                          quantities=[1, 1], rho=5.775),
    IridiumOxide=dict(formula='IrO2', elements=['Ir', 'O'],
                      quantities=[1, 2.0], rho=11.66),
    GalliumArsenide=dict(formula='GaAs', elements=['Ga', 'As'],
                         quantities=[1, 1], rho=5.316),
    GalliumNitride=dict(formula='GaN', elements=['Ga', 'N'],
                        quantities=[1, 1], rho=6.1),
    GalliumPhosphide=dict(formula='GaP', elements=['Ga', 'P'],
                          quantities=[1, 1], rho=4.13),
    HafniumOxide=dict(formula='HfO2', elements=['Hf', 'O'],
                      quantities=[1, 2.0], rho=9.68),
    LithiumFluoride=dict(formula='LiF', elements=['Li', 'F'],
                         quantities=[1, 1], rho=2.635),
    LithiumHydride=dict(formula='LiH', elements=['Li', 'H'],
                        quantities=[1, 1], rho=0.783),
    LithiumHydroxide=dict(formula='LiOH', elements=['Li', 'O', 'H'],
                          quantities=[1, 1, 1], rho=1.43),
    MagnesiumFluoride=dict(formula='MgF2', elements=['Mg', 'F'],
                           quantities=[1, 2.0], rho=3.18),
    MagnesiumOxide=dict(formula='MgO', elements=['Mg', 'O'],
                        quantities=[1, 1], rho=3.58),
    MagnesiumSilicide=dict(formula='Mg2Si', elements=['Mg', 'Si'],
                           quantities=[2.0, 1], rho=1.94),
    Mica=dict(formula='KAl3Si3O12H2', elements=['K', 'Al', 'Si', 'O', 'H'],
              quantities=[1, 3.0, 3.0, 12.0, 2.0], rho=2.83),
    Manganese2Oxide=dict(formula='MnO', elements=['Mn', 'O'],
                         quantities=[1, 1], rho=5.44),
    Manganese4Oxide=dict(formula='MnO2', elements=['Mn', 'O'],
                         quantities=[1, 2.0], rho=5.03),
    Molybdenum4Oxide=dict(formula='MoO2', elements=['Mo', 'O'],
                          quantities=[1, 2.0], rho=6.47),
    Molybdenum6Oxide=dict(formula='MoO3', elements=['Mo', 'O'],
                          quantities=[1, 3.0], rho=4.69),
    MolybdenumSilicide=dict(formula='MoSi2', elements=['Mo', 'Si'],
                            quantities=[1, 2.0], rho=6.31),
    RockSalt=dict(formula='NaCl', elements=['Na', 'Cl'],
                  quantities=[1, 1], rho=2.165),
    NiobiumSilicide=dict(formula='NbSi2', elements=['Nb', 'Si'],
                         quantities=[1, 2.0], rho=5.37),
    NiobiumNitride=dict(formula='NbN', elements=['Nb', 'N'],
                        quantities=[1, 1], rho=8.47),
    NickelOxide=dict(formula='NiO', elements=['Ni', 'O'],
                     quantities=[1, 1], rho=6.67),
    NickelSilicide=dict(formula='Ni2Si', elements=['Ni', 'Si'],
                        quantities=[2.0, 1], rho=7.2),
    RutheniumSilicide=dict(formula='Ru2Si3', elements=['Ru', 'Si'],
                           quantities=[2.0, 3.0], rho=6.96),
    Ruthenium4Oxide=dict(formula='RuO2', elements=['Ru', 'O'],
                         quantities=[1, 2.0], rho=6.97),
    SiliconCarbide=dict(formula='SiC', elements=['Si', 'C'],
                        quantities=[1, 1], rho=3.217),
    SiliconNitride=dict(formula='Si3N4', elements=['Si', 'N'],
                        quantities=[3.0, 4.0], rho=3.44),
    Silica=dict(formula='SiO2', elements=['Si', 'O'],
                quantities=[1, 2.0], rho=2.2),
    Quartz=dict(formula='SiO2', elements=['Si', 'O'],
                quantities=[1, 2.0], rho=2.65),
    TantalumNitride=dict(formula='TaN', elements=['Ta', 'N'],
                         quantities=[1, 1], rho=16.3),
    TantalumOxide=dict(formula='Ta2O5', elements=['Ta', 'O'],
                       quantities=[2.0, 5.0], rho=8.2),
    TitaniumNitride=dict(formula='TiN', elements=['Ti', 'N'],
                         quantities=[1, 1], rho=5.22),
    TitaniumSilicide=dict(formula='TiSi2', elements=['Ti', 'Si'],
                          quantities=[1, 2.0], rho=4.02),
    TantalumSilicide=dict(formula='Ta2Si', elements=['Ta', 'Si'],
                          quantities=[2.0, 1], rho=14),
    Rutile=dict(formula='TiO2', elements=['Ti', 'O'],
                quantities=[1, 2.0], rho=4.26),
    ULEGlass=dict(formula='Si.925Ti.075O2', elements=['Si', 'Ti', 'O'],
                  quantities=[0.925, 0.075, 2.0], rho=2.205),
    Uranium4Oxide=dict(formula='UO2', elements=['U', 'O'],
                       quantities=[1, 2.0], rho=10.96),
    VanadiumNitride=dict(formula='VN', elements=['V', 'N'],
                         quantities=[1, 1], rho=6.13),
    Water=dict(formula='H2O', elements=['H', 'O'],
               quantities=[2.0, 1], rho=1),
    TungstenCarbide=dict(formula='WC', elements=['W', 'C'],
                         quantities=[1, 1], rho=15.63),
    Zerodur=dict(formula='Si.56Al.5P.16Li.04Ti.02Zr.02Zn.03O2.46',
                 elements=['Si', 'Al', 'P', 'Li', 'Ti', 'Zr', 'Zn', 'O'],
                 quantities=[0.56, 0.5, 0.16, 0.04, 0.02, 0.02, 0.03, 2.46],
                 rho=2.53),
    ZinсOxide=dict(formula='ZnO', elements=['Zn', 'O'],
                   quantities=[1, 1], rho=5.675),
    ZinсSulfide=dict(formula='ZnS', elements=['Zn', 'S'],
                     quantities=[1, 1], rho=4.079),
    ZirconiumNitride=dict(formula='ZrN', elements=['Zr', 'N'],
                          quantities=[1, 1], rho=7.09),
    Zirconia=dict(formula='ZrO2', elements=['Zr', 'O'],
                  quantities=[1, 2.0], rho=5.68),
    ZirconiumSilicide=dict(formula='ZrSi2', elements=['Zr', 'Si'],
                           quantities=[1, 2.0], rho=4.88),
    Air=dict(formula='N0.781O0.209Ar0.009', elements=['N', 'O', 'Ar'],
             quantities=[0.781, 0.209, 0.009], rho=1.20E-06),
    CVDDiamond=dict(formula='C', elements=['C'],
                    quantities=[1], rho=3.52),
    )


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
            ef1f2 = (np.array(res[self.name+'_E'], dtype=np.float64),
                     np.array(res[self.name+'_f1'], dtype=np.float64),
                     np.array(res[self.name+f2key], dtype=np.float64))
        return ef1f2

    def get_f1f2(self, E):
        """Calculates (interpolates) f1 and f2 for the given array *E*."""
        if (E < self.E[0]) or (E > self.E[-1]):
            return (r'E={0} is out of the data table range [{1}, {2}]! ' +
                    r'Use another table.').format(E, self.E[0], self.E[-1])
        f1 = np.interp(E, self.E, self.f1)
        f2 = np.interp(E, self.E, self.f2)
        return f1 + 1j*f2
