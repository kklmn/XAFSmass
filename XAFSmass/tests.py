# -*- coding: utf-8 -*-
__author__ = "Konstantin Klementiev, Roman Chernikov"
__date__ = "9 Apr 2026"

import sys
import os
sys.path.append('..')  # analysis:ignore
from XAFSmass import XAFSmassCalc as xc

toSuper = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")


def test_formula():
    tests = """
        H
        NaCl
        HO
        H2O
        HOH
        (H2O)2
        (H2O)2OH
        ((H2O)2OH)12
        C6H5OH
        CuSO4
        CuSO3.8
        Fe%5SiO2
        (Fe2)%5SiO2
        (Al2O3)%10MgO
        Ni%1((Al2O3)%10MgO)
        Fe%2((Al2O3)%10SiO2)
        Fe%0.02((Al2O3)%10SiO2)
        FexSiO2
        """.splitlines()
    for t in tests:
        if t.strip():
            results = xc.formula.parseString(t, parseAll=True)
            print(t, '->', results.asList())
            print(t, '->', results)
            print(xc.reconstruct(results))


def test_calculate_element_dict(table='Chantler'):
    tests = [
        ['Na', 1100],
        ['Pd', 24522],
        ]
    for comp, E in tests:
        results = xc.formula.parseString(comp)
        res = xc.calculate_element_dict(results.asList(), E, table)
        print(res)


def test_powder_foil():
    tests = [
        ['CuSO4', 9100],
        ['Fe%5SiO2', 7200],
        ]
    for comp, E in tests:
        results = xc.formula.parseString(comp)
        nu, m, th, eDict = xc.calculate_powder(
            results.asList(), E, 2.5, 1.33, rho=5)
        print(u'{0} at {1}eV: nu={2}mmol, mass={3}mg, th={4}µm'.format(
              comp, E, xc.round_to_n(nu), xc.round_to_n(m), xc.round_to_n(th)))
        eDict2 = dict([k, [v[3], v[-2], v[-1]]] for k, v in eDict.items())
        print(eDict2)

        th, eDict = xc.calculate_foil(results.asList(), E, 2.5, rho=5)
        print(u'{0} at {1}eV: th={2}µm'.format(comp, E, xc.round_to_n(th)))

        eDict2 = dict([k, [v[-2], v[-1]]] for k, v in eDict.items())
        print(eDict2)


def test_gas():
    tests = [
        ['Ar', 9200],
        ]
    for comp, E in tests:
        results = xc.formula.parseString(comp)
        P = xc.calculate_gas(results.asList(), E, 0.1, 25)
        print(u'{0} at {1}eV: p={2}mbar'.format(comp, E, xc.round_to_n(P)))


def test_x():
    tests = [
        ['CuxSiO2', 9200],
        ]
    for comp, E in tests:
        results = xc.formula.parseString(comp)
        n, wt, wt1 = xc.calculate_x(results.asList(), E, 1, 0.1, 0.01)
        print(n, wt, wt1)


def test_parse_compound():
    tests = """
        MN
        Mn
        NaCl
        H2O
        Fe
        Fe%5SiO2
        Fe%5Si%O2
        (Fe2)%5SiO2
        (Al2O3)%10MgO
        Ni%1((Al2O3)%10MgO)
        Fe%2((Al2O3)%10SiO2)
        Fe%0.02((Al2O3)%10SiO2)
        FexSiO2
        """.splitlines()
    for t in tests:
        if t.strip():
            res = xc.parse_compound(t)
            print(t, '->', res)


def test_flux():
    # mix = [('N₂', 1013)]
    mix = [('N', 1013)]
    res, att = xc.calculate_flux(mix, energy=11208, length=40)
    fluxList = '{0:.1e}'.format(res).split('e+')
    fluxList[1] = fluxList[1].translate(toSuper)
    fluxStr = 'flux (ph/s/µA) = {0}, attenuation = {1:.1f} %'.format(
        "·10".join(fluxList), att*100)
    print(fluxStr)


def test_edge_jumps(edge='K', table='Chantler'):
    edges = ("K", "L1", "L2", "L3", "M1", "M2", "M3", "M4", "M5")

    selfDir = os.path.dirname(__file__)
    efname = os.path.join(selfDir, 'data', 'Energies.txt')
    Es = dict()
    iedge = edges.index(edge)
    with open(efname, 'r') as f:
        f.readline()
        f.readline()
        for line in f.readlines():
            cs = line.strip().split()
            try:
                Es[cs[1]] = eval(cs[2+iedge])
            except IndexError:
                pass

    for el, E in Es.items():
        formula = xc.formula.parseString(el)
        elementsDict = xc.calculate_element_dict(formula.asList(), E, table)[0]
        dSigma2 = elementsDict[el][2]
        print(el, dSigma2)


if __name__ == '__main__':
    # test_formula()
    # test_calculate_element_dict()
    # test_powder_foil()
    # test_gas()
    # test_x()
    # test_parse_compound()
    # test_flux()
    test_edge_jumps('K')
