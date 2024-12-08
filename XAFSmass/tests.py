# -*- coding: utf-8 -*-
__author__ = "Konstantin Klementiev, Roman Chernikov"
__date__ = "25 Nov 2024"

import sys
sys.path.append('..')  # analysis:ignore
from XAFSmass import XAFSmassCalc as xc


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


if __name__ == '__main__':
    # test_formula()
    # test_powder_foil()
    # test_gas()
    # test_x()
    test_parse_compound()
