# -*- coding: utf-8 -*-
__author__ = "Konstantin Klementiev, Roman Chernikov"
__date__ = "20 Feb 2023"

from math import log10, floor
import itertools
from collections import defaultdict, OrderedDict
import numpy as np
from pyparsing import (Suppress, Word, nums, Forward, Group,
                       Optional, OneOrMore, oneOf, ParseResults)
from . import materials_simple as rm
from .materials_simple import read_atomic_data

crossSection = 4.208e7  # 2*r_0*c*h*N_A [eV*cm^2/mol]
R = 8.314510  # J/K/mol={m^3Pa/K/mol}
T = 295  # K

LPAR, RPAR, PER = map(Suppress, "()#")
nreal = Word(nums + '.' + '%' + 'x')

element = oneOf(' '.join(rm.elementsList)[5:])

# forward declare 'formula' so it can be used in definition of 'term'
formula = Forward()
term = Group((element | Group(LPAR + formula + RPAR)("subgroup")) +
             Optional(nreal, default=1)("mult"))
formula << OneOrMore(term)

subGroupCount = itertools.count(1)


# add parse actions for parse-time processing:
def make_int(tokens):
    t = tokens[0]
    if t[0] == '%':
        return 1 + 1j*float(t[1:])
    elif t[0] == 'x':
        return -1
#    elif t[0] == 'y':
#        return -2
    ft = float(t)
    it = int(ft)
    return it if it == ft else ft
nreal.setParseAction(make_int)


def multiply_contents(tokens):
    t = tokens[0]
    if t.subgroup:
        mult = t.mult
        sg = next(subGroupCount)  # is needed to calculate wt%
        for term in t.subgroup:
            term[1] *= mult
            term.append(sg)
        return t.subgroup
term.setParseAction(multiply_contents)


def sum_by_element(tokens):
    nt = tokens
    dopes = []
    for t in nt:
        if isinstance(t[1], complex):
            sg = t[2] if len(t) > 2 else next(subGroupCount)
            dopes.append([t, sg])
    if dopes:
        matrix = [[t[0], t[1]] for t in nt if isinstance(t[1], (float, int))]
        matrixMass = sum([read_atomic_data(e[0])*e[1] for e in matrix])
        dopeGroups = []
        for e in dopes:
            if not (e[1] in dopeGroups):
                dopeGroups.append(e[1])
        dopeUnitMasses = []
        dopewtPercents = []
        for dg in dopeGroups:
            m = sum([read_atomic_data(e[0][0])*e[0][1].real for e in dopes
                     if e[1] == dg])
            dopeUnitMasses.append(m)
            for e in dopes:
                if e[1] == dg:
                    break
            wt = e[0][1].imag / e[0][1].real
            dopewtPercents.append(wt)
        dopewtPercentsSum = sum(dopewtPercents)
        for dg, m, wt in zip(dopeGroups, dopeUnitMasses, dopewtPercents):
            if wt > 0:  # without 'x' or 'y'
                # Ni%1Rh%1((SiO2)%10Al2O3)
                # Ni%1Rh%1(Al2O3)
                # (Si0.8Er0.2O2)%10Al2O3
                dopeMass = matrixMass * wt / (100. - dopewtPercentsSum)
                for e in dopes:
                    if e[1] == dg:
                        e[0][1] = e[0][1].real * dopeMass / m
    else:
        elementsList = [t[0] for t in tokens]
        duplicates = len(elementsList) > len(set(elementsList))
        if duplicates:
            dd = defaultdict(int)
            for t in tokens:
                dd[t[0]] += t[1]
            nt = ParseResults([ParseResults([k, v]) for k, v in dd.items()])
formula.setParseAction(sum_by_element)


def round_to_n(x, n=3):
    res = x
    try:
        res = round(x, -int(floor(log10(x))) + n-1) if \
            isinstance(x, (float, np.float32, np.float64)) else x
#        res = round(x, -int(floor(log10(x))) + n-1)
    except (ValueError, OverflowError):
        pass
    return res


def reconstruct(parsed):
    outStr = ''
    p = [[t[0], t[1]] for t in parsed]
    p2 = list(itertools.chain(*p))
    for s in p2:
        if isinstance(s, str):
            outStr += s
        else:
            ts = str(round_to_n(s))
            outStr += ts if ts != '1' else ''
    return outStr


def _simple_line(x1, x2, y1, y2):
    a = (y2 - y1) / (x2 - x1)
    b = -(y2*x1 - y1*x2) / (x2 - x1)
    return a, b


def find_edge_step(E, element):
    lenArr, dE, backStep, maxStep = 0, 10, 250, 100
    step = 0
    while lenArr < 3:
        mask = np.abs(E - element.E) < backStep
        f2 = element.f2[mask]
        ef2 = element.E[mask]
        backStep += dE  # eV
        lenArr = len(ef2)
        step += 1
        if step > maxStep:
            break
    df2 = np.diff(f2)
    dSigma2 = 0
    f2jump = 0
    sigma2x = 0
    try:
        iEdge = np.where(df2 > 0)[0][-1]
#        f2jump = df2[iEdge]  # simple
        a, b = _simple_line(ef2[iEdge-1], ef2[iEdge], f2[iEdge-1], f2[iEdge])
        f2jump = f2[iEdge+1] - a*ef2[iEdge+1] - b
        ef2jump = ef2[iEdge+1]
        dSigma2 = f2jump * crossSection / ef2jump
        sigma2x = f2[iEdge+1] * crossSection / ef2jump
    except IndexError:
        # print(e)
        pass
    return dSigma2, f2jump, sigma2x


def calculate_element_dict(formulaList, E, table):
    elementsDict = {}
    for t in formulaList:
        if t[0] not in elementsDict:
            el = rm.Element(t[0], table=table)
            f1f2 = el.get_f1f2(E)
            if isinstance(f1f2, str):
                return f1f2
            sigma2 = f1f2.imag * crossSection / E
            dSigma2, f2jump, sigma2x = find_edge_step(E, el)
            if sigma2x > 0:
                sigma2 = sigma2x
            elementsDict[t[0]] = [el, sigma2, dSigma2, 0, f1f2.imag, f2jump]

    sumSigma2 = 0.
    sumMass = 0.
    for t in formulaList:
        el, sigma2, dSigma2 = elementsDict[t[0]][0:3]
        elementsDict[t[0]][3] += t[1]
        if t[1] > 0:
            sumSigma2 += sigma2 * t[1]
            sumMass += el.mass * t[1]

    elementsDict = OrderedDict(sorted(elementsDict.items(),
                                      key=lambda t: t[1][0].Z))
    return elementsDict, sumSigma2, sumMass


def calculate_powder(formulaList, E, muTd, area=None, rho=None,
                     table='Chantler'):
    if isinstance(formulaList[0], dict):
        elementsDict, sumSigma2, sumMass = formulaList
    else:
        res = calculate_element_dict(formulaList, E, table)
        if isinstance(res, str):
            return res
        elementsDict, sumSigma2, sumMass = res

    if area:
        nu = muTd * area / sumSigma2
        mass = nu * sumMass * 1e3

    if rho:
        thickness = muTd / rho * sumMass / sumSigma2 * 1e4
    else:
        thickness = 0

    for elName, el in elementsDict.items():
        if area:
            el.append(nu * el[0].mass * el[3] * 1e3)  # mass mg
        el.append(muTd / sumSigma2 * el[2] * el[3])  # jump

    if area:
        return nu*1e3, mass, thickness, elementsDict
    else:
        return thickness, elementsDict


def calculate_foil(formulaList, E, muTd, rho, table='Chantler'):
    return calculate_powder(formulaList, E, muTd, rho=rho, table=table)


def calculate_gas(formulaList, E, attenuation, length, table='Chantler'):
    if isinstance(formulaList[0], dict):
        sumSigma2 = formulaList[1]
    else:
        res = calculate_element_dict(formulaList, E, table)
        if isinstance(res, str):
            return res
        sumSigma2 = res[1]

    nu = -np.log(1-attenuation) / sumSigma2
    P = nu * R * T / length * 1e4
    return P


def calculate_x(formulaList, E, muTd, Deltamud, deltamud=None,
                table='Chantler'):
    if isinstance(formulaList[0], dict):
        elementsDict, sumSigma2, sumMass = formulaList
    else:
        res = calculate_element_dict(formulaList, E, table)
        if isinstance(res, str):
            return res
        elementsDict, sumSigma2, sumMass = res

    nu = Deltamud / muTd
    xElement = None
    for k, e in elementsDict.items():
        if e[3] == -1:  # element with x
            e[3] = max(nu * sumSigma2 / (e[2] - nu*e[1]), 0)
            xElement = e
    wt = xElement[0].mass * xElement[3] / (sumMass + e[3]*e[0].mass) * 100

    wt1 = 0
    if deltamud:
        nu += deltamud * (1-nu)
        wt1 = abs((muTd-Deltamud)/(muTd)**2 * deltamud * wt)

    return xElement[3], wt, wt1


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
            results = formula.parseString(t, parseAll=True)
            print(t, '->', results.asList())
            print(t, '->', results)
            print(reconstruct(results))


def test_powder_foil():
    tests = [
        ['CuSO4', 9100],
        ['Fe%5SiO2', 7200],
        ]
    for comp, E in tests:
        results = formula.parseString(comp)
        nu, m, th, eDict = calculate_powder(results.asList(), E, 2.5, 1.33,
                                            rho=5)
        print(u'{0} at {1}eV: nu={2}mmol, mass={3}mg, th={4}µm'.format(
              comp, E, round_to_n(nu), round_to_n(m), round_to_n(th)))
        eDict2 = dict([k, [v[3], v[-2], v[-1]]] for k, v in eDict.items())
        print(eDict2)

        th, eDict = calculate_foil(results.asList(), E, 2.5, rho=5)
        print(u'{0} at {1}eV: th={2}µm'.format(comp, E, round_to_n(th)))

        eDict2 = dict([k, [v[-2], v[-1]]] for k, v in eDict.items())
        print(eDict2)


def test_gas():
    tests = [
        ['Ar', 9200],
        ]
    for comp, E in tests:
        results = formula.parseString(comp)
        P = calculate_gas(results.asList(), E, 0.1, 25)
        print(u'{0} at {1}eV: p={2}mbar'.format(comp, E, round_to_n(P)))


def test_x():
    tests = [
        ['CuxSiO2', 9200],
        ]
    for comp, E in tests:
        results = formula.parseString(comp)
        n, wt, wt1 = calculate_x(results.asList(), E, 1, 0.1, 0.01)
        print(n, wt, wt1)


if __name__ == '__main__':
    test_formula()
#    test_powder_foil()
#    test_gas()
#    test_x()
