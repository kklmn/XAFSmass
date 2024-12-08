# -*- coding: utf-8 -*-
__author__ = "Konstantin Klementiev, Roman Chernikov"
__date__ = "25 Nov 2024"

from math import log10, floor
import itertools
from collections import defaultdict, OrderedDict
import numpy as np
from pyparsing import (
    Suppress, Word, nums, Forward, Group, Optional, OneOrMore, oneOf,
    ParseResults, ParseBaseException)

from . import materials_simple as rm

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
        matrixMass = sum([rm.read_atomic_data(e[0])*e[1] for e in matrix])
        dopeGroups = []
        for e in dopes:
            if not (e[1] in dopeGroups):
                dopeGroups.append(e[1])
        dopeUnitMasses = []
        dopewtPercents = []
        for dg in dopeGroups:
            m = sum([rm.read_atomic_data(e[0][0])*e[0][1].real for e in dopes
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


def find_victoreen_f2(istart, iend, element):
    ef2 = element.E[istart:iend+1]
    f2 = element.f2[istart:iend+1]
    b = np.log(ef2)
    a = f2 * crossSection / ef2
    sum6 = np.exp(-6*b).sum()
    sum7 = np.exp(-7*b).sum()
    sum8 = np.exp(-8*b).sum()
    ysum1 = (a*np.exp(-3*b)).sum()
    ysum2 = (a*np.exp(-4*b)).sum()

    det = sum6*sum8 - sum7*sum7
    c = (sum8*ysum1 - sum7*ysum2) / det
    d = (-sum7*ysum1 + sum6*ysum2) / det
    return c, d


def find_edge_step(E, element):
    dE, backE, forwardE, maxStep = 50, 250, 50, 100
    dSigma2 = 0
    f2jump = 0
    sigma2x = 0
    iEdge = 0
    istep = 0
    while iEdge < 2:
        mask = (E - backE < element.E) & (element.E < E + forwardE)
        f2 = element.f2[mask]
        ef2 = element.E[mask]
        df2 = np.diff(f2)
        try:
            iEdge = np.where(df2 > 0)[0][-1]
            a, b = _simple_line(
                ef2[iEdge-1], ef2[iEdge], f2[iEdge-1], f2[iEdge])
            f2bknd = a*ef2[iEdge+1] + b
            f2jump = f2[iEdge+1] - f2bknd
            toSigma2 = crossSection / ef2[iEdge+1]
            dSigma2 = f2jump * toSigma2
            sigma2x = f2[iEdge+1] * toSigma2
        except IndexError:
            break
        backE += dE  # eV
        istep += 1
        if istep > maxStep:
            break

    if f2jump > 0:
        vicc, vicd = find_victoreen_f2(iEdge-2, iEdge, element)
    else:
        vicc, vicd = None, None
    return dSigma2, f2jump, sigma2x, vicc, vicd


def calculate_element_dict(formulaList, E, table):
    """
    For *formulaList* (a list of [element, mole_amount]), energy *E* (float or
    numpy array) and *table* (str) of scattering factors returns a tuple
    (elementsDict, sumSigma2, sumMass) that serves as input for
    :func:`calculate_powder`, :func:`calculate_foil` and :func:`calculate_gas`.

    *formulaList* is obtained by a previous call of `formula.parseString` or
    :func:`parse_compound`.

    *elementsDict* is a dict of element symbols with dict values
    [el, sigma2, dSigma2, n, f1f2.imag, f2jump, vicc, vicd], where *el* is
    rm.Element object, *sigma2* is absorption σ² interpolated at energies *E*,
    *dSigma2* is Δσ² (edje jump), *n* the total molar amount of this element in
    the compound, *f1f2.imag* is the factor f'' at energies *E*, *f2jump* its
    edge jump, and *vicc* and *vicd* are the Victoreen polynomial coefficients
    if an absorption edge was detected or Nones otherwise.
    """
    elementsDict = {}
    for t in formulaList:
        if t[0] not in elementsDict:
            el = rm.Element(t[0], table=table)
            f1f2 = el.get_f1f2(E)
            if isinstance(f1f2, str):
                return f1f2
            sigma2 = f1f2.imag * crossSection / E
            dSigma2, f2jump, sigma2x, vicc, vicd = find_edge_step(E, el)
            if sigma2x > 0:
                sigma2 = sigma2x
            elementsDict[t[0]] = [el, sigma2, dSigma2, 0, f1f2.imag, f2jump,
                                  vicc, vicd]

    sumSigma2 = 0.
    sumMass = 0.
    for t in formulaList:
        el, sigma2, dSigma2 = elementsDict[t[0]][0:3]
        elementsDict[t[0]][3] += t[1]  # total of element t[0] in the compound
        if t[1] > 0:
            sumSigma2 += sigma2 * t[1]
            sumMass += el.mass * t[1]

    elementsDict = OrderedDict(
        sorted(elementsDict.items(), key=lambda t: t[1][0].Z))
    return elementsDict, sumSigma2, sumMass


def calculate_absorption_background(elementsDict, E):
    """
    For the dict *elementsDict* returned by :func:`calculate_element_dict()`
    calculated absorption background (as σ²) at energies *E*. The background
    consists of σ² of all elements that do not have an absorption edge detected
    in *elementsDict* and an extrapolated pre0edge σ² for the element with a
    detected absorption edge.
    """
    bknd = np.zeros_like(E)
    for element in elementsDict.values():
        if element[2] > 0:  # has edge
            vicc, vicd = element[6:8]
            vic = vicc*np.exp(-3*np.log(E)) + vicd*np.exp(-4*np.log(E))
            bknd += vic * element[3]
        else:
            f1f2 = element[0].get_f1f2(E)
            bknd += f1f2.imag * element[3] * crossSection / E
    return bknd


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


def parse_compound(compound, mass_digit=5):
    """
    If successful, returns a tuple (parsed_result, mass_str), where
    *parsed_result* is a list of lists [element, mole_amount] and
    *mass_str* is a str representation of the compound mass.

    If unsuccessful, returns a str of the error statement.
    """

    if len(compound) == 0:
        return "no compound formula"
    try:
        res = formula.parseString(compound, parseAll=True)
        cMass = 0.
        for e in res.asList():
            if e[1] < 0:
                return "compound formula with x"
            else:
                cMass += rm.read_atomic_data(e[0])*e[1]
        if '%' in compound:
            return parse_compound(reconstruct(res))
        else:
            mass_str = '{0}'.format(round_to_n(cMass, mass_digit))
    except (ParseBaseException, ValueError, ZeroDivisionError):
        return "wrong compound formula"
    return res.asList(), mass_str
