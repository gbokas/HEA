import pandas as pd
from time import time
import collections
from pymatgen.entries.computed_entries import ComputedEntry
import os
import numpy as np
from pymatgen.analysis.phase_diagram import PhaseDiagram
from itertools import combinations, permutations, chain


#@numba.jit
def entropy_mix(n):
    '''entropy per system'''
    kb = 8.6173303e-5
    s = kb*np.log(n)
    return s


def computeValuesBinary(elems, nameDatabase):
    df_a = pd.read_hdf(nameDatabase, key='a')
    df_ab = pd.read_hdf(nameDatabase, key='ab')
    collect = combinations(elems, 2)
    numberofElements = []
    final_energy = []
    counter = 0
    for element in collect:
        try:
            if not totalDf.empty:
                pass
        except:
            totalDf = pd.DataFrame()
            totalDf_stable = pd.DataFrame()
        entriesInside = []
        counter += 1
        computedEntries = ComputedEntry(element[0],df_a[df_a['reduced']==element[0]].etotal.values[0],attribute=1)
        entriesInside.append(computedEntries)
        computedEntries = ComputedEntry(element[1],df_a[df_a['reduced']==element[1]].etotal.values[0],attribute=1)
        entriesInside.append(computedEntries)
        for index, row in df_ab.iterrows():
             if set(row['composition']).issubset(element):
                computedEntries = ComputedEntry(row['formula'], row['etotal'],attribute=2)
                entriesInside.append(computedEntries)
                break
        df, df_stable = convexHull(entriesInside)
        totalDf_stable = totalDf_stable.append(df_stable, ignore_index = True)
        totalDf = totalDf.append(df, ignore_index = True) 
    return totalDf, totalDf_stable


def computeValuesTernary(elems, nameDatabase):
    df_a = pd.read_hdf(nameDatabase, key='a')
    df_abc = pd.read_hdf(nameDatabase, key='abc_with_total')
    collect = combinations(elems, 3)
    numberofElements = []
    final_energy = []
    counter = 0
    for element in collect:
        df_stable=pd.read_hdf('data_FCC.h5', key="stable_binaries_Hull")
        counter += 1
        print(counter, element)
        try:
            if not totalDf.empty:
                pass
        except:
            totalDf = pd.DataFrame()
            totalDf_stable = pd.DataFrame()
        entriesInside = []
        counter += 1
        computedEntries = ComputedEntry(element[0],df_a[df_a['reduced']==element[0]].etotal.values[0],attribute=1)
        entriesInside.append(computedEntries)
        computedEntries = ComputedEntry(element[1],df_a[df_a['reduced']==element[1]].etotal.values[0],attribute=1)
        entriesInside.append(computedEntries)
        computedEntries = ComputedEntry(element[2],df_a[df_a['reduced']==element[2]].etotal.values[0],attribute=1)
        entriesInside.append(computedEntries)
        for index, row in df_stable.iterrows():
            if set(SplitUpperCase(row['Composition'])).issubset(element):
                computedEntries = ComputedEntry(row['Formula'], row['energy'],attribute=row['NElements'])
                entriesInside.append(computedEntries)
        for index, row in df_abc.iterrows():
             if set(row['decomposition']).issubset(element):
                computedEntries = ComputedEntry(row['decomposition'][0]+"12"+row['decomposition'][1]+"12"+row['decomposition'][2]+"12", row['Total_energy'],attribute=3)
                entriesInside.append(computedEntries)
                break
        df, df_stable = convexHull(entriesInside)
        totalDf_stable = totalDf_stable.append(df_stable, ignore_index = True)
        totalDf = totalDf.append(df, ignore_index = True) 
    return totalDf, totalDf_stable


def computeValuesQuaternary(elems, nameDatabase):
    df_a = pd.read_hdf(nameDatabase, key='a')
    df_abcd = pd.read_hdf(nameDatabase, key='abcd_with_total')
    collect = combinations(elems, 4)
    numberofElements = []
    final_energy = []
    counter = 0
    for element in collect:
        df_stable=pd.read_hdf('data_FCC.h5', key="stable_ternaries_Hull")
        counter += 1
        try:
            if not totalDf.empty:
                pass
        except:
            totalDf = pd.DataFrame()
            totalDf_stable = pd.DataFrame()
        entriesInside = []
        computedEntries = ComputedEntry(element[0],df_a[df_a['reduced']==element[0]].etotal.values[0],attribute=1)
        entriesInside.append(computedEntries)
        computedEntries = ComputedEntry(element[1],df_a[df_a['reduced']==element[1]].etotal.values[0],attribute=1)
        entriesInside.append(computedEntries)
        computedEntries = ComputedEntry(element[2],df_a[df_a['reduced']==element[2]].etotal.values[0],attribute=1)
        entriesInside.append(computedEntries)
        computedEntries = ComputedEntry(element[3],df_a[df_a['reduced']==element[3]].etotal.values[0],attribute=1)
        entriesInside.append(computedEntries)
        for index, row in df_stable.iterrows():
            if set(SplitUpperCase(row['Composition'])).issubset(element):
                computedEntries = ComputedEntry(row['Formula'], row['energy'],attribute=row['NElements'])
                entriesInside.append(computedEntries)
        for index, row in df_abcd.iterrows():
             if set(row['decomposition']).issubset(element):
                computedEntries = ComputedEntry(row['decomposition'][0]+"12"+row['decomposition'][1]+"12"+row['decomposition'][2]+"12"+row['decomposition'][3]+"12", row['Total_energy'],attribute=4)
                entriesInside.append(computedEntries)
                break
        df, df_stable = convexHull(entriesInside)
        totalDf_stable = totalDf_stable.append(df_stable, ignore_index = True)
        totalDf = totalDf.append(df, ignore_index = True) 
        print(counter, element)
    return totalDf, totalDf_stable

def computeValuesQuinary(elems, nameDatabase):
    df_a = pd.read_hdf(nameDatabase, key='a')
    df_abcde = pd.read_hdf(nameDatabase, key='abcde_with_total')
    collect = combinations(elems, 5)
    numberofElements = []
    final_energy = []
    counter = 0
    for element in collect:
        df_stable=pd.read_hdf('data_FCC.h5', key="stable_quaternaries_Hull")
        try:
            if not totalDf.empty:
                pass
        except:
            totalDf = pd.DataFrame()
            totalDf_stable = pd.DataFrame()
        entriesInside = []
        counter += 1
        computedEntries = ComputedEntry(element[0],df_a[df_a['reduced']==element[0]].etotal.values[0],attribute=1)
        entriesInside.append(computedEntries)
        computedEntries = ComputedEntry(element[1],df_a[df_a['reduced']==element[1]].etotal.values[0],attribute=1)
        entriesInside.append(computedEntries)
        computedEntries = ComputedEntry(element[2],df_a[df_a['reduced']==element[2]].etotal.values[0],attribute=1)
        entriesInside.append(computedEntries)
        computedEntries = ComputedEntry(element[3],df_a[df_a['reduced']==element[3]].etotal.values[0],attribute=1)
        entriesInside.append(computedEntries)
        computedEntries = ComputedEntry(element[4],df_a[df_a['reduced']==element[4]].etotal.values[0],attribute=1)
        entriesInside.append(computedEntries)
        for index, row in df_stable.iterrows():
            if set(SplitUpperCase(row['Composition'])).issubset(element):
                computedEntries = ComputedEntry(row['Formula'], row['energy'],attribute=row['NElements'])
                entriesInside.append(computedEntries)
        for index, row in df_abcde.iterrows():
             if set(row['decomposition']).issubset(element):
                computedEntries = ComputedEntry(row['decomposition'][0]+"12"+row['decomposition'][1]+"12"+row['decomposition'][2]+"12"+row['decomposition'][3]+"12"
				+row['decomposition'][4]+"12", row['Total_energy'],attribute=5)
                entriesInside.append(computedEntries)
                break
        df, df_stable = convexHull(entriesInside)
        totalDf_stable = totalDf_stable.append(df_stable, ignore_index = True)
        totalDf = totalDf.append(df, ignore_index = True) 
        print(counter, element)
    return totalDf, totalDf_stable


def computeValues (elems, nameDatabase):
    df_a = pd.read_hdf(nameDatabase, key='a')
    df_ab = pd.read_hdf(nameDatabase, key='ab')
    df_abc = pd.read_hdf(nameDatabase, key='abc_with_total')
    df_abcd = pd.read_hdf(nameDatabase, key='abcd_with_total')
    df_abcde = pd.read_hdf(nameDatabase, key='abcde_with_total')
    collect = list(combinations(elems, 5))[45000:47000]
    numberofElements = []
    final_energy = []
    counter = 0
    for element in collect:
        print(counter, element)
        try:
            if not totalDf.empty:
                pass
        except:
            totalDf = pd.DataFrame()
        entriesInside = []
        counter += 1
        computedEntries = ComputedEntry(str(element[0]),df_a[df_a['reduced']==element[0]].etotal.values[0],attribute=1)
        entriesInside.append(computedEntries)
        computedEntries = ComputedEntry(element[1],df_a[df_a['reduced']==element[1]].etotal.values[0],attribute=1)
        entriesInside.append(computedEntries)
        computedEntries = ComputedEntry(element[2],df_a[df_a['reduced']==element[2]].etotal.values[0],attribute=1)
        entriesInside.append(computedEntries)
        computedEntries = ComputedEntry(element[3],df_a[df_a['reduced']==element[3]].etotal.values[0],attribute=1)
        entriesInside.append(computedEntries)
        computedEntries = ComputedEntry(element[4],df_a[df_a['reduced']==element[4]].etotal.values[0],attribute=1)
        entriesInside.append(computedEntries)
        counterstop = 0
        for index, row in df_ab.iterrows():
             if set(row['composition']).issubset(element):
                counterstop += 1
                computedEntries = ComputedEntry(row['formula'], row['etotal'],attribute=2)
                entriesInside.append(computedEntries)
                if counterstop == 10:
                    break
        counterstop = 0
        for index, row in df_abc.iterrows():
             if set(row['decomposition']).issubset(element):
                counterstop += 1
                computedEntries = ComputedEntry(row['decomposition'][0]+"12"+row['decomposition'][1]+"12"+row['decomposition'][2]+"12", row['Total_energy'],attribute=3)
                entriesInside.append(computedEntries)
                if counterstop == 10:
                    break
        counterstop = 0
        for index, row in df_abcd.iterrows():
             if set(row['decomposition']).issubset(element):
                counterstop += 1
                computedEntries = ComputedEntry(row['decomposition'][0]+"12"+row['decomposition'][1]+"12"+row['decomposition'][2]+"12"+row['decomposition'][3]+"12", row['Total_energy'],attribute=4)
                entriesInside.append(computedEntries)
                if counterstop == 5:
                    break
        counterstop = 0
        for index, row in df_abcde.iterrows():
             if set(row['decomposition']).issubset(element):
                counterstop += 1
                computedEntries = ComputedEntry(row['decomposition'][0]+"12"+row['decomposition'][1]+"12"+row['decomposition'][2]+"12"+row['decomposition'][3]+"12"
				+row['decomposition'][4]+"12", row['Total_energy'],attribute=5)
                entriesInside.append(computedEntries)
                break
        df = convexHull(entriesInside)
        totalDf = totalDf.append(df, ignore_index = True) 
    return totalDf 
    
        

#                computedEntries = ComputedEntry(a + '12' + b + '12' + c + '12',
#                                                etotal - entropy_mix(36 + 0.00000000000001) * temperature * 36)



def convexHull(entries, tolerance=0, printing = False):
    # This initializes the REST adaptor. You may need to put your own API key in as an arg. If we need MP DB
    #with open("key.txt", "r") as myfile:
    #    data = myfile.readlines()
    #a = MPRester(data[0].rstrip())
    #    entries = a.get_entries_in_chemsys(elems)
    pd2 = PhaseDiagram(entries)
    data = collections.defaultdict(list)
    for e in pd2.stable_entries:
        ehull = pd2.get_equilibrium_reaction_energy(e)
        data["energy"].append(e.energy)
        data['Formula'].append(e.composition.alphabetical_formula.replace(" ",""))
        data["Composition"].append(e.composition.reduced_formula)
        data['NElements'].append(e.attribute)
        data["Ehull"].append(ehull)
        #data['Entropy'].append(e.data)
        #data['Gibbs'].append(ehull)
        data["Decomposition"].append("self")
        df_stable = pd.DataFrame(data, columns=["NElements", "Composition", "Formula", "Ehull", 'energy', "Decomposition"])
    for e in pd2.unstable_entries:
        decomp, ehull = pd2.get_decomp_and_e_above_hull(e, allow_negative=True)
#        decomp, ehull = pd2.get_decomp_and_e_above_hull(e)
        data["Composition"].append(e.composition.reduced_formula)
        data['NElements'].append(e.attribute)
        data["Ehull"].append(ehull)
        #data['Entropy'].append(e.data)
        #data['Gibbs'].append(ehull)
        data["Decomposition"].append(" + ".join(["%.2f %s" % (v, k.composition.formula) for k, v in decomp.items()])[0:70])
    df = pd.DataFrame(data, columns=["NElements", "Composition", "Ehull", "Decomposition"])
    if printing == True:
        #print(df[df.Ehull<=tolerance])
        print(df_stable)
    return df, df_stable


def SplitUpperCase(text):
    result = ''
    for char in text:
        if char.isupper():
            result += " " + char
        else:
            result += char
    return result.split()


if __name__ == "__main__":
    #df, df_stable = computeValuesBinary(['Al', 'Co', 'Cr', 'Cu', 'Fe', 'Hf', 'Mn','Mo', 'Nb', 'Ni', 'Ta', 'Ti', 'W', 'Zr','V', 'Mg', 'Re', 'Os', 'Rh', 'Ir', 'Pd','Pt', 'Ag', 'Au', 'Zn', 'Cd', 'Hg'], "data_FCC.h5")
    #df.drop_duplicates(inplace=True)
    #df_stable.drop_duplicates(inplace=True)
    #df.to_hdf('data_FCC.h5',key='binaries_Hull', format='table', min_itemsize={'Decomposition' : 70})
    #df_stable.to_hdf('data_FCC.h5', key="stable_binaries_Hull", format='table', min_itemsize={'Decomposition':70})
    #df, df_stable = computeValuesTernary(['Al', 'Co', 'Cr', 'Cu', 'Fe', 'Hf', 'Mn','Mo', 'Nb', 'Ni', 'Ta', 'Ti', 'W', 'Zr','V', 'Mg', 'Re', 'Os', 'Rh', 'Ir', 'Pd','Pt', 'Ag', 'Au', 'Zn', 'Cd', 'Hg'], "data_FCC.h5")
    #elements = ["Co", "Cr", "Fe"]

    #for compo, energy in zip(df_stable["Composition"], df_stable['Ehull']):
    #    if set(SplitUpperCase(compo)).issubset(elements):
    #        print(compo, energy)
    #df.drop_duplicates(inplace=True)
    #df_stable.drop_duplicates(inplace=True)
    #for compo, energy in zip(df_stable["Composition"], df_stable['Ehull']):
    #    if set(SplitUpperCase(compo)).issubset(elements):
    #        print(compo, energy)
    #df.to_hdf('data_FCC.h5',key='ternaries_Hull', format='table', min_itemsize={'Decomposition' : 70})
    #df_stable.to_hdf('data_FCC.h5', key="stable_ternaries_Hull", format='table', min_itemsize={'Decomposition':70})
    df, df_stable = computeValuesQuaternary(['Al', 'Co', 'Cr', 'Cu', 'Fe', 'Hf', 'Mn','Mo', 'Nb', 'Ni', 'Ta', 'Ti', 'W', 'Zr','V', 'Mg', 'Re', 'Os', 'Rh', 'Ir', 'Pd','Pt', 'Ag', 'Au', 'Zn', 'Cd', 'Hg'], "data_FCC.h5")
    df.drop_duplicates(inplace=True)
    df_stable.drop_duplicates(inplace=True)
    df.to_hdf('data_FCC.h5',key='quaternaries_Hull', format='table', min_itemsize={'Decomposition' : 70})
    df_stable.to_hdf('data_FCC.h5', key="stable_quaternaries_Hull", format='table', min_itemsize={'Decomposition':70})
    df, df_stable = computeValuesQuinary(['Al', 'Co', 'Cr', 'Cu', 'Fe', 'Hf', 'Mn','Mo', 'Nb', 'Ni', 'Ta', 'Ti', 'W', 'Zr','V', 'Mg', 'Re', 'Os', 'Rh', 'Ir', 'Pd','Pt', 'Ag', 'Au', 'Zn', 'Cd', 'Hg'], "data_FCC.h5")
    df.drop_duplicates(inplace=True)
    df_stable.drop_duplicates(inplace=True)
    df.to_hdf('data_FCC.h5',key='quinaries_Hull', format='table', min_itemsize={'Decomposition' : 70})
    df_stable.to_hdf('data_FCC.h5', key="stable_quinaries_Hull", format='table', min_itemsize={'Decomposition':70})
    #df.to_hdf('data_FCC.h5',key='all_together',format='table',min_itemsize={'Decomposition' : 70})
#    df_abc = CalculateEnergy("data_FCC.h5", 3, ['Al', 'Co', 'Cr', 'Cu', 'Fe', 'Hf', 'Mn','Mo', 'Nb', 'Ni', 'Ta', 'Ti', 'W', 'Zr','V', 'Mg', 'Re', 'Os', 'Rh', 'Ir', 'Pd','Pt', 'Ag', 'Au', 'Zn', 'Cd', 'Hg'])
#    print("finish 3")
#    df_abcd = CalculateEnergy("data_FCC.h5", 4, ['Al', 'Co', 'Cr', 'Cu', 'Fe', 'Hf', 'Mn','Mo', 'Nb', 'Ni', 'Ta', 'Ti', 'W', 'Zr','V', 'Mg', 'Re', 'Os', 'Rh', 'Ir', 'Pd','Pt', 'Ag', 'Au', 'Zn', 'Cd', 'Hg'])
#    print("finish 4")
#    df_abcde = CalculateEnergy("data_FCC.h5", 5, ['Al', 'Co', 'Cr', 'Cu', 'Fe', 'Hf', 'Mn','Mo', 'Nb', 'Ni', 'Ta', 'Ti', 'W', 'Zr','V', 'Mg', 'Re', 'Os', 'Rh', 'Ir', 'Pd','Pt', 'Ag', 'Au', 'Zn', 'Cd', 'Hg'])
#    print("finish 5")
    
#    df_a, df_ab = CalculateOmegas('totaldoublerelax_a.h5', 'totaldoublerelax_ab.h5', 'a', 'ab') 
#    df_abc = CalculateDH(df_a, df_ab, 3, ['Al', 'Co', 'Cr', 'Cu', 'Fe', 'Hf', 'Mn','Mo', 'Nb', 'Ni', 'Ta', 'Ti', 'W', 'Zr','V', 'Mg', 'Re', 'Os', 'Rh', 'Ir', 'Pd','Pt', 'Ag', 'Au', 'Zn', 'Cd', 'Hg'])
#    df_abcd = CalculateDH(df_a, df_ab, 4, ['Al', 'Co', 'Cr', 'Cu', 'Fe', 'Hf', 'Mn','Mo', 'Nb', 'Ni', 'Ta', 'Ti', 'W', 'Zr','V', 'Mg', 'Re', 'Os', 'Rh', 'Ir', 'Pd','Pt', 'Ag', 'Au', 'Zn', 'Cd', 'Hg'])
#    df_abcde = CalculateDH(df_a, df_ab, 5, ['Al', 'Co', 'Cr', 'Cu', 'Fe', 'Hf', 'Mn','Mo', 'Nb', 'Ni', 'Ta', 'Ti', 'W', 'Zr','V', 'Mg', 'Re', 'Os', 'Rh', 'Ir', 'Pd','Pt', 'Ag', 'Au', 'Zn', 'Cd', 'Hg'])
#    df_a.to_hdf('data_FCC.h5', key='a')
#    df_ab.to_hdf('data_FCC.h5', key='ab')
#    df_abc.to_hdf('data_FCC.h5', key='abc_with_total')
#    df_abcd.to_hdf('data_FCC.h5', key='abcd_with_total')
#    df_abcde.to_hdf('data_FCC.h5', key='abcde_with_total')

