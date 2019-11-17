import pandas as pd
from time import time
import collections
from pymatgen.entries.computed_entries import ComputedEntry
import os
import numpy as np
from pymatgen.analysis.phase_diagram import PhaseDiagram
from itertools import combinations, permutations, chain

def CalculateOmegas(hdf_a, hdf_ab, key_a, key_ab, TotalNumberAtomsStructure=16):
    '''
    Give two hdf data files containing the elemental 
    energies and the energies of the binary systems,
    in the following format:
       formula series  natoms composition reduced        magtot     etotal  etotal_per_atom
    0      Al1      a       1        [Al]      Al  2.297300e-03  -3.644295  -3.644295

       formula series  natoms composition reduced    magtot     etotal  etotal_per_atom 
    0   Al8Ag8     ab      16    [Ag, Al]    AgAl -0.000002 -52.260589   -3.266287 
    and also you have to give as an input the keys for the elemental and the binary database.
    In short the inputs are hdf_a, key_a, hdf_ab, key_ab'''
    df_a = pd.read_hdf(hdf_a, key_a)
    df_ab = pd.read_hdf(hdf_ab, key_ab)
    summary = []
    counter = -1
    for j in range(0, len(df_ab)):
        counter = counter + 1
        for i in range(0, len(df_a)):
            if df_ab.composition[j][0] == df_a.composition[i][0]:
                summary.append((df_ab.etotal[j] - 8*df_a.etotal[i]))
                for k in range(0, len(df_a)):
                    if df_ab.composition[j][1] == df_a.composition[k][0]:
                        summary[counter] = summary[counter] - 8*df_a.etotal[k]
                        break
                else:
                    summary[counter] = summary[counter] - 8*df_a.etotal[i] 
    for i in range(len(summary)):
        summary[i] = summary[i]/16
    df_ab['H'] = summary
    df_ab['Omega'] = df_ab['H']/0.25
    aggregation_functions = {'formula': 'first', 'series': 'first', 'natoms': 'first', 
                         'composition':'first', 'magtot':'mean', 'etotal':'mean',
                        'etotal_per_atom':'mean', 'H':'mean', 'Omega':'mean', 
                         'Omega_std':'first','etotal_std':'first'}
    df_ab['Omega_std'] = df_ab['Omega'].groupby(df_ab['reduced']).transform('std')
    df_ab['etotal_std'] = df_ab['etotal'].groupby(df_ab['reduced']).transform('std')
    df_ab = df_ab.groupby(df_ab['reduced']).aggregate(aggregation_functions)
    return df_a, df_ab


def CalculateDH(df_a, df_ab, numberElements, elems):
    """
    Give two databases containing the elemental 
    energies and the energies of the binary systems,
    the number of elements inside the alloy, and a list
    of the elements that you want to calculate the DH
    in the following format:
       formula series  natoms composition reduced        magtot     etotal	etotal_per_atom
    0      Al1      a       1        [Al]      Al  2.297300e-03  -3.644295  -3.644295

       formula series  natoms composition reduced    magtot     etotal  etotal_per_atom	H     Omega     Omega_std    etotal_std 
    0   Al8Ag8     ab      16    [Ag, Al]    AgAl -0.000002 -52.260589   -3.266287	-0.027820 -0.111281  7.996541e-02  3.198617e-01 
    and also you have to give as an input the keys for the elemental and the binary database.
    you can create the binary database by calling CalculateOmega module.
    """
    values = []
    reduced = []
    collect = combinations(elems, numberElements)  
    first = []
    second = []
    third = []
    forth = []
    fifth = []
    decomp = []
    if numberElements == 3:
        print("inside 3")
        for abc in collect:
            reduced.append("".join(sorted([abc[0], abc[1], abc[2]])))
            decomp.append([abc[0], abc[1], abc[2]])
            first.append(abc[0])
            second.append(abc[1])
            third.append(abc[2])
            sum = 0
            counterInside = 0
            for index, row in df_ab.iterrows():
                if set(row['composition']).issubset(abc):
                    sum += row['Omega']*1/9
            values.append(sum)
        abc_results = pd.DataFrame(data=first, columns=['first'])
        abc_results['second'] = second
        abc_results['third'] = third
        abc_results['reduced'] = reduced
        abc_results['decomposition'] = decomp      
        abc_results['H_predicted'] = values
        return abc_results
    elif numberElements == 4:
        print("inside 4")
        for abc in collect:
            reduced.append("".join(sorted([abc[0], abc[1], abc[2], abc[3]])))
            decomp.append([abc[0], abc[1], abc[2], abc[3]])
            first.append(abc[0])
            second.append(abc[1])
            third.append(abc[2])
            forth.append(abc[3])
            sum = 0
            for index, row in df_ab.iterrows():
                if set(row['composition']).issubset(abc):
                    sum += row['Omega']*1/16
            values.append(sum)
        abcd_results = pd.DataFrame(data=first, columns=['first'])
        abcd_results['second'] = second
        abcd_results['third'] = third
        abcd_results['forth'] = forth
        abcd_results['reduced'] = reduced
        abcd_results['decomposition'] = decomp      
        abcd_results['H_predicted'] = values
        return abcd_results
    elif numberElements == 5:
        print("inside 5")
        for abc in collect:
            reduced.append("".join(sorted([abc[0], abc[1], abc[2], abc[3],abc[4]])))
            decomp.append([abc[0], abc[1], abc[2], abc[3],abc[4]])
            first.append(abc[0])
            second.append(abc[1])
            third.append(abc[2])
            forth.append(abc[3])
            fifth.append(abc[4])
            sum = 0
            for index, row in df_ab.iterrows():
                if set(row['composition']).issubset(abc):
                    sum += row['Omega']*1/25
            values.append(sum)
        abcde_results = pd.DataFrame(data=first, columns=['first'])
        abcde_results['second'] = second
        abcde_results['third'] = third
        abcde_results['forth'] = forth
        abcde_results['fifth'] = fifth
        abcde_results['reduced'] = reduced
        abcde_results['decomposition'] = decomp      
        abcde_results['H_predicted'] = values
        return abcde_results

def CalculateEnergy(nameDatabase, numberElements, elems):
    df_a = pd.read_hdf(nameDatabase, key='a')
    df_ab = pd.read_hdf(nameDatabase, key='ab')
    df_abc = pd.read_hdf(nameDatabase, key='abc')
    df_abcd = pd.read_hdf(nameDatabase, key='abcd')
    df_abcde = pd.read_hdf(nameDatabase, key='abcde')
    values = []
    collect = combinations(elems, numberElements)  
    totalEnergy = []
    if numberElements == 3:
        for abc in collect:
            sum1 = 12*df_a[df_a['reduced']==abc[0]].etotal.values
            sum2 = 12*df_a[df_a['reduced']==abc[1]].etotal.values
            sum3 = 12*df_a[df_a['reduced']==abc[2]].etotal.values
            values = sum1+sum2+sum3
            totalEnergy.append(26*df_abc[df_abc['reduced'] == "".join(sorted([abc[0], abc[1], abc[2]]))]['H_predicted'].values[0]+values[0])
        df_abc['Total_energy'] = totalEnergy
        return df_abc
    if numberElements == 4:
        for abc in collect:
            sum1 = 12*df_a[df_a['reduced']==abc[0]].etotal.values
            sum2 = 12*df_a[df_a['reduced']==abc[1]].etotal.values
            sum3 = 12*df_a[df_a['reduced']==abc[2]].etotal.values
            sum4 = 12*df_a[df_a['reduced']==abc[3]].etotal.values
            values = sum1+sum2+sum3+sum4
            totalEnergy.append(48*df_abcd[df_abcd['reduced'] == "".join(sorted([abc[0], abc[1], abc[2],abc[3]]))]['H_predicted'].values[0]+values[0])
        df_abcd['Total_energy'] = totalEnergy
        return df_abcd
    if numberElements == 5:
        for abc in collect:
            sum1 = 12*df_a[df_a['reduced']==abc[0]].etotal.values
            sum2 = 12*df_a[df_a['reduced']==abc[1]].etotal.values
            sum3 = 12*df_a[df_a['reduced']==abc[2]].etotal.values
            sum4 = 12*df_a[df_a['reduced']==abc[3]].etotal.values
            sum5 = 12*df_a[df_a['reduced']==abc[4]].etotal.values
            values = sum1+sum2+sum3+sum4+sum5
            totalEnergy.append(60*df_abcde[df_abcde['reduced'] == "".join(sorted([abc[0], abc[1], abc[2],abc[3],abc[4]]))]['H_predicted'].values[0]+values[0])
        df_abcde['Total_energy'] = totalEnergy
        return df_abcde


#def addcomputedentries(namedb, howManyElements = 5, structurename='FCC', temperature=0):
#    entriesInside = []
#    NumberofAtoms = [1, 16, 36, 48, 60]
#    keys = ['a', 'ab', 'abc_with_total', 'abcd_with_total', 'abcde_with_total']
#    for i in range(1, howManyElements+1):
#        Datafr = pd.read_hdf(namedb, key=keys[i-1]))
#        for j in range(0, len(Datafr)):
#            energy = Datafr.etotal[j]
#            composition = Datafr.formula[j]
#            computedEntries = ComputedEntry(str(composition), energy - entropy_mix(i+0.00000000000001)*temperature*NumberofAtoms[i-1], attribute=structurename, data=entropy_mix(i+0.00000000000001)*temperature*NumberofAtoms[i-1])
#            entriesInside.append(computedEntries)
#    return entriesInside


#@numba.jit
def entropy_mix(n):
    '''entropy per system'''
    kb = 8.6173303e-5
    s = kb*np.log(n)
    return s

def computeValues (elems, nameDatabase):
    df_a = pd.read_hdf(nameDatabase, key='a')
    df_ab = pd.read_hdf(nameDatabase, key='ab')
    df_abc = pd.read_hdf(nameDatabase, key='abc_with_total')
    df_abcd = pd.read_hdf(nameDatabase, key='abcd_with_total')
    df_abcde = pd.read_hdf(nameDatabase, key='abcde_with_total')
    collect = list(combinations(elems, 5))[30000:40000]
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
        data["Composition"].append(e.composition.reduced_formula)
        data['NElements'].append(e.attribute)
        data["Ehull"].append(ehull)
        #data['Entropy'].append(e.data)
        #data['Gibbs'].append(ehull)
        data["Decomposition"].append("self")
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
        print(df)
    return(df)


def SplitUpperCase(text):
    result = ''
    for char in text:
        if char.isupper():
            result += " " + char
        else:
            result += char
    return result.split()


if __name__ == "__main__":
    df = computeValues(['Al', 'Co', 'Cr', 'Cu', 'Fe', 'Hf', 'Mn','Mo', 'Nb', 'Ni', 'Ta', 'Ti', 'W', 'Zr','V', 'Mg', 'Re', 'Os', 'Rh', 'Ir', 'Pd','Pt', 'Ag', 'Au', 'Zn', 'Cd', 'Hg'], "data_FCC.h5")
    df.to_hdf('data_FCC.h5',key='all_together', format='table', append=True, min_itemsize={'Decomposition' : 70})
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

