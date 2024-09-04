import pandas as pd
import re
import collections

#formats the cysteine and SS columns according to the format required for the rest of the code
def format_columns():
    c_list, d_list = [],[]
    for _, element in df.iterrows():
        cys = list(map(int,re.findall('\d+',element["Sequence C position"])))
        c_list.append([i+1 for i in cys])
        dsb = list(map(int,re.findall('\d+',element["Disulfide bond"])))
        d_list.append([[dsb[i],dsb[i+1]] for i in range(0, len(dsb)-1,2)])
    
    df['Cysteine positions'] = c_list
    df['Disulfide positions'] = d_list

#calculates the distance between the starting and ending positions of each disulfide bond
def calculate_intraSS_dist():
    dist = []
    for _, element in df.iterrows():
        value = [element['Disulfide positions'][i][1] - element['Disulfide positions'][i][0] for i in range(0, len(element['Disulfide positions']))]
        dist.append(value)
    df['Intra-SS distance'] = dist

#calculates the distance between each disulfide bond
def calculate_interSS_dist():
    dist = []
    for _, element in df.iterrows():
        value = [element['Disulfide positions'][i+1][0] - element['Disulfide positions'][i][1] for i in range(0, len(element['Disulfide positions'])-1)]
        dist.append(value)
    df['Inter-SS distance'] = dist

#calculates the distance between each cysteine residue in the sequence
def calculate_interCys_dist():
    dist = []
    for _, element in df.iterrows():
        value = [element['Cysteine positions'][i+1] - element['Cysteine positions'][i] for i in range(0, len(element['Cysteine positions'])-1)]
        dist.append(value)
    df['Inter-Cys distance'] = dist

#calculates the frequency of each distance between the cysteine residues
def calculate_interCys_frequency():
    freq = []
    for _, element in df.iterrows():
        counter = collections.Counter(element['Inter-Cys distance'])
        freq.append([str(x) + '-' + str(counter[x]) for x in sorted(counter.keys())])
    df['Inter-Cys frequency'] = freq

#cacluate the frequency of the distance between the starting and ending positions of a disulfide bond
def calculate_intraSS_frequency():
    freq = []
    for _, element in df.iterrows():
        counter = collections.Counter(element['Intra-SS distance'])
        freq.append([str(x) + '-' + str(counter[x]) for x in sorted(counter.keys())])
    df['Intra-SS frequency'] = freq

#calculate the freqeuncy of the distances between each disulfide bond
def calculate_interSS_frequency():
    freq = []
    for _, element in df.iterrows():
        counter = collections.Counter(element['Inter-SS distance'])
        freq.append([str(x) + '-' + str(counter[x]) for x in sorted(counter.keys())])
    df['Inter-SS frequency'] = freq

df = pd.read_csv('C_and_SS_input.csv')
format_columns()
calculate_interCys_dist()
calculate_interCys_frequency()
calculate_intraSS_dist()
calculate_intraSS_frequency()
calculate_interSS_dist()
calculate_interSS_frequency()

output_df = df.drop(['Sequence C position', 'Disulfide positions'],axis=1)
output_df.to_csv('C_and_SS_output.csv')