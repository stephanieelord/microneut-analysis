import pandas as pd
import math
import numpy as np
import warnings
from scipy.optimize import minimize
import matplotlib.pyplot as plt

import numpy as np

from scipy.optimize import fmin_l_bfgs_b


#get file path
#file_path =

#get sheet name
#sheet_name = input("Enter sheet name:")
sheet_name = r"C:\Users\steph\Downloads\XV-304 RSV A2 Microneuts D0.xlsx"
#make an alias --config file
#add a directory browser popup

#reading the excel
df = pd.read_excel(sheet_name, "day 0 Raw Data", header = None)

#constants
initial_cell = input("Enter first plate cell:")
initial_split = list(initial_cell)
initial_column = ord(initial_split[0])-65
print("initial column: ", initial_column)
initial_row = int(initial_split[1])-1
print("initial row: ", initial_row)
plate_data_rows = 8
plate_data_cols = 12
sample_label_rows = 4
sample_label_cols = 3
sample_label_offset = 15 #column R
plate_spacing = 12

#count plates
max_row = df.shape[0]
num_plates = (max_row + 4 - initial_row) // plate_spacing
#print("Max row: ", max_row)
#print(max_row," - ",initial_row," / ",plate_spacing," = ")
#print("Number of plates: ", num_plates)

#making array to store sample ID and IP
storage_array = []
rows = 0
cols = 2
storage_array = [[0 for _ in range(cols)] for _ in range(rows)]

skipped = []
storage_array = [[0 for _ in range(cols)] for _ in range(rows)]

minimum = 4.32
for plate in range(num_plates):

    #loop through plates
    plate_index = plate+1 #one for now
    data_row_start = initial_row + (plate_index-1) * plate_spacing
    label_row_start = data_row_start + 1

    #extraxt data

    #slice data
    plate_data = df.iloc[data_row_start : data_row_start + plate_data_rows, initial_column : initial_column+plate_data_cols] 
    #same thing
    sample_names = df.iloc[label_row_start : label_row_start + sample_label_rows,
                           sample_label_offset : sample_label_offset + sample_label_cols]

    #slice sample layout
    sample_layout = df.iloc[initial_row+1:initial_row+1+sample_label_rows, sample_label_offset: sample_label_offset+sample_label_cols]

        #print("data row start: ", data_row_start)
        #print("label row start: ", label_row_start)
    #print("Sample names: \n",sample_names,"\n")

    #creating df
    #sample labels creation
    #samples = [
     #   817, 818, 819,
      #  820, 821, 822,
       # 823, 824, 825,
        #826, 827, "Control"]

    samples = sample_names.values.flatten().tolist()
    #make it customizebale
    duplicate_suffixes = ["a", "b", "c", "d"]
    duplicate_IDs = [1, 2]
    labels = []

    #putting it in the correct order
    for i in range(0, len(samples), 3):
        sample_chunk = samples[i:i+3]
        for ID in duplicate_IDs:
            for sample in sample_chunk:
                for suffix in duplicate_suffixes:
                    labels.append(f"{sample}_{ID}{suffix}")
    #print("Labels: \n",labels)

    #dataframe with data
    dilution_map = {
        "1": 20,
        "2": 80,
        "3": 320,
        "4": 1280}
    samples = []
    dilution_ids = []
    dilutions = []
    duplicates = []


    #loop through to make dataframe
    for label in labels:
        sample, dil_duplicate = label.split("_")
        dil_id = dil_duplicate[0]
        duplicate = dil_duplicate[1]


        #map ID to value
        dilution = dilution_map[dil_id]

        #check if sample is control
        if sample == 'Control' or sample == 'Cntrls':  #or whatever control is
            if dil_id in ['1', '2']:
                control_type = 'Positive'
            elif dil_id in ['3', '4']:
                control_type = 'Negative'
            else:
                control_type = 'Unknown'
        else:
            control_type = 'Sample'

        #append to list
        samples.append(sample)
        dilution_ids.append(dil_id)
        dilutions.append(dilution)
        duplicates.append(duplicate)

    df_labels = pd.DataFrame({
        "Label": labels,
        "Sample": samples,
        "Dilution_ID": dilution_ids,
        "Dilution": dilutions,
        "Duplicate": duplicates })

    #flatten the OD readings and add to DF
    flat_od_values = plate_data.to_numpy().flatten()
    df_labels["OD_value"] = flat_od_values
    #print("Labels:\n",df_labels)

    #sort and average sample1_1a with sample1_2a etc
    grouped_df = df_labels.groupby(['Sample', 'Duplicate'])['OD_value'].mean().reset_index()
    #print("Grouped DF:\n",grouped_df)

    #dilution map and translate
    dilution_map = {
        "a": 20,
        "b": 80,
        "c": 320,
        "d": 1280}
    grouped_df['Dilution_values'] = grouped_df['Duplicate'].map(dilution_map)
    #print("Grouped DF2:\n",grouped_df)

    #get the log2 of all the dilutions
    dilution_logs = {
        20: math.log2(20),
        80: math.log2(80),
        320: math.log2(320),
        1280: math.log2(1280)}
    grouped_df['Log2'] = grouped_df['Dilution_values'].map(dilution_logs)

    #set up control averages
    controls_df = grouped_df[grouped_df['Sample'] == 'Cntrls']
        # positive control wells 
    positive = controls_df[controls_df['Duplicate'].str.startswith(('a', 'b'))]

        # negative control wells
    negative = controls_df[controls_df['Duplicate'].str.startswith(('c', 'd'))]
        #get averages of controls
    SP = positive['OD_value'].mean()
    upwards = negative['OD_value'].mean()
    print("\nSP: ", SP,"\nUpwards: ", upwards)


    #loop analyzing each sample on plate
    for sample_id, sample_group in grouped_df.groupby('Sample'):
        print(f"\nProcessing Sample {sample_id}")
        if sample_id == "Control" or sample_id == 'Cntrls':
            print("Skipping ",sample_id)
            break
        if pd.isna(sample_id):
            print("Found NaN sample at the end – skipping")
            break

        
        # Get the 4 OD values for a, b, c, d
        od_values = sample_group['OD_value'].values
        print(f"OD values: {od_values}")

        #analysis
        #fake test data
        dilutions = np.array([20, 80, 320, 1280])
        log2_dilutions = np.log2(dilutions)
        ydata = od_values
        #ydata = np.array([0.63, 1.17, 1.30, 1.21])

        #function
        def logistic_model(x, IP, slope, SP, upwards):
            #IP, slope, SP, upwards = params
            return np.exp((slope*(np.log2(x)-IP)))/(1+np.exp((slope*(np.log2(x)-IP))))*SP+upwards

        #model it
        from lmfit import Model
        Lmodel = Model(logistic_model)
        
        #unhardcode SP and upwards dolater

        positivecontrol = SP
        negativecontrol = upwards
        params = Lmodel.make_params(IP=1, slope=1, SP=positivecontrol, upwards=negativecontrol)


        #fix SP and Upwards
        params['SP'].set(vary=False)
        params['upwards'].set(vary=False)

        #fit model
        result = Lmodel.fit(ydata, params, x=dilutions)
        predicted = result.best_fit

        #set SP and Upwards
        print("\nSP: ", result.params['SP'].value,"\nUpwards: ", result.params['upwards'].value)

        print("IP: ", result.params['IP'].value)

        if(result.params['IP'].value < minimum):
            skipped.append([sample_id, result.params['IP'].value])

        if(result.params['IP'].value < minimum):
            storage_array.append([sample_id, 4.32])
        else:
            storage_array.append([sample_id, result.params['IP'].value])

df_storage = pd.DataFrame(storage_array)

df_skipped = pd.DataFrame(skipped)
pd.set_option('display.max_rows', None)
print(df_storage)
print("Skipped:")
print(df_skipped)




#end of stuff im using

#analysis

#function
def logistic_model(x, IP, slope, SP, upwards):
    #IP, slope, SP, upwards = params
    return np.exp((slope*(np.log2(x)-IP)))/(1+np.exp((slope*(np.log2(x)-IP))))*SP+upwards

#model it
from lmfit import Model
Lmodel = Model(logistic_model)
params = Lmodel.make_params(IP=1, slope=1, SP=1.268, upwards=0.14)

#fix SP and Upwards
params['SP'].set(vary=False)
params['upwards'].set(vary=False)

#fit model
result = Lmodel.fit(ydata, params, x=dilutions)
predicted = result.best_fit

#print("\nIP: ", result.params['IP'].value)
#print(result.fit_report())


#plot
import matplotlib.pyplot as plt
plt.scatter(dilutions, ydata, color='blue', label='Observed')
plt.scatter(dilutions, predicted, color='red', label='Fit')
plt.xscale('log', base=2)
plt.xlabel("Dilution (log2 scale)")
plt.ylabel("OD")
plt.legend()
plt.title("Excel-compatible logistic fit")
plt.grid(True)
#plt.show()

#print stuff with max columns

pd.set_option('display.max_columns', None)
print("Plate Data:")
print(plate_data)
print("Sample Layout:")
print(sample_layout)
