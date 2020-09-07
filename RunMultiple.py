import shutil
import glob
import os
import pandas as pd
import numpy as np
from pymc3 import SamplingError
from datetime import datetime
startTime = datetime.now()

print("Obtaining appropriate volcanoes")
GVP_volcanoes = pd.read_csv("Volcano_list.csv")
Whelley = pd.read_csv("Whelley_volcanoes.csv")
GVP_df = pd.read_csv("GVPDB2019.csv")
Volcanoes_all = pd.merge(GVP_volcanoes, Whelley, how='inner', on=['Volcano Number'])
Volcanoes_all_num = Volcanoes_all['Volcano Number']
np.savetxt('Volcano_numbers.txt', Volcanoes_all_num.values, fmt='%d')
count = 0

with open('Volcano_numbers.txt', 'r') as Volcano_numbers:
    print("Beginning Bayesian analysis...this is going to take a while.")
    for line in Volcano_numbers:
        print("Running analysis for volcano number:", line)
        f2 = open('Volc_num.txt', 'w')
        num = f2.writelines(line)
        f2.close()
        try:
            exec(open("FreqMagSEA.py").read())
        except (ValueError, KeyError, SamplingError, ZeroDivisionError):
            print("Error occurred with volcano:", line)
            with open('Volcano_not_run.txt', 'a+') as Volc_error:
                # Move read cursor to the start of file.
                Volc_error.seek(0)
                # If file is not empty then append '\n'
                #data = Volc_error.read(100)
                #if len(data) > 0:
                #    Volc_error.write("\n")
                # Append text at the end of file
                Volc_error.write(line)
            pass

        print("Analysis completed for volcano number:", line)
        count += 1
    Volcano_numbers.close()

dir = os.path.dirname(__file__)
filename = os.path.join(dir, 'Probabilities')

allFiles = glob.glob(filename + "/*.csv")
allFiles.sort()
with open('Probability_results.csv', 'wb') as outfile:
    for i, fname in enumerate(allFiles):
        with open(fname, 'rb') as infile:
            if i != 0:
                infile.readline()  # Throw away header on all but first file
            # Block copy rest of file from input to output without parsing
            shutil.copyfileobj(infile, outfile)
            print(fname + " has been imported.")

Duration = datetime.now() - startTime
print("Analysis completed in:", Duration)

