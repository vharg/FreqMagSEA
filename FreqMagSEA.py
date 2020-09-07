import pymc3 as pm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from scipy.stats import beta

with open("Volc_num.txt", 'r') as f:
  Volc_Num = f.readline().rstrip()
Volcano_Number = int(Volc_Num)
######
cores = 1
######


GVP_df = pd.read_csv("GVPDB2019.csv")
GVP_volcanoes = pd.read_csv("Volcano_list.csv")
GVP_df2 = GVP_df[
    ['Volcano Number', 'Volcano Name', 'Eruption Category', 'VEI', 'Start Year', 'Start Year Uncertainty']].copy()
GVP_df2a = pd.merge(GVP_df2, GVP_volcanoes, on=['Volcano Number'], how='left')
GVP_df3 = GVP_df2a[
    ['Volcano Number', 'Volcano Name_x', 'Eruption Category', 'Primary Volcano Type', 'VEI', 'Start Year',
     'Start Year Uncertainty', 'Region', 'Country', 'Subregion']].copy()

GVP_df3.loc[GVP_df3['Primary Volcano Type'] == 'Caldera', 'Volcano type'] = 'Caldera'
GVP_df3.loc[GVP_df3['Primary Volcano Type'] == 'Calderas', 'Volcano type'] = 'Caldera'
GVP_df3.loc[GVP_df3['Primary Volcano Type'] == 'Pyroclastic shield', 'Volcano type'] = 'Caldera'
GVP_df3.loc[GVP_df3['Primary Volcano Type'] == 'Caldera(s)', 'Volcano type'] = 'Caldera'

GVP_df3.loc[GVP_df3['Primary Volcano Type'] == 'Stratovolcano', 'Volcano type'] = 'Large cone'
GVP_df3.loc[GVP_df3['Primary Volcano Type'] == 'Complex', 'Volcano type'] = 'Large cone'
GVP_df3.loc[GVP_df3['Primary Volcano Type'] == 'Stratovolcano(es)', 'Volcano type'] = 'Large cone'
GVP_df3.loc[GVP_df3['Primary Volcano Type'] == 'Compound', 'Volcano type'] = 'Large cone'
GVP_df3.loc[GVP_df3['Primary Volcano Type'] == 'Stratovolcano?', 'Volcano type'] = 'Large cone'
GVP_df3.loc[GVP_df3['Primary Volcano Type'] == 'Complex(es)', 'Volcano type'] = 'Large cone'

GVP_df3.loc[GVP_df3['Primary Volcano Type'] == 'Shield', 'Volcano type'] = 'Shield'
GVP_df3.loc[GVP_df3['Primary Volcano Type'] == 'Shield(s)', 'Volcano type'] = 'Shield'

GVP_df3.loc[GVP_df3['Primary Volcano Type'] == 'Lava dome(s)', 'Volcano type'] = 'Lava dome'
GVP_df3.loc[GVP_df3['Primary Volcano Type'] == 'Lava dome', 'Volcano type'] = 'Lava dome'

GVP_df3.loc[GVP_df3['Primary Volcano Type'] == 'Pyroclastic cone(s)', 'Volcano type'] = 'Small cone'
GVP_df3.loc[GVP_df3['Primary Volcano Type'] == 'Pyroclastic cone', 'Volcano type'] = 'Small cone'
GVP_df3.loc[GVP_df3['Primary Volcano Type'] == 'Lava cone', 'Volcano type'] = 'Small cone'
GVP_df3.loc[GVP_df3['Primary Volcano Type'] == 'Explosion crater(s)', 'Volcano type'] = 'Small cone'
GVP_df3.loc[GVP_df3['Primary Volcano Type'] == 'Fissure vent(s)', 'Volcano type'] = 'Small cone'
GVP_df3.loc[GVP_df3['Primary Volcano Type'] == 'Volcanic field', 'Volcano type'] = 'Small cone'
GVP_df3.loc[GVP_df3['Primary Volcano Type'] == 'Maar(s)', 'Volcano type'] = 'Small cone'
GVP_df3.loc[GVP_df3['Primary Volcano Type'] == 'Cone(s)', 'Volcano type'] = 'Small cone'
GVP_df3.loc[GVP_df3['Primary Volcano Type'] == 'Fissure vent', 'Volcano type'] = 'Small cone'
GVP_df3.loc[GVP_df3['Primary Volcano Type'] == 'Crater rows', 'Volcano type'] = 'Small cone'
GVP_df3.loc[GVP_df3['Primary Volcano Type'] == 'Maar', 'Volcano type'] = 'Small cone'
GVP_df3.loc[GVP_df3['Primary Volcano Type'] == 'Tuff cone(s)', 'Volcano type'] = 'Small cone'
GVP_df3.loc[GVP_df3['Primary Volcano Type'] == 'Lava cone(s)', 'Volcano type'] = 'Small cone'
GVP_df3.loc[GVP_df3['Primary Volcano Type'] == 'Tuff cone', 'Volcano type'] = 'Small cone'
GVP_df3.loc[GVP_df3['Primary Volcano Type'] == 'Tuff ring(s)', 'Volcano type'] = 'Small cone'
GVP_df3.loc[GVP_df3['Primary Volcano Type'] == 'Lava cone(es)', 'Volcano type'] = 'Small cone'

# For volcanoes with no recorded eruptions
GVP_volcanoes.loc[GVP_volcanoes['Primary Volcano Type'] == 'Caldera', 'Volcano type'] = 'Caldera'
GVP_volcanoes.loc[GVP_volcanoes['Primary Volcano Type'] == 'Calderas', 'Volcano type'] = 'Caldera'
GVP_volcanoes.loc[GVP_volcanoes['Primary Volcano Type'] == 'Pyroclastic shield', 'Volcano type'] = 'Caldera'
GVP_volcanoes.loc[GVP_volcanoes['Primary Volcano Type'] == 'Caldera(s)', 'Volcano type'] = 'Caldera'

GVP_volcanoes.loc[GVP_volcanoes['Primary Volcano Type'] == 'Stratovolcano', 'Volcano type'] = 'Large cone'
GVP_volcanoes.loc[GVP_volcanoes['Primary Volcano Type'] == 'Complex', 'Volcano type'] = 'Large cone'
GVP_volcanoes.loc[GVP_volcanoes['Primary Volcano Type'] == 'Stratovolcano(es)', 'Volcano type'] = 'Large cone'
GVP_volcanoes.loc[GVP_volcanoes['Primary Volcano Type'] == 'Compound', 'Volcano type'] = 'Large cone'
GVP_volcanoes.loc[GVP_volcanoes['Primary Volcano Type'] == 'Stratovolcano?', 'Volcano type'] = 'Large cone'
GVP_volcanoes.loc[GVP_volcanoes['Primary Volcano Type'] == 'Complex(es)', 'Volcano type'] = 'Large cone'

GVP_volcanoes.loc[GVP_volcanoes['Primary Volcano Type'] == 'Shield', 'Volcano type'] = 'Shield'
GVP_volcanoes.loc[GVP_volcanoes['Primary Volcano Type'] == 'Shield(s)', 'Volcano type'] = 'Shield'

GVP_volcanoes.loc[GVP_volcanoes['Primary Volcano Type'] == 'Lava dome(s)', 'Volcano type'] = 'Lava dome'
GVP_volcanoes.loc[GVP_volcanoes['Primary Volcano Type'] == 'Lava dome', 'Volcano type'] = 'Lava dome'

GVP_volcanoes.loc[GVP_volcanoes['Primary Volcano Type'] == 'Pyroclastic cone(s)', 'Volcano type'] = 'Small cone'
GVP_volcanoes.loc[GVP_volcanoes['Primary Volcano Type'] == 'Pyroclastic cone', 'Volcano type'] = 'Small cone'
GVP_volcanoes.loc[GVP_volcanoes['Primary Volcano Type'] == 'Lava cone', 'Volcano type'] = 'Small cone'
GVP_volcanoes.loc[GVP_volcanoes['Primary Volcano Type'] == 'Explosion crater(s)', 'Volcano type'] = 'Small cone'
GVP_volcanoes.loc[GVP_volcanoes['Primary Volcano Type'] == 'Fissure vent(s)', 'Volcano type'] = 'Small cone'
GVP_volcanoes.loc[GVP_volcanoes['Primary Volcano Type'] == 'Volcanic field', 'Volcano type'] = 'Small cone'
GVP_volcanoes.loc[GVP_volcanoes['Primary Volcano Type'] == 'Maar(s)', 'Volcano type'] = 'Small cone'
GVP_volcanoes.loc[GVP_volcanoes['Primary Volcano Type'] == 'Cone(s)', 'Volcano type'] = 'Small cone'
GVP_volcanoes.loc[GVP_volcanoes['Primary Volcano Type'] == 'Fissure vent', 'Volcano type'] = 'Small cone'
GVP_volcanoes.loc[GVP_volcanoes['Primary Volcano Type'] == 'Crater rows', 'Volcano type'] = 'Small cone'
GVP_volcanoes.loc[GVP_volcanoes['Primary Volcano Type'] == 'Maar', 'Volcano type'] = 'Small cone'
GVP_volcanoes.loc[GVP_volcanoes['Primary Volcano Type'] == 'Tuff cone(s)', 'Volcano type'] = 'Small cone'
GVP_volcanoes.loc[GVP_volcanoes['Primary Volcano Type'] == 'Lava cone(s)', 'Volcano type'] = 'Small cone'
GVP_volcanoes.loc[GVP_volcanoes['Primary Volcano Type'] == 'Tuff cone', 'Volcano type'] = 'Small cone'
GVP_volcanoes.loc[GVP_volcanoes['Primary Volcano Type'] == 'Tuff ring(s)', 'Volcano type'] = 'Small cone'
GVP_volcanoes.loc[GVP_volcanoes['Primary Volcano Type'] == 'Lava cone(es)', 'Volcano type'] = 'Small cone'

Jenkins_VEI_all = GVP_df3[(GVP_df3['Eruption Category'] == 'Confirmed Eruption')].groupby('Volcano type').count()
Jenkins_VEI_Less_than_3 = GVP_df3[(GVP_df3.VEI <= 3) & (GVP_df3['Eruption Category'] == "Confirmed Eruption")].groupby(
    'Volcano type').count()
Jenkins_VEI_4 = GVP_df3[(GVP_df3.VEI == 4) & (GVP_df3['Eruption Category'] == "Confirmed Eruption")].groupby(
    'Volcano type').count()
Jenkins_VEI_5 = GVP_df3[(GVP_df3.VEI == 5) & (GVP_df3['Eruption Category'] == "Confirmed Eruption")].groupby(
    'Volcano type').count()
Jenkins_VEI_6 = GVP_df3[(GVP_df3.VEI == 6) & (GVP_df3['Eruption Category'] == "Confirmed Eruption")].groupby(
    'Volcano type').count()
Jenkins_VEI_7 = GVP_df3[(GVP_df3.VEI == 7) & (GVP_df3['Eruption Category'] == "Confirmed Eruption")].groupby(
    'Volcano type').count()

Jenkins_VEI_Less_than_3_df = pd.DataFrame(Jenkins_VEI_Less_than_3, columns=['Volcano Number'])
Jenkins_VEI_Less_than_3_df.columns = ['VEI <= 3']
Jenkins_VEI_Less_than_3_df_T = Jenkins_VEI_Less_than_3_df.T

Jenkins_VEI_4_df = pd.DataFrame(Jenkins_VEI_4, columns=['Volcano Number'])
Jenkins_VEI_4_df.columns = ['VEI 4']
Jenkins_VEI_4_df_T = Jenkins_VEI_4_df.T

Jenkins_VEI_5_df = pd.DataFrame(Jenkins_VEI_5, columns=['Volcano Number'])
Jenkins_VEI_5_df.columns = ['VEI 5']
Jenkins_VEI_5_df_T = Jenkins_VEI_5_df.T

Jenkins_VEI_6_df = pd.DataFrame(Jenkins_VEI_6, columns=['Volcano Number'])
Jenkins_VEI_6_df.columns = ['VEI 6']
Jenkins_VEI_6_df_T = Jenkins_VEI_6_df.T

Jenkins_VEI_7_df = pd.DataFrame(Jenkins_VEI_7, columns=['Volcano Number'])
Jenkins_VEI_7_df.columns = ['VEI 7']
Jenkins_VEI_7_df_T = Jenkins_VEI_7_df.T

Jenkins_VEI_frames = [Jenkins_VEI_Less_than_3_df_T, Jenkins_VEI_4_df_T, Jenkins_VEI_5_df_T, Jenkins_VEI_6_df_T,
                      Jenkins_VEI_7_df_T]
Jenkins_VEI_all_eruptions = pd.concat(Jenkins_VEI_frames, sort=True)

Jenkins_VEI_all_eruptions.fillna(0, inplace=True)

Jenkins_VEI_all_eruptions_T = Jenkins_VEI_all_eruptions.T

Jenkins_VEI_all_eruptions_T['Total'] = Jenkins_VEI_all_eruptions_T['VEI <= 3'] + Jenkins_VEI_all_eruptions_T['VEI 4'] + \
                                       Jenkins_VEI_all_eruptions_T['VEI 5'] + Jenkins_VEI_all_eruptions_T['VEI 6'] + \
                                       Jenkins_VEI_all_eruptions_T['VEI 7']

Jenkins_Total_eruptions_caldera = Jenkins_VEI_all_eruptions_T.iloc[0]['Total']
Jenkins_Total_eruptions_LC = Jenkins_VEI_all_eruptions_T.iloc[1]['Total']
Jenkins_Total_eruptions_LD = Jenkins_VEI_all_eruptions_T.iloc[2]['Total']
Jenkins_Total_eruptions_Shield = Jenkins_VEI_all_eruptions_T.iloc[3]['Total']
Jenkins_Total_eruptions_SC = Jenkins_VEI_all_eruptions_T.iloc[4]['Total']

# Caldera
Jenkins_VEI_3_eruptions_caldera = Jenkins_VEI_all_eruptions_T.iloc[0]['VEI <= 3']
Jenkins_VEI_4_eruptions_caldera = Jenkins_VEI_all_eruptions_T.iloc[0]['VEI 4']
Jenkins_VEI_5_eruptions_caldera = Jenkins_VEI_all_eruptions_T.iloc[0]['VEI 5']
Jenkins_VEI_6_eruptions_caldera = Jenkins_VEI_all_eruptions_T.iloc[0]['VEI 6']
Jenkins_VEI_7_eruptions_caldera = Jenkins_VEI_all_eruptions_T.iloc[0]['VEI 7']

# Large cone
Jenkins_VEI_3_eruptions_LC = Jenkins_VEI_all_eruptions_T.iloc[1]['VEI <= 3']
Jenkins_VEI_4_eruptions_LC = Jenkins_VEI_all_eruptions_T.iloc[1]['VEI 4']
Jenkins_VEI_5_eruptions_LC = Jenkins_VEI_all_eruptions_T.iloc[1]['VEI 5']
Jenkins_VEI_6_eruptions_LC = Jenkins_VEI_all_eruptions_T.iloc[1]['VEI 6']
Jenkins_VEI_7_eruptions_LC = Jenkins_VEI_all_eruptions_T.iloc[1]['VEI 7']

# Lava dome
Jenkins_VEI_3_eruptions_LD = Jenkins_VEI_all_eruptions_T.iloc[2]['VEI <= 3']
Jenkins_VEI_4_eruptions_LD = Jenkins_VEI_all_eruptions_T.iloc[2]['VEI 4']
Jenkins_VEI_5_eruptions_LD = Jenkins_VEI_all_eruptions_T.iloc[2]['VEI 5']
Jenkins_VEI_6_eruptions_LD = Jenkins_VEI_all_eruptions_T.iloc[2]['VEI 6']
Jenkins_VEI_7_eruptions_LD = Jenkins_VEI_all_eruptions_T.iloc[2]['VEI 7']

# Shield
Jenkins_VEI_3_eruptions_Shield = Jenkins_VEI_all_eruptions_T.iloc[3]['VEI <= 3']
Jenkins_VEI_4_eruptions_Shield = Jenkins_VEI_all_eruptions_T.iloc[3]['VEI 4']
Jenkins_VEI_5_eruptions_Shield = Jenkins_VEI_all_eruptions_T.iloc[3]['VEI 5']
Jenkins_VEI_6_eruptions_Shield = Jenkins_VEI_all_eruptions_T.iloc[3]['VEI 6']
Jenkins_VEI_7_eruptions_Shield = Jenkins_VEI_all_eruptions_T.iloc[3]['VEI 7']

# Small cone
Jenkins_VEI_3_eruptions_SC = Jenkins_VEI_all_eruptions_T.iloc[4]['VEI <= 3']
Jenkins_VEI_4_eruptions_SC = Jenkins_VEI_all_eruptions_T.iloc[4]['VEI 4']
Jenkins_VEI_5_eruptions_SC = Jenkins_VEI_all_eruptions_T.iloc[4]['VEI 5']
Jenkins_VEI_6_eruptions_SC = Jenkins_VEI_all_eruptions_T.iloc[4]['VEI 6']
Jenkins_VEI_7_eruptions_SC = Jenkins_VEI_all_eruptions_T.iloc[4]['VEI 7']

# determining conditional probabilities given an eruption of any magnitude
Jenkins_Cond_Prob_VEI_3_Caldera = Jenkins_VEI_3_eruptions_caldera / Jenkins_Total_eruptions_caldera
Jenkins_Cond_Prob_VEI_4_Caldera = Jenkins_VEI_4_eruptions_caldera / Jenkins_Total_eruptions_caldera
Jenkins_Cond_Prob_VEI_5_Caldera = Jenkins_VEI_5_eruptions_caldera / Jenkins_Total_eruptions_caldera
Jenkins_Cond_Prob_VEI_6_Caldera = Jenkins_VEI_6_eruptions_caldera / Jenkins_Total_eruptions_caldera
Jenkins_Cond_Prob_VEI_7_Caldera = Jenkins_VEI_7_eruptions_caldera / Jenkins_Total_eruptions_caldera

Jenkins_Cond_Prob_VEI_3_LC = Jenkins_VEI_3_eruptions_LC / Jenkins_Total_eruptions_LC
Jenkins_Cond_Prob_VEI_4_LC = Jenkins_VEI_4_eruptions_LC / Jenkins_Total_eruptions_LC
Jenkins_Cond_Prob_VEI_5_LC = Jenkins_VEI_5_eruptions_LC / Jenkins_Total_eruptions_LC
Jenkins_Cond_Prob_VEI_6_LC = Jenkins_VEI_6_eruptions_LC / Jenkins_Total_eruptions_LC
Jenkins_Cond_Prob_VEI_7_LC = Jenkins_VEI_7_eruptions_LC / Jenkins_Total_eruptions_LC

Jenkins_Cond_Prob_VEI_3_LD = Jenkins_VEI_3_eruptions_LD / Jenkins_Total_eruptions_LD
Jenkins_Cond_Prob_VEI_4_LD = Jenkins_VEI_4_eruptions_LD / Jenkins_Total_eruptions_LD
Jenkins_Cond_Prob_VEI_5_LD = Jenkins_VEI_5_eruptions_LD / Jenkins_Total_eruptions_LD
Jenkins_Cond_Prob_VEI_6_LD = Jenkins_VEI_6_eruptions_LD / Jenkins_Total_eruptions_LD
Jenkins_Cond_Prob_VEI_7_LD = Jenkins_VEI_7_eruptions_LD / Jenkins_Total_eruptions_LD

Jenkins_Cond_Prob_VEI_3_Shield = Jenkins_VEI_3_eruptions_Shield / Jenkins_Total_eruptions_Shield
Jenkins_Cond_Prob_VEI_4_Shield = Jenkins_VEI_4_eruptions_Shield / Jenkins_Total_eruptions_Shield
Jenkins_Cond_Prob_VEI_5_Shield = Jenkins_VEI_5_eruptions_Shield / Jenkins_Total_eruptions_Shield
Jenkins_Cond_Prob_VEI_6_Shield = Jenkins_VEI_6_eruptions_Shield / Jenkins_Total_eruptions_Shield
Jenkins_Cond_Prob_VEI_7_Shield = Jenkins_VEI_7_eruptions_Shield / Jenkins_Total_eruptions_Shield

Jenkins_Cond_Prob_VEI_3_SC = Jenkins_VEI_3_eruptions_SC / Jenkins_Total_eruptions_SC
Jenkins_Cond_Prob_VEI_4_SC = Jenkins_VEI_4_eruptions_SC / Jenkins_Total_eruptions_SC
Jenkins_Cond_Prob_VEI_5_SC = Jenkins_VEI_5_eruptions_SC / Jenkins_Total_eruptions_SC
Jenkins_Cond_Prob_VEI_6_SC = Jenkins_VEI_6_eruptions_SC / Jenkins_Total_eruptions_SC
Jenkins_Cond_Prob_VEI_7_SC = Jenkins_VEI_7_eruptions_SC / Jenkins_Total_eruptions_SC

Jenkins_Conditional_probabilities_all = pd.DataFrame(np.array([[Jenkins_Cond_Prob_VEI_3_Caldera,
                                                                Jenkins_Cond_Prob_VEI_3_LC, Jenkins_Cond_Prob_VEI_3_LD,
                                                                Jenkins_Cond_Prob_VEI_3_Shield,
                                                                Jenkins_Cond_Prob_VEI_3_SC],
                                                               [Jenkins_Cond_Prob_VEI_4_Caldera,
                                                                Jenkins_Cond_Prob_VEI_4_LC, Jenkins_Cond_Prob_VEI_4_LD,
                                                                Jenkins_Cond_Prob_VEI_4_Shield,
                                                                Jenkins_Cond_Prob_VEI_4_SC],
                                                               [Jenkins_Cond_Prob_VEI_5_Caldera,
                                                                Jenkins_Cond_Prob_VEI_5_LC, Jenkins_Cond_Prob_VEI_5_LD,
                                                                Jenkins_Cond_Prob_VEI_5_Shield,
                                                                Jenkins_Cond_Prob_VEI_5_SC],
                                                               [Jenkins_Cond_Prob_VEI_6_Caldera,
                                                                Jenkins_Cond_Prob_VEI_6_LC, Jenkins_Cond_Prob_VEI_6_LD,
                                                                Jenkins_Cond_Prob_VEI_6_Shield,
                                                                Jenkins_Cond_Prob_VEI_6_SC],
                                                               [Jenkins_Cond_Prob_VEI_7_Caldera,
                                                                Jenkins_Cond_Prob_VEI_7_LC, Jenkins_Cond_Prob_VEI_7_LD,
                                                                Jenkins_Cond_Prob_VEI_7_Shield,
                                                                Jenkins_Cond_Prob_VEI_7_SC]]),
                                                     columns=['Caldera', 'Large cone', 'Lava dome', 'Shield',
                                                              'Small cone'])

######

######
# Model 2: Based upon Whelley volcano classification within SEA

Whelley = pd.read_csv("Whelley_SEA.csv")
Whelley_GVP = pd.merge(Whelley, GVP_df2, on=['Volcano Number'], how='left')
Whelley_GVP_df = Whelley_GVP[['Volcano Number',
                              'Volcano Name_x',
                              'Eruption Category',
                              'Whelley classification',
                              'VEI',
                              'Start Year',
                              'Start Year Uncertainty']].copy()

######

######

Whelley_VEI_all = Whelley_GVP_df[(Whelley_GVP_df['Eruption Category'] == 'Confirmed Eruption')].groupby(
    'Whelley classification').count()
Whelley_VEI_Less_than_3 = Whelley_GVP_df[
    (Whelley_GVP_df.VEI <= 3) & (Whelley_GVP_df['Eruption Category'] == "Confirmed Eruption")].groupby(
    'Whelley classification').count()
Whelley_VEI_4 = Whelley_GVP_df[
    (Whelley_GVP_df.VEI == 4) & (Whelley_GVP_df['Eruption Category'] == "Confirmed Eruption")].groupby(
    'Whelley classification').count()
Whelley_VEI_5 = Whelley_GVP_df[
    (Whelley_GVP_df.VEI == 5) & (Whelley_GVP_df['Eruption Category'] == "Confirmed Eruption")].groupby(
    'Whelley classification').count()
Whelley_VEI_6 = Whelley_GVP_df[
    (Whelley_GVP_df.VEI == 6) & (Whelley_GVP_df['Eruption Category'] == "Confirmed Eruption")].groupby(
    'Whelley classification').count()
Whelley_VEI_7 = Whelley_GVP_df[
    (Whelley_GVP_df.VEI == 7) & (Whelley_GVP_df['Eruption Category'] == "Confirmed Eruption")].groupby(
    'Whelley classification').count()

Whelley_VEI_Less_than_3_df = pd.DataFrame(Whelley_VEI_Less_than_3, columns=['Volcano Number'])
Whelley_VEI_Less_than_3_df.columns = ['VEI <= 3']
Whelley_VEI_Less_than_3_df_T = Whelley_VEI_Less_than_3_df.T

Whelley_VEI_4_df = pd.DataFrame(Whelley_VEI_4, columns=['Volcano Number'])
Whelley_VEI_4_df.columns = ['VEI 4']
Whelley_VEI_4_df_T = Whelley_VEI_4_df.T

Whelley_VEI_5_df = pd.DataFrame(Whelley_VEI_5, columns=['Volcano Number'])
Whelley_VEI_5_df.columns = ['VEI 5']
Whelley_VEI_5_df_T = Whelley_VEI_5_df.T

Whelley_VEI_6_df = pd.DataFrame(Whelley_VEI_6, columns=['Volcano Number'])
Whelley_VEI_6_df.columns = ['VEI 6']
Whelley_VEI_6_df_T = Whelley_VEI_6_df.T

Whelley_VEI_7_df = pd.DataFrame(Whelley_VEI_7, columns=['Volcano Number'])
Whelley_VEI_7_df.columns = ['VEI 7']
Whelley_VEI_7_df_T = Whelley_VEI_7_df.T

Whelley_VEI_frames = [Whelley_VEI_Less_than_3_df_T,
                      Whelley_VEI_4_df_T,
                      Whelley_VEI_5_df_T,
                      Whelley_VEI_6_df_T,
                      Whelley_VEI_7_df_T]
Whelley_VEI_all_eruptions = pd.concat(Whelley_VEI_frames, sort=True)

Whelley_VEI_all_eruptions.fillna(0, inplace=True)

Whelley_VEI_all_eruptions_T = Whelley_VEI_all_eruptions.T

Whelley_VEI_all_eruptions_T['Total'] = Whelley_VEI_all_eruptions_T['VEI <= 3'] + Whelley_VEI_all_eruptions_T['VEI 4'] + \
                                       Whelley_VEI_all_eruptions_T['VEI 5'] + Whelley_VEI_all_eruptions_T['VEI 6'] + \
                                       Whelley_VEI_all_eruptions_T['VEI 7']

Whelley_Total_eruptions_SC = Whelley_VEI_all_eruptions_T.iloc[0]['Total']
Whelley_Total_eruptions_LC = Whelley_VEI_all_eruptions_T.iloc[1]['Total']
Whelley_Total_eruptions_OS = Whelley_VEI_all_eruptions_T.iloc[2]['Total']
Whelley_Total_eruptions_SPS = Whelley_VEI_all_eruptions_T.iloc[3]['Total']
Whelley_Total_eruptions_WPS = Whelley_VEI_all_eruptions_T.iloc[4]['Total']

# Distributed cones and fields
Whelley_VEI_3_eruptions_SC = Whelley_VEI_all_eruptions_T.iloc[0]['VEI <= 3']
Whelley_VEI_4_eruptions_SC = Whelley_VEI_all_eruptions_T.iloc[0]['VEI 4']
Whelley_VEI_5_eruptions_SC = Whelley_VEI_all_eruptions_T.iloc[0]['VEI 5']
Whelley_VEI_6_eruptions_SC = Whelley_VEI_all_eruptions_T.iloc[0]['VEI 6']
Whelley_VEI_7_eruptions_SC = Whelley_VEI_all_eruptions_T.iloc[0]['VEI 7']

# Large caldera
Whelley_VEI_3_eruptions_LC = Whelley_VEI_all_eruptions_T.iloc[1]['VEI <= 3']
Whelley_VEI_4_eruptions_LC = Whelley_VEI_all_eruptions_T.iloc[1]['VEI 4']
Whelley_VEI_5_eruptions_LC = Whelley_VEI_all_eruptions_T.iloc[1]['VEI 5']
Whelley_VEI_6_eruptions_LC = Whelley_VEI_all_eruptions_T.iloc[1]['VEI 6']
Whelley_VEI_7_eruptions_LC = Whelley_VEI_all_eruptions_T.iloc[1]['VEI 7']

# Open-vent stratocone
Whelley_VEI_3_eruptions_OS = Whelley_VEI_all_eruptions_T.iloc[2]['VEI <= 3']
Whelley_VEI_4_eruptions_OS = Whelley_VEI_all_eruptions_T.iloc[2]['VEI 4']
Whelley_VEI_5_eruptions_OS = Whelley_VEI_all_eruptions_T.iloc[2]['VEI 5']
Whelley_VEI_6_eruptions_OS = Whelley_VEI_all_eruptions_T.iloc[2]['VEI 6']
Whelley_VEI_7_eruptions_OS = Whelley_VEI_all_eruptions_T.iloc[2]['VEI 7']

# Semi-plugged stratocone
Whelley_VEI_3_eruptions_SPS = Whelley_VEI_all_eruptions_T.iloc[3]['VEI <= 3']
Whelley_VEI_4_eruptions_SPS = Whelley_VEI_all_eruptions_T.iloc[3]['VEI 4']
Whelley_VEI_5_eruptions_SPS = Whelley_VEI_all_eruptions_T.iloc[3]['VEI 5']
Whelley_VEI_6_eruptions_SPS = Whelley_VEI_all_eruptions_T.iloc[3]['VEI 6']
Whelley_VEI_7_eruptions_SPS = Whelley_VEI_all_eruptions_T.iloc[3]['VEI 7']

# Well-plugged stratocone
Whelley_VEI_3_eruptions_WPS = Whelley_VEI_all_eruptions_T.iloc[4]['VEI <= 3']
Whelley_VEI_4_eruptions_WPS = Whelley_VEI_all_eruptions_T.iloc[4]['VEI 4']
Whelley_VEI_5_eruptions_WPS = Whelley_VEI_all_eruptions_T.iloc[4]['VEI 5']
Whelley_VEI_6_eruptions_WPS = Whelley_VEI_all_eruptions_T.iloc[4]['VEI 6']
Whelley_VEI_7_eruptions_WPS = Whelley_VEI_all_eruptions_T.iloc[4]['VEI 7']

# determining conditional probabilities given an eruption of any magnitude
Whelley_Cond_Prob_VEI_3_SC = Whelley_VEI_3_eruptions_SC / Whelley_Total_eruptions_SC
Whelley_Cond_Prob_VEI_4_SC = Whelley_VEI_4_eruptions_SC / Whelley_Total_eruptions_SC
Whelley_Cond_Prob_VEI_5_SC = Whelley_VEI_5_eruptions_SC / Whelley_Total_eruptions_SC
Whelley_Cond_Prob_VEI_6_SC = Whelley_VEI_6_eruptions_SC / Whelley_Total_eruptions_SC
Whelley_Cond_Prob_VEI_7_SC = Whelley_VEI_7_eruptions_SC / Whelley_Total_eruptions_SC

Whelley_Cond_Prob_VEI_3_LC = Whelley_VEI_3_eruptions_LC / Whelley_Total_eruptions_LC
Whelley_Cond_Prob_VEI_4_LC = Whelley_VEI_4_eruptions_LC / Whelley_Total_eruptions_LC
Whelley_Cond_Prob_VEI_5_LC = Whelley_VEI_5_eruptions_LC / Whelley_Total_eruptions_LC
Whelley_Cond_Prob_VEI_6_LC = Whelley_VEI_6_eruptions_LC / Whelley_Total_eruptions_LC
Whelley_Cond_Prob_VEI_7_LC = Whelley_VEI_7_eruptions_LC / Whelley_Total_eruptions_LC

Whelley_Cond_Prob_VEI_3_OS = Whelley_VEI_3_eruptions_OS / Whelley_Total_eruptions_OS
Whelley_Cond_Prob_VEI_4_OS = Whelley_VEI_4_eruptions_OS / Whelley_Total_eruptions_OS
Whelley_Cond_Prob_VEI_5_OS = Whelley_VEI_5_eruptions_OS / Whelley_Total_eruptions_OS
Whelley_Cond_Prob_VEI_6_OS = Whelley_VEI_6_eruptions_OS / Whelley_Total_eruptions_OS
Whelley_Cond_Prob_VEI_7_OS = Whelley_VEI_7_eruptions_OS / Whelley_Total_eruptions_OS

Whelley_Cond_Prob_VEI_3_SPS = Whelley_VEI_3_eruptions_SPS / Whelley_Total_eruptions_SPS
Whelley_Cond_Prob_VEI_4_SPS = Whelley_VEI_4_eruptions_SPS / Whelley_Total_eruptions_SPS
Whelley_Cond_Prob_VEI_5_SPS = Whelley_VEI_5_eruptions_SPS / Whelley_Total_eruptions_SPS
Whelley_Cond_Prob_VEI_6_SPS = Whelley_VEI_6_eruptions_SPS / Whelley_Total_eruptions_SPS
Whelley_Cond_Prob_VEI_7_SPS = Whelley_VEI_7_eruptions_SPS / Whelley_Total_eruptions_SPS

Whelley_Cond_Prob_VEI_3_WPS = Whelley_VEI_3_eruptions_WPS / Whelley_Total_eruptions_WPS
Whelley_Cond_Prob_VEI_4_WPS = Whelley_VEI_4_eruptions_WPS / Whelley_Total_eruptions_WPS
Whelley_Cond_Prob_VEI_5_WPS = Whelley_VEI_5_eruptions_WPS / Whelley_Total_eruptions_WPS
Whelley_Cond_Prob_VEI_6_WPS = Whelley_VEI_6_eruptions_WPS / Whelley_Total_eruptions_WPS
Whelley_Cond_Prob_VEI_7_WPS = Whelley_VEI_7_eruptions_WPS / Whelley_Total_eruptions_WPS

Whelley_Conditional_probabilities_all = pd.DataFrame(np.array([[Whelley_Cond_Prob_VEI_3_SC, Whelley_Cond_Prob_VEI_3_LC,
                                                                Whelley_Cond_Prob_VEI_3_OS, Whelley_Cond_Prob_VEI_3_SPS,
                                                                Whelley_Cond_Prob_VEI_3_WPS],
                                                               [Whelley_Cond_Prob_VEI_4_SC, Whelley_Cond_Prob_VEI_4_LC,
                                                                Whelley_Cond_Prob_VEI_4_OS, Whelley_Cond_Prob_VEI_4_SPS,
                                                                Whelley_Cond_Prob_VEI_4_WPS],
                                                               [Whelley_Cond_Prob_VEI_5_SC, Whelley_Cond_Prob_VEI_5_LC,
                                                                Whelley_Cond_Prob_VEI_5_OS, Whelley_Cond_Prob_VEI_5_SPS,
                                                                Whelley_Cond_Prob_VEI_5_WPS],
                                                               [Whelley_Cond_Prob_VEI_6_SC, Whelley_Cond_Prob_VEI_6_LC,
                                                                Whelley_Cond_Prob_VEI_6_OS, Whelley_Cond_Prob_VEI_6_SPS,
                                                                Whelley_Cond_Prob_VEI_6_WPS],
                                                               [Whelley_Cond_Prob_VEI_7_SC, Whelley_Cond_Prob_VEI_7_LC,
                                                                Whelley_Cond_Prob_VEI_7_OS, Whelley_Cond_Prob_VEI_7_SPS,
                                                                Whelley_Cond_Prob_VEI_7_WPS]]),
                                                     columns=['Distributed cones and fields',
                                                              'Large caldera',
                                                              'Open-vent stratocone',
                                                              'Semi-plugged stratocone',
                                                              'Well-plugged stratocone'])

#####
# Volcano specifics
######
Change_points = pd.read_csv("Change_points.csv")
Change_points.set_index("Region", inplace=True)
Region = GVP_volcanoes[(GVP_volcanoes['Volcano Number'] == Volcano_Number)].iloc[0][6]
Complete_record_small = Change_points.loc[Region, 'All eruptions']
Complete_record_large = Change_points.loc[Region, 'Large eruptions']
Complete_record = 2019 - Complete_record_small
Volcano_name = GVP_volcanoes[(GVP_volcanoes['Volcano Number'] == Volcano_Number)].iloc[0][1]

######
# Setting analogue annual frequency of eruption based on GVP classification
######

# Caldera
Eruptions_caldera = GVP_df3[(GVP_df3['Eruption Category'] == 'Confirmed Eruption')
                            & (GVP_df3['Volcano type'] == 'Caldera')
                            & (GVP_df3['Start Year'] > 1960)].groupby('Volcano Number').count()
Number_of_eruptions_caldera = Eruptions_caldera['Volcano Name_x']
Rate_of_eruptions_caldera = Number_of_eruptions_caldera / 1960
Average_rate_caldera = Rate_of_eruptions_caldera.mean()
STD_rate_caldera = Rate_of_eruptions_caldera.std()

# Large cone
Eruptions_LC = GVP_df3[(GVP_df3['Eruption Category'] == 'Confirmed Eruption')
                       & (GVP_df3['Volcano type'] == 'Large cone')
                       & (GVP_df3['Start Year'] > 1960)].groupby('Volcano Number').count()
Number_of_eruptions_LC = Eruptions_LC['Volcano Name_x']
Rate_of_eruptions_LC = Number_of_eruptions_LC / 1960
Average_rate_LC = Rate_of_eruptions_LC.mean()
STD_rate_LC = Rate_of_eruptions_LC.std()

# Lava dome
Eruptions_LD = GVP_df3[(GVP_df3['Eruption Category'] == 'Confirmed Eruption')
                       & (GVP_df3['Volcano type'] == 'Lava dome')
                       & (GVP_df3['Start Year'] > 1960)].groupby('Volcano Number').count()
Number_of_eruptions_LD = Eruptions_LD['Volcano Name_x']
Rate_of_eruptions_LD = Number_of_eruptions_LD / 1960
Average_rate_LD = Rate_of_eruptions_LD.mean()
STD_rate_LD = Rate_of_eruptions_LD.std()

# Shield
Eruptions_Shield = GVP_df3[(GVP_df3['Eruption Category'] == 'Confirmed Eruption')
                           & (GVP_df3['Volcano type'] == 'Shield')
                           & (GVP_df3['Start Year'] > 1960)].groupby('Volcano Number').count()
Number_of_eruptions_Shield = Eruptions_Shield['Volcano Name_x']
Rate_of_eruptions_Shield = Number_of_eruptions_Shield / 1960
Average_rate_Shield = Rate_of_eruptions_Shield.mean()
STD_rate_Shield = Rate_of_eruptions_Shield.std()

# Small cone
Eruptions_SC = GVP_df3[(GVP_df3['Eruption Category'] == 'Confirmed Eruption')
                       & (GVP_df3['Volcano type'] == 'Small cone')
                       & (GVP_df3['Start Year'] > 1960)].groupby('Volcano Number').count()
Number_of_eruptions_SC = Eruptions_SC['Volcano Name_x']
Rate_of_eruptions_SC = Number_of_eruptions_SC / 1960
Average_rate_SC = Rate_of_eruptions_SC.mean()
STD_rate_SC = Rate_of_eruptions_SC.std()

######

######
# Setting analogue annual frequency of eruption based on Whelley classification

# Distributed cones and fields
Whelley_Eruptions_SC = Whelley_GVP_df[(Whelley_GVP_df['Eruption Category'] == 'Confirmed Eruption')
                                      & (Whelley_GVP_df['Whelley classification'] == 'Distributed cones and fields')
                                      & (Whelley_GVP_df['Start Year'] > 1960)].groupby('Volcano Number').count()
Whelley_Number_of_eruptions_SC = Whelley_Eruptions_SC['Volcano Name_x']
Whelley_Rate_of_eruptions_SC = Whelley_Number_of_eruptions_SC / 1960
Whelley_Average_rate_SC = Whelley_Rate_of_eruptions_SC.mean()
Whelley_STD_rate_SC = Whelley_Rate_of_eruptions_SC.std()

# Large caldera
Whelley_Eruptions_LC = Whelley_GVP_df[(Whelley_GVP_df['Eruption Category'] == 'Confirmed Eruption')
                                      & (Whelley_GVP_df['Whelley classification'] == 'Large caldera')
                                      & (Whelley_GVP_df['Start Year'] > 1960)].groupby('Volcano Number').count()
Whelley_Number_of_eruptions_LC = Whelley_Eruptions_LC['Volcano Name_x']
Whelley_Rate_of_eruptions_LC = Whelley_Number_of_eruptions_LC / 1960
Whelley_Average_rate_LC = Whelley_Rate_of_eruptions_LC.mean()
Whelley_STD_rate_LC = Whelley_Rate_of_eruptions_LC.std()

# Open-vent stratocone
Whelley_Eruptions_OS = Whelley_GVP_df[(Whelley_GVP_df['Eruption Category'] == 'Confirmed Eruption')
                                      & (Whelley_GVP_df['Whelley classification'] == 'Open-vent stratocone')
                                      & (Whelley_GVP_df['Start Year'] > 1960)].groupby('Volcano Number').count()
Whelley_Number_of_eruptions_OS = Whelley_Eruptions_OS['Volcano Name_x']
Whelley_Rate_of_eruptions_OS = Whelley_Number_of_eruptions_OS / 1960
Whelley_Average_rate_OS = Whelley_Rate_of_eruptions_OS.mean()
Whelley_STD_rate_OS = Whelley_Rate_of_eruptions_OS.std()

# Semi-plugged stratocone
Whelley_Eruptions_SPS = Whelley_GVP_df[(Whelley_GVP_df['Eruption Category'] == 'Confirmed Eruption')
                                       & (Whelley_GVP_df['Whelley classification'] == 'Semi-plugged stratocone')
                                       & (Whelley_GVP_df['Start Year'] > 1960)].groupby('Volcano Number').count()
Whelley_Number_of_eruptions_SPS = Whelley_Eruptions_SPS['Volcano Name_x']
Whelley_Rate_of_eruptions_SPS = Whelley_Number_of_eruptions_SPS / 1960
Whelley_Average_rate_SPS = Whelley_Rate_of_eruptions_SPS.mean()
Whelley_STD_rate_SPS = Whelley_Rate_of_eruptions_SPS.std()

# Well-plugged stratocone
Whelley_Eruptions_WPS = Whelley_GVP_df[(Whelley_GVP_df['Eruption Category'] == 'Confirmed Eruption')
                                       & (Whelley_GVP_df['Whelley classification'] == 'Well-plugged stratocone')
                                       & (Whelley_GVP_df['Start Year'] > 1960)].groupby('Volcano Number').count()
Whelley_Number_of_eruptions_WPS = Whelley_Eruptions_WPS['Volcano Name_x']
Whelley_Rate_of_eruptions_WPS = Whelley_Number_of_eruptions_WPS / 1960
Whelley_Average_rate_WPS = Whelley_Rate_of_eruptions_WPS.mean()
Whelley_STD_rate_WPS = Whelley_Rate_of_eruptions_WPS.std()

Change_points = pd.read_csv("Change_points.csv")
Change_points.set_index("Region", inplace=True)
######

######
# Observed rate of eruptions for selected volcano based on Jenkins approach

Small_eruptions = GVP_df3[(GVP_df3['Volcano Number'] == Volcano_Number)
                          & (GVP_df3['Eruption Category'] == 'Confirmed Eruption')
                          & (GVP_df3['VEI'] <= 3)
                          & (GVP_df3['Start Year'] > Complete_record_small)]
Number_of_small_eruptions = len(Small_eruptions)
Rate_of_small_eruptions = Number_of_small_eruptions / (2019 - Complete_record_small)

Large_eruptions = GVP_df3[(GVP_df3['Volcano Number'] == Volcano_Number)
                          & (GVP_df3['Eruption Category'] == 'Confirmed Eruption')
                          & (GVP_df3['VEI'] > 3)
                          & (GVP_df3['Start Year'] > Complete_record_large)]
Number_of_large_eruptions = len(Large_eruptions)
Rate_of_large_eruptions = Number_of_large_eruptions / (2019 - Complete_record_large)

Unknown_eruptions = GVP_df3[(GVP_df3['Volcano Number'] == Volcano_Number)
                            & (GVP_df3['Eruption Category'] == 'Confirmed Eruption')
                            & (GVP_df3.VEI.isnull())
                            & (GVP_df3['Start Year'] > Complete_record_large)]
Number_of_unknown_eruptions = len(Unknown_eruptions)
Rate_of_unknown_eruptions = Number_of_unknown_eruptions / (2019 - Complete_record_large)

Total_eruptions_per_year = Rate_of_unknown_eruptions + Rate_of_small_eruptions + Rate_of_large_eruptions
Total_eruptions_per_10000_year = Total_eruptions_per_year * 10000
Number_of_confirmed_VEI_eruptions = Number_of_small_eruptions + Number_of_large_eruptions

######

######
# Observed frequency-magnitude for selected volcano
VEI_3_or_less_eruptions = GVP_df3[(GVP_df3['Volcano Number'] == Volcano_Number) &
                                  (GVP_df3['Eruption Category'] == 'Confirmed Eruption') &
                                  (GVP_df3['Start Year'] > Complete_record_small) &
                                  (GVP_df3['VEI'] <= 3)]
Number_VEI_3_or_less_eruptions = len(VEI_3_or_less_eruptions)

VEI_4_eruptions = GVP_df3[(GVP_df3['Volcano Number'] == Volcano_Number) &
                          (GVP_df3['Eruption Category'] == 'Confirmed Eruption') &
                          (GVP_df3['Start Year'] > Complete_record_large) &
                          (GVP_df3['VEI'] == 4)]
Number_VEI_4_eruptions = len(VEI_4_eruptions)

VEI_5_eruptions = GVP_df3[(GVP_df3['Volcano Number'] == Volcano_Number) &
                          (GVP_df3['Eruption Category'] == 'Confirmed Eruption') &
                          (GVP_df3['Start Year'] > Complete_record_large) &
                          (GVP_df3['VEI'] == 5)]
Number_VEI_5_eruptions = len(VEI_5_eruptions)

VEI_6_eruptions = GVP_df3[(GVP_df3['Volcano Number'] == Volcano_Number) &
                          (GVP_df3['Eruption Category'] == 'Confirmed Eruption') &
                          (GVP_df3['Start Year'] > Complete_record_large) &
                          (GVP_df3['VEI'] == 6)]
Number_VEI_6_eruptions = len(VEI_6_eruptions)

VEI_7_eruptions = GVP_df3[(GVP_df3['Volcano Number'] == Volcano_Number) &
                          (GVP_df3['Eruption Category'] == 'Confirmed Eruption') &
                          (GVP_df3['Start Year'] > Complete_record_large) &
                          (GVP_df3['VEI'] == 7)]
Number_VEI_7_eruptions = len(VEI_7_eruptions)

#####
######
# Setting analogues frequency-magnitudes for selected volcano
try:
    Volcano_type_Jenkins = GVP_df3[(GVP_df3['Volcano Number'] == Volcano_Number)].values[0]
    Volcano_type_Jenkins = Volcano_type_Jenkins[10]
except IndexError:
    Volcano_type_Jenkins = GVP_volcanoes[(GVP_volcanoes['Volcano Number'] == Volcano_Number)].values[0]
    Volcano_type_Jenkins = Volcano_type_Jenkins[13]
Jenkins_Frequency_VEI_3 = Jenkins_Conditional_probabilities_all[Volcano_type_Jenkins][0]
Jenkins_Frequency_VEI_4 = Jenkins_Conditional_probabilities_all[Volcano_type_Jenkins][1]
Jenkins_Frequency_VEI_5 = Jenkins_Conditional_probabilities_all[Volcano_type_Jenkins][2]
Jenkins_Frequency_VEI_6 = Jenkins_Conditional_probabilities_all[Volcano_type_Jenkins][3]
Jenkins_Frequency_VEI_7 = Jenkins_Conditional_probabilities_all[Volcano_type_Jenkins][4]

Volcano_type_Whelley = Whelley_GVP_df[(Whelley_GVP_df['Volcano Number'] == Volcano_Number)].values[0]
Volcano_type_Whelley = Volcano_type_Whelley[3]
Whelley_Frequency_VEI_3 = Whelley_Conditional_probabilities_all[Volcano_type_Whelley][0]
Whelley_Frequency_VEI_4 = Whelley_Conditional_probabilities_all[Volcano_type_Whelley][1]
Whelley_Frequency_VEI_5 = Whelley_Conditional_probabilities_all[Volcano_type_Whelley][2]
Whelley_Frequency_VEI_6 = Whelley_Conditional_probabilities_all[Volcano_type_Whelley][3]
Whelley_Frequency_VEI_7 = Whelley_Conditional_probabilities_all[Volcano_type_Whelley][4]

# Prior for rate of eruptions

# Model 1 Jenkins

# Determining appropriate mean rate of eruptions
if (Volcano_type_Jenkins == "Caldera"):
    Average_eruptions_per_year = Average_rate_caldera
elif (Volcano_type_Jenkins == "Large cone"):
    Average_eruptions_per_year = Average_rate_LC
elif (Volcano_type_Jenkins == "Lava dome"):
    Average_eruptions_per_year = Average_rate_LD
elif (Volcano_type_Jenkins == "Shield"):
    Average_eruptions_per_year = Average_rate_Shield
elif (Volcano_type_Jenkins == "Small cone"):
    Average_eruptions_per_year = Average_rate_SC
else:
    print("ERROR")

# Determining appropriate standard deviation
if Volcano_type_Jenkins == "Caldera":
    STD_eruptions_per_year = STD_rate_caldera
elif Volcano_type_Jenkins == "Large cone":
    STD_eruptions_per_year = STD_rate_LC
elif Volcano_type_Jenkins == "Lava dome":
    STD_eruptions_per_year = STD_rate_LD
elif Volcano_type_Jenkins == "Shield":
    STD_eruptions_per_year = STD_rate_Shield
elif Volcano_type_Jenkins == "Small cone":
    STD_eruptions_per_year = STD_rate_SC
else:
    print("ERROR")

# Model 2 Whelley
# Determining appropriate mean rate of eruptions
if Volcano_type_Whelley == "Distributed cones and fields":
    Whelley_Average_eruptions_per_year = Whelley_Average_rate_SC
elif Volcano_type_Whelley == "Large caldera":
    Whelley_Average_eruptions_per_year = Whelley_Average_rate_LC
elif Volcano_type_Whelley == "Open-vent stratocone":
    Whelley_Average_eruptions_per_year = Whelley_Average_rate_OS
elif Volcano_type_Whelley == "Semi-plugged stratocone":
    Whelley_Average_eruptions_per_year = Whelley_Average_rate_SPS
elif Volcano_type_Whelley == "Well-plugged stratocone":
    Whelley_Average_eruptions_per_year = Whelley_Average_rate_WPS
else:
    print("ERROR")

# Determining appropriate standard deviation of eruptions
if Volcano_type_Whelley == "Distributed cones and fields":
    Whelley_STD_eruptions_per_year = Whelley_STD_rate_SC
elif Volcano_type_Whelley == "Large caldera":
    Whelley_STD_eruptions_per_year = Whelley_STD_rate_LC
elif Volcano_type_Whelley == "Open-vent stratocone":
    Whelley_STD_eruptions_per_year = Whelley_STD_rate_OS
elif Volcano_type_Whelley == "Semi-plugged stratocone":
    Whelley_STD_eruptions_per_year = Whelley_STD_rate_SPS
elif Volcano_type_Whelley == "Well-plugged stratocone":
    Whelley_STD_eruptions_per_year = Whelley_STD_rate_WPS
else:
    print("ERROR")

######

######
# Assessing frequency of eruption (any VEI/magnitude) - Bayesian model
# Jenkins Model

y = np.array([Total_eruptions_per_10000_year])
n = np.array([10000])
N = len(n)

Jenkins_1 = np.array([Jenkins_Frequency_VEI_3*10000,
                      Jenkins_Frequency_VEI_4*10000,
                      Jenkins_Frequency_VEI_5*10000,
                      Jenkins_Frequency_VEI_6*10000,
                      Jenkins_Frequency_VEI_7*10000])

Jenkins_2 = np.array([Jenkins_Frequency_VEI_3*10000,
                      Jenkins_Frequency_VEI_4*10000,
                      Jenkins_Frequency_VEI_5*10000,
                      Jenkins_Frequency_VEI_5*10000*0.1,
                      Jenkins_Frequency_VEI_5*10000*0.01])

Jenkins_3 = np.array([Jenkins_Frequency_VEI_3*10000,
                      Jenkins_Frequency_VEI_4*10000,
                      Jenkins_Frequency_VEI_4*10000*0.1,
                      Jenkins_Frequency_VEI_4*10000*0.01,
                      Jenkins_Frequency_VEI_4*10000*0.001])

if Jenkins_Frequency_VEI_7 > 0:
    Jenkins_c = Jenkins_1
elif Jenkins_Frequency_VEI_5 > 0 and Jenkins_Frequency_VEI_7 == 0:
    Jenkins_c = Jenkins_2
elif Jenkins_Frequency_VEI_5 == 0 and Jenkins_Frequency_VEI_7 == 0:
    Jenkins_c = Jenkins_3

# Jenkins_c = Jenkins_1
Jenkins_shape = Jenkins_c.size

Std_eruptions_per_year = STD_eruptions_per_year
alpha_analogue_jenkins = Average_eruptions_per_year ** 2 * (
            (1 - Average_eruptions_per_year) / Std_eruptions_per_year ** 2 - 1 / Average_eruptions_per_year)
beta_analogue_jenkins = alpha_analogue_jenkins * (1 / Average_eruptions_per_year - 1)

Whelley_Std_eruptions_per_year = Whelley_STD_eruptions_per_year
alpha_analogue_Whelley = Whelley_Average_eruptions_per_year ** 2 * ((1 - Whelley_Average_eruptions_per_year) / Whelley_Std_eruptions_per_year ** 2 - 1 / Whelley_Average_eruptions_per_year)
beta_analogue_Whelley = alpha_analogue_Whelley * (1 / Whelley_Average_eruptions_per_year - 1)

Freq_mag_analogue = stats.dirichlet.rvs(Jenkins_c, size=10000, random_state=None)

Freq_mag_samples = pd.DataFrame({'VEI3': Freq_mag_analogue[:, 0],
                        'VEI4': Freq_mag_analogue[:, 1],
                        'VEI5': Freq_mag_analogue[:, 2],
                        'VEI6': Freq_mag_analogue[:, 3],
                        'VEI7': Freq_mag_analogue[:, 4]})

Jenkins_freq_mag_means = Freq_mag_samples.mean()
Jenkins_Frequency_VEI_3 = Jenkins_freq_mag_means[0]
Jenkins_Frequency_VEI_4 = Jenkins_freq_mag_means[1]
Jenkins_Frequency_VEI_5 = Jenkins_freq_mag_means[2]
Jenkins_Frequency_VEI_6 = Jenkins_freq_mag_means[3]
Jenkins_Frequency_VEI_7 = Jenkins_freq_mag_means[4]

Jenkins_freq_mag_std = Freq_mag_samples.std()
Jenkins_Frequency_VEI_3_std = Jenkins_freq_mag_std[0]
Jenkins_Frequency_VEI_4_std = Jenkins_freq_mag_std[1]
Jenkins_Frequency_VEI_5_std = Jenkins_freq_mag_std[2]
Jenkins_Frequency_VEI_6_std = Jenkins_freq_mag_std[3]
Jenkins_Frequency_VEI_7_std = Jenkins_freq_mag_std[4]

#####
if Number_of_confirmed_VEI_eruptions == 0:
    Method = "Analogue"
    a_ave = alpha_analogue_jenkins #Average_eruptions_per_year ** 2 * ((1 - Average_eruptions_per_year) / STD_eruptions_per_year ** 2 - 1 / Average_eruptions_per_year)
    b_ave = beta_analogue_jenkins #a_ave * (1 / Average_eruptions_per_year - 1)

    # Combining probability of eruption and relative probability of different VEI
    AP_VEI_3 = []
    Number = 10000
    lower, upper = 0, 1
    for i in range(Number):
        """Monte Carlo simulation"""
        AP_VEI3 = (np.random.beta(alpha_analogue_jenkins, beta_analogue_jenkins) * stats.truncnorm.rvs((lower - Jenkins_Frequency_VEI_3) / Jenkins_Frequency_VEI_3_std,
                                                                      (upper - Jenkins_Frequency_VEI_3) / Jenkins_Frequency_VEI_3_std, loc=Jenkins_Frequency_VEI_3,
                                                                      scale=Jenkins_Frequency_VEI_3_std))
        AP_VEI_3 = AP_VEI_3 + [AP_VEI3]

    AP_VEI_4 = []
    Number = 10000
    lower, upper = 0, 1
    for i in range(Number):
        """Monte Carlo simulation"""
        AP_VEI4 = (np.random.beta(a_ave, b_ave) * stats.truncnorm.rvs((lower - Jenkins_Frequency_VEI_4) / Jenkins_Frequency_VEI_4_std,
                                                                      (upper - Jenkins_Frequency_VEI_4) / Jenkins_Frequency_VEI_4_std, loc=Jenkins_Frequency_VEI_4,
                                                                      scale=Jenkins_Frequency_VEI_4_std))
        AP_VEI_4 = AP_VEI_4 + [AP_VEI4]

    AP_VEI_5 = []
    Number = 10000
    lower, upper = 0, 1
    for i in range(Number):
        """Monte Carlo simulation"""
        AP_VEI5 = (np.random.beta(a_ave, b_ave) * stats.truncnorm.rvs((lower - Jenkins_Frequency_VEI_5) / Jenkins_Frequency_VEI_5_std,
                                                                      (upper - Jenkins_Frequency_VEI_5) / Jenkins_Frequency_VEI_5_std, loc=Jenkins_Frequency_VEI_5,
                                                                      scale=Jenkins_Frequency_VEI_5_std))
        AP_VEI_5 = AP_VEI_5 + [AP_VEI5]

    AP_VEI_6 = []
    Number = 10000
    lower, upper = 0, 1
    for i in range(Number):
        """Monte Carlo simulation"""
        AP_VEI6 = (np.random.beta(a_ave, b_ave) * stats.truncnorm.rvs((lower - Jenkins_Frequency_VEI_6) / Jenkins_Frequency_VEI_6_std,
                                                                      (upper - Jenkins_Frequency_VEI_6) / Jenkins_Frequency_VEI_6_std, loc=Jenkins_Frequency_VEI_6,
                                                                      scale=Jenkins_Frequency_VEI_6_std))
        AP_VEI_6 = AP_VEI_6 + [AP_VEI6]

    AP_VEI_7 = []
    Number = 10000
    lower, upper = 0, 1
    for i in range(Number):
        """Monte Carlo simulation"""
        AP_VEI7 = (np.random.beta(a_ave, b_ave) * stats.truncnorm.rvs((lower - Jenkins_Frequency_VEI_7) / Jenkins_Frequency_VEI_7_std,
                                                                      (upper - Jenkins_Frequency_VEI_7) / Jenkins_Frequency_VEI_7_std, loc=Jenkins_Frequency_VEI_7,
                                                                      scale=Jenkins_Frequency_VEI_7_std))
        AP_VEI_7 = AP_VEI_7 + [AP_VEI7]

else:
    Method = "Bayesian update"

    ######
    # Jenkins Global Analogue Probabilities
    if __name__ == '__main__':
        with pm.Model() as AP_model_GVP:
            theta = pm.Beta('theta', alpha=alpha_analogue_jenkins, beta=beta_analogue_jenkins)
            p = pm.Binomial('y', p=theta, observed=y, n=n)
            trace_AP = pm.sample(20000, tune=10000, chains=2, cores=cores, target_accept=0.99)

        #pm.traceplot(trace_AP);
        Alpha_Beta = pm.summary(trace_AP)
        print(Alpha_Beta)

    # Whelley Global Analogue Probabilities
        y = np.array([Total_eruptions_per_10000_year])
        n = np.array([10000])
        N = len(n)

        with pm.Model() as AP_model_Whelley:
            theta_Whelley = pm.Beta('theta', alpha=alpha_analogue_Whelley, beta=beta_analogue_Whelley)
            p_whelley = pm.Binomial('y', p=theta_Whelley, observed=y, n=n)
            trace_AP_Whelley = pm.sample(20000, tune=10000, chains=2, cores=cores, target_accept=0.99)

        #pm.traceplot(trace_AP_Whelley);
        Alpha_Beta_Whelley = pm.summary(trace_AP_Whelley)
        print(Alpha_Beta_Whelley)

        ######

        ######
        # Model averaging for frequency of eruption
        ##Comparing the models and assigning a weight.
        traces = [trace_AP, trace_AP_Whelley]
        model_dict = dict(zip([AP_model_GVP, AP_model_Whelley], traces))
        comp = pm.compare(model_dict, method='stacking')
        comp['index'] = [0, 1]
        comp.set_index("index", inplace = True)

        #Combining the models

        ppc_ave = pm.sample_posterior_predictive_w(traces, 100000, [AP_model_GVP, AP_model_Whelley],
                                weights=comp.weight.sort_index(ascending=True),
                                progressbar=False)

        ppc_Jenkins = pm.sample_posterior_predictive(trace_AP, 100000, AP_model_GVP,
                             progressbar=False)

        ppc_Whelley = pm.sample_posterior_predictive(trace_AP_Whelley, 100000, AP_model_Whelley,
                             progressbar=False)

        ppc_ave_adjust = ppc_ave['y']/10000
        ppc_Jenkins_adjust = ppc_Jenkins['y']/10000
        ppc_Whelley_adjust = ppc_Whelley['y']/10000

        #Visualising the combined models.

        Ave = ppc_ave_adjust.mean()
        Std = ppc_ave_adjust.std()

        Ave_Jenkins = ppc_Jenkins_adjust.mean()
        Std_Jenkins = ppc_Jenkins_adjust.std()
        Ave_Whelley = ppc_Whelley_adjust.mean()
        Std_Whelley = ppc_Whelley_adjust.std()


        a_ave = Ave**2 * ((1-Ave)/Std**2 -1/Ave)
        b_ave = a_ave * (1/Ave-1)

        a_Jenks = Ave_Jenkins**2 * ((1-Ave_Jenkins)/Std_Jenkins**2 -1/Ave_Jenkins)
        b_Jenks = a_ave * (1/Ave_Jenkins-1)
        a_Whelley = Ave_Whelley**2 * ((1-Ave_Whelley)/Std_Whelley**2 -1/Ave_Whelley)
        b_Whelley = a_ave * (1/Ave_Whelley-1)

        data_beta = beta.rvs(a=a_ave, b=b_ave, size=1000000)
        data_Jenkins = beta.rvs(a=a_Jenks, b=b_Jenks, size=1000000)
        data_Whelley = beta.rvs(a=a_Whelley, b=b_Whelley, size=1000000)

        ######
        #Jenkins F-M model
        ######

        Observed_1 = np.array([Number_VEI_3_or_less_eruptions,
                               Number_VEI_4_eruptions,
                               Number_VEI_5_eruptions,
                               Number_VEI_6_eruptions,
                               Number_VEI_7_eruptions])

        # Observed_2 = np.array([Number_VEI_3_or_less_eruptions,
        #                        Number_VEI_4_eruptions,
        #                        Number_VEI_5_eruptions,
        #                        Number_VEI_5_eruptions,
        #                        Number_VEI_5_eruptions])
        #
        # Observed_3 = np.array([Number_VEI_3_or_less_eruptions,
        #                        Number_VEI_4_eruptions,
        #                        Number_VEI_4_eruptions,
        #                        Number_VEI_4_eruptions,
        #                        Number_VEI_4_eruptions])
        #
        # if Jenkins_Frequency_VEI_7 > 0:
        #     Observed = Observed_1
        # elif Jenkins_Frequency_VEI_5 > 0 and Jenkins_Frequency_VEI_7 == 0:
        #     Observed = Observed_2
        # elif Jenkins_Frequency_VEI_5 == 0 and Jenkins_Frequency_VEI_7 == 0:
        #     Observed = Observed_3

        Observed = Observed_1
        Observed_sum = np.sum(Observed)

        # Create model
        with pm.Model() as VEI_model_GVP:
            parameters = pm.Dirichlet('parameters', a=Jenkins_c, shape=Jenkins_shape)
            # Observed data is from a Multinomial distribution
            observed_data = pm.Multinomial(
                'observed_data', n=Observed_sum, p=parameters, shape=Jenkins_shape, observed=Observed)

        with VEI_model_GVP:
            # Sample from the posterior
            trace_VEI_GVP = pm.sample(draws=10000, chains=2, cores=cores, tune=5000,
                                      discard_tuned_samples=True, target_accept=0.99)

        #pm.densityplot(trace_VEI_GVP)
        #pm.summary(trace_VEI_GVP)

        ######
        #Whelley F-M model
        ######

        Whelley_1 = np.array([Whelley_Frequency_VEI_3*10000,
                              Whelley_Frequency_VEI_4*10000,
                              Whelley_Frequency_VEI_5*10000,
                              Whelley_Frequency_VEI_6*10000,
                              Whelley_Frequency_VEI_7*10000])

        Whelley_2 = np.array([Whelley_Frequency_VEI_3*10000,
                              Whelley_Frequency_VEI_4*10000,
                              Whelley_Frequency_VEI_5*10000,
                              Whelley_Frequency_VEI_5*10000*0.1,
                              Whelley_Frequency_VEI_5*10000*0.01])

        Whelley_3 = np.array([Whelley_Frequency_VEI_3*10000,
                              Whelley_Frequency_VEI_4*10000,
                              Whelley_Frequency_VEI_4*10000*0.1,
                              Whelley_Frequency_VEI_4*10000*0.01,
                              Whelley_Frequency_VEI_4*10000*0.001])

        if Whelley_Frequency_VEI_7 > 0:
            Whelley_c = Whelley_1
        elif Whelley_Frequency_VEI_5 > 0 and Whelley_Frequency_VEI_7 == 0:
            Whelley_c = Whelley_2
        elif Whelley_Frequency_VEI_5 == 0 and Whelley_Frequency_VEI_7 == 0:
            Whelley_c = Whelley_3
        else:
           print("ERROR")

        #Whelley_c = Whelley_1
        Whelley_shape = Whelley_c.size

        Observed_1 = np.array([Number_VEI_3_or_less_eruptions,
                               Number_VEI_4_eruptions,
                               Number_VEI_5_eruptions,
                               Number_VEI_6_eruptions,
                               Number_VEI_7_eruptions])
#
#         # Observed_2 = np.array([Number_VEI_3_or_less_eruptions * 1000,
#         #                        Number_VEI_4_eruptions * 1000,
#         #                        Number_VEI_5_eruptions * 1000])
#         #
#         # Observed_3 = np.array([Number_VEI_3_or_less_eruptions * 1000,
#         #                        Number_VEI_4_eruptions * 1000])
#         #
#         # if Whelley_Frequency_VEI_7 > 0:
#         #     Observed = Observed_1
#         # elif Whelley_Frequency_VEI_5 > 0 and Whelley_Frequency_VEI_7 == 0:
#         #     Observed = Observed_2
#         # elif Whelley_Frequency_VEI_5 == 0 and Whelley_Frequency_VEI_7 == 0:
#         #     Observed = Observed_3
#
        Observed = Observed_1
        Observed_sum = np.sum(Observed)

        # Create model
        with pm.Model() as VEI_model_Whelley:
            parameters = pm.Dirichlet('parameters', a=Whelley_c, shape=Whelley_shape)
            # Observed data is from a Multinomial distribution
            observed_data = pm.Multinomial(
                'observed_data', n=Observed_sum, p=parameters, shape=Whelley_shape, observed=Observed)

        with VEI_model_Whelley:
            # Sample from the posterior
            trace_VEI_Whelley = pm.sample(draws=10000, chains=2, cores=cores, tune=5000,
                                          discard_tuned_samples=True, target_accept=0.99)

        #pm.densityplot(trace_VEI_Whelley)

        #####

        #####
        #Comparing the F-M models

        traces_VEI = [trace_VEI_GVP, trace_VEI_Whelley]
        model_dict_VEI = dict(zip([VEI_model_GVP, VEI_model_Whelley], traces_VEI))
        comp_VEI = pm.compare(model_dict_VEI, method='stacking') #BB-pseudo-BMA
        comp_VEI['index'] = [0, 1]
        comp_VEI.set_index("index", inplace = True)

        # Combining the models based on the weighting

        ppc_ave_VEI = pm.sample_posterior_predictive_w(traces_VEI, 1000, [VEI_model_GVP, VEI_model_Whelley],
                                weights=comp.weight.sort_index(ascending=True),
                                progressbar=False)


        try:
            ppc_ave_VEI_df = pd.DataFrame(ppc_ave_VEI)
            ppc_ave_VEI_df2 = pd.DataFrame(ppc_ave_VEI_df['observed_data'].to_list(),
                                       columns=['VEI3', 'VEI4', 'VEI5', 'VEI6', 'VEI7'])
            ppc_ave_VEI_df2 = ppc_ave_VEI_df2.fillna(0)
        except ValueError:
            ppc_ave_VEI_df = pd.DataFrame(ppc_ave_VEI['observed_data'])
            ppc_ave_VEI_df2 = ppc_ave_VEI_df
            ppc_ave_VEI_df2.columns = ['VEI3', 'VEI4', 'VEI5', 'VEI6', 'VEI7']

        ppc_ave_VEI_df2['total'] = ppc_ave_VEI_df2['VEI3'] + ppc_ave_VEI_df2['VEI4'] + ppc_ave_VEI_df2['VEI5'] + \
                                   ppc_ave_VEI_df2['VEI6'] + ppc_ave_VEI_df2['VEI7']
        ppc_ave_VEI_df2['pVEI3'] = ppc_ave_VEI_df2['VEI3'] / ppc_ave_VEI_df2['total']
        ppc_ave_VEI_df2['pVEI4'] = ppc_ave_VEI_df2['VEI4'] / ppc_ave_VEI_df2['total']
        ppc_ave_VEI_df2['pVEI5'] = ppc_ave_VEI_df2['VEI5'] / ppc_ave_VEI_df2['total']
        ppc_ave_VEI_df2['pVEI6'] = ppc_ave_VEI_df2['VEI6'] / ppc_ave_VEI_df2['total']
        ppc_ave_VEI_df2['pVEI7'] = ppc_ave_VEI_df2['VEI7'] / ppc_ave_VEI_df2['total']


        Ave_VEI3 = ppc_ave_VEI_df2['pVEI3'].mean()
        Std_VEI3 = ppc_ave_VEI_df2['pVEI3'].std()

        Ave_VEI4 = ppc_ave_VEI_df2['pVEI4'].mean()
        Std_VEI4 = ppc_ave_VEI_df2['pVEI4'].std()

        Ave_VEI5 = ppc_ave_VEI_df2['pVEI5'].mean()
        Std_VEI5 = ppc_ave_VEI_df2['pVEI5'].std()

        Ave_VEI6 = ppc_ave_VEI_df2['pVEI6'].mean()
        Std_VEI6 = ppc_ave_VEI_df2['pVEI6'].std()

        Ave_VEI7 = ppc_ave_VEI_df2['pVEI7'].mean()
        Std_VEI7 = ppc_ave_VEI_df2['pVEI7'].std()

        # if Ave_VEI3 > 0:
        #     a_ave_VEI3 = Ave_VEI3**2 * ((1-Ave_VEI3)/Std_VEI3**2 -1/Ave_VEI3)
        #     b_ave_VEI3 = a_ave_VEI3 * (1/Ave_VEI3-1)
        #     data_beta_VEI3 = beta.rvs(a=a_ave_VEI3, b=b_ave_VEI3, size=1000000)
        # else:
        #     a_ave_VEI3 = 0
        #     b_ave_VEI3 = 0
        #     data_beta_VEI3 = 0
        #
        # if Ave_VEI4 > 0:
        #     a_ave_VEI4 = Ave_VEI4**2 * ((1-Ave_VEI4)/Std_VEI4**2 -1/Ave_VEI4)
        #     b_ave_VEI4 = a_ave_VEI4 * (1/Ave_VEI4-1)
        #     data_beta_VEI4 = beta.rvs(a=a_ave_VEI4, b=b_ave_VEI4, size=1000000)
        # else:
        #     a_ave_VEI4 = 0
        #     b_ave_VEI4 = 0
        #     data_beta_VEI4 = 0
        #
        # if Ave_VEI5 > 0:
        #     a_ave_VEI5 = Ave_VEI5**2 * ((1-Ave_VEI5)/Std_VEI5**2 -1/Ave_VEI5)
        #     b_ave_VEI5 = a_ave_VEI5 * (1/Ave_VEI5-1)
        #     data_beta_VEI5 = beta.rvs(a=a_ave_VEI5, b=b_ave_VEI5, size=1000000)
        # else:
        #     a_ave_VEI5 = 0
        #     b_ave_VEI5 = 0
        #     data_beta_VEI5 = 0
        #
        # if Ave_VEI6 > 0:
        #     a_ave_VEI6 = Ave_VEI6**2 * ((1-Ave_VEI6)/Std_VEI6**2 -1/Ave_VEI6)
        #     b_ave_VEI6 = a_ave_VEI6 * (1/Ave_VEI6-1)
        #     data_beta_VEI6 = beta.rvs(a=a_ave_VEI6, b=b_ave_VEI6, size=1000000)
        # else:
        #     a_ave_VEI6 = 0
        #     b_ave_VEI6 = 0
        #     data_beta_VEI6 = 0
        #
        # if Ave_VEI6 > 0:
        #     a_ave_VEI7 = Ave_VEI7**2 * ((1-Ave_VEI7)/Std_VEI7**2 -1/Ave_VEI7)
        #     b_ave_VEI7 = a_ave_VEI7 * (1/Ave_VEI7-1)
        #     data_beta_VEI7 = beta.rvs(a=a_ave_VEI7, b=b_ave_VEI7, size=1000000)
        # else:
        #     a_ave_VEI6 = 0
        #     b_ave_VEI6 = 0
        #     data_beta_VEI6 = 0

        VEI_3_mean = Ave_VEI3
        VEI_4_mean = Ave_VEI4
        VEI_5_mean = Ave_VEI5
        VEI_6_mean = Ave_VEI6
        VEI_7_mean = Ave_VEI7

        VEI_3_sd = Std_VEI3
        VEI_4_sd = Std_VEI4
        VEI_5_sd = Std_VEI5
        VEI_6_sd = Std_VEI6
        VEI_7_sd = Std_VEI7

        ######
        # Combining probability of eruption and relative probability of different VEI
        ######
        if VEI_3_mean > 0:
            AP_VEI_3 = []
            Number = 10000
            lower, upper = 0, 1
            for i in range(Number):
                """Monte Carlo simulation"""
                AP_VEI3 = (np.random.beta(a_ave, b_ave) * stats.truncnorm.rvs((lower - VEI_3_mean) / VEI_3_sd,
                                                                              (upper - VEI_3_mean) / VEI_3_sd, loc=VEI_3_mean,
                                                                              scale=VEI_3_sd))
                AP_VEI_3 = AP_VEI_3 + [AP_VEI3]
        else:
            AP_VEI_3 = 0

        if VEI_4_mean > 0:
            AP_VEI_4 = []
            Number = 10000
            lower, upper = 0, 1
            for i in range(Number):
                """Monte Carlo simulation"""
                AP_VEI4 = (np.random.beta(a_ave, b_ave) * stats.truncnorm.rvs((lower - VEI_4_mean) / VEI_4_sd,
                                                                              (upper - VEI_4_mean) / VEI_4_sd, loc=VEI_4_mean,
                                                                              scale=VEI_4_sd))
                AP_VEI_4 = AP_VEI_4 + [AP_VEI4]
        else:
            AP_VEI_4 = 0

        if VEI_5_mean > 0:
            AP_VEI_5 = []
            Number = 10000
            lower, upper = 0, 1
            for i in range(Number):
                """Monte Carlo simulation"""
                AP_VEI5 = (np.random.beta(a_ave, b_ave) * stats.truncnorm.rvs((lower - VEI_5_mean) / VEI_5_sd,
                                                                              (upper - VEI_5_mean) / VEI_5_sd, loc=VEI_5_mean,
                                                                              scale=VEI_5_sd))
                AP_VEI_5 = AP_VEI_5 + [AP_VEI5]
        else:
            AP_VEI_5 = 0

        if VEI_6_mean > 0:
            AP_VEI_6 = []
            Number = 10000
            lower, upper = 0, 1
            for i in range(Number):
                """Monte Carlo simulation"""
                AP_VEI6 = (np.random.beta(a_ave, b_ave) * stats.truncnorm.rvs((lower - VEI_6_mean) / VEI_6_sd,
                                                                              (upper - VEI_6_mean) / VEI_6_sd, loc=VEI_6_mean,
                                                                              scale=VEI_6_sd))
                AP_VEI_6 = AP_VEI_6 + [AP_VEI6]
        else:
            AP_VEI_5 = 0

        if VEI_7_mean > 0:
            AP_VEI_7 = []
            Number = 10000
            lower, upper = 0, 1
            for i in range(Number):
                """Monte Carlo simulation"""
                AP_VEI7 = (np.random.beta(a_ave, b_ave) * stats.truncnorm.rvs((lower - VEI_7_mean) / VEI_7_sd,
                                                                              (upper - VEI_7_mean) / VEI_7_sd, loc=VEI_7_mean,
                                                                              scale=VEI_7_sd))
                AP_VEI_7 = AP_VEI_7 + [AP_VEI7]
        else:
            AP_VEI_7 = 0

#####
# Visualisation of median, 5th, and 95th percentile.
#####

VEI_3_Median = stats.scoreatpercentile(AP_VEI_3, 50)
VEI_3_Percentile_05 = stats.scoreatpercentile(AP_VEI_3, 5)
VEI_3_Percentile_95 = stats.scoreatpercentile(AP_VEI_3, 95)

VEI_4_Median = stats.scoreatpercentile(AP_VEI_4, 50)
VEI_4_Percentile_05 = stats.scoreatpercentile(AP_VEI_4, 5)
VEI_4_Percentile_95 = stats.scoreatpercentile(AP_VEI_4, 95)

VEI_5_Median = stats.scoreatpercentile(AP_VEI_5, 50)
VEI_5_Percentile_05 = stats.scoreatpercentile(AP_VEI_5, 5)
VEI_5_Percentile_95 = stats.scoreatpercentile(AP_VEI_5, 95)

VEI_6_Median = stats.scoreatpercentile(AP_VEI_6, 50)
VEI_6_Percentile_05 = stats.scoreatpercentile(AP_VEI_6, 5)
VEI_6_Percentile_95 = stats.scoreatpercentile(AP_VEI_6, 95)

VEI_7_Median = stats.scoreatpercentile(AP_VEI_7, 50)
VEI_7_Percentile_05 = stats.scoreatpercentile(AP_VEI_7, 5)
VEI_7_Percentile_95 = stats.scoreatpercentile(AP_VEI_7, 95)

Prob_erupt_Median = VEI_3_Median + VEI_4_Median + VEI_5_Median + VEI_6_Median + VEI_7_Median
Prob_erupt_Percentile_05 = VEI_3_Percentile_05 + VEI_4_Percentile_05 + VEI_5_Percentile_05 + VEI_6_Percentile_05 + VEI_7_Percentile_05
Prob_erupt_Percentile_95 = VEI_3_Percentile_95 + VEI_4_Percentile_95 + VEI_5_Percentile_95 + VEI_6_Percentile_95 + VEI_7_Percentile_95

Median = "50th percentile"
Percentile_05 = "5th percentile"
Percentile_95 = "95th percentile"

Percentile_05s = [[3, 4, 5, 6, 7], [VEI_3_Percentile_05, VEI_4_Percentile_05, VEI_5_Percentile_05, VEI_6_Percentile_05, VEI_7_Percentile_05]]
Percentile_95s = [[3, 4, 5, 6, 7], [VEI_3_Percentile_95, VEI_4_Percentile_95, VEI_5_Percentile_95, VEI_6_Percentile_95, VEI_7_Percentile_95]]

Probabilities = pd.DataFrame([[3, VEI_3_Median],
                        [4, VEI_4_Median],
                        [5, VEI_5_Median],
                        [6, VEI_6_Median],
                        [7, VEI_7_Median],
                        [3, VEI_3_Percentile_05],
                        [4, VEI_4_Percentile_05],
                        [5, VEI_5_Percentile_05],
                        [6, VEI_6_Percentile_05],
                        [7, VEI_7_Percentile_05],
                        [3, VEI_3_Percentile_95],
                        [4, VEI_4_Percentile_95],
                        [5, VEI_5_Percentile_95],
                        [6, VEI_6_Percentile_95],
                        [7, VEI_7_Percentile_95]],
                       columns=["VEI", "Annual probability"])

Volcano_name = Volcano_name
textstr = Volcano_name

Probabilities_df = pd.DataFrame([[Volcano_Number, Volcano_name, Region, Method,
                                  Number_of_small_eruptions, Number_of_large_eruptions, Number_of_unknown_eruptions,
                                  Volcano_type_Jenkins, Volcano_type_Whelley,
                                  Prob_erupt_Percentile_05, Prob_erupt_Median, Prob_erupt_Percentile_95,
                                  VEI_3_Percentile_05, VEI_3_Median, VEI_3_Percentile_95,
                                  VEI_4_Percentile_05, VEI_4_Median, VEI_4_Percentile_95,
                                  VEI_5_Percentile_05, VEI_5_Median, VEI_5_Percentile_95,
                                  VEI_6_Percentile_05, VEI_6_Median, VEI_6_Percentile_95,
                                  VEI_7_Percentile_05, VEI_7_Median, VEI_7_Percentile_95]],
                                  columns=["Volcano Number",
                                           "Volcano Name",
                                           "Region",
                                           "Estimate method",
                                           "Recorded small eruptions",
                                           "Recorded large eruptions",
                                           "Recorded unknown eruptions",
                                           "Jenkins classification",
                                           "Whelley classification",
                                           "Probability of eruption 5th percentile",
                                           "Probability of eruption 50th percentile",
                                           "Probability of eruption 95th percentile",
                                           "<= VEI 3 5th percentile",
                                           "<= VEI 3 50th percentile",
                                           "<= VEI 3 95th percentile",
                                           "VEI 4 5th percentile",
                                           "VEI 4 50th percentile",
                                           "VEI 4 95th percentile",
                                           "VEI 5 5th percentile",
                                           "VEI 5 50th percentile",
                                           "VEI 5 95th percentile",
                                           "VEI 6 5th percentile",
                                           "VEI 6 50th percentile",
                                           "VEI 6 95th percentile",
                                           "VEI 7 5th percentile",
                                           "VEI 7 50th percentile",
                                           "VEI 7 95th percentile"])

path_csv = "Probabilities/" + Volcano_name + ".csv"
Probabilities_df.to_csv(path_csv,  index=False)

ax = sns.lineplot(x="VEI", y="Annual probability", markers=True, data=Probabilities)
ax.set(yscale="log")
ax.set_xticks(range(3, 8)) # <--- set the ticks first
ax.tick_params(labelsize=10)
ax.set_xticklabels(['<=3','4','5','6','7'])
ax.set_xlabel("VEI",fontsize=12)
ax.set_ylabel("Annual probability",fontsize=12)
props = dict(boxstyle='round', facecolor='white', alpha=0.5)
# place a text box in upper left in axes coords
ax.text(0.80, 0.95, textstr, transform=ax.transAxes, fontsize=12,
        verticalalignment='top', bbox=props)
path_fig = "Figures/" + Volcano_name + ".png"
plt.savefig(path_fig)
plt.close()