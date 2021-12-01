import pymc3 as pm
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from scipy.stats import beta
import arviz


import csv

def get_change_point (GVP_df, GVP_volcanoes, csv):
    """
This function obtains the change points in the global volcanic record for different regions of the world.

The function obtains change points by rotating the cumulative eruptions for a region through time, rotating the graph
and finding the index of the maximum value. This is then used to find the year of the change point.

    :param GVP_df (csv): GVP global eruptions database.
    :param GVP_volcanoes (csv): GVP global Holocene volcanoes database
    :param csv (Bool): Defines whether a CSV file is save (True) or not (False)
    :return:
    """
    import kneebow
    from kneebow.rotor import Rotor
    # Combing data from the GVP eruption catalogue with the volcano list to obtain regions and countries
    GVP_df2 = GVP_df[['Volcano Number', 'Volcano Name', 'Eruption Category', 'VEI', 'Start Year', 'Start Year Uncertainty']].copy()
    GVP_df2a = pd.merge(GVP_df2, GVP_volcanoes, on=['Volcano Number'], how='left')
    GVP_df3 = GVP_df2a[
        ['Volcano Number', 'Volcano Name_x', 'Eruption Category', 'Primary Volcano Type', 'VEI', 'Start Year',
         'Start Year Uncertainty', 'Region', 'Country', 'Subregion']].copy()
    GVP_df3 = GVP_df3[GVP_df3.VEI.notnull()]
    GVP_df3['quantity'] = 1
    Africa_large = GVP_df3.loc[(GVP_df3['Region'] =='Africa and Red Sea') & (GVP_df3['VEI'] > 3)].sort_values(by="Start Year")
    Alaska_large = GVP_df3.loc[(GVP_df3['Region'] == 'Alaska') & (GVP_df3['VEI'] > 3)].sort_values(by="Start Year")
    AtlanticOcean_large = GVP_df3.loc[(GVP_df3['Region'] == 'Atlantic Ocean') & (GVP_df3['VEI'] > 3)].sort_values(by="Start Year")
    CanadaWesternUSA_large = GVP_df3.loc[(GVP_df3['Region'] == 'Canada and Western USA') & (GVP_df3['VEI'] > 3)].sort_values(by="Start Year")
    HawaiiPacificOcean_large = GVP_df3.loc[(GVP_df3['Region'] == 'Hawaii and Pacific Ocean') & (GVP_df3['VEI'] > 3)].sort_values(by="Start Year")
    IcelandArcticOcean_large = GVP_df3.loc[(GVP_df3['Region'] == 'Iceland and Arctic Ocean') & (GVP_df3['VEI'] > 3)].sort_values(by="Start Year")
    Indonesia_large = GVP_df3.loc[(GVP_df3['Region'] == 'Indonesia') & (GVP_df3['VEI'] > 3)].sort_values(by="Start Year")
    JapanTaiwanMarianas_large = GVP_df3.loc[(GVP_df3['Region'] == 'Japan, Taiwan, Marianas') & (GVP_df3['VEI'] > 3)].sort_values(by="Start Year")
    KamchatkaMainlandAsia_large = GVP_df3.loc[(GVP_df3['Region'] == 'Kamchatka and Mainland Asia') & (GVP_df3['VEI'] > 3)].sort_values(by="Start Year")
    KurilIslands_large = GVP_df3.loc[(GVP_df3['Region'] == 'Kuril Islands') & (GVP_df3['VEI'] > 3)].sort_values(by="Start Year")
    MediterraneanWesternAsia_large = GVP_df3.loc[(GVP_df3['Region'] == 'Mediterranean and Western Asia') & (GVP_df3['VEI'] > 3)].sort_values(by="Start Year")
    MelanesiaAustralia_large = GVP_df3.loc[(GVP_df3['Region'] == 'Melanesia and Australia') & (GVP_df3['VEI'] > 3)].sort_values(by="Start Year")
    MexicoCentralAmerica_large = GVP_df3.loc[(GVP_df3['Region'] == 'México and Central America') & (GVP_df3['VEI'] > 3)].sort_values(by="Start Year")
    MiddleEastIndianOcean_large = GVP_df3.loc[(GVP_df3['Region'] == 'Middle East and Indian Ocean') & (GVP_df3['VEI'] > 3)].sort_values(by="Start Year")
    NewZealandtoFiji_large = GVP_df3.loc[(GVP_df3['Region'] == 'New Zealand to Fiji') & (GVP_df3['VEI'] > 3)].sort_values(by="Start Year")
    PhilippinesSEAsia_large = GVP_df3.loc[(GVP_df3['Region'] == 'Philippines and SE Asia') & (GVP_df3['VEI'] > 3)].sort_values(by="Start Year")
    SouthAmerica_large = GVP_df3.loc[(GVP_df3['Region'] == 'South America') & (GVP_df3['VEI'] > 3)].sort_values(by="Start Year")
    WestIndies_large = GVP_df3.loc[(GVP_df3['Region'] == 'West Indies') & (GVP_df3['VEI'] > 3)].sort_values(by="Start Year")

    Africa_large['CumSum'] = Africa_large['quantity'].cumsum()
    Alaska_large['CumSum'] = Alaska_large['quantity'].cumsum()
    AtlanticOcean_large['CumSum'] = AtlanticOcean_large['quantity'].cumsum()
    CanadaWesternUSA_large['CumSum'] = CanadaWesternUSA_large['quantity'].cumsum()
    HawaiiPacificOcean_large['CumSum'] = HawaiiPacificOcean_large['quantity'].cumsum()
    IcelandArcticOcean_large['CumSum'] = IcelandArcticOcean_large['quantity'].cumsum()
    Indonesia_large['CumSum'] = Indonesia_large['quantity'].cumsum()
    JapanTaiwanMarianas_large['CumSum'] = JapanTaiwanMarianas_large['quantity'].cumsum()
    KamchatkaMainlandAsia_large['CumSum'] = KamchatkaMainlandAsia_large['quantity'].cumsum()
    KurilIslands_large['CumSum'] = KurilIslands_large['quantity'].cumsum()
    MediterraneanWesternAsia_large['CumSum'] = MediterraneanWesternAsia_large['quantity'].cumsum()
    MelanesiaAustralia_large['CumSum'] = MelanesiaAustralia_large['quantity'].cumsum()
    MexicoCentralAmerica_large['CumSum'] = MexicoCentralAmerica_large['quantity'].cumsum()
    MiddleEastIndianOcean_large['CumSum'] = MiddleEastIndianOcean_large['quantity'].cumsum()
    NewZealandtoFiji_large['CumSum'] = NewZealandtoFiji_large['quantity'].cumsum()
    PhilippinesSEAsia_large['CumSum'] = PhilippinesSEAsia_large['quantity'].cumsum()
    SouthAmerica_large['CumSum'] = SouthAmerica_large['quantity'].cumsum()
    WestIndies_large['CumSum'] = WestIndies_large['quantity'].cumsum()

    Africa_array_large = Africa_large[['Start Year', 'CumSum']].to_numpy()
    Alaska_array_large = Alaska_large[['Start Year', 'CumSum']].to_numpy()
    AtlanticOcean_array_large = AtlanticOcean_large[['Start Year', 'CumSum']].to_numpy()
    CanadaWesternUSA_array_large = CanadaWesternUSA_large[['Start Year', 'CumSum']].to_numpy()
    HawaiiPacificOcean_array_large = HawaiiPacificOcean_large[['Start Year', 'CumSum']].to_numpy()
    IcelandArcticOcean_array_large = IcelandArcticOcean_large[['Start Year', 'CumSum']].to_numpy()
    Indonesia_array_large = Indonesia_large[['Start Year', 'CumSum']].to_numpy()
    JapanTaiwanMarianas_array_large = JapanTaiwanMarianas_large[['Start Year', 'CumSum']].to_numpy()
    KamchatkaMainlandAsia_array_large = KamchatkaMainlandAsia_large[['Start Year', 'CumSum']].to_numpy()
    KurilIslands_array_large = KurilIslands_large[['Start Year', 'CumSum']].to_numpy()
    MediterraneanWesternAsia_array_large = MediterraneanWesternAsia_large[['Start Year', 'CumSum']].to_numpy()
    MelanesiaAustralia_array_large = MelanesiaAustralia_large[['Start Year', 'CumSum']].to_numpy()
    MexicoCentralAmerica_array_large = MexicoCentralAmerica_large[['Start Year', 'CumSum']].to_numpy()
    MiddleEastIndianOcean_array_large = MiddleEastIndianOcean_large[['Start Year', 'CumSum']].to_numpy()
    NewZealandtoFiji_array_large = NewZealandtoFiji_large[['Start Year', 'CumSum']].to_numpy()
    PhilippinesSEAsia_array_large = PhilippinesSEAsia_large[['Start Year', 'CumSum']].to_numpy()
    SouthAmerica_array_large = SouthAmerica_large[['Start Year', 'CumSum']].to_numpy()
    WestIndies_array_large = WestIndies_large[['Start Year', 'CumSum']].to_numpy()

    rotor_Africa_large = Rotor()
    rotor_Alaska_large = Rotor()
    rotor_AtlanticOcean_large = Rotor()
    rotor_CanadaWesternUSA_large = Rotor()
    rotor_HawaiiPacificOcean_large = Rotor()
    rotor_IcelandArcticOcean_large = Rotor()
    rotor_Indonesia_large = Rotor()
    rotor_JapanTaiwanMarianas_large = Rotor()
    rotor_KamchatkaMainlandAsia_large = Rotor()
    rotor_KurilIslands_large = Rotor()
    rotor_MediterraneanWesternAsia_large = Rotor()
    rotor_MelanesiaAustralia_large = Rotor()
    rotor_MexicoCentralAmerica_large = Rotor()
    rotor_MiddleEastIndianOcean_large = Rotor()
    rotor_NewZealandtoFiji_large = Rotor()
    rotor_PhilippinesSEAsia_large = Rotor()
    rotor_SouthAmerica_large = Rotor()
    rotor_WestIndies_large = Rotor()

    rotor_Africa_large.fit_rotate(Africa_array_large)
    rotor_Alaska_large.fit_rotate(Alaska_array_large)
    rotor_AtlanticOcean_large.fit_rotate(AtlanticOcean_array_large)
    rotor_CanadaWesternUSA_large.fit_rotate(CanadaWesternUSA_array_large)
    rotor_HawaiiPacificOcean_large.fit_rotate(HawaiiPacificOcean_array_large)
    rotor_IcelandArcticOcean_large.fit_rotate(IcelandArcticOcean_array_large)
    rotor_Indonesia_large.fit_rotate(Indonesia_array_large)
    rotor_JapanTaiwanMarianas_large.fit_rotate(JapanTaiwanMarianas_array_large)
    rotor_KamchatkaMainlandAsia_large.fit_rotate(KamchatkaMainlandAsia_array_large)
    rotor_KurilIslands_large.fit_rotate(KurilIslands_array_large)
    rotor_MediterraneanWesternAsia_large.fit_rotate(MediterraneanWesternAsia_array_large)
    rotor_MelanesiaAustralia_large.fit_rotate(MelanesiaAustralia_array_large)
    rotor_MexicoCentralAmerica_large.fit_rotate(MexicoCentralAmerica_array_large)
    rotor_MiddleEastIndianOcean_large.fit_rotate(MiddleEastIndianOcean_array_large)
    rotor_NewZealandtoFiji_large.fit_rotate(NewZealandtoFiji_array_large)
    rotor_PhilippinesSEAsia_large.fit_rotate(PhilippinesSEAsia_array_large)
    rotor_SouthAmerica_large.fit_rotate(SouthAmerica_array_large)
    rotor_WestIndies_large.fit_rotate(WestIndies_array_large)

    elbow_idx_Africa_large = rotor_Africa_large.get_elbow_index()
    elbow_idx_Alaska_large = rotor_Alaska_large.get_elbow_index()
    elbow_idx_AtlanticOcean_large = rotor_AtlanticOcean_large.get_elbow_index()
    elbow_idx_CanadaWesternUSA_large = rotor_CanadaWesternUSA_large.get_elbow_index()
    elbow_idx_HawaiiPacificOcean_large = rotor_HawaiiPacificOcean_large.get_elbow_index()
    elbow_idx_IcelandArcticOcean_large = rotor_IcelandArcticOcean_large.get_elbow_index()
    elbow_idx_Indonesia_large = rotor_Indonesia_large.get_elbow_index()
    elbow_idx_JapanTaiwanMarianas_large = rotor_JapanTaiwanMarianas_large.get_elbow_index()
    elbow_idx_KamchatkaMainlandAsia_large = rotor_KamchatkaMainlandAsia_large.get_elbow_index()
    elbow_idx_KurilIslands_large = rotor_KurilIslands_large.get_elbow_index()
    elbow_idx_MediterraneanWesternAsia_large = rotor_MediterraneanWesternAsia_large.get_elbow_index()
    elbow_idx_MelanesiaAustralia_large = rotor_MelanesiaAustralia_large.get_elbow_index()
    elbow_idx_MexicoCentralAmerica_large = rotor_MexicoCentralAmerica_large.get_elbow_index()
    elbow_idx_MiddleEastIndianOcean_large = rotor_MiddleEastIndianOcean_large.get_elbow_index()
    elbow_idx_NewZealandtoFiji_large = rotor_NewZealandtoFiji_large.get_elbow_index()
    elbow_idx_PhilippinesSEAsia_large = rotor_PhilippinesSEAsia_large.get_elbow_index()
    elbow_idx_SouthAmerica_large = rotor_SouthAmerica_large.get_elbow_index()
    elbow_idx_WestIndies_large = rotor_WestIndies_large.get_elbow_index()

    change_point_Africa_large = Africa_large['Start Year'].iloc[elbow_idx_Africa_large]
    change_point_Alaska_large = Alaska_large['Start Year'].iloc[elbow_idx_Alaska_large]
    change_point_AtlanticOcean_large = AtlanticOcean_large['Start Year'].iloc[elbow_idx_AtlanticOcean_large]
    change_point_CanadaWesternUSA_large = CanadaWesternUSA_large['Start Year'].iloc[elbow_idx_CanadaWesternUSA_large]
    change_point_HawaiiPacificOcean_large = HawaiiPacificOcean_large['Start Year'].iloc[elbow_idx_HawaiiPacificOcean_large]
    change_point_IcelandArcticOcean_large = IcelandArcticOcean_large['Start Year'].iloc[elbow_idx_IcelandArcticOcean_large]
    change_point_Indonesia_large = Indonesia_large['Start Year'].iloc[elbow_idx_Indonesia_large]
    change_point_JapanTaiwanMarianas_large = JapanTaiwanMarianas_large['Start Year'].iloc[elbow_idx_JapanTaiwanMarianas_large]
    change_point_KamchatkaMainlandAsia_large = KamchatkaMainlandAsia_large['Start Year'].iloc[elbow_idx_KamchatkaMainlandAsia_large]
    change_point_KurilIslands_large = KurilIslands_large['Start Year'].iloc[elbow_idx_KurilIslands_large]
    change_point_MediterraneanWesternAsia_large = MediterraneanWesternAsia_large['Start Year'].iloc[
        elbow_idx_MediterraneanWesternAsia_large]
    change_point_MelanesiaAustralia_large = MelanesiaAustralia_large['Start Year'].iloc[elbow_idx_MelanesiaAustralia_large]
    change_point_MexicoCentralAmerica_large = MexicoCentralAmerica_large['Start Year'].iloc[elbow_idx_MexicoCentralAmerica_large]
    change_point_MiddleEastIndianOcean_large = MiddleEastIndianOcean_large['Start Year'].iloc[elbow_idx_MiddleEastIndianOcean_large]
    change_point_NewZealandtoFiji_large = NewZealandtoFiji_large['Start Year'].iloc[elbow_idx_NewZealandtoFiji_large]
    change_point_PhilippinesSEAsia_large = PhilippinesSEAsia_large['Start Year'].iloc[elbow_idx_PhilippinesSEAsia_large]
    change_point_SouthAmerica_large = SouthAmerica_large['Start Year'].iloc[elbow_idx_SouthAmerica_large]
    change_point_WestIndies_large = WestIndies_large['Start Year'].iloc[elbow_idx_WestIndies_large]

    Africa_small = GVP_df3.loc[(GVP_df3['Region'] =='Africa and Red Sea') & (GVP_df3['VEI'] < 4)].sort_values(by="Start Year")
    Alaska_small = GVP_df3.loc[(GVP_df3['Region'] == 'Alaska') & (GVP_df3['VEI'] < 4)].sort_values(by="Start Year")
    AtlanticOcean_small = GVP_df3.loc[(GVP_df3['Region'] == 'Atlantic Ocean') & (GVP_df3['VEI'] < 4)].sort_values(by="Start Year")
    CanadaWesternUSA_small = GVP_df3.loc[(GVP_df3['Region'] == 'Canada and Western USA') & (GVP_df3['VEI'] < 4)].sort_values(by="Start Year")
    HawaiiPacificOcean_small = GVP_df3.loc[(GVP_df3['Region'] == 'Hawaii and Pacific Ocean') & (GVP_df3['VEI'] < 4)].sort_values(by="Start Year")
    IcelandArcticOcean_small = GVP_df3.loc[(GVP_df3['Region'] == 'Iceland and Arctic Ocean') & (GVP_df3['VEI'] < 4)].sort_values(by="Start Year")
    Indonesia_small = GVP_df3.loc[(GVP_df3['Region'] == 'Indonesia') & (GVP_df3['VEI'] < 4)].sort_values(by="Start Year")
    JapanTaiwanMarianas_small = GVP_df3.loc[(GVP_df3['Region'] == 'Japan, Taiwan, Marianas') & (GVP_df3['VEI'] < 4)].sort_values(by="Start Year")
    KamchatkaMainlandAsia_small = GVP_df3.loc[(GVP_df3['Region'] == 'Kamchatka and Mainland Asia') & (GVP_df3['VEI'] < 4)].sort_values(by="Start Year")
    KurilIslands_small = GVP_df3.loc[(GVP_df3['Region'] == 'Kuril Islands') & (GVP_df3['VEI'] < 4)].sort_values(by="Start Year")
    MediterraneanWesternAsia_small = GVP_df3.loc[(GVP_df3['Region'] == 'Mediterranean and Western Asia') & (GVP_df3['VEI'] < 4)].sort_values(by="Start Year")
    MelanesiaAustralia_small = GVP_df3.loc[(GVP_df3['Region'] == 'Melanesia and Australia') & (GVP_df3['VEI'] < 4)].sort_values(by="Start Year")
    MexicoCentralAmerica_small = GVP_df3.loc[(GVP_df3['Region'] == 'México and Central America') & (GVP_df3['VEI'] < 4)].sort_values(by="Start Year")
    MiddleEastIndianOcean_small = GVP_df3.loc[(GVP_df3['Region'] == 'Middle East and Indian Ocean') & (GVP_df3['VEI'] < 4)].sort_values(by="Start Year")
    NewZealandtoFiji_small = GVP_df3.loc[(GVP_df3['Region'] == 'New Zealand to Fiji') & (GVP_df3['VEI'] < 4)].sort_values(by="Start Year")
    PhilippinesSEAsia_small = GVP_df3.loc[(GVP_df3['Region'] == 'Philippines and SE Asia') & (GVP_df3['VEI'] < 4)].sort_values(by="Start Year")
    SouthAmerica_small = GVP_df3.loc[(GVP_df3['Region'] == 'South America') & (GVP_df3['VEI'] < 4)].sort_values(by="Start Year")
    WestIndies_small = GVP_df3.loc[(GVP_df3['Region'] == 'West Indies') & (GVP_df3['VEI'] < 4)].sort_values(by="Start Year")

    Africa_small['CumSum'] = Africa_small['quantity'].cumsum()
    Alaska_small['CumSum'] = Alaska_small['quantity'].cumsum()
    AtlanticOcean_small['CumSum'] = AtlanticOcean_small['quantity'].cumsum()
    CanadaWesternUSA_small['CumSum'] = CanadaWesternUSA_small['quantity'].cumsum()
    HawaiiPacificOcean_small['CumSum'] = HawaiiPacificOcean_small['quantity'].cumsum()
    IcelandArcticOcean_small['CumSum'] = IcelandArcticOcean_small['quantity'].cumsum()
    Indonesia_small['CumSum'] = Indonesia_small['quantity'].cumsum()
    JapanTaiwanMarianas_small['CumSum'] = JapanTaiwanMarianas_small['quantity'].cumsum()
    KamchatkaMainlandAsia_small['CumSum'] = KamchatkaMainlandAsia_small['quantity'].cumsum()
    KurilIslands_small['CumSum'] = KurilIslands_small['quantity'].cumsum()
    MediterraneanWesternAsia_small['CumSum'] = MediterraneanWesternAsia_small['quantity'].cumsum()
    MelanesiaAustralia_small['CumSum'] = MelanesiaAustralia_small['quantity'].cumsum()
    MexicoCentralAmerica_small['CumSum'] = MexicoCentralAmerica_small['quantity'].cumsum()
    MiddleEastIndianOcean_small['CumSum'] = MiddleEastIndianOcean_small['quantity'].cumsum()
    NewZealandtoFiji_small['CumSum'] = NewZealandtoFiji_small['quantity'].cumsum()
    PhilippinesSEAsia_small['CumSum'] = PhilippinesSEAsia_small['quantity'].cumsum()
    SouthAmerica_small['CumSum'] = SouthAmerica_small['quantity'].cumsum()
    WestIndies_small['CumSum'] = WestIndies_small['quantity'].cumsum()

    Africa_array_small = Africa_small[['Start Year', 'CumSum']].to_numpy()
    Alaska_array_small = Alaska_small[['Start Year', 'CumSum']].to_numpy()
    AtlanticOcean_array_small = AtlanticOcean_small[['Start Year', 'CumSum']].to_numpy()
    CanadaWesternUSA_array_small = CanadaWesternUSA_small[['Start Year', 'CumSum']].to_numpy()
    HawaiiPacificOcean_array_small = HawaiiPacificOcean_small[['Start Year', 'CumSum']].to_numpy()
    IcelandArcticOcean_array_small = IcelandArcticOcean_small[['Start Year', 'CumSum']].to_numpy()
    Indonesia_array_small = Indonesia_small[['Start Year', 'CumSum']].to_numpy()
    JapanTaiwanMarianas_array_small = JapanTaiwanMarianas_small[['Start Year', 'CumSum']].to_numpy()
    KamchatkaMainlandAsia_array_small = KamchatkaMainlandAsia_small[['Start Year', 'CumSum']].to_numpy()
    KurilIslands_array_small = KurilIslands_small[['Start Year', 'CumSum']].to_numpy()
    MediterraneanWesternAsia_array_small = MediterraneanWesternAsia_small[['Start Year', 'CumSum']].to_numpy()
    MelanesiaAustralia_array_small = MelanesiaAustralia_small[['Start Year', 'CumSum']].to_numpy()
    MexicoCentralAmerica_array_small = MexicoCentralAmerica_small[['Start Year', 'CumSum']].to_numpy()
    MiddleEastIndianOcean_array_small = MiddleEastIndianOcean_small[['Start Year', 'CumSum']].to_numpy()
    NewZealandtoFiji_array_small = NewZealandtoFiji_small[['Start Year', 'CumSum']].to_numpy()
    PhilippinesSEAsia_array_small = PhilippinesSEAsia_small[['Start Year', 'CumSum']].to_numpy()
    SouthAmerica_array_small = SouthAmerica_small[['Start Year', 'CumSum']].to_numpy()
    WestIndies_array_small = WestIndies_small[['Start Year', 'CumSum']].to_numpy()

    rotor_Africa_small = Rotor()
    rotor_Alaska_small = Rotor()
    rotor_AtlanticOcean_small = Rotor()
    rotor_CanadaWesternUSA_small = Rotor()
    rotor_HawaiiPacificOcean_small = Rotor()
    rotor_IcelandArcticOcean_small = Rotor()
    rotor_Indonesia_small = Rotor()
    rotor_JapanTaiwanMarianas_small = Rotor()
    rotor_KamchatkaMainlandAsia_small = Rotor()
    rotor_KurilIslands_small = Rotor()
    rotor_MediterraneanWesternAsia_small = Rotor()
    rotor_MelanesiaAustralia_small = Rotor()
    rotor_MexicoCentralAmerica_small = Rotor()
    rotor_MiddleEastIndianOcean_small = Rotor()
    rotor_NewZealandtoFiji_small = Rotor()
    rotor_PhilippinesSEAsia_small = Rotor()
    rotor_SouthAmerica_small = Rotor()
    rotor_WestIndies_small = Rotor()

    rotor_Africa_small.fit_rotate(Africa_array_small)
    rotor_Alaska_small.fit_rotate(Alaska_array_small)
    rotor_AtlanticOcean_small.fit_rotate(AtlanticOcean_array_small)
    rotor_CanadaWesternUSA_small.fit_rotate(CanadaWesternUSA_array_small)
    rotor_HawaiiPacificOcean_small.fit_rotate(HawaiiPacificOcean_array_small)
    rotor_IcelandArcticOcean_small.fit_rotate(IcelandArcticOcean_array_small)
    rotor_Indonesia_small.fit_rotate(Indonesia_array_small)
    rotor_JapanTaiwanMarianas_small.fit_rotate(JapanTaiwanMarianas_array_small)
    rotor_KamchatkaMainlandAsia_small.fit_rotate(KamchatkaMainlandAsia_array_small)
    rotor_KurilIslands_small.fit_rotate(KurilIslands_array_small)
    rotor_MediterraneanWesternAsia_small.fit_rotate(MediterraneanWesternAsia_array_small)
    rotor_MelanesiaAustralia_small.fit_rotate(MelanesiaAustralia_array_small)
    rotor_MexicoCentralAmerica_small.fit_rotate(MexicoCentralAmerica_array_small)
    rotor_MiddleEastIndianOcean_small.fit_rotate(MiddleEastIndianOcean_array_small)
    rotor_NewZealandtoFiji_small.fit_rotate(NewZealandtoFiji_array_small)
    rotor_PhilippinesSEAsia_small.fit_rotate(PhilippinesSEAsia_array_small)
    rotor_SouthAmerica_small.fit_rotate(SouthAmerica_array_small)
    rotor_WestIndies_small.fit_rotate(WestIndies_array_small)

    elbow_idx_Africa_small = rotor_Africa_small.get_elbow_index()
    elbow_idx_Alaska_small = rotor_Alaska_small.get_elbow_index()
    elbow_idx_AtlanticOcean_small = rotor_AtlanticOcean_small.get_elbow_index()
    elbow_idx_CanadaWesternUSA_small = rotor_CanadaWesternUSA_small.get_elbow_index()
    elbow_idx_HawaiiPacificOcean_small = rotor_HawaiiPacificOcean_small.get_elbow_index()
    elbow_idx_IcelandArcticOcean_small = rotor_IcelandArcticOcean_small.get_elbow_index()
    elbow_idx_Indonesia_small = rotor_Indonesia_small.get_elbow_index()
    elbow_idx_JapanTaiwanMarianas_small = rotor_JapanTaiwanMarianas_small.get_elbow_index()
    elbow_idx_KamchatkaMainlandAsia_small = rotor_KamchatkaMainlandAsia_small.get_elbow_index()
    elbow_idx_KurilIslands_small = rotor_KurilIslands_small.get_elbow_index()
    elbow_idx_MediterraneanWesternAsia_small = rotor_MediterraneanWesternAsia_small.get_elbow_index()
    elbow_idx_MelanesiaAustralia_small = rotor_MelanesiaAustralia_small.get_elbow_index()
    elbow_idx_MexicoCentralAmerica_small = rotor_MexicoCentralAmerica_small.get_elbow_index()
    elbow_idx_MiddleEastIndianOcean_small = rotor_MiddleEastIndianOcean_small.get_elbow_index()
    elbow_idx_NewZealandtoFiji_small = rotor_NewZealandtoFiji_small.get_elbow_index()
    elbow_idx_PhilippinesSEAsia_small = rotor_PhilippinesSEAsia_small.get_elbow_index()
    elbow_idx_SouthAmerica_small = rotor_SouthAmerica_small.get_elbow_index()
    elbow_idx_WestIndies_small = rotor_WestIndies_small.get_elbow_index()

    change_point_Africa_small = Africa_small['Start Year'].iloc[elbow_idx_Africa_small]
    change_point_Alaska_small = Alaska_small['Start Year'].iloc[elbow_idx_Alaska_small]
    change_point_AtlanticOcean_small = AtlanticOcean_small['Start Year'].iloc[elbow_idx_AtlanticOcean_small]
    change_point_CanadaWesternUSA_small = CanadaWesternUSA_small['Start Year'].iloc[elbow_idx_CanadaWesternUSA_small]
    change_point_HawaiiPacificOcean_small = HawaiiPacificOcean_small['Start Year'].iloc[elbow_idx_HawaiiPacificOcean_small]
    change_point_IcelandArcticOcean_small = IcelandArcticOcean_small['Start Year'].iloc[elbow_idx_IcelandArcticOcean_small]
    change_point_Indonesia_small = Indonesia_small['Start Year'].iloc[elbow_idx_Indonesia_small]
    change_point_JapanTaiwanMarianas_small = JapanTaiwanMarianas_small['Start Year'].iloc[elbow_idx_JapanTaiwanMarianas_small]
    change_point_KamchatkaMainlandAsia_small = KamchatkaMainlandAsia_small['Start Year'].iloc[elbow_idx_KamchatkaMainlandAsia_small]
    change_point_KurilIslands_small = KurilIslands_small['Start Year'].iloc[elbow_idx_KurilIslands_small]
    change_point_MediterraneanWesternAsia_small = MediterraneanWesternAsia_small['Start Year'].iloc[elbow_idx_MediterraneanWesternAsia_small]
    change_point_MelanesiaAustralia_small = MelanesiaAustralia_small['Start Year'].iloc[elbow_idx_MelanesiaAustralia_small]
    change_point_MexicoCentralAmerica_small = MexicoCentralAmerica_small['Start Year'].iloc[elbow_idx_MexicoCentralAmerica_small]
    change_point_MiddleEastIndianOcean_small = MiddleEastIndianOcean_small['Start Year'].iloc[elbow_idx_MiddleEastIndianOcean_small]
    change_point_NewZealandtoFiji_small = NewZealandtoFiji_small['Start Year'].iloc[elbow_idx_NewZealandtoFiji_small]
    change_point_PhilippinesSEAsia_small = PhilippinesSEAsia_small['Start Year'].iloc[elbow_idx_PhilippinesSEAsia_small]
    change_point_SouthAmerica_small = SouthAmerica_small['Start Year'].iloc[elbow_idx_SouthAmerica_small]
    change_point_WestIndies_small = WestIndies_small['Start Year'].iloc[elbow_idx_WestIndies_small]

    change_point_dict = {'Africa': [change_point_Africa_small,
                                    change_point_Africa_large],
                           'Alaska' : [change_point_Alaska_small,
                                       change_point_Alaska_large],
                           'Atlantic Ocean' : [change_point_AtlanticOcean_small,
                                               change_point_AtlanticOcean_large],
                           'Canada and Western USA' : [change_point_CanadaWesternUSA_small,
                                                       change_point_CanadaWesternUSA_large],
                           'Hawaii and Pacific Ocean' : [change_point_HawaiiPacificOcean_small,
                                                         change_point_HawaiiPacificOcean_large],
                           'Iceland and Arctic Ocean' : [change_point_IcelandArcticOcean_small,
                                                         change_point_IcelandArcticOcean_large],
                           'Indonesia' : [change_point_Indonesia_small,
                                          change_point_Indonesia_large],
                           'Japan, Tawain, Marianas' : [change_point_JapanTaiwanMarianas_small,
                                                        change_point_JapanTaiwanMarianas_large],
                           'Kamchatka and Mainland Asia' : [change_point_KamchatkaMainlandAsia_small,
                                                            change_point_KamchatkaMainlandAsia_large],
                           'Kuril Islands' : [change_point_KurilIslands_small,
                                              change_point_KurilIslands_large],
                           'Mediterranean and Western Asia' : [change_point_MediterraneanWesternAsia_small,
                                                               change_point_MediterraneanWesternAsia_large],
                           'Melanesia and Australia' : [change_point_MelanesiaAustralia_small,
                                                        change_point_MelanesiaAustralia_large],
                           'Mexico and Central America' : [change_point_MexicoCentralAmerica_small,
                                                           change_point_MexicoCentralAmerica_large],
                           'Middle East and Indian Ocean' : [change_point_MiddleEastIndianOcean_small,
                                                             change_point_MiddleEastIndianOcean_large],
                           'New Zealand to Fiji' : [change_point_NewZealandtoFiji_small,
                                                    change_point_NewZealandtoFiji_large],
                           'Philippines and SE Asia' : [change_point_PhilippinesSEAsia_small,
                                                        change_point_PhilippinesSEAsia_large],
                           'South America' : [change_point_SouthAmerica_small,
                                              change_point_SouthAmerica_large],
                           'West Indies' : [change_point_WestIndies_small,
                                            change_point_WestIndies_large]}
    change_point_df = pd.DataFrame.from_dict(change_point_dict, orient='index', columns=['small', 'large'])

    if csv==True:
        change_point_df = pd.DataFrame.from_dict(change_point_dict, orient='index', columns=['small','large'])
        change_point_df.to_csv('change_point_auto.csv')

    return(change_point_df)

def get_complete_record (GVP_df, GVP_volcanoes, change_points, confirmed, method):
    """
    GVP_df: Dataframe of all eruptions in the eruption record
    GVP_volcanoes: all volcanoes in the eruption record
    change_points: Dataframe of change points for each region/country
    confirmed: Defines whether only confirmed (True) or confirmed and uncertain (False) eruptions are included
    Returns:
        Complete_record: Returns three dataframes containing all eruptions within the record considered complete at 5/50/95th percentile.
    """
    # Combing data from the GVP eruption catalogue with the volcano list to obtain regions and countries
    Eruptions = None
    Eruptions_5 = None
    Eruptions_50 = None
    Eruptions_95 = None
    GVP_df2 = GVP_df[['Volcano Number', 'Volcano Name', 'Eruption Category', 'VEI', 'Start Year', 'Start Year Uncertainty']].copy()
    GVP_df2a = pd.merge(GVP_df2, GVP_volcanoes, on=['Volcano Number'], how='left')
    GVP_df3 = GVP_df2a[
        ['Volcano Number', 'Volcano Name_x', 'Eruption Category', 'Primary Volcano Type', 'VEI', 'Start Year',
         'Start Year Uncertainty', 'Region', 'Country', 'Subregion']].copy()
    GVP_df3a = pd.merge(GVP_df3, change_points, on=['Region'], how='left')
    if method == 'Auto':
        if confirmed==True:
            Eruptions_small = GVP_df3a[(GVP_df3a['Eruption Category'] == 'Confirmed Eruption')
                                         & (GVP_df3a['VEI'] < 4)
                                         & (GVP_df3a['Start Year'] > (GVP_df3a['Small']))]
            Eruptions_large = GVP_df3a[(GVP_df3a['Eruption Category'] == 'Confirmed Eruption')
                                         & (GVP_df3a['VEI'] > 3)
                                         & (GVP_df3a['Start Year'] > (GVP_df3a['Large']))]
        else:
            Eruptions_small = GVP_df3a[(GVP_df3a['VEI'] < 4)
                                         & (GVP_df3a['Start Year'] > (GVP_df3a['Small']))]
            Eruptions_large = GVP_df3a[(GVP_df3a['VEI'] > 3)
                                         & (GVP_df3a['Start Year'] > (GVP_df3a['Large']))]
        Eruptions_frames = [Eruptions_small, Eruptions_large]
        Eruptions = pd.concat(Eruptions_frames)
        print("Number of eruptions in complete global record: ", len(Eruptions))
    else:
        if confirmed==True:
            Eruptions_small_50 = GVP_df3a[(GVP_df3a['Eruption Category'] == 'Confirmed Eruption')
                                         & (GVP_df3a['VEI'] < 4)
                                         & (GVP_df3a['Start Year'] > (GVP_df3a['Small 50']))]
            Eruptions_large_50 = GVP_df3a[(GVP_df3a['Eruption Category'] == 'Confirmed Eruption')
                                         & (GVP_df3a['VEI'] > 3)
                                         & (GVP_df3a['Start Year'] > (GVP_df3a['Large 50']))]
            Eruptions_small_5 = GVP_df3a[(GVP_df3a['Eruption Category'] == 'Confirmed Eruption')
                                         & (GVP_df3a['VEI'] < 4)
                                         & (GVP_df3a['Start Year'] > (GVP_df3a['Small 5']))]
            Eruptions_large_5 = GVP_df3a[(GVP_df3a['Eruption Category'] == 'Confirmed Eruption')
                                         & (GVP_df3a['VEI'] > 3)
                                         & (GVP_df3a['Start Year'] > (GVP_df3a['Large 5']))]
            Eruptions_small_95 = GVP_df3a[(GVP_df3a['Eruption Category'] == 'Confirmed Eruption')
                                         & (GVP_df3a['VEI'] < 4)
                                         & (GVP_df3a['Start Year'] > (GVP_df3a['Small 95']))]
            Eruptions_large_95 = GVP_df3a[(GVP_df3a['Eruption Category'] == 'Confirmed Eruption')
                                         & (GVP_df3a['VEI'] > 3)
                                         & (GVP_df3a['Start Year'] > (GVP_df3a['Large 95']))]
            Eruptions_5_frames = [Eruptions_small_5, Eruptions_large_5]
            Eruptions_5 = pd.concat(Eruptions_5_frames)
            Eruptions_50_frames = [Eruptions_small_50, Eruptions_large_50]
            Eruptions_50 = pd.concat(Eruptions_50_frames)
            Eruptions_95_frames = [Eruptions_small_95, Eruptions_large_95]
            Eruptions_95 = pd.concat(Eruptions_95_frames)

        else:
            Eruptions_small_50 = GVP_df3a[(GVP_df3a['VEI'] < 4)
                                         & (GVP_df3a['Start Year'] > (GVP_df3a['Small 50']))]
            Eruptions_large_50 = GVP_df3a[(GVP_df3a['VEI'] > 3)
                                         & (GVP_df3a['Start Year'] > (GVP_df3a['Large 50']))]
            Eruptions_small_5 = GVP_df3a[(GVP_df3a['VEI'] < 4)
                                         & (GVP_df3a['Start Year'] > (GVP_df3a['Small 5']))]
            Eruptions_large_5 = GVP_df3a[(GVP_df3a['VEI'] > 3)
                                         & (GVP_df3a['Start Year'] > (GVP_df3a['Large 5']))]
            Eruptions_small_95 = GVP_df3a[(GVP_df3a['VEI'] < 4)
                                         & (GVP_df3a['Start Year'] > (GVP_df3a['Small 95']))]
            Eruptions_large_95 = GVP_df3a[(GVP_df3a['VEI'] > 3)
                                         & (GVP_df3a['Start Year'] > (GVP_df3a['Large 95']))]
            Eruptions_5_frames = [Eruptions_small_5, Eruptions_large_5]
            Eruptions_5 = pd.concat(Eruptions_5_frames)
            Eruptions_50_frames = [Eruptions_small_50, Eruptions_large_50]
            Eruptions_50 = pd.concat(Eruptions_50_frames)
            Eruptions_95_frames = [Eruptions_small_95, Eruptions_large_95]
            Eruptions_95 = pd.concat(Eruptions_95_frames)

    if method == "MeadMagill5":
        Eruptions = Eruptions_5
        print("Complete global record using 5th percentile change point: ", len(Eruptions_5))

    elif method == "MeadMagill50":
        Eruptions = Eruptions_50
        print("Complete global record using 50th percentile change point: ", len(Eruptions_50))

    elif method == "MeadMagill95":
        Eruptions = Eruptions_95
        print("Complete global record using 95th percentile change point: ", len(Eruptions_95))

    return (Eruptions)

def get_eruption_count (complete_record, Volcano_Number):
    """
    This function will count the number of eruptions for a paticular volcano within an eruption catalogue
    :param complete_record: Dataframe of the complete eruption record
    :param Volcano_Number: Int of the GVP volcano number for the volcano of interest
    :return: Dictionary of the eruption counts
    """
    Small_eruptions = complete_record[(complete_record['Volcano Number'] == Volcano_Number)
                              & (complete_record['VEI'] < 4)]
    Number_of_small_eruptions = len(Small_eruptions)

    Large_eruptions = complete_record[(complete_record['Volcano Number'] == Volcano_Number)
                              & (complete_record['VEI'] > 3)]
    Number_of_large_eruptions = len(Large_eruptions)

    Unknown_eruptions = complete_record[(complete_record['Volcano Number'] == Volcano_Number)
                                & (complete_record.VEI.isnull())]
    Number_of_unknown_eruptions = len(Unknown_eruptions)
    Total_number_of_eruptions = Number_of_small_eruptions+Number_of_large_eruptions+Number_of_unknown_eruptions
    count = {'Number of small eruptions' : Number_of_small_eruptions,
             'Number of large eruptions' : Number_of_large_eruptions,
             'Number of unknown VEI eruption' : Number_of_unknown_eruptions,
             'Total number of eruptions' : Total_number_of_eruptions}
    return (count)

def assign_A1_classification (GVP_volcanoes):
    """
    This function assigns analogue classes to each volcano within the GVP Holocene Volcano List.
    The analogue classes use the Jenkins et al. (2012) criteria for assigning analogue classes
    :param GVP_volcanoes: Dataframe of the GVP Holocene volcano list
    :return: A dataframe with analogue classes added to the GVP Holocene Volcano list.
    """
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

    A1_volcano_classified = GVP_volcanoes[['Volcano Number', 'Volcano Name', 'Volcano type']]
    nan_value = float("NaN")
    A1_volcano_classified.replace("", nan_value)
    A1_volcano_classified.dropna(subset=["Volcano type"])

    return(A1_volcano_classified)

def calc_mean_annual_frequency_analogue (Analogue_volcanoes, complete_record, confirmed, method, year):
    """
    This function obtains the mean annual frequency of eruption for each analogue class using the global eruption record.

    :param Analogue_volcanoes: Dataframe of the volcanoes with analogue system assigned to them (see assign_A1_classification)
    :param complete_record:
    :param confirmed:
    :param method:
    :return:
    """
    print("calculating mean annual frequency")
    classified_eruptions = pd.merge(complete_record, Analogue_volcanoes, on=['Volcano Number'], how='inner')
    if method == "MeadMagill50":
        small = "Small 50"
        large = "Large 50"
    elif method == "MeadMagill5":
        small = "Small 5"
        large = "Large 5"
    elif method == "MeadMagill95":
        small = "Small 95"
        large = "Large 95"
    else:
        print("Error with method classification: check MeadMagill50, MeadMagill95, MeadMagill5")
    if confirmed == "True":
        Freq_small = classified_eruptions[(classified_eruptions['Eruption Category'] == 'Confirmed Eruption') &
                                          (classified_eruptions['Start Year'] >= classified_eruptions[small]) &
                                          (classified_eruptions['VEI'] < 4)].groupby(['Volcano Number']).count()
        Freq_small = Freq_small.rename(columns={'Volcano type': 'Frequency'})

        Freq_small_volcanoes = classified_eruptions[(classified_eruptions['Eruption Category'] == 'Confirmed Eruption') &
                                                    (classified_eruptions['VEI'] < 4)].groupby(['Volcano Number']).max(small)
        Freq_small_volcanoes.reset_index(level=0, inplace=True)

        Freq_small = pd.merge(Analogue_volcanoes, Freq_small[['Frequency']], on=['Volcano Number'], how='inner')

        Freq_small = pd.merge(Freq_small, Freq_small_volcanoes[['Volcano Number', small]], on=['Volcano Number'],
                              how='inner')
        Freq_small['Small_years'] = year - Freq_small[small]
        Freq_small['Small_rate'] = Freq_small['Frequency'] / Freq_small['Small_years']

        Freq_large = classified_eruptions[(classified_eruptions['Eruption Category'] == 'Confirmed Eruption') &
                                          (classified_eruptions['Start Year'] >= classified_eruptions[large]) &
                                          (classified_eruptions['VEI'] > 3)].groupby(['Volcano Number']).count()
        Freq_large = Freq_large.rename(columns={'Volcano type': 'Frequency'})
        Freq_large_volcanoes = classified_eruptions[
            (classified_eruptions['Eruption Category'] == 'Confirmed Eruption')
            & (classified_eruptions['VEI'] > 3)].groupby(['Volcano Number']).max(large)
        Freq_large_volcanoes.reset_index(level=0, inplace=True)

        Freq_large = pd.merge(Analogue_volcanoes, Freq_large[['Frequency']], on=['Volcano Number'], how='inner')
        Freq_large = pd.merge(Freq_large, Freq_large_volcanoes[['Volcano Number', large]], on=['Volcano Number'],
                              how='inner')
        Freq_large['Large_years'] = year - Freq_large[large]
        Freq_large['Large_rate'] = Freq_large['Frequency'] / Freq_large['Large_years']

        Frequency = pd.merge(Freq_large, Freq_small, on=['Volcano Number'], how='outer')
        Frequency['Small_rate'] = Frequency['Small_rate'].fillna(0)
        Frequency['Large_rate'] = Frequency['Large_rate'].fillna(0)
        Frequency['Eruption_rate'] = Frequency['Large_rate'] + Frequency['Small_rate']
        Frequency_eruptions = Frequency[['Volcano Number', 'Eruption_rate']]
        Frequency_eruptions_type = pd.merge(Frequency_eruptions, Analogue_volcanoes, on=['Volcano Number'],
                                            how='left')
        Frequency_eruptions_type = Frequency_eruptions_type.dropna(subset=['Volcano type'])

        Frequency_eruptions_type.to_csv("frequency_eruptions_test_globe.csv")
        print("saved frequency file")

        Freq_rate = Frequency_eruptions_type.groupby('Volcano type').mean('Eruption_rate')
        Freq_rate = Freq_rate['Eruption_rate']

    elif confirmed == "False":
        print("eruptions are not confirmed in the mean annual frequency")

        Freq_small = classified_eruptions[(classified_eruptions['Start Year'] >= classified_eruptions[small]) &
                                          (classified_eruptions['VEI'] < 4)].groupby(['Volcano Number']).count()
        Freq_small = Freq_small.rename(columns={'Volcano type': 'Frequency'})
        Freq_small_volcanoes = classified_eruptions[(classified_eruptions['VEI'] < 4)].groupby(
            ['Volcano Number']).max(small)
        Freq_small_volcanoes.reset_index(level=0, inplace=True)

        Freq_small = pd.merge(Analogue_volcanoes, Freq_small[['Frequency']], on=['Volcano Number'], how='inner')
        Freq_small = pd.merge(Freq_small, Freq_small_volcanoes[['Volcano Number', small]], on=['Volcano Number'],
                              how='inner')
        Freq_small['Small_years'] = year - Freq_small[small]
        Freq_small['Small_rate'] = Freq_small['Frequency'] / Freq_small['Small_years']

        Freq_large = classified_eruptions[(classified_eruptions['Start Year'] >= classified_eruptions[large]) &
                                          (classified_eruptions['VEI'] > 3)].groupby(['Volcano Number']).count()
        Freq_large = Freq_large.rename(columns={'Volcano type': 'Frequency'})
        Freq_large_volcanoes = classified_eruptions[(classified_eruptions['VEI'] > 3)].groupby(
            ['Volcano Number']).max(large)
        Freq_large_volcanoes.reset_index(level=0, inplace=True)

        Freq_large = pd.merge(Analogue_volcanoes, Freq_large[['Frequency']], on=['Volcano Number'], how='inner')
        Freq_large = pd.merge(Freq_large, Freq_large_volcanoes[['Volcano Number', large]], on=['Volcano Number'],
                              how='inner')
        Freq_large['Large_years'] = year - Freq_large[large]
        Freq_large['Large_rate'] = Freq_large['Frequency'] / Freq_large['Large_years']

        Frequency = pd.merge(Freq_large, Freq_small, on=['Volcano Number'], how='outer')
        Frequency['Small_rate'] = Frequency['Small_rate'].fillna(0)
        Frequency['Large_rate'] = Frequency['Large_rate'].fillna(0)
        Frequency['Eruption_rate'] = Frequency['Large_rate'] + Frequency['Small_rate']
        Frequency_eruptions = Frequency[['Volcano Number', 'Eruption_rate']]
        Frequency_eruptions_type = pd.merge(Frequency_eruptions, Analogue_volcanoes, on=['Volcano Number'],
                                            how='left')
        Frequency_eruptions_type = Frequency_eruptions_type.dropna(subset=['Volcano type'])

        Freq_rate = Frequency_eruptions_type.groupby('Volcano type').mean('Eruption_rate')
        Freq_rate = Freq_rate['Eruption_rate']

    else:
        print("There is a problem")

    return(Freq_rate)

def calc_std_annual_frequency_analogue (Analogue_volcanoes, complete_record, confirmed, method, year):
    classified_eruptions = pd.merge(complete_record, Analogue_volcanoes, on=['Volcano Number'], how='inner')
    saving_eruptions = classified_eruptions.to_csv("classified_eruptions_test2.csv")

    if method == "MeadMagill50":
        small = "Small 50"
        large = "Large 50"
    elif method == "MeadMagill5":
        small = "Small 5"
        large = "Large 5"
    elif method == "MeadMagill95":
        small = "Small 95"
        large = "Large 95"
    else:
        print("Error with method classification: check MeadMagill50, MeadMagill95, MeadMagill5")
    if confirmed == "True":

        Freq_small = classified_eruptions[(classified_eruptions['Eruption Category'] == 'Confirmed Eruption') &
                                          (classified_eruptions['Start Year'] >= classified_eruptions[small]) &
                                          (classified_eruptions['VEI'] < 4)].groupby(['Volcano Number']).count()
        Freq_small = Freq_small.rename(columns={'Volcano type': 'Frequency'})
        Freq_small_volcanoes = classified_eruptions[
            (classified_eruptions['Eruption Category'] == 'Confirmed Eruption')
            & (classified_eruptions['VEI'] < 4)].groupby(['Volcano Number']).max(small)
        Freq_small_volcanoes.reset_index(level=0, inplace=True)

        Freq_small = pd.merge(Analogue_volcanoes, Freq_small[['Frequency']], on=['Volcano Number'], how='inner')
        Freq_small = pd.merge(Freq_small, Freq_small_volcanoes[['Volcano Number', small]], on=['Volcano Number'],
                              how='inner')
        Freq_small['Small_years'] = year - Freq_small[small]
        Freq_small['Small_rate'] = Freq_small['Frequency'] / Freq_small['Small_years']

        Freq_large = classified_eruptions[(classified_eruptions['Eruption Category'] == 'Confirmed Eruption') &
                                          (classified_eruptions['Start Year'] >= classified_eruptions[large]) &
                                          (classified_eruptions['VEI'] > 3)].groupby(['Volcano Number']).count()
        Freq_large = Freq_large.rename(columns={'Volcano type': 'Frequency'})
        Freq_large_volcanoes = classified_eruptions[
            (classified_eruptions['Eruption Category'] == 'Confirmed Eruption')
            & (classified_eruptions['VEI'] > 3)].groupby(['Volcano Number']).max(large)
        Freq_large_volcanoes.reset_index(level=0, inplace=True)

        Freq_large = pd.merge(Analogue_volcanoes, Freq_large[['Frequency']], on=['Volcano Number'], how='inner')
        Freq_large = pd.merge(Freq_large, Freq_large_volcanoes[['Volcano Number', large]], on=['Volcano Number'],
                              how='inner')
        Freq_large['Large_years'] = year - Freq_large[large]
        Freq_large['Large_rate'] = Freq_large['Frequency'] / Freq_large['Large_years']

        Frequency = pd.merge(Freq_large, Freq_small, on=['Volcano Number'], how='outer')
        Frequency['Small_rate'] = Frequency['Small_rate'].fillna(0)
        Frequency['Large_rate'] = Frequency['Large_rate'].fillna(0)
        Frequency['Eruption_rate'] = Frequency['Large_rate'] + Frequency['Small_rate']
        Frequency_eruptions = Frequency[['Volcano Number', 'Eruption_rate']]
        Frequency_eruptions_type = pd.merge(Frequency_eruptions, Analogue_volcanoes, on=['Volcano Number'],
                                            how='left')
        Frequency_eruptions_type = Frequency_eruptions_type.dropna(subset=['Volcano type'])

        Freq_rate_std = Frequency_eruptions_type.groupby('Volcano type')['Eruption_rate'].std()

    else:
        Freq_small = classified_eruptions[(classified_eruptions['Start Year'] >= classified_eruptions[small]) &
                                          (classified_eruptions['VEI'] < 4)].groupby(['Volcano Number']).count()
        Freq_small = Freq_small.rename(columns={'Volcano type': 'Frequency'})
        Freq_small_volcanoes = classified_eruptions[(classified_eruptions['VEI'] < 4)].groupby(
            ['Volcano Number']).max(small)
        Freq_small_volcanoes.reset_index(level=0, inplace=True)

        Freq_small = pd.merge(Analogue_volcanoes, Freq_small[['Frequency']], on=['Volcano Number'], how='inner')
        Freq_small = pd.merge(Freq_small, Freq_small_volcanoes[['Volcano Number', small]], on=['Volcano Number'],
                              how='inner')
        Freq_small['Small_years'] = year - Freq_small[small]
        Freq_small['Small_rate'] = Freq_small['Frequency'] / Freq_small['Small_years']

        Freq_large = classified_eruptions[(classified_eruptions['Start Year'] >= classified_eruptions[large]) &
                                          (classified_eruptions['VEI'] > 3)].groupby(['Volcano Number']).count()
        Freq_large = Freq_large.rename(columns={'Volcano type': 'Frequency'})
        Freq_large_volcanoes = classified_eruptions[(classified_eruptions['VEI'] > 3)].groupby(
            ['Volcano Number']).max(large)
        Freq_large_volcanoes.reset_index(level=0, inplace=True)

        Freq_large = pd.merge(Analogue_volcanoes, Freq_large[['Frequency']], on=['Volcano Number'], how='inner')
        Freq_large = pd.merge(Freq_large, Freq_large_volcanoes[['Volcano Number', large]], on=['Volcano Number'],
                              how='inner')
        Freq_large['Large_years'] = year - Freq_large[large]
        Freq_large['Large_rate'] = Freq_large['Frequency'] / Freq_large['Large_years']

        Frequency = pd.merge(Freq_large, Freq_small, on=['Volcano Number'], how='outer')
        Frequency['Small_rate'] = Frequency['Small_rate'].fillna(0)
        Frequency['Large_rate'] = Frequency['Large_rate'].fillna(0)
        Frequency['Eruption_rate'] = Frequency['Large_rate'] + Frequency['Small_rate']
        Frequency_eruptions = Frequency[['Volcano Number', 'Eruption_rate']]
        Frequency_eruptions_type = pd.merge(Frequency_eruptions, Analogue_volcanoes, on=['Volcano Number'],
                                            how='left')
        Frequency_eruptions_type = Frequency_eruptions_type.dropna(subset=['Volcano type'])

        Freq_rate_std = Frequency_eruptions_type.groupby('Volcano type')['Eruption_rate'].std()

    Freq_rate_std = Freq_rate_std.fillna(0)

    return(Freq_rate_std)

def calc_relative_probability (Analogue_volcanoes, Full_record, confirmed, VEI_schema):
    """
    This function takes a volcano classification system (output from the 'assign_A1/A2_classification' function) and
    the complete volcano record to determine the relative probability for each VEI for each analogue class.
    :param Analogue_volcanoes:
    :param complete_record:
    :param confirmed:
    :return:
    """
    classified_eruptions = pd.merge(Full_record, Analogue_volcanoes, on=['Volcano Number'], how='inner')
    classified_eruptions.to_csv("classified_eruptions_A2.csv")
    print("classified eruptions saved")

    if VEI_schema == 1:
        # VEI <=3, 4, 5, 6, 7
        if confirmed==True:
            All_VEI = classified_eruptions[(classified_eruptions['Eruption Category'] == 'Confirmed Eruption')].groupby('Volcano type').count()
            All_VEI = All_VEI.rename(columns={'Volcano Number': 'All VEI'})
            All_VEI = All_VEI[['All VEI']]
            All_VEI_T = All_VEI.T


            VEI_3 = classified_eruptions[(classified_eruptions.VEI < 4) & (classified_eruptions['Eruption Category'] == 'Confirmed Eruption')].groupby('Volcano type').count()
            VEI_3 = VEI_3.rename(columns={'Volcano Number': 'VEI 3'})
            VEI_3 = VEI_3[['VEI 3']]
            VEI_3_T = VEI_3.T
            VEI_4 = classified_eruptions[(classified_eruptions.VEI == 4) & (classified_eruptions['Eruption Category'] == 'Confirmed Eruption')].groupby('Volcano type').count()
            VEI_4 = VEI_4.rename(columns={'Volcano Number': 'VEI 4'})
            VEI_4 = VEI_4[['VEI 4']]
            VEI_4_T = VEI_4.T
            VEI_5 = classified_eruptions[(classified_eruptions.VEI == 5) & (classified_eruptions['Eruption Category'] == 'Confirmed Eruption')].groupby('Volcano type').count()
            VEI_5 = VEI_5.rename(columns={'Volcano Number': 'VEI 5'})
            VEI_5 = VEI_5[['VEI 5']]
            VEI_5_T = VEI_5.T
            VEI_6 = classified_eruptions[(classified_eruptions.VEI == 6) & (classified_eruptions['Eruption Category'] == 'Confirmed Eruption')].groupby('Volcano type').count()
            VEI_6 = VEI_6.rename(columns={'Volcano Number': 'VEI 6'})
            VEI_6 = VEI_6[['VEI 6']]
            VEI_6_T = VEI_6.T
            VEI_7 = classified_eruptions[(classified_eruptions.VEI == 7) & (classified_eruptions['Eruption Category'] == 'Confirmed Eruption')].groupby('Volcano type').count()
            VEI_7 = VEI_7.rename(columns={'Volcano Number': 'VEI 7'})
            VEI_7 = VEI_7[['VEI 7']]
            VEI_7_T = VEI_7.T
        else:
            All_VEI = classified_eruptions.groupby('Volcano type').count()
            All_VEI = All_VEI.rename(columns={'Volcano Number': 'All VEI'})
            All_VEI = All_VEI[['All VEI']]
            All_VEI_T = All_VEI.T

            VEI_3 = classified_eruptions[(classified_eruptions.VEI < 4)].groupby('Volcano type').count()
            VEI_3 = VEI_3.rename(columns={'Volcano Number': 'VEI 3'})
            VEI_3 = VEI_3[['VEI 3']]
            VEI_3_T = VEI_3.T
            VEI_4 = classified_eruptions[(classified_eruptions.VEI == 4)].groupby('Volcano type').count()
            VEI_4 = VEI_4.rename(columns={'Volcano Number': 'VEI 4'})
            VEI_4 = VEI_4[['VEI 4']]
            VEI_4_T = VEI_4.T
            VEI_5 = classified_eruptions[(classified_eruptions.VEI == 5)].groupby('Volcano type').count()
            VEI_5 = VEI_5.rename(columns={'Volcano Number': 'VEI 5'})
            VEI_5 = VEI_5[['VEI 5']]
            VEI_5_T = VEI_5.T
            VEI_6 = classified_eruptions[(classified_eruptions.VEI == 6)].groupby('Volcano type').count()
            VEI_6 = VEI_6.rename(columns={'Volcano Number': 'VEI 6'})
            VEI_6 = VEI_6[['VEI 6']]
            VEI_6_T = VEI_6.T
            VEI_7 = classified_eruptions[(classified_eruptions.VEI == 7)].groupby('Volcano type').count()
            VEI_7 = VEI_7.rename(columns={'Volcano Number': 'VEI 7'})
            VEI_7 = VEI_7[['VEI 7']]
            VEI_7_T = VEI_7.T

        VEI_frames = [VEI_3_T, VEI_4_T, VEI_5_T, VEI_6_T,VEI_7_T, All_VEI_T]
        VEI_counts = pd.concat(VEI_frames, sort=True)

        VEI_counts.fillna(0, inplace=True)

        VEI_counts_T = VEI_counts.T
        typology_list = list(VEI_counts_T.index.values)
        Total_A = VEI_counts_T.iloc[0]['All VEI']
        Total_B = VEI_counts_T.iloc[1]['All VEI']
        Total_C = VEI_counts_T.iloc[2]['All VEI']
        Total_D = VEI_counts_T.iloc[3]['All VEI']
        Total_E = VEI_counts_T.iloc[4]['All VEI']

        Type_A_VEI_3 = VEI_counts_T.iloc[0]['VEI 3']
        Type_A_VEI_4 = VEI_counts_T.iloc[0]['VEI 4']
        Type_A_VEI_5 = VEI_counts_T.iloc[0]['VEI 5']
        Type_A_VEI_6 = VEI_counts_T.iloc[0]['VEI 6']
        Type_A_VEI_7 = VEI_counts_T.iloc[0]['VEI 7']

        Type_B_VEI_3 = VEI_counts_T.iloc[1]['VEI 3']
        Type_B_VEI_4 = VEI_counts_T.iloc[1]['VEI 4']
        Type_B_VEI_5 = VEI_counts_T.iloc[1]['VEI 5']
        Type_B_VEI_6 = VEI_counts_T.iloc[1]['VEI 6']
        Type_B_VEI_7 = VEI_counts_T.iloc[1]['VEI 7']

        Type_C_VEI_3 = VEI_counts_T.iloc[2]['VEI 3']
        Type_C_VEI_4 = VEI_counts_T.iloc[2]['VEI 4']
        Type_C_VEI_5 = VEI_counts_T.iloc[2]['VEI 5']
        Type_C_VEI_6 = VEI_counts_T.iloc[2]['VEI 6']
        Type_C_VEI_7 = VEI_counts_T.iloc[2]['VEI 7']

        Type_D_VEI_3 = VEI_counts_T.iloc[3]['VEI 3']
        Type_D_VEI_4 = VEI_counts_T.iloc[3]['VEI 4']
        Type_D_VEI_5 = VEI_counts_T.iloc[3]['VEI 5']
        Type_D_VEI_6 = VEI_counts_T.iloc[3]['VEI 6']
        Type_D_VEI_7 = VEI_counts_T.iloc[3]['VEI 7']

        Type_E_VEI_3 = VEI_counts_T.iloc[4]['VEI 3']
        Type_E_VEI_4 = VEI_counts_T.iloc[4]['VEI 4']
        Type_E_VEI_5 = VEI_counts_T.iloc[4]['VEI 5']
        Type_E_VEI_6 = VEI_counts_T.iloc[4]['VEI 6']
        Type_E_VEI_7 = VEI_counts_T.iloc[4]['VEI 7']

        Type_A_VEI_3_cond_prob = Type_A_VEI_3 / Total_A
        Type_A_VEI_4_cond_prob = Type_A_VEI_4 / Total_A
        Type_A_VEI_5_cond_prob = Type_A_VEI_5 / Total_A
        Type_A_VEI_6_cond_prob = Type_A_VEI_6 / Total_A
        Type_A_VEI_7_cond_prob = Type_A_VEI_7 / Total_A

        Type_B_VEI_3_cond_prob = Type_B_VEI_3 / Total_B
        Type_B_VEI_4_cond_prob = Type_B_VEI_4 / Total_B
        Type_B_VEI_5_cond_prob = Type_B_VEI_5 / Total_B
        Type_B_VEI_6_cond_prob = Type_B_VEI_6 / Total_B
        Type_B_VEI_7_cond_prob = Type_B_VEI_7 / Total_B

        Type_C_VEI_3_cond_prob = Type_C_VEI_3 / Total_C
        Type_C_VEI_4_cond_prob = Type_C_VEI_4 / Total_C
        Type_C_VEI_5_cond_prob = Type_C_VEI_5 / Total_C
        Type_C_VEI_6_cond_prob = Type_C_VEI_6 / Total_C
        Type_C_VEI_7_cond_prob = Type_C_VEI_7 / Total_C

        Type_D_VEI_3_cond_prob = Type_D_VEI_3 / Total_D
        Type_D_VEI_4_cond_prob = Type_D_VEI_4 / Total_D
        Type_D_VEI_5_cond_prob = Type_D_VEI_5 / Total_D
        Type_D_VEI_6_cond_prob = Type_D_VEI_6 / Total_D
        Type_D_VEI_7_cond_prob = Type_D_VEI_7 / Total_D

        Type_E_VEI_3_cond_prob = Type_E_VEI_3 / Total_E
        Type_E_VEI_4_cond_prob = Type_E_VEI_4 / Total_E
        Type_E_VEI_5_cond_prob = Type_E_VEI_5 / Total_E
        Type_E_VEI_6_cond_prob = Type_E_VEI_6 / Total_E
        Type_E_VEI_7_cond_prob = Type_E_VEI_7 / Total_E

        Conditional_probabilities = pd.DataFrame(np.array([[Type_A_VEI_3_cond_prob,
                                                            Type_B_VEI_3_cond_prob,
                                                            Type_C_VEI_3_cond_prob,
                                                            Type_D_VEI_3_cond_prob,
                                                            Type_E_VEI_3_cond_prob],
                                                           [Type_A_VEI_4_cond_prob,
                                                            Type_B_VEI_4_cond_prob,
                                                            Type_C_VEI_4_cond_prob,
                                                            Type_D_VEI_4_cond_prob,
                                                            Type_E_VEI_4_cond_prob],
                                                           [Type_A_VEI_5_cond_prob,
                                                            Type_B_VEI_5_cond_prob,
                                                            Type_C_VEI_5_cond_prob,
                                                            Type_D_VEI_5_cond_prob,
                                                            Type_E_VEI_5_cond_prob],
                                                           [Type_A_VEI_6_cond_prob,
                                                            Type_B_VEI_6_cond_prob,
                                                            Type_C_VEI_6_cond_prob,
                                                            Type_D_VEI_6_cond_prob,
                                                            Type_E_VEI_6_cond_prob],
                                                           [Type_A_VEI_7_cond_prob,
                                                            Type_B_VEI_7_cond_prob,
                                                            Type_C_VEI_7_cond_prob,
                                                            Type_D_VEI_7_cond_prob,
                                                            Type_E_VEI_7_cond_prob]]),
                                                 columns=typology_list)

    elif VEI_schema == 2:
        # VEI <= 2, 3, 4, 5, 6, 7
        if confirmed == True:
            All_VEI = classified_eruptions[(classified_eruptions['Eruption Category'] == 'Confirmed Eruption')].groupby(
                'Volcano type').count()
            All_VEI = All_VEI.rename(columns={'Volcano Number': 'All VEI'})
            All_VEI = All_VEI[['All VEI']]
            All_VEI_T = All_VEI.T

            VEI_2 = classified_eruptions[(classified_eruptions.VEI < 3) & (
                        classified_eruptions['Eruption Category'] == 'Confirmed Eruption')].groupby(
                'Volcano type').count()
            VEI_2 = VEI_2.rename(columns={'Volcano Number': 'VEI 2'})
            VEI_2 = VEI_2[['VEI 2']]
            VEI_2_T = VEI_2.T

            VEI_3 = classified_eruptions[(classified_eruptions.VEI == 3) & (
                        classified_eruptions['Eruption Category'] == 'Confirmed Eruption')].groupby(
                'Volcano type').count()
            VEI_3 = VEI_3.rename(columns={'Volcano Number': 'VEI 3'})
            VEI_3 = VEI_3[['VEI 3']]
            VEI_3_T = VEI_3.T
            VEI_4 = classified_eruptions[(classified_eruptions.VEI == 4) & (
                        classified_eruptions['Eruption Category'] == 'Confirmed Eruption')].groupby(
                'Volcano type').count()
            VEI_4 = VEI_4.rename(columns={'Volcano Number': 'VEI 4'})
            VEI_4 = VEI_4[['VEI 4']]
            VEI_4_T = VEI_4.T
            VEI_5 = classified_eruptions[(classified_eruptions.VEI == 5) & (
                        classified_eruptions['Eruption Category'] == 'Confirmed Eruption')].groupby(
                'Volcano type').count()
            VEI_5 = VEI_5.rename(columns={'Volcano Number': 'VEI 5'})
            VEI_5 = VEI_5[['VEI 5']]
            VEI_5_T = VEI_5.T
            VEI_6 = classified_eruptions[(classified_eruptions.VEI == 6) & (
                        classified_eruptions['Eruption Category'] == 'Confirmed Eruption')].groupby(
                'Volcano type').count()
            VEI_6 = VEI_6.rename(columns={'Volcano Number': 'VEI 6'})
            VEI_6 = VEI_6[['VEI 6']]
            VEI_6_T = VEI_6.T
            VEI_7 = classified_eruptions[(classified_eruptions.VEI == 7) & (
                        classified_eruptions['Eruption Category'] == 'Confirmed Eruption')].groupby(
                'Volcano type').count()
            VEI_7 = VEI_7.rename(columns={'Volcano Number': 'VEI 7'})
            VEI_7 = VEI_7[['VEI 7']]
            VEI_7_T = VEI_7.T
        else:
            All_VEI = classified_eruptions.groupby('Volcano type').count()
            All_VEI = All_VEI.rename(columns={'Volcano Number': 'All VEI'})
            All_VEI = All_VEI[['All VEI']]
            All_VEI_T = All_VEI.T

            VEI_2 = classified_eruptions[(classified_eruptions.VEI < 3)].groupby('Volcano type').count()
            VEI_2 = VEI_2.rename(columns={'Volcano Number': 'VEI 2'})
            VEI_2 = VEI_2[['VEI 2']]
            VEI_2_T = VEI_2.T

            VEI_3 = classified_eruptions[(classified_eruptions.VEI == 3)].groupby('Volcano type').count()
            VEI_3 = VEI_3.rename(columns={'Volcano Number': 'VEI 3'})
            VEI_3 = VEI_3[['VEI 3']]
            VEI_3_T = VEI_3.T
            VEI_4 = classified_eruptions[(classified_eruptions.VEI == 4)].groupby('Volcano type').count()
            VEI_4 = VEI_4.rename(columns={'Volcano Number': 'VEI 4'})
            VEI_4 = VEI_4[['VEI 4']]
            VEI_4_T = VEI_4.T
            VEI_5 = classified_eruptions[(classified_eruptions.VEI == 5)].groupby('Volcano type').count()
            VEI_5 = VEI_5.rename(columns={'Volcano Number': 'VEI 5'})
            VEI_5 = VEI_5[['VEI 5']]
            VEI_5_T = VEI_5.T
            VEI_6 = classified_eruptions[(classified_eruptions.VEI == 6)].groupby('Volcano type').count()
            VEI_6 = VEI_6.rename(columns={'Volcano Number': 'VEI 6'})
            VEI_6 = VEI_6[['VEI 6']]
            VEI_6_T = VEI_6.T
            VEI_7 = classified_eruptions[(classified_eruptions.VEI == 7)].groupby('Volcano type').count()
            VEI_7 = VEI_7.rename(columns={'Volcano Number': 'VEI 7'})
            VEI_7 = VEI_7[['VEI 7']]
            VEI_7_T = VEI_7.T

        VEI_frames = [VEI_2_T, VEI_3_T, VEI_4_T, VEI_5_T, VEI_6_T, VEI_7_T, All_VEI_T]
        VEI_counts = pd.concat(VEI_frames, sort=True)

        VEI_counts.fillna(0, inplace=True)

        VEI_counts_T = VEI_counts.T
        typology_list = list(VEI_counts_T.index.values)
        Total_A = VEI_counts_T.iloc[0]['All VEI']
        Total_B = VEI_counts_T.iloc[1]['All VEI']
        Total_C = VEI_counts_T.iloc[2]['All VEI']
        Total_D = VEI_counts_T.iloc[3]['All VEI']
        Total_E = VEI_counts_T.iloc[4]['All VEI']

        Type_A_VEI_2 = VEI_counts_T.iloc[0]['VEI 2']
        Type_A_VEI_3 = VEI_counts_T.iloc[0]['VEI 3']
        Type_A_VEI_4 = VEI_counts_T.iloc[0]['VEI 4']
        Type_A_VEI_5 = VEI_counts_T.iloc[0]['VEI 5']
        Type_A_VEI_6 = VEI_counts_T.iloc[0]['VEI 6']
        Type_A_VEI_7 = VEI_counts_T.iloc[0]['VEI 7']

        Type_B_VEI_2 = VEI_counts_T.iloc[1]['VEI 2']
        Type_B_VEI_3 = VEI_counts_T.iloc[1]['VEI 3']
        Type_B_VEI_4 = VEI_counts_T.iloc[1]['VEI 4']
        Type_B_VEI_5 = VEI_counts_T.iloc[1]['VEI 5']
        Type_B_VEI_6 = VEI_counts_T.iloc[1]['VEI 6']
        Type_B_VEI_7 = VEI_counts_T.iloc[1]['VEI 7']

        Type_C_VEI_2 = VEI_counts_T.iloc[2]['VEI 2']
        Type_C_VEI_3 = VEI_counts_T.iloc[2]['VEI 3']
        Type_C_VEI_4 = VEI_counts_T.iloc[2]['VEI 4']
        Type_C_VEI_5 = VEI_counts_T.iloc[2]['VEI 5']
        Type_C_VEI_6 = VEI_counts_T.iloc[2]['VEI 6']
        Type_C_VEI_7 = VEI_counts_T.iloc[2]['VEI 7']

        Type_D_VEI_2 = VEI_counts_T.iloc[3]['VEI 2']
        Type_D_VEI_3 = VEI_counts_T.iloc[3]['VEI 3']
        Type_D_VEI_4 = VEI_counts_T.iloc[3]['VEI 4']
        Type_D_VEI_5 = VEI_counts_T.iloc[3]['VEI 5']
        Type_D_VEI_6 = VEI_counts_T.iloc[3]['VEI 6']
        Type_D_VEI_7 = VEI_counts_T.iloc[3]['VEI 7']

        Type_E_VEI_2 = VEI_counts_T.iloc[4]['VEI 2']
        Type_E_VEI_3 = VEI_counts_T.iloc[4]['VEI 3']
        Type_E_VEI_4 = VEI_counts_T.iloc[4]['VEI 4']
        Type_E_VEI_5 = VEI_counts_T.iloc[4]['VEI 5']
        Type_E_VEI_6 = VEI_counts_T.iloc[4]['VEI 6']
        Type_E_VEI_7 = VEI_counts_T.iloc[4]['VEI 7']

        Type_A_VEI_2_cond_prob = Type_A_VEI_2 / Total_A
        Type_A_VEI_3_cond_prob = Type_A_VEI_3 / Total_A
        Type_A_VEI_4_cond_prob = Type_A_VEI_4 / Total_A
        Type_A_VEI_5_cond_prob = Type_A_VEI_5 / Total_A
        Type_A_VEI_6_cond_prob = Type_A_VEI_6 / Total_A
        Type_A_VEI_7_cond_prob = Type_A_VEI_7 / Total_A

        Type_B_VEI_2_cond_prob = Type_B_VEI_2 / Total_B
        Type_B_VEI_3_cond_prob = Type_B_VEI_3 / Total_B
        Type_B_VEI_4_cond_prob = Type_B_VEI_4 / Total_B
        Type_B_VEI_5_cond_prob = Type_B_VEI_5 / Total_B
        Type_B_VEI_6_cond_prob = Type_B_VEI_6 / Total_B
        Type_B_VEI_7_cond_prob = Type_B_VEI_7 / Total_B

        Type_C_VEI_2_cond_prob = Type_C_VEI_2 / Total_C
        Type_C_VEI_3_cond_prob = Type_C_VEI_3 / Total_C
        Type_C_VEI_4_cond_prob = Type_C_VEI_4 / Total_C
        Type_C_VEI_5_cond_prob = Type_C_VEI_5 / Total_C
        Type_C_VEI_6_cond_prob = Type_C_VEI_6 / Total_C
        Type_C_VEI_7_cond_prob = Type_C_VEI_7 / Total_C

        Type_D_VEI_2_cond_prob = Type_D_VEI_2 / Total_D
        Type_D_VEI_3_cond_prob = Type_D_VEI_3 / Total_D
        Type_D_VEI_4_cond_prob = Type_D_VEI_4 / Total_D
        Type_D_VEI_5_cond_prob = Type_D_VEI_5 / Total_D
        Type_D_VEI_6_cond_prob = Type_D_VEI_6 / Total_D
        Type_D_VEI_7_cond_prob = Type_D_VEI_7 / Total_D

        Type_E_VEI_2_cond_prob = Type_E_VEI_2 / Total_E
        Type_E_VEI_3_cond_prob = Type_E_VEI_3 / Total_E
        Type_E_VEI_4_cond_prob = Type_E_VEI_4 / Total_E
        Type_E_VEI_5_cond_prob = Type_E_VEI_5 / Total_E
        Type_E_VEI_6_cond_prob = Type_E_VEI_6 / Total_E
        Type_E_VEI_7_cond_prob = Type_E_VEI_7 / Total_E

        Conditional_probabilities = pd.DataFrame(np.array([[Type_A_VEI_2_cond_prob,
                                                            Type_B_VEI_2_cond_prob,
                                                            Type_C_VEI_2_cond_prob,
                                                            Type_D_VEI_2_cond_prob,
                                                            Type_E_VEI_2_cond_prob],
                                                           [Type_A_VEI_3_cond_prob,
                                                            Type_B_VEI_3_cond_prob,
                                                            Type_C_VEI_3_cond_prob,
                                                            Type_D_VEI_3_cond_prob,
                                                            Type_E_VEI_3_cond_prob],
                                                           [Type_A_VEI_4_cond_prob,
                                                            Type_B_VEI_4_cond_prob,
                                                            Type_C_VEI_4_cond_prob,
                                                            Type_D_VEI_4_cond_prob,
                                                            Type_E_VEI_4_cond_prob],
                                                           [Type_A_VEI_5_cond_prob,
                                                            Type_B_VEI_5_cond_prob,
                                                            Type_C_VEI_5_cond_prob,
                                                            Type_D_VEI_5_cond_prob,
                                                            Type_E_VEI_5_cond_prob],
                                                           [Type_A_VEI_6_cond_prob,
                                                            Type_B_VEI_6_cond_prob,
                                                            Type_C_VEI_6_cond_prob,
                                                            Type_D_VEI_6_cond_prob,
                                                            Type_E_VEI_6_cond_prob],
                                                           [Type_A_VEI_7_cond_prob,
                                                            Type_B_VEI_7_cond_prob,
                                                            Type_C_VEI_7_cond_prob,
                                                            Type_D_VEI_7_cond_prob,
                                                            Type_E_VEI_7_cond_prob]]),
                                                 columns=typology_list)

    elif VEI_schema == 3:
        # VEI <=1, 2, 3, 4, 5, 6, 7
        if confirmed == True:
            All_VEI = classified_eruptions[(classified_eruptions['Eruption Category'] == 'Confirmed Eruption')].groupby(
                'Volcano type').count()
            All_VEI = All_VEI.rename(columns={'Volcano Number': 'All VEI'})
            All_VEI = All_VEI[['All VEI']]
            All_VEI_T = All_VEI.T

            VEI_1 = classified_eruptions[(classified_eruptions.VEI < 2) & (
                    classified_eruptions['Eruption Category'] == 'Confirmed Eruption')].groupby(
                'Volcano type').count()
            VEI_1 = VEI_1.rename(columns={'Volcano Number': 'VEI 1'})
            VEI_1 = VEI_1[['VEI 1']]
            VEI_1_T = VEI_1.T

            VEI_2 = classified_eruptions[(classified_eruptions.VEI == 2) & (
                    classified_eruptions['Eruption Category'] == 'Confirmed Eruption')].groupby(
                'Volcano type').count()
            VEI_2 = VEI_2.rename(columns={'Volcano Number': 'VEI 2'})
            VEI_2 = VEI_2[['VEI 2']]
            VEI_2_T = VEI_2.T

            VEI_3 = classified_eruptions[(classified_eruptions.VEI == 3) & (
                    classified_eruptions['Eruption Category'] == 'Confirmed Eruption')].groupby(
                'Volcano type').count()
            VEI_3 = VEI_3.rename(columns={'Volcano Number': 'VEI 3'})
            VEI_3 = VEI_3[['VEI 3']]
            VEI_3_T = VEI_3.T
            VEI_4 = classified_eruptions[(classified_eruptions.VEI == 4) & (
                    classified_eruptions['Eruption Category'] == 'Confirmed Eruption')].groupby(
                'Volcano type').count()
            VEI_4 = VEI_4.rename(columns={'Volcano Number': 'VEI 4'})
            VEI_4 = VEI_4[['VEI 4']]
            VEI_4_T = VEI_4.T
            VEI_5 = classified_eruptions[(classified_eruptions.VEI == 5) & (
                    classified_eruptions['Eruption Category'] == 'Confirmed Eruption')].groupby(
                'Volcano type').count()
            VEI_5 = VEI_5.rename(columns={'Volcano Number': 'VEI 5'})
            VEI_5 = VEI_5[['VEI 5']]
            VEI_5_T = VEI_5.T
            VEI_6 = classified_eruptions[(classified_eruptions.VEI == 6) & (
                    classified_eruptions['Eruption Category'] == 'Confirmed Eruption')].groupby(
                'Volcano type').count()
            VEI_6 = VEI_6.rename(columns={'Volcano Number': 'VEI 6'})
            VEI_6 = VEI_6[['VEI 6']]
            VEI_6_T = VEI_6.T
            VEI_7 = classified_eruptions[(classified_eruptions.VEI == 7) & (
                    classified_eruptions['Eruption Category'] == 'Confirmed Eruption')].groupby(
                'Volcano type').count()
            VEI_7 = VEI_7.rename(columns={'Volcano Number': 'VEI 7'})
            VEI_7 = VEI_7[['VEI 7']]
            VEI_7_T = VEI_7.T
        else:
            All_VEI = classified_eruptions.groupby('Volcano type').count()
            All_VEI = All_VEI.rename(columns={'Volcano Number': 'All VEI'})
            All_VEI = All_VEI[['All VEI']]
            All_VEI_T = All_VEI.T

            VEI_1 = classified_eruptions[(classified_eruptions.VEI < 2)].groupby('Volcano type').count()
            VEI_1 = VEI_1.rename(columns={'Volcano Number': 'VEI 1'})
            VEI_1 = VEI_1[['VEI 1']]
            VEI_1_T = VEI_1.T

            VEI_2 = classified_eruptions[(classified_eruptions.VEI == 2)].groupby('Volcano type').count()
            VEI_2 = VEI_2.rename(columns={'Volcano Number': 'VEI 2'})
            VEI_2 = VEI_2[['VEI 2']]
            VEI_2_T = VEI_2.T

            VEI_3 = classified_eruptions[(classified_eruptions.VEI == 3)].groupby('Volcano type').count()
            VEI_3 = VEI_3.rename(columns={'Volcano Number': 'VEI 3'})
            VEI_3 = VEI_3[['VEI 3']]
            VEI_3_T = VEI_3.T
            VEI_4 = classified_eruptions[(classified_eruptions.VEI == 4)].groupby('Volcano type').count()
            VEI_4 = VEI_4.rename(columns={'Volcano Number': 'VEI 4'})
            VEI_4 = VEI_4[['VEI 4']]
            VEI_4_T = VEI_4.T
            VEI_5 = classified_eruptions[(classified_eruptions.VEI == 5)].groupby('Volcano type').count()
            VEI_5 = VEI_5.rename(columns={'Volcano Number': 'VEI 5'})
            VEI_5 = VEI_5[['VEI 5']]
            VEI_5_T = VEI_5.T
            VEI_6 = classified_eruptions[(classified_eruptions.VEI == 6)].groupby('Volcano type').count()
            VEI_6 = VEI_6.rename(columns={'Volcano Number': 'VEI 6'})
            VEI_6 = VEI_6[['VEI 6']]
            VEI_6_T = VEI_6.T
            VEI_7 = classified_eruptions[(classified_eruptions.VEI == 7)].groupby('Volcano type').count()
            VEI_7 = VEI_7.rename(columns={'Volcano Number': 'VEI 7'})
            VEI_7 = VEI_7[['VEI 7']]
            VEI_7_T = VEI_7.T

        VEI_frames = [VEI_1_T, VEI_2_T, VEI_3_T, VEI_4_T, VEI_5_T, VEI_6_T, VEI_7_T, All_VEI_T]
        VEI_counts = pd.concat(VEI_frames, sort=True)

        VEI_counts.fillna(0, inplace=True)

        VEI_counts_T = VEI_counts.T
        typology_list = list(VEI_counts_T.index.values)
        Total_A = VEI_counts_T.iloc[0]['All VEI']
        Total_B = VEI_counts_T.iloc[1]['All VEI']
        Total_C = VEI_counts_T.iloc[2]['All VEI']
        Total_D = VEI_counts_T.iloc[3]['All VEI']
        Total_E = VEI_counts_T.iloc[4]['All VEI']

        Type_A_VEI_1 = VEI_counts_T.iloc[0]['VEI 1']
        Type_A_VEI_2 = VEI_counts_T.iloc[0]['VEI 2']
        Type_A_VEI_3 = VEI_counts_T.iloc[0]['VEI 3']
        Type_A_VEI_4 = VEI_counts_T.iloc[0]['VEI 4']
        Type_A_VEI_5 = VEI_counts_T.iloc[0]['VEI 5']
        Type_A_VEI_6 = VEI_counts_T.iloc[0]['VEI 6']
        Type_A_VEI_7 = VEI_counts_T.iloc[0]['VEI 7']

        Type_B_VEI_1 = VEI_counts_T.iloc[1]['VEI 1']
        Type_B_VEI_2 = VEI_counts_T.iloc[1]['VEI 2']
        Type_B_VEI_3 = VEI_counts_T.iloc[1]['VEI 3']
        Type_B_VEI_4 = VEI_counts_T.iloc[1]['VEI 4']
        Type_B_VEI_5 = VEI_counts_T.iloc[1]['VEI 5']
        Type_B_VEI_6 = VEI_counts_T.iloc[1]['VEI 6']
        Type_B_VEI_7 = VEI_counts_T.iloc[1]['VEI 7']

        Type_C_VEI_1 = VEI_counts_T.iloc[2]['VEI 1']
        Type_C_VEI_2 = VEI_counts_T.iloc[2]['VEI 2']
        Type_C_VEI_3 = VEI_counts_T.iloc[2]['VEI 3']
        Type_C_VEI_4 = VEI_counts_T.iloc[2]['VEI 4']
        Type_C_VEI_5 = VEI_counts_T.iloc[2]['VEI 5']
        Type_C_VEI_6 = VEI_counts_T.iloc[2]['VEI 6']
        Type_C_VEI_7 = VEI_counts_T.iloc[2]['VEI 7']

        Type_D_VEI_1 = VEI_counts_T.iloc[3]['VEI 1']
        Type_D_VEI_2 = VEI_counts_T.iloc[3]['VEI 2']
        Type_D_VEI_3 = VEI_counts_T.iloc[3]['VEI 3']
        Type_D_VEI_4 = VEI_counts_T.iloc[3]['VEI 4']
        Type_D_VEI_5 = VEI_counts_T.iloc[3]['VEI 5']
        Type_D_VEI_6 = VEI_counts_T.iloc[3]['VEI 6']
        Type_D_VEI_7 = VEI_counts_T.iloc[3]['VEI 7']

        Type_E_VEI_1 = VEI_counts_T.iloc[4]['VEI 1']
        Type_E_VEI_2 = VEI_counts_T.iloc[4]['VEI 2']
        Type_E_VEI_3 = VEI_counts_T.iloc[4]['VEI 3']
        Type_E_VEI_4 = VEI_counts_T.iloc[4]['VEI 4']
        Type_E_VEI_5 = VEI_counts_T.iloc[4]['VEI 5']
        Type_E_VEI_6 = VEI_counts_T.iloc[4]['VEI 6']
        Type_E_VEI_7 = VEI_counts_T.iloc[4]['VEI 7']

        Type_A_VEI_1_cond_prob = Type_A_VEI_1 / Total_A
        Type_A_VEI_2_cond_prob = Type_A_VEI_2 / Total_A
        Type_A_VEI_3_cond_prob = Type_A_VEI_3 / Total_A
        Type_A_VEI_4_cond_prob = Type_A_VEI_4 / Total_A
        Type_A_VEI_5_cond_prob = Type_A_VEI_5 / Total_A
        Type_A_VEI_6_cond_prob = Type_A_VEI_6 / Total_A
        Type_A_VEI_7_cond_prob = Type_A_VEI_7 / Total_A

        Type_B_VEI_1_cond_prob = Type_B_VEI_1 / Total_B
        Type_B_VEI_2_cond_prob = Type_B_VEI_2 / Total_B
        Type_B_VEI_3_cond_prob = Type_B_VEI_3 / Total_B
        Type_B_VEI_4_cond_prob = Type_B_VEI_4 / Total_B
        Type_B_VEI_5_cond_prob = Type_B_VEI_5 / Total_B
        Type_B_VEI_6_cond_prob = Type_B_VEI_6 / Total_B
        Type_B_VEI_7_cond_prob = Type_B_VEI_7 / Total_B

        Type_C_VEI_1_cond_prob = Type_C_VEI_1 / Total_C
        Type_C_VEI_2_cond_prob = Type_C_VEI_2 / Total_C
        Type_C_VEI_3_cond_prob = Type_C_VEI_3 / Total_C
        Type_C_VEI_4_cond_prob = Type_C_VEI_4 / Total_C
        Type_C_VEI_5_cond_prob = Type_C_VEI_5 / Total_C
        Type_C_VEI_6_cond_prob = Type_C_VEI_6 / Total_C
        Type_C_VEI_7_cond_prob = Type_C_VEI_7 / Total_C

        Type_D_VEI_1_cond_prob = Type_D_VEI_1 / Total_D
        Type_D_VEI_2_cond_prob = Type_D_VEI_2 / Total_D
        Type_D_VEI_3_cond_prob = Type_D_VEI_3 / Total_D
        Type_D_VEI_4_cond_prob = Type_D_VEI_4 / Total_D
        Type_D_VEI_5_cond_prob = Type_D_VEI_5 / Total_D
        Type_D_VEI_6_cond_prob = Type_D_VEI_6 / Total_D
        Type_D_VEI_7_cond_prob = Type_D_VEI_7 / Total_D

        Type_E_VEI_1_cond_prob = Type_E_VEI_1 / Total_E
        Type_E_VEI_2_cond_prob = Type_E_VEI_2 / Total_E
        Type_E_VEI_3_cond_prob = Type_E_VEI_3 / Total_E
        Type_E_VEI_4_cond_prob = Type_E_VEI_4 / Total_E
        Type_E_VEI_5_cond_prob = Type_E_VEI_5 / Total_E
        Type_E_VEI_6_cond_prob = Type_E_VEI_6 / Total_E
        Type_E_VEI_7_cond_prob = Type_E_VEI_7 / Total_E

        Conditional_probabilities = pd.DataFrame(np.array([[Type_A_VEI_1_cond_prob,
                                                            Type_B_VEI_1_cond_prob,
                                                            Type_C_VEI_1_cond_prob,
                                                            Type_D_VEI_1_cond_prob,
                                                            Type_E_VEI_1_cond_prob],
                                                           [Type_A_VEI_2_cond_prob,
                                                            Type_B_VEI_2_cond_prob,
                                                            Type_C_VEI_2_cond_prob,
                                                            Type_D_VEI_2_cond_prob,
                                                            Type_E_VEI_2_cond_prob],
                                                           [Type_A_VEI_3_cond_prob,
                                                            Type_B_VEI_3_cond_prob,
                                                            Type_C_VEI_3_cond_prob,
                                                            Type_D_VEI_3_cond_prob,
                                                            Type_E_VEI_3_cond_prob],
                                                           [Type_A_VEI_4_cond_prob,
                                                            Type_B_VEI_4_cond_prob,
                                                            Type_C_VEI_4_cond_prob,
                                                            Type_D_VEI_4_cond_prob,
                                                            Type_E_VEI_4_cond_prob],
                                                           [Type_A_VEI_5_cond_prob,
                                                            Type_B_VEI_5_cond_prob,
                                                            Type_C_VEI_5_cond_prob,
                                                            Type_D_VEI_5_cond_prob,
                                                            Type_E_VEI_5_cond_prob],
                                                           [Type_A_VEI_6_cond_prob,
                                                            Type_B_VEI_6_cond_prob,
                                                            Type_C_VEI_6_cond_prob,
                                                            Type_D_VEI_6_cond_prob,
                                                            Type_E_VEI_6_cond_prob],
                                                           [Type_A_VEI_7_cond_prob,
                                                            Type_B_VEI_7_cond_prob,
                                                            Type_C_VEI_7_cond_prob,
                                                            Type_D_VEI_7_cond_prob,
                                                            Type_E_VEI_7_cond_prob]]),
                                                 columns=typology_list)



    return(Conditional_probabilities)

def get_VEI_eruptions_count (volc_num, complete_record):
    VEI_eruptions = complete_record[(complete_record['Volcano Number'] == volc_num)
                                    & (complete_record['Eruption Category'] == 'Confirmed Eruption')
                                    & (complete_record['VEI'] >=1)]
    VEI_eruptions_count = len(VEI_eruptions)

    return(VEI_eruptions_count)

def get_observed_eruption_rate (volc_num, complete_record, confirmed, year, Full_record, method):
    GVP = pd.read_csv("GVPDB2019.csv")
    if method == "MeadMagill50":
        small = "Small 50"
        large = "Large 50"
    elif method == "MeadMagill5":
        small = "Small 5"
        large = "Large 5"
    elif method == "MeadMagill95":
        small = "Small 95"
        large = "Large 95"
    else:
        print("Error with method classification: check MeadMagill50, MeadMagill95, MeadMagill5")
    if confirmed == "True":

        complete_record_length_small = next(iter(complete_record.loc[complete_record['Volcano Number'] == volc_num, small]), 'no match for volc_num')
        small_eruptions = complete_record[(complete_record['Volcano Number'] == volc_num)
                                          & (complete_record['Eruption Category'] == 'Confirmed Eruption')
                                          & (complete_record['VEI'] < 4)]
        number_of_small_eruptions = len(small_eruptions)
        if number_of_small_eruptions >0:
            rate_small_eruptions = number_of_small_eruptions / (year - complete_record_length_small)
        else:
            rate_small_eruptions = 0

        complete_record_length_large = next(iter(complete_record.loc[complete_record['Volcano Number'] == volc_num, large]), 'no match for volc_num')
        large_eruptions = complete_record[(complete_record['Volcano Number'] == volc_num)
                                          & (complete_record['Eruption Category'] == 'Confirmed Eruption')
                                          & (complete_record['VEI'] > 3)]
        number_of_large_eruptions = len(large_eruptions)
        if number_of_large_eruptions >0:
            rate_large_eruptions = number_of_large_eruptions / (year - complete_record_length_large)
        else:
            rate_large_eruptions = 0

        unknown_eruptions = GVP[(GVP['Volcano Number'] == volc_num)
                                          & (GVP['Eruption Category'] == 'Confirmed Eruption')
                                          & (GVP['VEI'].isna())]
        number_of_unknown_eruptions = len(unknown_eruptions)
        if number_of_unknown_eruptions > 0:
            unknown_record_length = unknown_eruptions['Start Year'].min()
            rate_unknown_eruptions = number_of_unknown_eruptions / (year - unknown_record_length)
        else:
            rate_unknown_eruptions = 0

        rate_of_eruptions = rate_small_eruptions + rate_large_eruptions + rate_unknown_eruptions

    else:
        complete_record_length_small = next(iter(complete_record.loc[complete_record['Volcano Number'] == volc_num, small]), 'no match for volc_num')
        small_eruptions = complete_record[(complete_record['Volcano Number'] == volc_num)
                                          & (complete_record['VEI'] < 4)]
        number_of_small_eruptions = len(small_eruptions)
        if number_of_small_eruptions >0:
            rate_small_eruptions = number_of_small_eruptions / (year - complete_record_length_small)
        else:
            rate_small_eruptions = 0

        complete_record_length_large = next(iter(complete_record.loc[complete_record['Volcano Number'] == volc_num, large]), 'no match for volc_num')
        large_eruptions = complete_record[(complete_record['Volcano Number'] == volc_num)
                                          & (complete_record['VEI'] > 3)]
        number_of_large_eruptions = len(large_eruptions)
        if number_of_large_eruptions >0:
            rate_large_eruptions = number_of_large_eruptions / (year - complete_record_length_large)
        else:
            rate_large_eruptions = 0

        unknown_eruptions = complete_record[(complete_record['Volcano Number'] == volc_num)
                                          & (complete_record['VEI'].isna())]
        number_of_unknown_eruptions = len(unknown_eruptions)
        if number_of_unknown_eruptions > 0:
            rate_unknown_eruptions = number_of_large_eruptions / (year - complete_record_length_large)
        else:
            rate_unknown_eruptions = 0

        rate_of_eruptions = rate_small_eruptions + rate_large_eruptions + rate_unknown_eruptions

    return(rate_of_eruptions)

def get_observed_eruption_rate_full_record (volc_num, year, confirmed, Full_record):
    Full_record = Full_record
    if confirmed == "True":
        Eruption_record = Full_record[(Full_record['Volcano Number'] == volc_num)
                                          & (Full_record['Eruption Category'] == 'Confirmed Eruption')]
        Eruption_record_length = Eruption_record['Start Year'].min()
        Number_of_eruptions = len(Eruption_record)
        if Number_of_eruptions > 0:
            rate_of_eruptions_full_record = Number_of_eruptions / (year - Eruption_record_length)
        else:
            rate_of_eruptions_full_record = 0

    else:
        Eruption_record = Full_record[(Full_record['Volcano Number'] == volc_num)]
        Eruption_record_length = Eruption_record['Start Year'].min()
        Number_of_eruptions = len(Eruption_record)
        if Number_of_eruptions > 0:
            rate_of_eruptions_full_record = Number_of_eruptions / (year - Eruption_record_length)
        else:
            rate_of_eruptions_full_record = 0

    return (rate_of_eruptions_full_record)

def get_observed_relative_frequency (volc_num, complete_record, confirmed, Full_record, VEI_schema):

    if VEI_schema == 1:

        if confirmed == "True":
           VEI_3_or_less = complete_record[(complete_record['Volcano Number'] == volc_num)
                                           & (complete_record['Eruption Category'] == 'Confirmed Eruption')
                                           & (complete_record['VEI'] < 4)]
           Number_VEI_3_or_less = len(VEI_3_or_less)

           VEI_4 = complete_record[(complete_record['Volcano Number'] == volc_num)
                                           & (complete_record['Eruption Category'] == 'Confirmed Eruption')
                                           & (complete_record['VEI'] == 4)]
           Number_VEI_4 = len(VEI_4)

           VEI_5 = complete_record[(complete_record['Volcano Number'] == volc_num)
                                           & (complete_record['Eruption Category'] == 'Confirmed Eruption')
                                           & (complete_record['VEI'] == 5)]
           Number_VEI_5 = len(VEI_5)

           VEI_6 = complete_record[(complete_record['Volcano Number'] == volc_num)
                                           & (complete_record['Eruption Category'] == 'Confirmed Eruption')
                                           & (complete_record['VEI'] == 6)]
           Number_VEI_6 = len(VEI_6)

           VEI_7 = complete_record[(complete_record['Volcano Number'] == volc_num)
                                           & (complete_record['Eruption Category'] == 'Confirmed Eruption')
                                           & (complete_record['VEI'] == 7)]
           Number_VEI_7 = len(VEI_7)
        else:
            VEI_3_or_less = complete_record[(complete_record['Volcano Number'] == volc_num)
                                            & (complete_record['VEI'] < 4)]
            Number_VEI_3_or_less = len(VEI_3_or_less)

            VEI_4 = complete_record[(complete_record['Volcano Number'] == volc_num)
                                    & (complete_record['VEI'] == 4)]
            Number_VEI_4 = len(VEI_4)

            VEI_5 = complete_record[(complete_record['Volcano Number'] == volc_num)
                                    & (complete_record['VEI'] == 5)]
            Number_VEI_5 = len(VEI_5)

            VEI_6 = complete_record[(complete_record['Volcano Number'] == volc_num)
                                    & (complete_record['VEI'] == 6)]
            Number_VEI_6 = len(VEI_6)

            VEI_7 = complete_record[(complete_record['Volcano Number'] == volc_num)
                                    & (complete_record['VEI'] == 7)]
            Number_VEI_7 = len(VEI_7)

        Observed_VEI = [Number_VEI_3_or_less, Number_VEI_4, Number_VEI_5, Number_VEI_6, Number_VEI_7]

    elif VEI_schema == 2:
        #VEI <=2, 3, 4, 5, 6, 7

        if confirmed == "True":
            VEI_2_or_less = complete_record[(complete_record['Volcano Number'] == volc_num)
                                            & (complete_record['Eruption Category'] == 'Confirmed Eruption')
                                            & (complete_record['VEI'] < 3)]
            Number_VEI_2_or_less = len(VEI_2_or_less)

            VEI_3 = complete_record[(complete_record['Volcano Number'] == volc_num)
                                            & (complete_record['Eruption Category'] == 'Confirmed Eruption')
                                            & (complete_record['VEI'] == 3)]
            Number_VEI_3 = len(VEI_3)

            VEI_4 = complete_record[(complete_record['Volcano Number'] == volc_num)
                                    & (complete_record['Eruption Category'] == 'Confirmed Eruption')
                                    & (complete_record['VEI'] == 4)]
            Number_VEI_4 = len(VEI_4)

            VEI_5 = complete_record[(complete_record['Volcano Number'] == volc_num)
                                    & (complete_record['Eruption Category'] == 'Confirmed Eruption')
                                    & (complete_record['VEI'] == 5)]
            Number_VEI_5 = len(VEI_5)

            VEI_6 = complete_record[(complete_record['Volcano Number'] == volc_num)
                                    & (complete_record['Eruption Category'] == 'Confirmed Eruption')
                                    & (complete_record['VEI'] == 6)]
            Number_VEI_6 = len(VEI_6)

            VEI_7 = complete_record[(complete_record['Volcano Number'] == volc_num)
                                    & (complete_record['Eruption Category'] == 'Confirmed Eruption')
                                    & (complete_record['VEI'] == 7)]
            Number_VEI_7 = len(VEI_7)
        else:
            VEI_2_or_less = complete_record[(complete_record['Volcano Number'] == volc_num)
                                            & (complete_record['VEI'] < 3)]
            Number_VEI_2_or_less = len(VEI_2_or_less)

            VEI_3 = complete_record[(complete_record['Volcano Number'] == volc_num)
                                            & (complete_record['VEI'] == 3)]
            Number_VEI_3 = len(VEI_3)

            VEI_4 = complete_record[(complete_record['Volcano Number'] == volc_num)
                                    & (complete_record['VEI'] == 4)]
            Number_VEI_4 = len(VEI_4)

            VEI_5 = complete_record[(complete_record['Volcano Number'] == volc_num)
                                    & (complete_record['VEI'] == 5)]
            Number_VEI_5 = len(VEI_5)

            VEI_6 = complete_record[(complete_record['Volcano Number'] == volc_num)
                                    & (complete_record['VEI'] == 6)]
            Number_VEI_6 = len(VEI_6)

            VEI_7 = complete_record[(complete_record['Volcano Number'] == volc_num)
                                    & (complete_record['VEI'] == 7)]
            Number_VEI_7 = len(VEI_7)

        Observed_VEI = [Number_VEI_2_or_less, Number_VEI_3, Number_VEI_4, Number_VEI_5, Number_VEI_6, Number_VEI_7]

    elif VEI_schema == 3:
        # VEI <=1, 2, 3, 4, 5, 6, 7

        if confirmed == "True":
            VEI_1_or_less = complete_record[(complete_record['Volcano Number'] == volc_num)
                                            & (complete_record['Eruption Category'] == 'Confirmed Eruption')
                                            & (complete_record['VEI'] < 2)]
            Number_VEI_1_or_less = len(VEI_1_or_less)

            VEI_2 = complete_record[(complete_record['Volcano Number'] == volc_num)
                                            & (complete_record['Eruption Category'] == 'Confirmed Eruption')
                                            & (complete_record['VEI'] == 2)]
            Number_VEI_2 = len(VEI_2)

            VEI_3 = complete_record[(complete_record['Volcano Number'] == volc_num)
                                    & (complete_record['Eruption Category'] == 'Confirmed Eruption')
                                    & (complete_record['VEI'] == 3)]
            Number_VEI_3 = len(VEI_3)

            VEI_4 = complete_record[(complete_record['Volcano Number'] == volc_num)
                                    & (complete_record['Eruption Category'] == 'Confirmed Eruption')
                                    & (complete_record['VEI'] == 4)]
            Number_VEI_4 = len(VEI_4)

            VEI_5 = complete_record[(complete_record['Volcano Number'] == volc_num)
                                    & (complete_record['Eruption Category'] == 'Confirmed Eruption')
                                    & (complete_record['VEI'] == 5)]
            Number_VEI_5 = len(VEI_5)

            VEI_6 = complete_record[(complete_record['Volcano Number'] == volc_num)
                                    & (complete_record['Eruption Category'] == 'Confirmed Eruption')
                                    & (complete_record['VEI'] == 6)]
            Number_VEI_6 = len(VEI_6)

            VEI_7 = complete_record[(complete_record['Volcano Number'] == volc_num)
                                    & (complete_record['Eruption Category'] == 'Confirmed Eruption')
                                    & (complete_record['VEI'] == 7)]
            Number_VEI_7 = len(VEI_7)
        else:
            VEI_1_or_less = complete_record[(complete_record['Volcano Number'] == volc_num)
                                            & (complete_record['VEI'] < 2)]
            Number_VEI_1_or_less = len(VEI_1_or_less)

            VEI_2 = complete_record[(complete_record['Volcano Number'] == volc_num)
                                            & (complete_record['VEI'] == 2)]
            Number_VEI_2 = len(VEI_2)

            VEI_3 = complete_record[(complete_record['Volcano Number'] == volc_num)
                                    & (complete_record['VEI'] == 3)]
            Number_VEI_3 = len(VEI_3)

            VEI_4 = complete_record[(complete_record['Volcano Number'] == volc_num)
                                    & (complete_record['VEI'] == 4)]
            Number_VEI_4 = len(VEI_4)

            VEI_5 = complete_record[(complete_record['Volcano Number'] == volc_num)
                                    & (complete_record['VEI'] == 5)]
            Number_VEI_5 = len(VEI_5)

            VEI_6 = complete_record[(complete_record['Volcano Number'] == volc_num)
                                    & (complete_record['VEI'] == 6)]
            Number_VEI_6 = len(VEI_6)

            VEI_7 = complete_record[(complete_record['Volcano Number'] == volc_num)
                                    & (complete_record['VEI'] == 7)]
            Number_VEI_7 = len(VEI_7)

        Observed_VEI = [Number_VEI_1_or_less, Number_VEI_2, Number_VEI_3, Number_VEI_4, Number_VEI_5, Number_VEI_6, Number_VEI_7]

    return(Observed_VEI)

#Fix this later - add in the VEI schema
def get_observed_relative_frequency_full_record (volc_num, Full_record, confirmed):
    if confirmed=="True":
       VEI_3_or_less = Full_record[(Full_record['Volcano Number'] == volc_num)
                                       & (Full_record['Eruption Category'] == 'Confirmed Eruption')
                                       & (Full_record['VEI'] < 4)]
       Number_VEI_3_or_less = len(VEI_3_or_less)

       VEI_4 = Full_record[(Full_record['Volcano Number'] == volc_num)
                                       & (Full_record['Eruption Category'] == 'Confirmed Eruption')
                                       & (Full_record['VEI'] == 4)]
       Number_VEI_4 = len(VEI_4)

       VEI_5 = Full_record[(Full_record['Volcano Number'] == volc_num)
                                       & (Full_record['Eruption Category'] == 'Confirmed Eruption')
                                       & (Full_record['VEI'] == 5)]
       Number_VEI_5 = len(VEI_5)

       VEI_6 = Full_record[(Full_record['Volcano Number'] == volc_num)
                                       & (Full_record['Eruption Category'] == 'Confirmed Eruption')
                                       & (Full_record['VEI'] == 6)]
       Number_VEI_6 = len(VEI_6)

       VEI_7 = Full_record[(Full_record['Volcano Number'] == volc_num)
                                       & (Full_record['Eruption Category'] == 'Confirmed Eruption')
                                       & (Full_record['VEI'] == 7)]
       Number_VEI_7 = len(VEI_7)
    else:
        VEI_3_or_less = Full_record[(Full_record['Volcano Number'] == volc_num)
                                        & (Full_record['VEI'] < 4)]
        Number_VEI_3_or_less = len(VEI_3_or_less)

        VEI_4 = Full_record[(Full_record['Volcano Number'] == volc_num)
                                & (Full_record['VEI'] == 4)]
        Number_VEI_4 = len(VEI_4)

        VEI_5 = Full_record[(Full_record['Volcano Number'] == volc_num)
                                & (Full_record['VEI'] == 5)]
        Number_VEI_5 = len(VEI_5)

        VEI_6 = Full_record[(Full_record['Volcano Number'] == volc_num)
                                & (Full_record['VEI'] == 6)]
        Number_VEI_6 = len(VEI_6)

        VEI_7 = Full_record[(Full_record['Volcano Number'] == volc_num)
                                & (Full_record['VEI'] == 7)]
        Number_VEI_7 = len(VEI_7)

    Observed_VEI_full_record = [Number_VEI_3_or_less, Number_VEI_4, Number_VEI_5, Number_VEI_6, Number_VEI_7]

    return(Observed_VEI_full_record)

def get_freq_mag_dirichlet (analogue_freq_mag, Analogue_volcanoes, volc_num, Power_law, VEI_schema):

    get_volcano_type = Analogue_volcanoes.loc[Analogue_volcanoes['Volcano Number'] == volc_num, 'Volcano type'].iloc[0]

    find_FM = analogue_freq_mag[get_volcano_type]*100

    if VEI_schema == 1:

        if find_FM[4] > 0:
            print("Dirichlet 1")
            freq_mag = find_FM.values
            freq_mag = freq_mag*100
            size = (int(find_FM.sum())*100)
            Freq_mag_samples = stats.dirichlet.rvs(freq_mag, size=size, random_state=None)
            Freq_mag_samples_df = pd.DataFrame({'3': Freq_mag_samples[:, 0],
                                    '4': Freq_mag_samples[:, 1],
                                    '5': Freq_mag_samples[:, 2],
                                    '6': Freq_mag_samples[:, 3],
                                    '7': Freq_mag_samples[:, 4]})

            # Frequency_VEI_3 = Freq_mag_samples_df['3'].mean()
            # Frequency_VEI_4 = Freq_mag_samples_df['4'].mean()
            # Frequency_VEI_5 = Freq_mag_samples_df['5'].mean()
            # Frequency_VEI_6 = Freq_mag_samples_df['6'].mean()
            # Frequency_VEI_7 = Freq_mag_samples_df['7'].mean()
            #
            # Frequency_VEI_3_std = Freq_mag_samples_df['3'].std()
            # Frequency_VEI_4_std = Freq_mag_samples_df['4'].std()
            # Frequency_VEI_5_std = Freq_mag_samples_df['5'].std()
            # Frequency_VEI_6_std = Freq_mag_samples_df['6'].std()
            # Frequency_VEI_7_std = Freq_mag_samples_df['7'].std()

        elif find_FM[3] > 0 and find_FM[4] == 0:
            print("Dirichlet 2")
            freq_mag_ = find_FM.iloc[0:4]
            freq_mag = freq_mag_.values
            freq_mag = freq_mag*10
            size = (int(find_FM.sum())*10)

            if Power_law == "True":
                Freq_mag_samples = stats.dirichlet.rvs(freq_mag, size=size, random_state=None)
                Freq_mag_samples_df = pd.DataFrame({'3': Freq_mag_samples[:, 0],
                                                    '4': Freq_mag_samples[:, 1],
                                                    '5': Freq_mag_samples[:, 2],
                                                    '6': Freq_mag_samples[:, 3]})
                Freq_mag_samples_df['7'] = Freq_mag_samples_df['6'] * 0.1

                # Frequency_VEI_3 = Freq_mag_samples_df['3'].mean()
                # Frequency_VEI_4 = Freq_mag_samples_df['4'].mean()
                # Frequency_VEI_5 = Freq_mag_samples_df['5'].mean()
                # Frequency_VEI_6 = Freq_mag_samples_df['6'].mean()
                # Frequency_VEI_7 = Frequency_VEI_6*0.1
                #
                # Frequency_VEI_3_std = Freq_mag_samples_df['3'].std()
                # Frequency_VEI_4_std = Freq_mag_samples_df['4'].std()
                # Frequency_VEI_5_std = Freq_mag_samples_df['5'].std()
                # Frequency_VEI_6_std = Freq_mag_samples_df['6'].std()
                # Frequency_VEI_7_std = Freq_mag_samples_df['6'].std()
            else:
                Freq_mag_samples = stats.dirichlet.rvs(freq_mag, size=size, random_state=None)
                Freq_mag_samples_df = pd.DataFrame({'3': Freq_mag_samples[:, 0],
                                                    '4': Freq_mag_samples[:, 1],
                                                    '5': Freq_mag_samples[:, 2],
                                                    '6': Freq_mag_samples[:, 3]})
                Freq_mag_samples_df['7'] = 0

                # Frequency_VEI_3 = Freq_mag_samples_df['3'].mean()
                # Frequency_VEI_4 = Freq_mag_samples_df['4'].mean()
                # Frequency_VEI_5 = Freq_mag_samples_df['5'].mean()
                # Frequency_VEI_6 = Freq_mag_samples_df['6'].mean()
                # Frequency_VEI_7 = 0
                #
                # Frequency_VEI_3_std = Freq_mag_samples_df['3'].std()
                # Frequency_VEI_4_std = Freq_mag_samples_df['4'].std()
                # Frequency_VEI_5_std = Freq_mag_samples_df['5'].std()
                # Frequency_VEI_6_std = Freq_mag_samples_df['6'].std()
                # Frequency_VEI_7_std = 0

        elif find_FM[2] > 0 and find_FM[3] == 0:
            print("Dirichlet 3")
            freq_mag_ = find_FM.iloc[0:3]
            freq_mag = freq_mag_.values
            freq_mag = freq_mag*10
            size = (int(find_FM.sum())*10)

            if Power_law == "True":
                Freq_mag_samples = stats.dirichlet.rvs(freq_mag, size=size, random_state=None)
                Freq_mag_samples_df = pd.DataFrame({'3': Freq_mag_samples[:, 0],
                                                    '4': Freq_mag_samples[:, 1],
                                                    '5': Freq_mag_samples[:, 2]})
                Freq_mag_samples_df['6'] = Freq_mag_samples_df['5'] * 0.1
                Freq_mag_samples_df['7'] = Freq_mag_samples_df['5'] * 0.01

                # Frequency_VEI_3 = Freq_mag_samples_df['3'].mean()
                # Frequency_VEI_4 = Freq_mag_samples_df['4'].mean()
                # Frequency_VEI_5 = Freq_mag_samples_df['5'].mean()
                # Frequency_VEI_6 = Frequency_VEI_5*0.1
                # Frequency_VEI_7 = Frequency_VEI_5*0.01
                #
                # Frequency_VEI_3_std = Freq_mag_samples_df['3'].std()
                # Frequency_VEI_4_std = Freq_mag_samples_df['4'].std()
                # Frequency_VEI_5_std = Freq_mag_samples_df['5'].std()
                # Frequency_VEI_6_std = Freq_mag_samples_df['5'].std()
                # Frequency_VEI_7_std = Freq_mag_samples_df['5'].std()

            else:
                Freq_mag_samples = stats.dirichlet.rvs(freq_mag, size=size, random_state=None)
                Freq_mag_samples_df = pd.DataFrame({'3': Freq_mag_samples[:, 0],
                                                    '4': Freq_mag_samples[:, 1],
                                                    '5': Freq_mag_samples[:, 2]})
                Freq_mag_samples_df['6'] = 0
                Freq_mag_samples_df['7'] = 0

                # Frequency_VEI_3 = Freq_mag_samples_df['3'].mean()
                # Frequency_VEI_4 = Freq_mag_samples_df['4'].mean()
                # Frequency_VEI_5 = Freq_mag_samples_df['5'].mean()
                # Frequency_VEI_6 = 0
                # Frequency_VEI_7 = 0
                #
                # Frequency_VEI_3_std = Freq_mag_samples_df['3'].std()
                # Frequency_VEI_4_std = Freq_mag_samples_df['4'].std()
                # Frequency_VEI_5_std = Freq_mag_samples_df['5'].std()
                # Frequency_VEI_6_std = 0
                # Frequency_VEI_7_std = 0

        elif find_FM[1] > 0 and find_FM[2] == 0:
            print("Dirichlet 4")
            freq_mag_ = find_FM.iloc[0:2]
            freq_mag = freq_mag_.values
            freq_mag = freq_mag*10
            size = (int(find_FM.sum())*10)
            Freq_mag_samples = stats.dirichlet.rvs(freq_mag, size=size, random_state=None)


            if Power_law == "True":

                Freq_mag_samples_df = pd.DataFrame({'3': Freq_mag_samples[:, 0],
                                                    '4': Freq_mag_samples[:, 1]})
                Freq_mag_samples_df['5'] = Freq_mag_samples_df['4'] * 0.1
                Freq_mag_samples_df['6'] = Freq_mag_samples_df['4'] * 0.01
                Freq_mag_samples_df['7'] = Freq_mag_samples_df['4'] * 0.001

                # Frequency_VEI_3 = Freq_mag_samples_df['3'].mean()
                # Frequency_VEI_4 = Freq_mag_samples_df['4'].mean()
                # Frequency_VEI_5 = Frequency_VEI_4*0.1
                # Frequency_VEI_6 = Frequency_VEI_4*0.01
                # Frequency_VEI_7 = Frequency_VEI_4*0.001
                #
                # Frequency_VEI_3_std = Freq_mag_samples_df['3'].std()
                # Frequency_VEI_4_std = Freq_mag_samples_df['4'].std()
                # Frequency_VEI_5_std = Freq_mag_samples_df['4'].std()
                # Frequency_VEI_6_std = Freq_mag_samples_df['4'].std()
                # Frequency_VEI_7_std = Freq_mag_samples_df['4'].std()
            else:
                Freq_mag_samples_df = pd.DataFrame({'3': Freq_mag_samples[:, 0],
                                                    '4': Freq_mag_samples[:, 1]})
                Freq_mag_samples_df['5'] = 0
                Freq_mag_samples_df['6'] = 0
                Freq_mag_samples_df['7'] = 0


                # Frequency_VEI_3 = Freq_mag_samples_df['3'].mean()
                # Frequency_VEI_4 = Freq_mag_samples_df['4'].mean()
                # Frequency_VEI_5 = 0
                # Frequency_VEI_6 = 0
                # Frequency_VEI_7 = 0
                #
                # Frequency_VEI_3_std = Freq_mag_samples_df['3'].std()
                # Frequency_VEI_4_std = Freq_mag_samples_df['4'].std()
                # Frequency_VEI_5_std = 0
                # Frequency_VEI_6_std = 0
                # Frequency_VEI_7_std = 0

        #return(Frequency_VEI_3, Frequency_VEI_4, Frequency_VEI_5, Frequency_VEI_6, Frequency_VEI_7,
         #      Frequency_VEI_3_std, Frequency_VEI_4_std, Frequency_VEI_5_std, Frequency_VEI_6_std, Frequency_VEI_7_std)

    elif VEI_schema == 2:

        if find_FM[5] > 0:
            print("Dirichlet 1")
            freq_mag = find_FM.values
            freq_mag = freq_mag * 100
            size = (int(find_FM.sum()) * 100)
            Freq_mag_samples = stats.dirichlet.rvs(freq_mag, size=size, random_state=None)
            Freq_mag_samples_df = pd.DataFrame({'2': Freq_mag_samples[:, 0],
                                                '3': Freq_mag_samples[:, 1],
                                                '4': Freq_mag_samples[:, 2],
                                                '5': Freq_mag_samples[:, 3],
                                                '6': Freq_mag_samples[:, 4],
                                                '7': Freq_mag_samples[:, 5]})

            # Frequency_VEI_3 = Freq_mag_samples_df['3'].mean()
            # Frequency_VEI_4 = Freq_mag_samples_df['4'].mean()
            # Frequency_VEI_5 = Freq_mag_samples_df['5'].mean()
            # Frequency_VEI_6 = Freq_mag_samples_df['6'].mean()
            # Frequency_VEI_7 = Freq_mag_samples_df['7'].mean()
            #
            # Frequency_VEI_3_std = Freq_mag_samples_df['3'].std()
            # Frequency_VEI_4_std = Freq_mag_samples_df['4'].std()
            # Frequency_VEI_5_std = Freq_mag_samples_df['5'].std()
            # Frequency_VEI_6_std = Freq_mag_samples_df['6'].std()
            # Frequency_VEI_7_std = Freq_mag_samples_df['7'].std()

        elif find_FM[4] > 0 and find_FM[5] == 0:
            print("Dirichlet 2")
            freq_mag_ = find_FM.iloc[0:5]
            freq_mag = freq_mag_.values
            freq_mag = freq_mag * 10
            size = (int(find_FM.sum()) * 10)

            if Power_law == "True":
                Freq_mag_samples = stats.dirichlet.rvs(freq_mag, size=size, random_state=None)
                Freq_mag_samples_df = pd.DataFrame({'2': Freq_mag_samples[:, 0],
                                                    '3': Freq_mag_samples[:, 1],
                                                    '4': Freq_mag_samples[:, 2],
                                                    '5': Freq_mag_samples[:, 3],
                                                    '6': Freq_mag_samples[:, 4]})
                Freq_mag_samples_df['7'] = Freq_mag_samples_df['6'] * 0.1

                # Frequency_VEI_3 = Freq_mag_samples_df['3'].mean()
                # Frequency_VEI_4 = Freq_mag_samples_df['4'].mean()
                # Frequency_VEI_5 = Freq_mag_samples_df['5'].mean()
                # Frequency_VEI_6 = Freq_mag_samples_df['6'].mean()
                # Frequency_VEI_7 = Frequency_VEI_6*0.1
                #
                # Frequency_VEI_3_std = Freq_mag_samples_df['3'].std()
                # Frequency_VEI_4_std = Freq_mag_samples_df['4'].std()
                # Frequency_VEI_5_std = Freq_mag_samples_df['5'].std()
                # Frequency_VEI_6_std = Freq_mag_samples_df['6'].std()
                # Frequency_VEI_7_std = Freq_mag_samples_df['6'].std()
            else:
                Freq_mag_samples = stats.dirichlet.rvs(freq_mag, size=size, random_state=None)
                Freq_mag_samples_df = pd.DataFrame({'2': Freq_mag_samples[:, 0],
                                                    '3': Freq_mag_samples[:, 1],
                                                    '4': Freq_mag_samples[:, 2],
                                                    '5': Freq_mag_samples[:, 3],
                                                    '6': Freq_mag_samples[:, 4]})
                Freq_mag_samples_df['7'] = 0

                # Frequency_VEI_3 = Freq_mag_samples_df['3'].mean()
                # Frequency_VEI_4 = Freq_mag_samples_df['4'].mean()
                # Frequency_VEI_5 = Freq_mag_samples_df['5'].mean()
                # Frequency_VEI_6 = Freq_mag_samples_df['6'].mean()
                # Frequency_VEI_7 = 0
                #
                # Frequency_VEI_3_std = Freq_mag_samples_df['3'].std()
                # Frequency_VEI_4_std = Freq_mag_samples_df['4'].std()
                # Frequency_VEI_5_std = Freq_mag_samples_df['5'].std()
                # Frequency_VEI_6_std = Freq_mag_samples_df['6'].std()
                # Frequency_VEI_7_std = 0

        elif find_FM[3] > 0 and find_FM[4] == 0:
            print("Dirichlet 3")
            freq_mag_ = find_FM.iloc[0:4]
            freq_mag = freq_mag_.values
            freq_mag = freq_mag * 10
            size = (int(find_FM.sum()) * 10)

            if Power_law == "True":
                Freq_mag_samples = stats.dirichlet.rvs(freq_mag, size=size, random_state=None)
                Freq_mag_samples_df = pd.DataFrame({'2': Freq_mag_samples[:, 0],
                                                    '3': Freq_mag_samples[:, 1],
                                                    '4': Freq_mag_samples[:, 2],
                                                    '5': Freq_mag_samples[:, 3]})
                Freq_mag_samples_df['6'] = Freq_mag_samples_df['5'] * 0.1
                Freq_mag_samples_df['7'] = Freq_mag_samples_df['5'] * 0.01

                # Frequency_VEI_3 = Freq_mag_samples_df['3'].mean()
                # Frequency_VEI_4 = Freq_mag_samples_df['4'].mean()
                # Frequency_VEI_5 = Freq_mag_samples_df['5'].mean()
                # Frequency_VEI_6 = Frequency_VEI_5*0.1
                # Frequency_VEI_7 = Frequency_VEI_5*0.01
                #
                # Frequency_VEI_3_std = Freq_mag_samples_df['3'].std()
                # Frequency_VEI_4_std = Freq_mag_samples_df['4'].std()
                # Frequency_VEI_5_std = Freq_mag_samples_df['5'].std()
                # Frequency_VEI_6_std = Freq_mag_samples_df['5'].std()
                # Frequency_VEI_7_std = Freq_mag_samples_df['5'].std()

            else:
                Freq_mag_samples = stats.dirichlet.rvs(freq_mag, size=size, random_state=None)
                Freq_mag_samples_df = pd.DataFrame({'2': Freq_mag_samples[:, 0],
                                                    '3': Freq_mag_samples[:, 1],
                                                    '4': Freq_mag_samples[:, 2],
                                                    '5': Freq_mag_samples[:, 3]})
                Freq_mag_samples_df['6'] = 0
                Freq_mag_samples_df['7'] = 0

                # Frequency_VEI_3 = Freq_mag_samples_df['3'].mean()
                # Frequency_VEI_4 = Freq_mag_samples_df['4'].mean()
                # Frequency_VEI_5 = Freq_mag_samples_df['5'].mean()
                # Frequency_VEI_6 = 0
                # Frequency_VEI_7 = 0
                #
                # Frequency_VEI_3_std = Freq_mag_samples_df['3'].std()
                # Frequency_VEI_4_std = Freq_mag_samples_df['4'].std()
                # Frequency_VEI_5_std = Freq_mag_samples_df['5'].std()
                # Frequency_VEI_6_std = 0
                # Frequency_VEI_7_std = 0

        elif find_FM[2] > 0 and find_FM[3] == 0:
            print("Dirichlet 4")
            freq_mag_ = find_FM.iloc[0:3]
            freq_mag = freq_mag_.values
            freq_mag = freq_mag * 10
            size = (int(find_FM.sum()) * 10)
            Freq_mag_samples = stats.dirichlet.rvs(freq_mag, size=size, random_state=None)

            if Power_law == "True":

                Freq_mag_samples_df = pd.DataFrame({'2': Freq_mag_samples[:, 0],
                                                    '3': Freq_mag_samples[:, 1],
                                                    '4': Freq_mag_samples[:, 2]})
                Freq_mag_samples_df['5'] = Freq_mag_samples_df['4'] * 0.1
                Freq_mag_samples_df['6'] = Freq_mag_samples_df['4'] * 0.01
                Freq_mag_samples_df['7'] = Freq_mag_samples_df['4'] * 0.001

                # Frequency_VEI_3 = Freq_mag_samples_df['3'].mean()
                # Frequency_VEI_4 = Freq_mag_samples_df['4'].mean()
                # Frequency_VEI_5 = Frequency_VEI_4*0.1
                # Frequency_VEI_6 = Frequency_VEI_4*0.01
                # Frequency_VEI_7 = Frequency_VEI_4*0.001
                #
                # Frequency_VEI_3_std = Freq_mag_samples_df['3'].std()
                # Frequency_VEI_4_std = Freq_mag_samples_df['4'].std()
                # Frequency_VEI_5_std = Freq_mag_samples_df['4'].std()
                # Frequency_VEI_6_std = Freq_mag_samples_df['4'].std()
                # Frequency_VEI_7_std = Freq_mag_samples_df['4'].std()
            else:
                Freq_mag_samples_df = pd.DataFrame({'2': Freq_mag_samples[:, 0],
                                                    '3': Freq_mag_samples[:, 1],
                                                    '4': Freq_mag_samples[:, 2]})
                Freq_mag_samples_df['5'] = 0
                Freq_mag_samples_df['6'] = 0
                Freq_mag_samples_df['7'] = 0

                # Frequency_VEI_3 = Freq_mag_samples_df['3'].mean()
                # Frequency_VEI_4 = Freq_mag_samples_df['4'].mean()
                # Frequency_VEI_5 = 0
                # Frequency_VEI_6 = 0
                # Frequency_VEI_7 = 0
                #
                # Frequency_VEI_3_std = Freq_mag_samples_df['3'].std()
                # Frequency_VEI_4_std = Freq_mag_samples_df['4'].std()
                # Frequency_VEI_5_std = 0
                # Frequency_VEI_6_std = 0
                # Frequency_VEI_7_std = 0

        # return(Frequency_VEI_3, Frequency_VEI_4, Frequency_VEI_5, Frequency_VEI_6, Frequency_VEI_7,
        #      Frequency_VEI_3_std, Frequency_VEI_4_std, Frequency_VEI_5_std, Frequency_VEI_6_std, Frequency_VEI_7_std)

    elif VEI_schema == 3:
        # VEI <=1, 2, 3, 4, 5, 6, 7
        if find_FM[6] > 0:
            print("Dirichlet 1")
            freq_mag = find_FM.values
            freq_mag = freq_mag * 100
            size = (int(find_FM.sum()) * 100)
            Freq_mag_samples = stats.dirichlet.rvs(freq_mag, size=size, random_state=None)
            Freq_mag_samples_df = pd.DataFrame({'1': Freq_mag_samples[:, 0],
                                                '2': Freq_mag_samples[:, 1],
                                                '3': Freq_mag_samples[:, 2],
                                                '4': Freq_mag_samples[:, 3],
                                                '5': Freq_mag_samples[:, 4],
                                                '6': Freq_mag_samples[:, 5],
                                                '7': Freq_mag_samples[:, 6]})

            # Frequency_VEI_3 = Freq_mag_samples_df['3'].mean()
            # Frequency_VEI_4 = Freq_mag_samples_df['4'].mean()
            # Frequency_VEI_5 = Freq_mag_samples_df['5'].mean()
            # Frequency_VEI_6 = Freq_mag_samples_df['6'].mean()
            # Frequency_VEI_7 = Freq_mag_samples_df['7'].mean()
            #
            # Frequency_VEI_3_std = Freq_mag_samples_df['3'].std()
            # Frequency_VEI_4_std = Freq_mag_samples_df['4'].std()
            # Frequency_VEI_5_std = Freq_mag_samples_df['5'].std()
            # Frequency_VEI_6_std = Freq_mag_samples_df['6'].std()
            # Frequency_VEI_7_std = Freq_mag_samples_df['7'].std()

        elif find_FM[5] > 0 and find_FM[6] == 0:
            print("Dirichlet 2")
            freq_mag_ = find_FM.iloc[0:6]
            freq_mag = freq_mag_.values
            freq_mag = freq_mag * 10
            size = (int(find_FM.sum()) * 10)

            if Power_law == "True":
                Freq_mag_samples = stats.dirichlet.rvs(freq_mag, size=size, random_state=None)
                Freq_mag_samples_df = pd.DataFrame({'1': Freq_mag_samples[:, 0],
                                                    '2': Freq_mag_samples[:, 1],
                                                    '3': Freq_mag_samples[:, 2],
                                                    '4': Freq_mag_samples[:, 3],
                                                    '5': Freq_mag_samples[:, 4],
                                                    '6': Freq_mag_samples[:, 5]})
                Freq_mag_samples_df['7'] = Freq_mag_samples_df['6'] * 0.1

                # Frequency_VEI_3 = Freq_mag_samples_df['3'].mean()
                # Frequency_VEI_4 = Freq_mag_samples_df['4'].mean()
                # Frequency_VEI_5 = Freq_mag_samples_df['5'].mean()
                # Frequency_VEI_6 = Freq_mag_samples_df['6'].mean()
                # Frequency_VEI_7 = Frequency_VEI_6*0.1
                #
                # Frequency_VEI_3_std = Freq_mag_samples_df['3'].std()
                # Frequency_VEI_4_std = Freq_mag_samples_df['4'].std()
                # Frequency_VEI_5_std = Freq_mag_samples_df['5'].std()
                # Frequency_VEI_6_std = Freq_mag_samples_df['6'].std()
                # Frequency_VEI_7_std = Freq_mag_samples_df['6'].std()
            else:
                Freq_mag_samples = stats.dirichlet.rvs(freq_mag, size=size, random_state=None)
                Freq_mag_samples_df = pd.DataFrame({'1': Freq_mag_samples[:, 0],
                                                    '2': Freq_mag_samples[:, 1],
                                                    '3': Freq_mag_samples[:, 2],
                                                    '4': Freq_mag_samples[:, 3],
                                                    '5': Freq_mag_samples[:, 4],
                                                    '6': Freq_mag_samples[:, 5]})
                Freq_mag_samples_df['7'] = 0

                # Frequency_VEI_3 = Freq_mag_samples_df['3'].mean()
                # Frequency_VEI_4 = Freq_mag_samples_df['4'].mean()
                # Frequency_VEI_5 = Freq_mag_samples_df['5'].mean()
                # Frequency_VEI_6 = Freq_mag_samples_df['6'].mean()
                # Frequency_VEI_7 = 0
                #
                # Frequency_VEI_3_std = Freq_mag_samples_df['3'].std()
                # Frequency_VEI_4_std = Freq_mag_samples_df['4'].std()
                # Frequency_VEI_5_std = Freq_mag_samples_df['5'].std()
                # Frequency_VEI_6_std = Freq_mag_samples_df['6'].std()
                # Frequency_VEI_7_std = 0

        elif find_FM[4] > 0 and find_FM[5] == 0:
            print("Dirichlet 3")
            freq_mag_ = find_FM.iloc[0:5]
            freq_mag = freq_mag_.values
            freq_mag = freq_mag * 10
            size = (int(find_FM.sum()) * 10)

            if Power_law == "True":
                Freq_mag_samples = stats.dirichlet.rvs(freq_mag, size=size, random_state=None)
                Freq_mag_samples_df = pd.DataFrame({'1': Freq_mag_samples[:, 0],
                                                    '2': Freq_mag_samples[:, 1],
                                                    '3': Freq_mag_samples[:, 2],
                                                    '4': Freq_mag_samples[:, 3],
                                                    '5': Freq_mag_samples[:, 4]})
                Freq_mag_samples_df['6'] = Freq_mag_samples_df['5'] * 0.1
                Freq_mag_samples_df['7'] = Freq_mag_samples_df['5'] * 0.01

                # Frequency_VEI_3 = Freq_mag_samples_df['3'].mean()
                # Frequency_VEI_4 = Freq_mag_samples_df['4'].mean()
                # Frequency_VEI_5 = Freq_mag_samples_df['5'].mean()
                # Frequency_VEI_6 = Frequency_VEI_5*0.1
                # Frequency_VEI_7 = Frequency_VEI_5*0.01
                #
                # Frequency_VEI_3_std = Freq_mag_samples_df['3'].std()
                # Frequency_VEI_4_std = Freq_mag_samples_df['4'].std()
                # Frequency_VEI_5_std = Freq_mag_samples_df['5'].std()
                # Frequency_VEI_6_std = Freq_mag_samples_df['5'].std()
                # Frequency_VEI_7_std = Freq_mag_samples_df['5'].std()

            else:
                Freq_mag_samples = stats.dirichlet.rvs(freq_mag, size=size, random_state=None)
                Freq_mag_samples_df = pd.DataFrame({'1': Freq_mag_samples[:, 0],
                                                    '2': Freq_mag_samples[:, 1],
                                                    '3': Freq_mag_samples[:, 2],
                                                    '4': Freq_mag_samples[:, 3],
                                                    '5': Freq_mag_samples[:, 4]})
                Freq_mag_samples_df['6'] = 0
                Freq_mag_samples_df['7'] = 0

                # Frequency_VEI_3 = Freq_mag_samples_df['3'].mean()
                # Frequency_VEI_4 = Freq_mag_samples_df['4'].mean()
                # Frequency_VEI_5 = Freq_mag_samples_df['5'].mean()
                # Frequency_VEI_6 = 0
                # Frequency_VEI_7 = 0
                #
                # Frequency_VEI_3_std = Freq_mag_samples_df['3'].std()
                # Frequency_VEI_4_std = Freq_mag_samples_df['4'].std()
                # Frequency_VEI_5_std = Freq_mag_samples_df['5'].std()
                # Frequency_VEI_6_std = 0
                # Frequency_VEI_7_std = 0

        elif find_FM[3] > 0 and find_FM[4] == 0:
            print("Dirichlet 4")
            freq_mag_ = find_FM.iloc[0:4]
            freq_mag = freq_mag_.values
            freq_mag = freq_mag * 10
            size = (int(find_FM.sum()) * 10)
            Freq_mag_samples = stats.dirichlet.rvs(freq_mag, size=size, random_state=None)

            if Power_law == "True":

                Freq_mag_samples_df = pd.DataFrame({'1': Freq_mag_samples[:, 0],
                                                    '2': Freq_mag_samples[:, 1],
                                                    '3': Freq_mag_samples[:, 2],
                                                    '4': Freq_mag_samples[:, 3]})
                Freq_mag_samples_df['5'] = Freq_mag_samples_df['4'] * 0.1
                Freq_mag_samples_df['6'] = Freq_mag_samples_df['4'] * 0.01
                Freq_mag_samples_df['7'] = Freq_mag_samples_df['4'] * 0.001

                # Frequency_VEI_3 = Freq_mag_samples_df['3'].mean()
                # Frequency_VEI_4 = Freq_mag_samples_df['4'].mean()
                # Frequency_VEI_5 = Frequency_VEI_4*0.1
                # Frequency_VEI_6 = Frequency_VEI_4*0.01
                # Frequency_VEI_7 = Frequency_VEI_4*0.001
                #
                # Frequency_VEI_3_std = Freq_mag_samples_df['3'].std()
                # Frequency_VEI_4_std = Freq_mag_samples_df['4'].std()
                # Frequency_VEI_5_std = Freq_mag_samples_df['4'].std()
                # Frequency_VEI_6_std = Freq_mag_samples_df['4'].std()
                # Frequency_VEI_7_std = Freq_mag_samples_df['4'].std()
            else:
                Freq_mag_samples_df = pd.DataFrame({'1': Freq_mag_samples[:, 0],
                                                    '2': Freq_mag_samples[:, 1],
                                                    '3': Freq_mag_samples[:, 2],
                                                    '4': Freq_mag_samples[:, 3]})
                Freq_mag_samples_df['5'] = 0
                Freq_mag_samples_df['6'] = 0
                Freq_mag_samples_df['7'] = 0

                # Frequency_VEI_3 = Freq_mag_samples_df['3'].mean()
                # Frequency_VEI_4 = Freq_mag_samples_df['4'].mean()
                # Frequency_VEI_5 = 0
                # Frequency_VEI_6 = 0
                # Frequency_VEI_7 = 0
                #
                # Frequency_VEI_3_std = Freq_mag_samples_df['3'].std()
                # Frequency_VEI_4_std = Freq_mag_samples_df['4'].std()
                # Frequency_VEI_5_std = 0
                # Frequency_VEI_6_std = 0
                # Frequency_VEI_7_std = 0

        # return(Frequency_VEI_3, Frequency_VEI_4, Frequency_VEI_5, Frequency_VEI_6, Frequency_VEI_7,
        #      Frequency_VEI_3_std, Frequency_VEI_4_std, Frequency_VEI_5_std, Frequency_VEI_6_std, Frequency_VEI_7_std)

    return(Freq_mag_samples_df)

def set_analogue_eruption_rate (volc_num, Analogue_volcanoes, Freq_rate, Freq_rate_std):
    average_eruptions_per_year = None
    std_eruptions_per_year = None
    volcano_type = Analogue_volcanoes[(Analogue_volcanoes['Volcano Number'] == volc_num)].values[0]
    volcano_type = volcano_type[2]
    Freq_rate_t = Freq_rate.T
    if (volcano_type == "Caldera"):
        average_eruptions_per_year = Freq_rate[0]
    elif (volcano_type == "Distributed cones and fields"):
        average_eruptions_per_year = Freq_rate[0]
    elif (volcano_type == "Large cone"):
        average_eruptions_per_year = Freq_rate[1]
    elif (volcano_type == "Large caldera"):
        average_eruptions_per_year = Freq_rate[1]
    elif (volcano_type == "Lava dome"):
        average_eruptions_per_year = Freq_rate[2]
    elif (volcano_type == "Open-vent stratocone"):
        average_eruptions_per_year = Freq_rate[2]
    elif (volcano_type == "Shield"):
        average_eruptions_per_year = Freq_rate[3]
    elif (volcano_type == "Semi-plugged stratocone"):
        average_eruptions_per_year = Freq_rate[3]
    elif (volcano_type == "Small cone"):
        average_eruptions_per_year = Freq_rate[4]
    elif (volcano_type == "Well-plugged stratocone"):
        average_eruptions_per_year = Freq_rate[4]
    if (volcano_type == "Caldera"):
        std_eruptions_per_year = Freq_rate_std[0]
    elif (volcano_type == "Distributed cones and fields"):
        std_eruptions_per_year = Freq_rate_std[0]
    elif (volcano_type == "Large cone"):
        std_eruptions_per_year = Freq_rate_std[1]
    elif (volcano_type == "Large caldera"):
        std_eruptions_per_year = Freq_rate_std[1]
    elif (volcano_type == "Lava dome"):
        std_eruptions_per_year = Freq_rate_std[2]
    elif (volcano_type == "Open-vent stratocone"):
        std_eruptions_per_year = Freq_rate_std[2]
    elif (volcano_type == "Shield"):
        std_eruptions_per_year = Freq_rate_std[3]
    elif (volcano_type == "Semi-plugged stratocone"):
        std_eruptions_per_year = Freq_rate_std[3]
    elif (volcano_type == "Small cone"):
        std_eruptions_per_year = Freq_rate_std[4]
    elif (volcano_type == "Well-plugged stratocone"):
        std_eruptions_per_year = Freq_rate_std[4]

    return(average_eruptions_per_year, std_eruptions_per_year)

#Fix this function -Remove?
def get_FM_relationship (analogue_freq_mag, volc_num, GVP_volcanoes, Volcano_type, Observed_rate, Observed_rate_full_record, Observed_VEI, Observed_VEI_full_record, confirmed_VEI_eruptions, Freq_rate, Freq_rate_std, VEI_freq, get_trace, show_trace):
    AP_VEI_3 = None
    AP_VEI_4 = None
    AP_VEI_5 = None
    AP_VEI_6 = None
    AP_VEI_7 = None
    Volcano_data = GVP_volcanoes[(GVP_volcanoes['Volcano Number'] == volc_num)]
    Volcano_name = Volcano_data.iloc[0]['Volcano Name']
    alpha_analogue = Freq_rate[Volcano_type] ** 2 * (
                (1 - Freq_rate[Volcano_type]) / Freq_rate_std[Volcano_type] ** 2 - 1 / Freq_rate[Volcano_type])
    beta_analogue = alpha_analogue * (1 / Freq_rate[Volcano_type] - 1)
    if confirmed_VEI_eruptions == 0:
        Method = "Bayesian update using full record"
        print("Initiating model using the: ", Method, " Method")
        # Probability of an eruption of any VEI
        cores=1
        y = np.array([Observed_rate_full_record*10000])
        n = np.array([10000])
        N = len(n)
        RANDOM_SEED = 42
        #if __name__ == '__main__':
        # Probability of eruption model
        with pm.Model() as AP_model_GVP:
            theta = pm.Beta('theta', alpha=alpha_analogue, beta=beta_analogue)
            p = pm.Binomial('y', p=theta, observed=y, n=n)
            trace_AP = pm.sample(20000, tune=10000, chains=2, cores=cores, target_accept=0.99)

            if get_trace == True:
                arviz.plot_trace(trace_AP)
                pthfigAP = "Japan/Figures/TracePlots/" + Volcano_name + "_traceplot.png"
                plt.savefig(pthfigAP)
                if show_trace == True:
                    plt.show()
                else:
                    plt.close()
            else:
                pass
            ppc_AP = pm.sample_posterior_predictive(trace_AP, var_names=['y'], random_seed=RANDOM_SEED)
            ppc_shape  = ppc_AP["y"].shape
            # arviz.plot_autocorr(trace_AP)
        ppc_ave_adjust = ppc_AP['y'] / 10000
        AP_ave = ppc_ave_adjust.mean()
        AP_std = ppc_ave_adjust.std()

        find_FM = analogue_freq_mag[Volcano_type]*100

        ######
        Jenkins_1 = np.array([find_FM[0],
                              find_FM[1],
                              find_FM[2],
                              find_FM[3],
                              find_FM[4]])

        Jenkins_2 = np.array([find_FM[0],
                              find_FM[1],
                              find_FM[2]])

        Jenkins_3 = np.array([find_FM[0],
                              find_FM[1]])

        if find_FM[4] > 0:
            Jenkins_c = Jenkins_1
            print("Using Jenkins_1")
        elif find_FM[2] > 0 and find_FM[4] == 0:
            Jenkins_c = Jenkins_2
            print("Using Jenkins_2")
        elif find_FM[2] == 0 and find_FM[4] == 0:
            Jenkins_c = Jenkins_3
            print("Using Jenkins_3")

        #Jenkins_c = Jenkins_1
        Jenkins_shape = Jenkins_c.size
        Observed_1 = np.array([Observed_VEI_full_record[0],
                               Observed_VEI_full_record[1],
                               Observed_VEI_full_record[2],
                               Observed_VEI_full_record[3],
                               Observed_VEI_full_record[4]])

        Observed_2 = np.array([Observed_VEI_full_record[0],
                               Observed_VEI_full_record[1],
                               Observed_VEI_full_record[2]])

        Observed_3 = np.array([Observed_VEI_full_record[0],
                               Observed_VEI_full_record[1]])

        if find_FM[4] > 0:
            Observed = Observed_1
        elif find_FM[2] > 0 and find_FM[4] == 0:
            Observed = Observed_2
        elif find_FM[2] == 0 and find_FM[4] == 0:
            Observed = Observed_3
        else:
            Observed = print("ERROR - Check observed version being used")

        #Observed = Observed_1
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

            if get_trace == True:
                arviz.plot_density(trace_VEI_GVP)
                pthfigVEI = "Japan/Figures/TracePlots/" + Volcano_name + "_densityplot.png"
                plt.savefig(pthfigVEI)
                if show_trace == True:
                    plt.show()
                else:
                    plt.close()
            else:
                pass
            # arviz.plot_autocorr(trace_VEI_GVP)
            # plt.show()
            VEI_trace_data = trace_VEI_GVP.get_values(parameters, burn=5000)
            VEI_trace_data_df = pd.DataFrame(VEI_trace_data)
            ppc_VEI = pm.sample_posterior_predictive(trace_VEI_GVP, var_names=['observed_data'], random_seed=RANDOM_SEED)


        if Jenkins_shape == 5:
            try:
                ppc_ave_VEI_df = pd.DataFrame(ppc_VEI)
                ppc_ave_VEI_df2 = pd.DataFrame(ppc_ave_VEI_df['observed_data'].to_list(),
                                               columns=['VEI3', 'VEI4', 'VEI5', 'VEI6', 'VEI7'])
                ppc_ave_VEI_df2 = ppc_ave_VEI_df2.fillna(0)
            except ValueError:
                ppc_ave_VEI_df = pd.DataFrame(ppc_VEI['observed_data'])
                ppc_ave_VEI_df2 = ppc_ave_VEI_df
                ppc_ave_VEI_df2.columns = ['VEI3', 'VEI4', 'VEI5', 'VEI6', 'VEI7']

            #arviz.plot_density(trace_VEI_GVP)
            #plt.show()

            Ave_VEI3 = VEI_trace_data_df.iloc[:, 0].mean()
            Std_VEI3 = VEI_trace_data_df.iloc[:, 0].std()

            Ave_VEI4 = VEI_trace_data_df.iloc[:, 1].mean()
            Std_VEI4 = VEI_trace_data_df.iloc[:, 1].std()

            Ave_VEI5 = VEI_trace_data_df.iloc[:, 2].mean()
            Std_VEI5 = VEI_trace_data_df.iloc[:, 2].std()

            Ave_VEI6 = VEI_trace_data_df.iloc[:, 3].mean()
            Std_VEI6 = VEI_trace_data_df.iloc[:, 3].std()

            Ave_VEI7 = VEI_trace_data_df.iloc[:, 4].mean()
            Std_VEI7 = VEI_trace_data_df.iloc[:, 4].std()

        elif Jenkins_shape == 3:
            try:
                ppc_ave_VEI_df = pd.DataFrame(ppc_VEI)
                ppc_ave_VEI_df2 = pd.DataFrame(ppc_ave_VEI_df['observed_data'].to_list(),
                                               columns=['VEI3', 'VEI4', 'VEI5'])
                ppc_ave_VEI_df2 = ppc_ave_VEI_df2.fillna(0)
            except ValueError:
                ppc_ave_VEI_df = pd.DataFrame(ppc_VEI['observed_data'])
                ppc_ave_VEI_df2 = ppc_ave_VEI_df
                ppc_ave_VEI_df2.columns = ['VEI3', 'VEI4', 'VEI5']

            # arviz.plot_density(trace_VEI_GVP)
            # plt.show()

            Ave_VEI3 = VEI_trace_data_df.iloc[:, 0].mean()
            Std_VEI3 = VEI_trace_data_df.iloc[:, 0].std()

            Ave_VEI4 = VEI_trace_data_df.iloc[:, 1].mean()
            Std_VEI4 = VEI_trace_data_df.iloc[:, 1].std()

            Ave_VEI5 = VEI_trace_data_df.iloc[:, 2].mean()
            Std_VEI5 = VEI_trace_data_df.iloc[:, 2].std()

            Ave_VEI6 = 0
            Std_VEI6 = 0

            Ave_VEI7 = 0
            Std_VEI7 = 0

        elif Jenkins_shape == 2:
            try:
                ppc_ave_VEI_df = pd.DataFrame(ppc_VEI)
                ppc_ave_VEI_df2 = pd.DataFrame(ppc_ave_VEI_df['observed_data'].to_list(),
                                               columns=['VEI3', 'VEI4'])
                ppc_ave_VEI_df2 = ppc_ave_VEI_df2.fillna(0)
            except ValueError:
                ppc_ave_VEI_df = pd.DataFrame(ppc_VEI['observed_data'])
                ppc_ave_VEI_df2 = ppc_ave_VEI_df
                ppc_ave_VEI_df2.columns = ['VEI3', 'VEI4']

            # arviz.plot_density(trace_VEI_GVP)
            # plt.show()

            Ave_VEI3 = VEI_trace_data_df.iloc[:, 0].mean()
            Std_VEI3 = VEI_trace_data_df.iloc[:, 0].std()

            Ave_VEI4 = VEI_trace_data_df.iloc[:, 1].mean()
            Std_VEI4 = VEI_trace_data_df.iloc[:, 1].std()

            Ave_VEI5 = 0
            Std_VEI5 = 0

            Ave_VEI6 = 0
            Std_VEI6 = 0

            Ave_VEI7 = 0
            Std_VEI7 = 0

        if Ave_VEI3 > 0:
            a_ave_VEI3 = Ave_VEI3 ** 2 * ((1 - Ave_VEI3) / Std_VEI3 ** 2 - 1 / Ave_VEI3)
            b_ave_VEI3 = a_ave_VEI3 * (1 / Ave_VEI3 - 1)
            data_beta_VEI3 = beta.rvs(a=a_ave_VEI3, b=b_ave_VEI3, size=1000000)
        else:
            a_ave_VEI3 = 0
            b_ave_VEI3 = 0
            data_beta_VEI3 = 0

        if Ave_VEI4 > 0:
            a_ave_VEI4 = Ave_VEI4 ** 2 * ((1 - Ave_VEI4) / Std_VEI4 ** 2 - 1 / Ave_VEI4)
            b_ave_VEI4 = a_ave_VEI4 * (1 / Ave_VEI4 - 1)
            data_beta_VEI4 = beta.rvs(a=a_ave_VEI4, b=b_ave_VEI4, size=1000000)
        else:
            a_ave_VEI4 = 0
            b_ave_VEI4 = 0
            data_beta_VEI4 = 0

        if Ave_VEI5 > 0:
            a_ave_VEI5 = Ave_VEI5 ** 2 * ((1 - Ave_VEI5) / Std_VEI5 ** 2 - 1 / Ave_VEI5)
            b_ave_VEI5 = a_ave_VEI5 * (1 / Ave_VEI5 - 1)
            data_beta_VEI5 = beta.rvs(a=a_ave_VEI5, b=b_ave_VEI5, size=1000000)
        else:
            a_ave_VEI5 = 0
            b_ave_VEI5 = 0
            data_beta_VEI5 = 0

        if Ave_VEI6 > 0:
            a_ave_VEI6 = Ave_VEI6 ** 2 * ((1 - Ave_VEI6) / Std_VEI6 ** 2 - 1 / Ave_VEI6)
            b_ave_VEI6 = a_ave_VEI6 * (1 / Ave_VEI6 - 1)
            data_beta_VEI6 = beta.rvs(a=a_ave_VEI6, b=b_ave_VEI6, size=1000000)
        else:
            a_ave_VEI6 = 0
            b_ave_VEI6 = 0
            data_beta_VEI6 = 0

        if Ave_VEI7 > 0:
            a_ave_VEI7 = Ave_VEI7 ** 2 * ((1 - Ave_VEI7) / Std_VEI7 ** 2 - 1 / Ave_VEI7)
            b_ave_VEI7 = a_ave_VEI7 * (1 / Ave_VEI7 - 1)
            data_beta_VEI7 = beta.rvs(a=a_ave_VEI7, b=b_ave_VEI7, size=1000000)
        else:
            a_ave_VEI7 = 0
            b_ave_VEI7 = 0
            data_beta_VEI7 = 0

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
        a_ave = AP_ave ** 2 * ((1 - AP_ave) / AP_std ** 2 - 1 / AP_ave)
        b_ave = a_ave * (1 / AP_ave - 1)

        if VEI_3_mean > 0:
            AP_VEI_3 = []
            Number = 10000
            lower, upper = 0, 1
            for i in range(Number):
                """Monte Carlo simulation"""
                AP_VEI3 = (np.random.beta(a_ave, b_ave) * stats.truncnorm.rvs((lower - VEI_3_mean) / VEI_3_sd,
                                                                              (upper - VEI_3_mean) / VEI_3_sd,
                                                                              loc=VEI_3_mean,
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
                                                                              (upper - VEI_4_mean) / VEI_4_sd,
                                                                              loc=VEI_4_mean,
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
                                                                              (upper - VEI_5_mean) / VEI_5_sd,
                                                                              loc=VEI_5_mean,
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
                                                                              (upper - VEI_6_mean) / VEI_6_sd,
                                                                              loc=VEI_6_mean,
                                                                              scale=VEI_6_sd))
                AP_VEI_6 = AP_VEI_6 + [AP_VEI6]
        else:
            AP_VEI_6 = 0

        if VEI_7_mean > 0:
            AP_VEI_7 = []
            Number = 10000
            lower, upper = 0, 1
            for i in range(Number):
                """Monte Carlo simulation"""
                AP_VEI7 = (np.random.beta(a_ave, b_ave) * stats.truncnorm.rvs((lower - VEI_7_mean) / VEI_7_sd,
                                                                              (upper - VEI_7_mean) / VEI_7_sd,
                                                                              loc=VEI_7_mean,
                                                                              scale=VEI_7_sd))
                AP_VEI_7 = AP_VEI_7 + [AP_VEI7]
        else:
            AP_VEI_7 = 0


    else:
        Method = "Bayesian update using complete record"
        print("Initiating model using the: ", Method, " Method")
        # Probability of an eruption of any VEI
        cores=1
        y = np.array([Observed_rate*10000])
        n = np.array([10000])
        N = len(n)
        RANDOM_SEED = 42
        #if __name__ == '__main__':
        # Probability of eruption model
        with pm.Model() as AP_model_GVP:
            theta = pm.Beta('theta', alpha=alpha_analogue, beta=beta_analogue)
            p = pm.Binomial('y', p=theta, observed=y, n=n)
            trace_AP = pm.sample(20000, tune=10000, chains=2, cores=cores, target_accept=0.99)
            if get_trace == True:
                arviz.plot_trace(trace_AP)
                pthfigAP = "Japan/Figures/TracePlots/" + Volcano_name + "_traceplot.png"
                plt.savefig(pthfigAP)
                if show_trace == True:
                    plt.show()
                else:
                    plt.close()
            else:
                pass
            ppc_AP = pm.sample_posterior_predictive(trace_AP, var_names=['y'], random_seed=RANDOM_SEED)

            ppc_shape  = ppc_AP["y"].shape
            # arviz.plot_autocorr(trace_AP)
            #plt.show()
        ppc_ave_adjust = ppc_AP['y'] / 10000
        AP_ave = ppc_ave_adjust.mean()
        AP_std = ppc_ave_adjust.std()
        # a_ave = Ave ** 2 * ((1 - Ave) / Std ** 2 - 1 / Ave)
        # b_ave = a_ave * (1 / Ave - 1)
        # data_beta = beta.rvs(a=a_ave, b=b_ave, size=1000000)


        # #Visualization and summary
        # Alpha_Beta = arviz.summary(trace_AP)
        # print(Alpha_Beta)
        # arviz.plot_trace(trace_AP)
        # plt.show()
        # AP_trace_data = trace_AP.get_values('theta', burn=10000)
        # AP_trace_data_df = pd.DataFrame(AP_trace_data)
        # AP_ave = AP_trace_data_df.mean()
        # AP_std = AP_trace_data_df.std()

        #Probability of VEI model
        #get_volcano_type = Analogue_volcanoes.loc[Analogue_volcanoes['Volcano Number'] == volc_num, 'Volcano type'].iloc[0]

        find_FM = analogue_freq_mag[Volcano_type]*100

        ######
        Jenkins_1 = np.array([find_FM[0],
                              find_FM[1],
                              find_FM[2],
                              find_FM[3],
                              find_FM[4]])

        Jenkins_2 = np.array([find_FM[0],
                              find_FM[1],
                              find_FM[2]])

        Jenkins_3 = np.array([find_FM[0],
                              find_FM[1]])

        if find_FM[4] > 0:
            Jenkins_c = Jenkins_1
            print("Using Jenkins_1")
        elif find_FM[2] > 0 and find_FM[4] == 0:
            Jenkins_c = Jenkins_2
            print("Using Jenkins_2")
        elif find_FM[2] == 0 and find_FM[4] == 0:
            Jenkins_c = Jenkins_3
            print("Using Jenkins_3")

        #Jenkins_c = Jenkins_1
        Jenkins_shape = Jenkins_c.size
        Observed_1 = np.array([Observed_VEI[0],
                               Observed_VEI[1],
                               Observed_VEI[2],
                               Observed_VEI[3],
                               Observed_VEI[4]])

        Observed_2 = np.array([Observed_VEI[0],
                               Observed_VEI[1],
                               Observed_VEI[2]])

        Observed_3 = np.array([Observed_VEI[0],
                               Observed_VEI[1]])

        if find_FM[4] > 0:
            Observed = Observed_1
        elif find_FM[2] > 0 and find_FM[4] == 0:
            Observed = Observed_2
        elif find_FM[2] == 0 and find_FM[4] == 0:
            Observed = Observed_3

        #Observed = Observed_1
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

            if get_trace == True:
                arviz.plot_density(trace_VEI_GVP)
                pthfigVEI = "Japan/Figures/TracePlots/" + Volcano_name + "_densityplot.png"
                plt.savefig(pthfigVEI)
                if show_trace == True:
                    plt.show()
                else:
                    plt.close()
            else:
                pass
            # arviz.plot_autocorr(trace_VEI_GVP)
            # plt.show()
            VEI_trace_data = trace_VEI_GVP.get_values(parameters, burn=5000)
            VEI_trace_data_df = pd.DataFrame(VEI_trace_data)
            ppc_VEI = pm.sample_posterior_predictive(trace_VEI_GVP, var_names=['observed_data'], random_seed=RANDOM_SEED)


        if Jenkins_shape == 5:
            try:
                ppc_ave_VEI_df = pd.DataFrame(ppc_VEI)
                ppc_ave_VEI_df2 = pd.DataFrame(ppc_ave_VEI_df['observed_data'].to_list(),
                                               columns=['VEI3', 'VEI4', 'VEI5', 'VEI6', 'VEI7'])
                ppc_ave_VEI_df2 = ppc_ave_VEI_df2.fillna(0)
            except ValueError:
                ppc_ave_VEI_df = pd.DataFrame(ppc_VEI['observed_data'])
                ppc_ave_VEI_df2 = ppc_ave_VEI_df
                ppc_ave_VEI_df2.columns = ['VEI3', 'VEI4', 'VEI5', 'VEI6', 'VEI7']

            #arviz.plot_density(trace_VEI_GVP)
            #plt.show()

            Ave_VEI3 = VEI_trace_data_df.iloc[:, 0].mean()
            Std_VEI3 = VEI_trace_data_df.iloc[:, 0].std()

            Ave_VEI4 = VEI_trace_data_df.iloc[:, 1].mean()
            Std_VEI4 = VEI_trace_data_df.iloc[:, 1].std()

            Ave_VEI5 = VEI_trace_data_df.iloc[:, 2].mean()
            Std_VEI5 = VEI_trace_data_df.iloc[:, 2].std()

            Ave_VEI6 = VEI_trace_data_df.iloc[:, 3].mean()
            Std_VEI6 = VEI_trace_data_df.iloc[:, 3].std()

            Ave_VEI7 = VEI_trace_data_df.iloc[:, 4].mean()
            Std_VEI7 = VEI_trace_data_df.iloc[:, 4].std()

        elif Jenkins_shape == 3:
            try:
                ppc_ave_VEI_df = pd.DataFrame(ppc_VEI)
                ppc_ave_VEI_df2 = pd.DataFrame(ppc_ave_VEI_df['observed_data'].to_list(),
                                               columns=['VEI3', 'VEI4', 'VEI5'])
                ppc_ave_VEI_df2 = ppc_ave_VEI_df2.fillna(0)
            except ValueError:
                ppc_ave_VEI_df = pd.DataFrame(ppc_VEI['observed_data'])
                ppc_ave_VEI_df2 = ppc_ave_VEI_df
                ppc_ave_VEI_df2.columns = ['VEI3', 'VEI4', 'VEI5']

            # arviz.plot_density(trace_VEI_GVP)
            # plt.show()

            Ave_VEI3 = VEI_trace_data_df.iloc[:, 0].mean()
            Std_VEI3 = VEI_trace_data_df.iloc[:, 0].std()

            Ave_VEI4 = VEI_trace_data_df.iloc[:, 1].mean()
            Std_VEI4 = VEI_trace_data_df.iloc[:, 1].std()

            Ave_VEI5 = VEI_trace_data_df.iloc[:, 2].mean()
            Std_VEI5 = VEI_trace_data_df.iloc[:, 2].std()

            Ave_VEI6 = 0
            Std_VEI6 = 0

            Ave_VEI7 = 0
            Std_VEI7 = 0

        elif Jenkins_shape == 2:
            try:
                ppc_ave_VEI_df = pd.DataFrame(ppc_VEI)
                ppc_ave_VEI_df2 = pd.DataFrame(ppc_ave_VEI_df['observed_data'].to_list(),
                                               columns=['VEI3', 'VEI4'])
                ppc_ave_VEI_df2 = ppc_ave_VEI_df2.fillna(0)
            except ValueError:
                ppc_ave_VEI_df = pd.DataFrame(ppc_VEI['observed_data'])
                ppc_ave_VEI_df2 = ppc_ave_VEI_df
                ppc_ave_VEI_df2.columns = ['VEI3', 'VEI4']

            print(ppc_ave_VEI_df2)
            # arviz.plot_density(trace_VEI_GVP)
            # plt.show()

            Ave_VEI3 = VEI_trace_data_df.iloc[:, 0].mean()
            Std_VEI3 = VEI_trace_data_df.iloc[:, 0].std()

            Ave_VEI4 = VEI_trace_data_df.iloc[:, 1].mean()
            Std_VEI4 = VEI_trace_data_df.iloc[:, 1].std()

            Ave_VEI5 = 0
            Std_VEI5 = 0

            Ave_VEI6 = 0
            Std_VEI6 = 0

            Ave_VEI7 = 0
            Std_VEI7 = 0

        if Ave_VEI3 > 0:
            a_ave_VEI3 = Ave_VEI3 ** 2 * ((1 - Ave_VEI3) / Std_VEI3 ** 2 - 1 / Ave_VEI3)
            b_ave_VEI3 = a_ave_VEI3 * (1 / Ave_VEI3 - 1)
            data_beta_VEI3 = beta.rvs(a=a_ave_VEI3, b=b_ave_VEI3, size=1000000)
        else:
            a_ave_VEI3 = 0
            b_ave_VEI3 = 0
            data_beta_VEI3 = 0

        if Ave_VEI4 > 0:
            a_ave_VEI4 = Ave_VEI4 ** 2 * ((1 - Ave_VEI4) / Std_VEI4 ** 2 - 1 / Ave_VEI4)
            b_ave_VEI4 = a_ave_VEI4 * (1 / Ave_VEI4 - 1)
            data_beta_VEI4 = beta.rvs(a=a_ave_VEI4, b=b_ave_VEI4, size=1000000)
        else:
            a_ave_VEI4 = 0
            b_ave_VEI4 = 0
            data_beta_VEI4 = 0

        if Ave_VEI5 > 0:
            a_ave_VEI5 = Ave_VEI5 ** 2 * ((1 - Ave_VEI5) / Std_VEI5 ** 2 - 1 / Ave_VEI5)
            b_ave_VEI5 = a_ave_VEI5 * (1 / Ave_VEI5 - 1)
            data_beta_VEI5 = beta.rvs(a=a_ave_VEI5, b=b_ave_VEI5, size=1000000)
        else:
            a_ave_VEI5 = 0
            b_ave_VEI5 = 0
            data_beta_VEI5 = 0

        if Ave_VEI6 > 0:
            a_ave_VEI6 = Ave_VEI6 ** 2 * ((1 - Ave_VEI6) / Std_VEI6 ** 2 - 1 / Ave_VEI6)
            b_ave_VEI6 = a_ave_VEI6 * (1 / Ave_VEI6 - 1)
            data_beta_VEI6 = beta.rvs(a=a_ave_VEI6, b=b_ave_VEI6, size=1000000)
        else:
            a_ave_VEI6 = 0
            b_ave_VEI6 = 0
            data_beta_VEI6 = 0

        if Ave_VEI7 > 0:
            a_ave_VEI7 = Ave_VEI7 ** 2 * ((1 - Ave_VEI7) / Std_VEI7 ** 2 - 1 / Ave_VEI7)
            b_ave_VEI7 = a_ave_VEI7 * (1 / Ave_VEI7 - 1)
            data_beta_VEI7 = beta.rvs(a=a_ave_VEI7, b=b_ave_VEI7, size=1000000)
        else:
            a_ave_VEI7 = 0
            b_ave_VEI7 = 0
            data_beta_VEI7 = 0

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
        a_ave = AP_ave ** 2 * ((1 - AP_ave) / AP_std ** 2 - 1 / AP_ave)
        b_ave = a_ave * (1 / AP_ave - 1)

        if VEI_3_mean > 0:
            AP_VEI_3 = []
            Number = 10000
            lower, upper = 0, 1
            for i in range(Number):
                """Monte Carlo simulation"""
                AP_VEI3 = (np.random.beta(a_ave, b_ave) * stats.truncnorm.rvs((lower - VEI_3_mean) / VEI_3_sd,
                                                                              (upper - VEI_3_mean) / VEI_3_sd,
                                                                              loc=VEI_3_mean,
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
                                                                              (upper - VEI_4_mean) / VEI_4_sd,
                                                                              loc=VEI_4_mean,
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
                                                                              (upper - VEI_5_mean) / VEI_5_sd,
                                                                              loc=VEI_5_mean,
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
                                                                              (upper - VEI_6_mean) / VEI_6_sd,
                                                                              loc=VEI_6_mean,
                                                                              scale=VEI_6_sd))
                AP_VEI_6 = AP_VEI_6 + [AP_VEI6]
        else:
            AP_VEI_6 = 0

        if VEI_7_mean > 0:
            AP_VEI_7 = []
            Number = 10000
            lower, upper = 0, 1
            for i in range(Number):
                """Monte Carlo simulation"""
                AP_VEI7 = (np.random.beta(a_ave, b_ave) * stats.truncnorm.rvs((lower - VEI_7_mean) / VEI_7_sd,
                                                                              (upper - VEI_7_mean) / VEI_7_sd,
                                                                              loc=VEI_7_mean,
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

    Volcano_data = GVP_volcanoes[(GVP_volcanoes['Volcano Number'] == volc_num)]
    Volcano_name = Volcano_data.iloc[0]['Volcano Name']
    Volcano_Number = volc_num
    textstr = Volcano_name

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

    Probabilities_df = pd.DataFrame([[Volcano_Number, Volcano_name, Method,
                                      Prob_erupt_Percentile_05, Prob_erupt_Median, Prob_erupt_Percentile_95,
                                      VEI_3_Percentile_05, VEI_3_Median, VEI_3_Percentile_95,
                                      VEI_4_Percentile_05, VEI_4_Median, VEI_4_Percentile_95,
                                      VEI_5_Percentile_05, VEI_5_Median, VEI_5_Percentile_95,
                                      VEI_6_Percentile_05, VEI_6_Median, VEI_6_Percentile_95,
                                      VEI_7_Percentile_05, VEI_7_Median, VEI_7_Percentile_95]],
                                    columns=["Volcano Number",
                                             "Volcano Name",
                                             "Estimate method",
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
    ID = str(volc_num)
    path_csv = "Japan/csv/Probabilities_" + Volcano_name + ".csv"
    Probabilities_df.to_csv(path_csv,  index=False)

    # ax = sns.lineplot(x="VEI", y="Annual probability", markers=True, data=Probabilities)
    # ax.set(yscale="log")
    # ax.set_xticks(range(3, 8)) # <--- set the ticks first
    # ax.tick_params(labelsize=10)
    # ax.set_xticklabels(['<=3','4','5','6','7'])
    # ax.set_xlabel("VEI",fontsize=12)
    # ax.set_ylabel("Annual probability",fontsize=12)
    # props = dict(boxstyle='round', facecolor='white', alpha=0.5)
    # # place a text box in upper left in axes coords
    # ax.text(0.80, 0.95, textstr, transform=ax.transAxes, fontsize=12,
    #         verticalalignment='top', bbox=props)
    # path_fig = "Japan/" + Volcano_name + ".png"
    # plt.savefig(path_fig)
    # plt.show()

    return(Probabilities)

def save_plot (Volcano_name, volc_num, FM, show_fig):
    print("Volcano Number: ", volc_num)
    Probabilities = FM
    textstr = Volcano_name
    #path_csv = "Japan/Probabilities_" + Volcano_name + ".csv"
    #Probabilities_df.to_csv(path_csv,  index=False)
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
    path_fig = "Japan/figures/" + Volcano_name + ".png"
    plt.savefig(path_fig)
    if show_fig == True:
        plt.show()
    else:
        plt.close()

def get_master_csv (path_csv, save_file_name):
    import glob
    all_files = glob.glob(path_csv + "*.csv")
    file_list = []
    for filename in all_files:
        df = pd.read_csv(filename, index_col=None, header=0)
        file_list.append(df)
    frame = pd.concat(file_list, axis=0, ignore_index=True)
    master_FM = frame.to_csv(save_file_name, index=False)

    return()


def get_model_average (GVP_volcanoes, volc_num, Volcano_type_1, Volcano_type_2, Observed_rate, Observed_VEI,
                       confirmed_VEI_eruptions, Freq_rate_1, Freq_rate_2, Freq_rate_std_1, Freq_rate_std_2, VEI_freq_1,
                       VEI_freq_2, Power_law, averaging_method, year, confirmed, method, project_folder, VEI_schema, Uncertainty_schema, show_trace):
    AP_VEI_3 = None
    AP_VEI_4 = None
    AP_VEI_5 = None
    AP_VEI_6 = None
    AP_VEI_7 = None
    yearstr = str(year)
    if confirmed == "True":
        eruption_certainty = "Confirmed"
    else:
        eruption_certainty = "Confirmed + Uncertain"
    alpha_analogue_1 = Freq_rate_1[Volcano_type_1] ** 2 * (
            (1 - Freq_rate_1[Volcano_type_1]) / Freq_rate_std_1[Volcano_type_1] ** 2 - 1 / Freq_rate_1[Volcano_type_1])
    beta_analogue_1 = alpha_analogue_1 * (1 / Freq_rate_1[Volcano_type_1] - 1)

    if confirmed_VEI_eruptions == 0:
        print("No confirmed VEI eruptions")

    else:
        Analysis_method = "Bayesian model average"
        print("Initiating model using: ", Analysis_method)

        # Model averaging for frequency of eruption of any magnitude
        ##Comparing the models and assigning a weight.
        ### Model 1
        # Probability of an eruption of any VEI
        cores = 1
        y = np.array([Observed_rate * 10000])
        n = np.array([10000])
        N = len(n)
        # if __name__ == '__main__':
        # Probability of eruption model
        with pm.Model() as AP_model_1:
            theta = pm.Beta('theta', alpha=alpha_analogue_1, beta=beta_analogue_1)
            p = pm.Binomial('y', p=theta, observed=y, n=n)
            AP_trace_1 = pm.sample(20000, tune=10000, chains=2, cores=cores, target_accept=0.99)
            arviz.plot_trace(AP_trace_1, legend=True)
            if show_trace == "True":
                plt.show()
            else:
                plt.close()

        ### Model 2
        alpha_analogue_2 = Freq_rate_2[Volcano_type_2] ** 2 * (
                (1 - Freq_rate_2[Volcano_type_2]) / Freq_rate_std_2[Volcano_type_2] ** 2 - 1 / Freq_rate_2[Volcano_type_2])
        beta_analogue_2 = alpha_analogue_2 * (1 / Freq_rate_2[Volcano_type_2] - 1)
        with pm.Model() as AP_model_2:
            theta = pm.Beta('theta', alpha=alpha_analogue_2, beta=beta_analogue_2)
            p = pm.Binomial('y', p=theta, observed=y, n=n)
            AP_trace_2 = pm.sample(20000, tune=10000, chains=2, cores=cores, target_accept=0.99)
            arviz.plot_trace(AP_trace_2, legend=True)
            if show_trace == "True":
                plt.show()
            else:
                plt.close()

        traces = [AP_trace_1, AP_trace_2]
        model_dict = dict(zip([AP_model_1, AP_model_2], traces))
        comp = arviz.compare(model_dict, method=averaging_method)
        comp['index'] = [0, 1]
        comp.set_index("index", inplace = True)
        ppc_ave = pm.sample_posterior_predictive_w(traces, 100000, [AP_model_1, AP_model_2],
                                        weights=comp.weight.sort_index(ascending=True),
                                        progressbar=False)


        ppc_ave_adjust = ppc_ave['y']/10000
        Ave = ppc_ave_adjust.mean()
        Std = ppc_ave_adjust.std()
        a_ave = Ave**2 * ((1-Ave)/Std**2 -1/Ave)
        b_ave = a_ave * (1/Ave-1)
        data_beta = beta.rvs(a=a_ave, b=b_ave, size=1000000)

        # # Comparing conditional probability of each VEI eruption

        # # # Model 1
        # print("VEI_freq_1 4", VEI_freq_1["4"])

        if VEI_schema == 1:

            # VEI <=3, 4, 5, 6, 7
            if VEI_freq_1["7"].mean() > 0:
                print("Using first concentration parameter for Model 1")
                concentration_1 = np.array([VEI_freq_1["3"].mean() * 10000,
                                            VEI_freq_1["4"].mean() * 10000,
                                            VEI_freq_1["5"].mean() * 10000,
                                            VEI_freq_1["6"].mean() * 10000,
                                            VEI_freq_1["7"].mean() * 10000])
                Observed = np.array([Observed_VEI[0],
                                     Observed_VEI[1],
                                     Observed_VEI[2],
                                     Observed_VEI[3],
                                     Observed_VEI[4]])

            elif VEI_freq_1["6"].mean() > 0 and VEI_freq_1["7"].mean() == 0:
                print("Using second concentration parameter for Model 1")
                concentration_1 = np.array([VEI_freq_1["3"] * 10000,
                                            VEI_freq_1["4"] * 10000,
                                            VEI_freq_1["5"] * 10000,
                                            VEI_freq_1["6"] * 10000])
                Observed = np.array([Observed_VEI[0],
                                     Observed_VEI[1],
                                     Observed_VEI[2],
                                     Observed_VEI[3]])

            elif VEI_freq_1["5"].mean() > 0 and VEI_freq_1["6"].mean() == 0:
                print("Using third concentration parameter for Model 1")
                concentration_1 = np.array([VEI_freq_1["3"].mean() * 10000,
                                            VEI_freq_1["4"].mean() * 10000,
                                            VEI_freq_1["5"].mean() * 10000])
                Observed = np.array([Observed_VEI[0],
                                     Observed_VEI[1],
                                     Observed_VEI[2]])

            elif VEI_freq_1["4"].mean() > 0 and VEI_freq_1["5"].mean() == 0:
                print("Using fourth concentration parameter for Model 1")
                concentration_1 = np.array([VEI_freq_1["3"].mean() * 10000,
                                            VEI_freq_1["4"].mean() * 10000])
                Observed = np.array([Observed_VEI[0],
                                     Observed_VEI[1]])

            else:
                print("Error with concentration parameter or observed parameter for Model 1")
            # Jenkins_c = Jenkins_1
            shape_1 = concentration_1.size


            Observed_sum = np.sum(Observed)

            # Create model 1
            #try:
            Create_graph = "Yes"
            with pm.Model() as VEI_model_1:
                parameters = pm.Dirichlet('parameters', a=concentration_1, shape=shape_1)
                # Observed data is from a Multinomial distribution
                observed_data = pm.Multinomial(
                    'observed_data', n=Observed_sum, p=parameters, shape=shape_1, observed=Observed)

            with VEI_model_1:
                # Sample from the posterior
                VEI_trace_1 = pm.sample(draws=10000, chains=2, cores=cores, tune=5000,
                                          discard_tuned_samples=True, target_accept=0.99)
            print(VEI_trace_1)
            arviz.plot_trace(VEI_trace_1, legend=True)
            if show_trace == "True":
                plt.show()
            else:
                plt.close()
            # Create model 2

            if VEI_freq_2["7"].mean() > 0:
                print("Using first concentration parameter for Model 2")

                concentration_2 = np.array([VEI_freq_2["3"].mean() * 10000,
                                            VEI_freq_2["4"].mean() * 10000,
                                            VEI_freq_2["5"].mean() * 10000,
                                            VEI_freq_2["6"].mean() * 10000,
                                            VEI_freq_2["7"].mean() * 10000])
                Observed = np.array([Observed_VEI[0],
                                     Observed_VEI[1],
                                     Observed_VEI[2],
                                     Observed_VEI[3],
                                     Observed_VEI[4]])

            elif VEI_freq_2["6"].mean() > 0 and VEI_freq_2["7"].mean() == 0:
                print("Using second concentration parameter for Model 2")
                concentration_2 = np.array([VEI_freq_2["3"].mean() * 10000,
                                            VEI_freq_2["4"].mean() * 10000,
                                            VEI_freq_2["5"].mean() * 10000,
                                            VEI_freq_2["6"].mean() * 10000])
                Observed = np.array([Observed_VEI[0],
                                     Observed_VEI[1],
                                     Observed_VEI[2],
                                     Observed_VEI[3]])

            elif VEI_freq_2["5"].mean() > 0 and VEI_freq_2["6"].mean() == 0:
                print("Using third concentration parameter for Model 2")
                concentration_2 = np.array([VEI_freq_2["3"].mean() * 10000,
                                            VEI_freq_2["4"].mean() * 10000,
                                            VEI_freq_2["5"].mean() * 10000])
                Observed = np.array([Observed_VEI[0],
                                     Observed_VEI[1],
                                     Observed_VEI[2]])

            elif VEI_freq_2["4"].mean() > 0 and VEI_freq_2["4"].mean() == 0:
                print("Using fourth concentration parameter for Model 2")
                concentration_2 = np.array([VEI_freq_2["3"].mean() * 10000,
                                            VEI_freq_2["4"].mean() * 10000])
                Observed = np.array([Observed_VEI[0],
                                     Observed_VEI[1]])

            else:
                print("Error with concentration parameter or observed parameter for Model 2")

            # Jenkins_c = Jenkins_1
            shape_2 = concentration_2.size

            with pm.Model() as VEI_model_2:
                parameters = pm.Dirichlet('parameters', a=concentration_2, shape=shape_2)
                # Observed data is from a Multinomial distribution
                observed_data = pm.Multinomial(
                    'observed_data', n=Observed_sum, p=parameters, shape=shape_2, observed=Observed)

            with VEI_model_2:
                # Sample from the posterior
                VEI_trace_2 = pm.sample(draws=10000, chains=2, cores=cores, tune=5000,
                                          discard_tuned_samples=True, target_accept=0.99)

            arviz.plot_trace(VEI_trace_2, legend=True)
            if show_trace == "True":
                plt.show()
            else:
                plt.close()
            traces_VEI = [VEI_trace_1, VEI_trace_2]
            model_dict_VEI = dict(zip([VEI_model_1, VEI_model_2], traces_VEI))
            comp_VEI = arviz.compare(model_dict_VEI, method=averaging_method) #BB-pseudo-BMA
            comp_VEI['index'] = [0, 1]
            comp_VEI.set_index("index", inplace = True)

            ppc_ave_VEI = pm.sample_posterior_predictive_w(traces_VEI, 1000, [VEI_model_1, VEI_model_2],
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
                if VEI_freq_2["7"].mean() > 0:
                    ppc_ave_VEI_df2.columns = ['VEI3', 'VEI4', 'VEI5', 'VEI6', 'VEI7']
                elif VEI_freq_2["5"].mean() > 0 and VEI_freq_2["7"].mean() == 0:
                    if len(ppc_ave_VEI_df2.columns) < 5:
                        print("Number of columns is:", len(ppc_ave_VEI_df2.columns))
                        ppc_ave_VEI_df2.columns = ['VEI3', 'VEI4', 'VEI5']
                        ppc_ave_VEI_df2['VEI6'] = 0
                        ppc_ave_VEI_df2['VEI7'] = 0
                    else:
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

            if Ave_VEI3 > 0:
                a_ave_VEI3 = Ave_VEI3**2 * ((1-Ave_VEI3)/Std_VEI3**2 -1/Ave_VEI3)
                b_ave_VEI3 = a_ave_VEI3 * (1/Ave_VEI3-1)
                data_beta_VEI3 = beta.rvs(a=a_ave_VEI3, b=b_ave_VEI3, size=1000000)
            else:
                a_ave_VEI3 = 0
                b_ave_VEI3 = 0
                data_beta_VEI3 = 0

            if Ave_VEI4 > 0:
                a_ave_VEI4 = Ave_VEI4**2 * ((1-Ave_VEI4)/Std_VEI4**2 -1/Ave_VEI4)
                b_ave_VEI4 = a_ave_VEI4 * (1/Ave_VEI4-1)
                data_beta_VEI4 = beta.rvs(a=a_ave_VEI4, b=b_ave_VEI4, size=1000000)
            else:
                a_ave_VEI4 = 0
                b_ave_VEI4 = 0
                data_beta_VEI4 = 0

            if Ave_VEI5 > 0:
                a_ave_VEI5 = Ave_VEI5**2 * ((1-Ave_VEI5)/Std_VEI5**2 -1/Ave_VEI5)
                b_ave_VEI5 = a_ave_VEI5 * (1/Ave_VEI5-1)
                data_beta_VEI5 = beta.rvs(a=a_ave_VEI5, b=b_ave_VEI5, size=1000000)
            else:
                a_ave_VEI5 = 0
                b_ave_VEI5 = 0
                data_beta_VEI5 = 0

            if Ave_VEI6 > 0:
                a_ave_VEI6 = Ave_VEI6**2 * ((1-Ave_VEI6)/Std_VEI6**2 -1/Ave_VEI6)
                b_ave_VEI6 = a_ave_VEI6 * (1/Ave_VEI6-1)
                data_beta_VEI6 = beta.rvs(a=a_ave_VEI6, b=b_ave_VEI6, size=1000000)
            else:
                a_ave_VEI6 = 0
                b_ave_VEI6 = 0
                data_beta_VEI6 = 0

            if Ave_VEI7 > 0:
                a_ave_VEI7 = Ave_VEI7**2 * ((1-Ave_VEI7)/Std_VEI7**2 -1/Ave_VEI7)
                b_ave_VEI7 = a_ave_VEI7 * (1/Ave_VEI7-1)
                data_beta_VEI7 = beta.rvs(a=a_ave_VEI7, b=b_ave_VEI7, size=1000000)
            else:
                a_ave_VEI7 = 0
                b_ave_VEI7 = 0
                data_beta_VEI7 = 0

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


            # Combining probability of eruption and relative probability of different VEI

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
                AP_VEI_6 = 0

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
            # Visualisation of data and production of csv files.
            #####

            if Uncertainty_schema == 1:
                Percentile_1 = 5
                p1_string = "5th percentile"
                Percentile_2 = 50
                p2_string = "50th percentile"
                Percentile_3 = 95
                p3_string = "95th percentile"
            elif Uncertainty_schema == 2:
                Percentile_1 = 10
                p1_string = "10th percentile"
                Percentile_2 = 50
                p2_string = "50th percentile"
                Percentile_3 = 90
                p3_string = "590th percentile"
            elif Uncertainty_schema == 3:
                Percentile_1 = 25
                p1_string = "25th percentile"
                Percentile_2 = 50
                p2_string = "50th percentile"
                Percentile_3 = 75
                p3_string = "75th percentile"

            VEI_3_Median = stats.scoreatpercentile(AP_VEI_3, Percentile_2)
            VEI_3_Percentile_1 = stats.scoreatpercentile(AP_VEI_3, Percentile_1)
            VEI_3_Percentile_3 = stats.scoreatpercentile(AP_VEI_3, Percentile_3)

            VEI_4_Median = stats.scoreatpercentile(AP_VEI_4, Percentile_2)
            VEI_4_Percentile_1 = stats.scoreatpercentile(AP_VEI_4, Percentile_1)
            VEI_4_Percentile_3 = stats.scoreatpercentile(AP_VEI_4, Percentile_3)

            VEI_5_Median = stats.scoreatpercentile(AP_VEI_5, Percentile_2)
            VEI_5_Percentile_1 = stats.scoreatpercentile(AP_VEI_5, Percentile_1)
            VEI_5_Percentile_3 = stats.scoreatpercentile(AP_VEI_5, Percentile_3)

            VEI_6_Median = stats.scoreatpercentile(AP_VEI_6, Percentile_2)
            VEI_6_Percentile_1 = stats.scoreatpercentile(AP_VEI_6, Percentile_1)
            VEI_6_Percentile_3 = stats.scoreatpercentile(AP_VEI_6, Percentile_3)

            VEI_7_Median = stats.scoreatpercentile(AP_VEI_7, Percentile_2)
            VEI_7_Percentile_1 = stats.scoreatpercentile(AP_VEI_7, Percentile_2)
            VEI_7_Percentile_3 = stats.scoreatpercentile(AP_VEI_7, Percentile_3)

            Prob_erupt_Median = VEI_3_Median + VEI_4_Median + VEI_5_Median + VEI_6_Median + VEI_7_Median
            Prob_erupt_Percentile_05 = VEI_3_Percentile_1 + VEI_4_Percentile_1 + VEI_5_Percentile_1 + VEI_6_Percentile_1 + VEI_7_Percentile_1
            Prob_erupt_Percentile_95 = VEI_3_Percentile_3 + VEI_4_Percentile_3 + VEI_5_Percentile_3 + VEI_6_Percentile_3 + VEI_7_Percentile_3

            Median = "50th percentile"
            Percentile_05 = "10th percentile"
            Percentile_95 = "90th percentile"

            Percentile_05s = [[3, 4, 5, 6, 7], [VEI_3_Percentile_1, VEI_4_Percentile_1, VEI_5_Percentile_1, VEI_6_Percentile_1, VEI_7_Percentile_1]]
            Percentile_95s = [[3, 4, 5, 6, 7], [VEI_3_Percentile_3, VEI_4_Percentile_3, VEI_5_Percentile_3, VEI_6_Percentile_3, VEI_7_Percentile_3]]

            Probabilities = pd.DataFrame([[3, VEI_3_Median],
                                    [4, VEI_4_Median],
                                    [5, VEI_5_Median],
                                    [6, VEI_6_Median],
                                    [7, VEI_7_Median],
                                    [3, VEI_3_Percentile_1],
                                    [4, VEI_4_Percentile_1],
                                    [5, VEI_5_Percentile_1],
                                    [6, VEI_6_Percentile_1],
                                    [7, VEI_7_Percentile_1],
                                    [3, VEI_3_Percentile_3],
                                    [4, VEI_4_Percentile_3],
                                    [5, VEI_5_Percentile_3],
                                    [6, VEI_6_Percentile_3],
                                    [7, VEI_7_Percentile_3]],
                                   columns=["VEI", "Annual probability"])

            Volcano_data = GVP_volcanoes[(GVP_volcanoes['Volcano Number'] == volc_num)]
            Volcano_name = Volcano_data.iloc[0]['Volcano Name']
            Volcano_Number = volc_num
            textstr = Volcano_name

            Probabilities_df = pd.DataFrame([[Volcano_Number, Volcano_name, year, eruption_certainty, Analysis_method, method, averaging_method, Power_law,
                                              Prob_erupt_Percentile_05, Prob_erupt_Median, Prob_erupt_Percentile_95,
                                              VEI_3_Percentile_1, VEI_3_Median, VEI_3_Percentile_3,
                                              VEI_4_Percentile_1, VEI_4_Median, VEI_4_Percentile_3,
                                              VEI_5_Percentile_1, VEI_5_Median, VEI_5_Percentile_3,
                                              VEI_6_Percentile_1, VEI_6_Median, VEI_6_Percentile_3,
                                              VEI_7_Percentile_1, VEI_7_Median, VEI_7_Percentile_3]],
                                            columns=["Volcano Number",
                                                     "Volcano Name",
                                                     "GVP DB Year",
                                                     "Included eruptions",
                                                     "Estimate method",
                                                     "Change point",
                                                     "Average method",
                                                     "Power law",
                                                     "Probability of eruption " + p1_string,
                                                     "Probability of eruption " + p2_string,
                                                     "Probability of eruption " + p3_string,
                                                     "<= VEI 3 " + p1_string,
                                                     "<= VEI 3 " + p2_string,
                                                     "<= VEI 3 " + p3_string,
                                                     "VEI 4 " + p1_string,
                                                     "VEI 4 " + p2_string,
                                                     "VEI 4 " + p3_string,
                                                     "VEI 5 " + p1_string,
                                                     "VEI 5 " + p2_string,
                                                     "VEI 5 " + p3_string,
                                                     "VEI 6 " + p1_string,
                                                     "VEI 6 " + p2_string,
                                                     "VEI 6 " + p3_string,
                                                     "VEI 7 " + p1_string,
                                                     "VEI 7 " + p2_string,
                                                     "VEI 7 " + p3_string])
            ID = str(volc_num)
            if confirmed == "True":
                certainty = "Certain"
            else:
                certainty = "Uncertain"
            if Power_law == "True":
                power = "PowerLaw"
            else:
                power = "NoPowerLaw"

            path_csv = "Probabilities/" + project_folder + "/" + Volcano_name + "_" + certainty + "_" + yearstr + "_" + method + "_" + averaging_method + "_" + power + ".csv"
            Probabilities_df.to_csv(path_csv,  index=False)

            if Create_graph == "No":
                print("No graph created due to sampling error")
            else:
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
                path_fig = "Figures/" + project_folder + "/" + Volcano_name + "_" + certainty + "_" + yearstr + "_" + method + "_" + "_" + averaging_method + "_" + power + ".png"
                plt.savefig(path_fig)
                plt.close()
                #plt.show()


        if VEI_schema == 2:
            # VEI <=2, 3, 4, 5, 6, 7
            if VEI_freq_1["7"].mean() > 0:
                print("Using first concentration parameter for Model 1")
                concentration_1 = np.array([VEI_freq_1["2"].mean() * 10000,
                                            VEI_freq_1["3"].mean() * 10000,
                                            VEI_freq_1["4"].mean() * 10000,
                                            VEI_freq_1["5"].mean() * 10000,
                                            VEI_freq_1["6"].mean() * 10000,
                                            VEI_freq_1["7"].mean() * 10000])
                Observed = np.array([Observed_VEI[0],
                                     Observed_VEI[1],
                                     Observed_VEI[2],
                                     Observed_VEI[3],
                                     Observed_VEI[4],
                                     Observed_VEI[5]])

            elif VEI_freq_1["6"].mean() > 0 and VEI_freq_1["7"].mean() == 0:
                print("Using second concentration parameter for Model 1")
                concentration_1 = np.array([VEI_freq_1["2"] * 10000,
                                            VEI_freq_1["3"] * 10000,
                                            VEI_freq_1["4"] * 10000,
                                            VEI_freq_1["5"] * 10000,
                                            VEI_freq_1["6"] * 10000])
                Observed = np.array([Observed_VEI[0],
                                     Observed_VEI[1],
                                     Observed_VEI[2],
                                     Observed_VEI[3],
                                     Observed_VEI[4]])

            elif VEI_freq_1["5"].mean() > 0 and VEI_freq_1["6"].mean() == 0:
                print("Using third concentration parameter for Model 1")
                concentration_1 = np.array([VEI_freq_1["2"].mean() * 10000,
                                            VEI_freq_1["3"].mean() * 10000,
                                            VEI_freq_1["4"].mean() * 10000,
                                            VEI_freq_1["5"].mean() * 10000])
                Observed = np.array([Observed_VEI[0],
                                     Observed_VEI[1],
                                     Observed_VEI[2],
                                     Observed_VEI[3]])

            elif VEI_freq_1["4"].mean() > 0 and VEI_freq_1["5"].mean() == 0:
                print("Using fourth concentration parameter for Model 1")
                concentration_1 = np.array([VEI_freq_1["2"].mean() * 10000,
                                            VEI_freq_1["3"].mean() * 10000,
                                            VEI_freq_1["4"].mean() * 10000])
                Observed = np.array([Observed_VEI[0],
                                     Observed_VEI[1],
                                     Observed_VEI[2]])

            else:
                print("Error with concentration parameter or observed parameter for Model 1")
            # Jenkins_c = Jenkins_1
            shape_1 = concentration_1.size

            Observed_sum = np.sum(Observed)

            # Create model 1
            # try:
            Create_graph = "Yes"
            with pm.Model() as VEI_model_1:
                parameters = pm.Dirichlet('parameters', a=concentration_1, shape=shape_1)
                # Observed data is from a Multinomial distribution
                observed_data = pm.Multinomial(
                    'observed_data', n=Observed_sum, p=parameters, shape=shape_1, observed=Observed)

            with VEI_model_1:
                # Sample from the posterior
                VEI_trace_1 = pm.sample(draws=10000, chains=2, cores=cores, tune=5000,
                                        discard_tuned_samples=True, target_accept=0.99)

            # Create model 2

            if VEI_freq_2["7"].mean() > 0:
                print("Using first concentration parameter for Model 2")

                concentration_2 = np.array([VEI_freq_2["2"].mean() * 10000,
                                            VEI_freq_2["3"].mean() * 10000,
                                            VEI_freq_2["4"].mean() * 10000,
                                            VEI_freq_2["5"].mean() * 10000,
                                            VEI_freq_2["6"].mean() * 10000,
                                            VEI_freq_2["7"].mean() * 10000])
                Observed = np.array([Observed_VEI[0],
                                     Observed_VEI[1],
                                     Observed_VEI[2],
                                     Observed_VEI[3],
                                     Observed_VEI[4],
                                     Observed_VEI[5]])

            elif VEI_freq_2["6"].mean() > 0 and VEI_freq_2["7"].mean() == 0:
                print("Using second concentration parameter for Model 2")
                concentration_2 = np.array([VEI_freq_2["2"].mean() * 10000,
                                            VEI_freq_2["3"].mean() * 10000,
                                            VEI_freq_2["4"].mean() * 10000,
                                            VEI_freq_2["5"].mean() * 10000,
                                            VEI_freq_2["6"].mean() * 10000])
                Observed = np.array([Observed_VEI[0],
                                     Observed_VEI[1],
                                     Observed_VEI[2],
                                     Observed_VEI[3],
                                     Observed_VEI[4]])

            elif VEI_freq_2["5"].mean() > 0 and VEI_freq_2["6"].mean() == 0:
                print("Using third concentration parameter for Model 2")
                concentration_2 = np.array([VEI_freq_2["2"].mean() * 10000,
                                            VEI_freq_2["3"].mean() * 10000,
                                            VEI_freq_2["4"].mean() * 10000,
                                            VEI_freq_2["5"].mean() * 10000])
                Observed = np.array([Observed_VEI[0],
                                     Observed_VEI[1],
                                     Observed_VEI[2],
                                     Observed_VEI[3]])

            elif VEI_freq_2["4"].mean() > 0 and VEI_freq_2["4"].mean() == 0:
                print("Using fourth concentration parameter for Model 2")
                concentration_2 = np.array([VEI_freq_2["2"].mean() * 10000,
                                            VEI_freq_2["3"].mean() * 10000,
                                            VEI_freq_2["4"].mean() * 10000])
                Observed = np.array([Observed_VEI[0],
                                     Observed_VEI[1],
                                     Observed_VEI[2]])

            else:
                print("Error with concentration parameter or observed parameter for Model 2")

            # Jenkins_c = Jenkins_1
            shape_2 = concentration_2.size

            with pm.Model() as VEI_model_2:
                parameters = pm.Dirichlet('parameters', a=concentration_2, shape=shape_2)
                # Observed data is from a Multinomial distribution
                observed_data = pm.Multinomial(
                    'observed_data', n=Observed_sum, p=parameters, shape=shape_2, observed=Observed)

            with VEI_model_2:
                # Sample from the posterior
                VEI_trace_2 = pm.sample(draws=10000, chains=2, cores=cores, tune=5000,
                                        discard_tuned_samples=True, target_accept=0.99)

            traces_VEI = [VEI_trace_1, VEI_trace_2]
            model_dict_VEI = dict(zip([VEI_model_1, VEI_model_2], traces_VEI))
            comp_VEI = arviz.compare(model_dict_VEI, method=averaging_method)  # BB-pseudo-BMA
            comp_VEI['index'] = [0, 1]
            comp_VEI.set_index("index", inplace=True)

            ppc_ave_VEI = pm.sample_posterior_predictive_w(traces_VEI, 1000, [VEI_model_1, VEI_model_2],
                                                           weights=comp.weight.sort_index(ascending=True),
                                                           progressbar=False)

            try:
                ppc_ave_VEI_df = pd.DataFrame(ppc_ave_VEI)
                ppc_ave_VEI_df2 = pd.DataFrame(ppc_ave_VEI_df['observed_data'].to_list(),
                                               columns=['VEI2', 'VEI3', 'VEI4', 'VEI5', 'VEI6', 'VEI7'])
                ppc_ave_VEI_df2 = ppc_ave_VEI_df2.fillna(0)
            except ValueError:
                ppc_ave_VEI_df = pd.DataFrame(ppc_ave_VEI['observed_data'])
                ppc_ave_VEI_df2 = ppc_ave_VEI_df
                if VEI_freq_2["7"].mean() > 0:
                    ppc_ave_VEI_df2.columns = ['VEI2', 'VEI3', 'VEI4', 'VEI5', 'VEI6', 'VEI7']
                elif VEI_freq_2["5"].mean() > 0 and VEI_freq_2["7"].mean() == 0:
                    if len(ppc_ave_VEI_df2.columns) < 5:
                        print("Number of columns is:", len(ppc_ave_VEI_df2.columns))
                        ppc_ave_VEI_df2.columns = ['VEI2', 'VEI3', 'VEI4', 'VEI5']
                        ppc_ave_VEI_df2['VEI6'] = 0
                        ppc_ave_VEI_df2['VEI7'] = 0
                    else:
                        ppc_ave_VEI_df2.columns = ['VEI2', 'VEI3', 'VEI4', 'VEI5', 'VEI6', 'VEI7']

            ppc_ave_VEI_df2['total'] = ppc_ave_VEI_df2['VEI2'] + ppc_ave_VEI_df2['VEI3'] + \
                                       ppc_ave_VEI_df2['VEI4'] + ppc_ave_VEI_df2['VEI5'] + \
                                       ppc_ave_VEI_df2['VEI6'] + ppc_ave_VEI_df2['VEI7']
            ppc_ave_VEI_df2['pVEI2'] = ppc_ave_VEI_df2['VEI2'] / ppc_ave_VEI_df2['total']
            ppc_ave_VEI_df2['pVEI3'] = ppc_ave_VEI_df2['VEI3'] / ppc_ave_VEI_df2['total']
            ppc_ave_VEI_df2['pVEI4'] = ppc_ave_VEI_df2['VEI4'] / ppc_ave_VEI_df2['total']
            ppc_ave_VEI_df2['pVEI5'] = ppc_ave_VEI_df2['VEI5'] / ppc_ave_VEI_df2['total']
            ppc_ave_VEI_df2['pVEI6'] = ppc_ave_VEI_df2['VEI6'] / ppc_ave_VEI_df2['total']
            ppc_ave_VEI_df2['pVEI7'] = ppc_ave_VEI_df2['VEI7'] / ppc_ave_VEI_df2['total']

            Ave_VEI2 = ppc_ave_VEI_df2['pVEI2'].mean()
            Std_VEI2 = ppc_ave_VEI_df2['pVEI2'].std()

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

            if Ave_VEI2 > 0:
                a_ave_VEI2 = Ave_VEI2 ** 2 * ((1 - Ave_VEI2) / Std_VEI2 ** 2 - 1 / Ave_VEI2)
                b_ave_VEI2 = a_ave_VEI2 * (1 / Ave_VEI2 - 1)
                data_beta_VEI2 = beta.rvs(a=a_ave_VEI2, b=b_ave_VEI2, size=1000000)
            else:
                a_ave_VEI2 = 0
                b_ave_VEI2 = 0
                data_beta_VEI2 = 0

            if Ave_VEI3 > 0:
                a_ave_VEI3 = Ave_VEI3 ** 2 * ((1 - Ave_VEI3) / Std_VEI3 ** 2 - 1 / Ave_VEI3)
                b_ave_VEI3 = a_ave_VEI3 * (1 / Ave_VEI3 - 1)
                data_beta_VEI3 = beta.rvs(a=a_ave_VEI3, b=b_ave_VEI3, size=1000000)
            else:
                a_ave_VEI3 = 0
                b_ave_VEI3 = 0
                data_beta_VEI3 = 0

            if Ave_VEI4 > 0:
                a_ave_VEI4 = Ave_VEI4 ** 2 * ((1 - Ave_VEI4) / Std_VEI4 ** 2 - 1 / Ave_VEI4)
                b_ave_VEI4 = a_ave_VEI4 * (1 / Ave_VEI4 - 1)
                data_beta_VEI4 = beta.rvs(a=a_ave_VEI4, b=b_ave_VEI4, size=1000000)
            else:
                a_ave_VEI4 = 0
                b_ave_VEI4 = 0
                data_beta_VEI4 = 0

            if Ave_VEI5 > 0:
                a_ave_VEI5 = Ave_VEI5 ** 2 * ((1 - Ave_VEI5) / Std_VEI5 ** 2 - 1 / Ave_VEI5)
                b_ave_VEI5 = a_ave_VEI5 * (1 / Ave_VEI5 - 1)
                data_beta_VEI5 = beta.rvs(a=a_ave_VEI5, b=b_ave_VEI5, size=1000000)
            else:
                a_ave_VEI5 = 0
                b_ave_VEI5 = 0
                data_beta_VEI5 = 0

            if Ave_VEI6 > 0:
                a_ave_VEI6 = Ave_VEI6 ** 2 * ((1 - Ave_VEI6) / Std_VEI6 ** 2 - 1 / Ave_VEI6)
                b_ave_VEI6 = a_ave_VEI6 * (1 / Ave_VEI6 - 1)
                data_beta_VEI6 = beta.rvs(a=a_ave_VEI6, b=b_ave_VEI6, size=1000000)
            else:
                a_ave_VEI6 = 0
                b_ave_VEI6 = 0
                data_beta_VEI6 = 0

            if Ave_VEI7 > 0:
                a_ave_VEI7 = Ave_VEI7 ** 2 * ((1 - Ave_VEI7) / Std_VEI7 ** 2 - 1 / Ave_VEI7)
                b_ave_VEI7 = a_ave_VEI7 * (1 / Ave_VEI7 - 1)
                data_beta_VEI7 = beta.rvs(a=a_ave_VEI7, b=b_ave_VEI7, size=1000000)
            else:
                a_ave_VEI7 = 0
                b_ave_VEI7 = 0
                data_beta_VEI7 = 0

            VEI_2_mean = Ave_VEI2
            VEI_3_mean = Ave_VEI3
            VEI_4_mean = Ave_VEI4
            VEI_5_mean = Ave_VEI5
            VEI_6_mean = Ave_VEI6
            VEI_7_mean = Ave_VEI7

            VEI_2_sd = Std_VEI2
            VEI_3_sd = Std_VEI3
            VEI_4_sd = Std_VEI4
            VEI_5_sd = Std_VEI5
            VEI_6_sd = Std_VEI6
            VEI_7_sd = Std_VEI7

            # Combining probability of eruption and relative probability of different VEI

            if VEI_2_mean > 0:
                AP_VEI_2 = []
                Number = 10000
                lower, upper = 0, 1
                for i in range(Number):
                    """Monte Carlo simulation"""
                    AP_VEI2 = (np.random.beta(a_ave, b_ave) * stats.truncnorm.rvs((lower - VEI_2_mean) / VEI_2_sd,
                                                                                  (upper - VEI_2_mean) / VEI_2_sd,
                                                                                  loc=VEI_2_mean,
                                                                                  scale=VEI_2_sd))
                    AP_VEI_2 = AP_VEI_2 + [AP_VEI2]
            else:
                AP_VEI_2 = 0

            if VEI_3_mean > 0:
                AP_VEI_3 = []
                Number = 10000
                lower, upper = 0, 1
                for i in range(Number):
                    """Monte Carlo simulation"""
                    AP_VEI3 = (np.random.beta(a_ave, b_ave) * stats.truncnorm.rvs((lower - VEI_3_mean) / VEI_3_sd,
                                                                                  (upper - VEI_3_mean) / VEI_3_sd,
                                                                                  loc=VEI_3_mean,
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
                                                                                  (upper - VEI_4_mean) / VEI_4_sd,
                                                                                  loc=VEI_4_mean,
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
                                                                                  (upper - VEI_5_mean) / VEI_5_sd,
                                                                                  loc=VEI_5_mean,
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
                                                                                  (upper - VEI_6_mean) / VEI_6_sd,
                                                                                  loc=VEI_6_mean,
                                                                                  scale=VEI_6_sd))
                    AP_VEI_6 = AP_VEI_6 + [AP_VEI6]
            else:
                AP_VEI_6 = 0

            if VEI_7_mean > 0:
                AP_VEI_7 = []
                Number = 10000
                lower, upper = 0, 1
                for i in range(Number):
                    """Monte Carlo simulation"""
                    AP_VEI7 = (np.random.beta(a_ave, b_ave) * stats.truncnorm.rvs((lower - VEI_7_mean) / VEI_7_sd,
                                                                                  (upper - VEI_7_mean) / VEI_7_sd,
                                                                                  loc=VEI_7_mean,
                                                                                  scale=VEI_7_sd))
                    AP_VEI_7 = AP_VEI_7 + [AP_VEI7]
            else:
                AP_VEI_7 = 0

            #####
            # Visualisation of percentiles and csv files.
            #####

            if Uncertainty_schema == 1:
                Percentile_1 = 5
                p1_string = "5th percentile"
                Percentile_2 = 50
                p2_string = "50th percentile"
                Percentile_3 = 95
                p3_string = "95th percentile"
            elif Uncertainty_schema == 2:
                Percentile_1 = 10
                p1_string = "10th percentile"
                Percentile_2 = 50
                p2_string = "50th percentile"
                Percentile_3 = 90
                p3_string = "590th percentile"
            elif Uncertainty_schema == 3:
                Percentile_1 = 25
                p1_string = "25th percentile"
                Percentile_2 = 50
                p2_string = "50th percentile"
                Percentile_3 = 75
                p3_string = "75th percentile"

            VEI_2_Median = stats.scoreatpercentile(AP_VEI_2, Percentile_2)
            VEI_2_Percentile_1 = stats.scoreatpercentile(AP_VEI_2, Percentile_1)
            VEI_2_Percentile_3 = stats.scoreatpercentile(AP_VEI_2, Percentile_3)

            VEI_3_Median = stats.scoreatpercentile(AP_VEI_3, Percentile_2)
            VEI_3_Percentile_1 = stats.scoreatpercentile(AP_VEI_3, Percentile_1)
            VEI_3_Percentile_3 = stats.scoreatpercentile(AP_VEI_3, Percentile_3)

            VEI_4_Median = stats.scoreatpercentile(AP_VEI_4, Percentile_2)
            VEI_4_Percentile_1 = stats.scoreatpercentile(AP_VEI_4, Percentile_1)
            VEI_4_Percentile_3 = stats.scoreatpercentile(AP_VEI_4, Percentile_3)

            VEI_5_Median = stats.scoreatpercentile(AP_VEI_5, Percentile_2)
            VEI_5_Percentile_1 = stats.scoreatpercentile(AP_VEI_5, Percentile_1)
            VEI_5_Percentile_3 = stats.scoreatpercentile(AP_VEI_5, Percentile_3)

            VEI_6_Median = stats.scoreatpercentile(AP_VEI_6, Percentile_2)
            VEI_6_Percentile_1 = stats.scoreatpercentile(AP_VEI_6, Percentile_1)
            VEI_6_Percentile_3 = stats.scoreatpercentile(AP_VEI_6, Percentile_3)

            VEI_7_Median = stats.scoreatpercentile(AP_VEI_7, Percentile_2)
            VEI_7_Percentile_1 = stats.scoreatpercentile(AP_VEI_7, Percentile_1)
            VEI_7_Percentile_3 = stats.scoreatpercentile(AP_VEI_7, Percentile_3)

            Prob_erupt_Median = VEI_2_Median + VEI_3_Median + VEI_4_Median + VEI_5_Median + VEI_6_Median + VEI_7_Median
            Prob_erupt_Percentile_1 = VEI_2_Percentile_1 + VEI_3_Percentile_1 + VEI_4_Percentile_1 + VEI_5_Percentile_1 + VEI_6_Percentile_1 + VEI_7_Percentile_1
            Prob_erupt_Percentile_3 = VEI_2_Percentile_3 + VEI_3_Percentile_3 + VEI_4_Percentile_3 + VEI_5_Percentile_3 + VEI_6_Percentile_3 + VEI_7_Percentile_3

            Median = "50th percentile"
            Percentile_05 = "10th percentile"
            Percentile_95 = "90th percentile"

            Percentile_05s = [[2, 3, 4, 5, 6, 7],
                              [VEI_2_Percentile_1,
                               VEI_3_Percentile_1,
                               VEI_4_Percentile_1,
                               VEI_5_Percentile_1,
                               VEI_6_Percentile_1,
                               VEI_7_Percentile_1]]
            Percentile_95s = [[2, 3, 4, 5, 6, 7],
                              [VEI_2_Percentile_3,
                               VEI_3_Percentile_3,
                               VEI_4_Percentile_3,
                               VEI_5_Percentile_3,
                               VEI_6_Percentile_3,
                               VEI_7_Percentile_3]]

            Probabilities = pd.DataFrame([[2, VEI_2_Median],
                                          [3, VEI_3_Median],
                                          [4, VEI_4_Median],
                                          [5, VEI_5_Median],
                                          [6, VEI_6_Median],
                                          [7, VEI_7_Median],
                                          [2, VEI_2_Percentile_1],
                                          [3, VEI_3_Percentile_1],
                                          [4, VEI_4_Percentile_1],
                                          [5, VEI_5_Percentile_1],
                                          [6, VEI_6_Percentile_1],
                                          [7, VEI_7_Percentile_1],
                                          [2, VEI_2_Percentile_3],
                                          [3, VEI_3_Percentile_3],
                                          [4, VEI_4_Percentile_3],
                                          [5, VEI_5_Percentile_3],
                                          [6, VEI_6_Percentile_3],
                                          [7, VEI_7_Percentile_3]],
                                         columns=["VEI", "Annual probability"])

            Volcano_data = GVP_volcanoes[(GVP_volcanoes['Volcano Number'] == volc_num)]
            Volcano_name = Volcano_data.iloc[0]['Volcano Name']
            Volcano_Number = volc_num
            textstr = Volcano_name

            Probabilities_df = pd.DataFrame([[Volcano_Number, Volcano_name, year, eruption_certainty, Analysis_method,
                                              method, averaging_method, Power_law,
                                              Prob_erupt_Percentile_1, Prob_erupt_Median, Prob_erupt_Percentile_3,
                                              VEI_2_Percentile_1, VEI_2_Median, VEI_2_Percentile_3,
                                              VEI_3_Percentile_1, VEI_3_Median, VEI_3_Percentile_3,
                                              VEI_4_Percentile_1, VEI_4_Median, VEI_4_Percentile_3,
                                              VEI_5_Percentile_1, VEI_5_Median, VEI_5_Percentile_3,
                                              VEI_6_Percentile_1, VEI_6_Median, VEI_6_Percentile_3,
                                              VEI_7_Percentile_1, VEI_7_Median, VEI_7_Percentile_3]],
                                            columns=["Volcano Number",
                                                     "Volcano Name",
                                                     "GVP DB Year",
                                                     "Included eruptions",
                                                     "Estimate method",
                                                     "Change point",
                                                     "Average method",
                                                     "Power law",
                                                     "Probability of eruption "+ p1_string,
                                                     "Probability of eruption " + p2_string,
                                                     "Probability of eruption " + p3_string,
                                                     "<= VEI 2 " + p1_string,
                                                     "<= VEI 2 " + p2_string,
                                                     "<= VEI 2 " + p3_string,
                                                     "VEI 3 " + p1_string,
                                                     "VEI 3 " + p2_string,
                                                     "VEI 3 " + p3_string,
                                                     "VEI 4 " + p1_string,
                                                     "VEI 4 " + p2_string,
                                                     "VEI 4 " + p3_string,
                                                     "VEI 5 " + p1_string,
                                                     "VEI 5 " + p2_string,
                                                     "VEI 5 " + p3_string,
                                                     "VEI 6 " + p1_string,
                                                     "VEI 6 " + p2_string,
                                                     "VEI 6 " + p3_string,
                                                     "VEI 7 " + p1_string,
                                                     "VEI 7 " + p2_string,
                                                     "VEI 7 " + p3_string])
            ID = str(volc_num)
            if confirmed == "True":
                certainty = "Certain"
            else:
                certainty = "Uncertain"
            if Power_law == "True":
                power = "PowerLaw"
            else:
                power = "NoPowerLaw"

            path_csv = "Probabilities/" + project_folder + "/" + Volcano_name + "_" + certainty + "_" + yearstr + "_" + method + "_" + averaging_method + "_" + power + ".csv"
            Probabilities_df.to_csv(path_csv, index=False)

            if Create_graph == "No":
                print("No graph created due to sampling error")
            else:
                ax = sns.lineplot(x="VEI", y="Annual probability", markers=True, data=Probabilities)
                ax.set(yscale="log")
                ax.set_xticks(range(2, 8))  # <--- set the ticks first
                ax.tick_params(labelsize=10)
                ax.set_xticklabels(['<=2','3', '4', '5', '6', '7'])
                ax.set_xlabel("VEI", fontsize=12)
                ax.set_ylabel("Annual probability", fontsize=12)
                props = dict(boxstyle='round', facecolor='white', alpha=0.5)
                # place a text box in upper left in axes coords
                ax.text(0.80, 0.95, textstr, transform=ax.transAxes, fontsize=12,
                        verticalalignment='top', bbox=props)
                path_fig = "Figures/" + project_folder + "/" + Volcano_name + "_" + certainty + "_" + yearstr + "_" + method + "_" + "_" + averaging_method + "_" + power + ".png"
                plt.savefig(path_fig)
                plt.close()
                # plt.show()

        if VEI_schema == 3:
            # VEI <=1, 2, 3, 4, 5, 6, 7
            if VEI_freq_1["7"].mean() > 0:
                print("Using first concentration parameter for Model 1")
                concentration_1 = np.array([VEI_freq_1["1"].mean() * 10000,
                                            VEI_freq_1["2"].mean() * 10000,
                                            VEI_freq_1["3"].mean() * 10000,
                                            VEI_freq_1["4"].mean() * 10000,
                                            VEI_freq_1["5"].mean() * 10000,
                                            VEI_freq_1["6"].mean() * 10000,
                                            VEI_freq_1["7"].mean() * 10000])
                Observed = np.array([Observed_VEI[0],
                                     Observed_VEI[1],
                                     Observed_VEI[2],
                                     Observed_VEI[3],
                                     Observed_VEI[4],
                                     Observed_VEI[5],
                                     Observed_VEI[6]])

            elif VEI_freq_1["6"].mean() > 0 and VEI_freq_1["7"].mean() == 0:
                print("Using second concentration parameter for Model 1")
                concentration_1 = np.array([VEI_freq_1["1"] * 10000,
                                            VEI_freq_1["2"] * 10000,
                                            VEI_freq_1["3"] * 10000,
                                            VEI_freq_1["4"] * 10000,
                                            VEI_freq_1["5"] * 10000,
                                            VEI_freq_1["6"] * 10000])
                Observed = np.array([Observed_VEI[0],
                                     Observed_VEI[1],
                                     Observed_VEI[2],
                                     Observed_VEI[3],
                                     Observed_VEI[4],
                                     Observed_VEI[5]])

            elif VEI_freq_1["5"].mean() > 0 and VEI_freq_1["6"].mean() == 0:
                print("Using third concentration parameter for Model 1")
                concentration_1 = np.array([VEI_freq_1["1"].mean() * 10000,
                                            VEI_freq_1["2"].mean() * 10000,
                                            VEI_freq_1["3"].mean() * 10000,
                                            VEI_freq_1["4"].mean() * 10000,
                                            VEI_freq_1["5"].mean() * 10000])
                Observed = np.array([Observed_VEI[0],
                                     Observed_VEI[1],
                                     Observed_VEI[2],
                                     Observed_VEI[3],
                                     Observed_VEI[4]])

            elif VEI_freq_1["4"].mean() > 0 and VEI_freq_1["5"].mean() == 0:
                print("Using fourth concentration parameter for Model 1")
                concentration_1 = np.array([VEI_freq_1["1"].mean() * 10000,
                                            VEI_freq_1["2"].mean() * 10000,
                                            VEI_freq_1["3"].mean() * 10000,
                                            VEI_freq_1["4"].mean() * 10000])
                Observed = np.array([Observed_VEI[0],
                                     Observed_VEI[1],
                                     Observed_VEI[2],
                                     Observed_VEI[3]])

            else:
                print("Error with concentration parameter or observed parameter for Model 1")
            # Jenkins_c = Jenkins_1
            shape_1 = concentration_1.size

            Observed_sum = np.sum(Observed)

            # Create model 1
            # try:
            Create_graph = "Yes"
            with pm.Model() as VEI_model_1:
                parameters = pm.Dirichlet('parameters', a=concentration_1, shape=shape_1)
                # Observed data is from a Multinomial distribution
                observed_data = pm.Multinomial(
                    'observed_data', n=Observed_sum, p=parameters, shape=shape_1, observed=Observed)

            with VEI_model_1:
                # Sample from the posterior
                VEI_trace_1 = pm.sample(draws=10000, chains=2, cores=cores, tune=5000,
                                        discard_tuned_samples=True, target_accept=0.99)

            # Create model 2

            if VEI_freq_2["7"].mean() > 0:
                print("Using first concentration parameter for Model 2")

                concentration_2 = np.array([VEI_freq_2["1"].mean() * 10000,
                                            VEI_freq_2["2"].mean() * 10000,
                                            VEI_freq_2["3"].mean() * 10000,
                                            VEI_freq_2["4"].mean() * 10000,
                                            VEI_freq_2["5"].mean() * 10000,
                                            VEI_freq_2["6"].mean() * 10000,
                                            VEI_freq_2["7"].mean() * 10000])
                Observed = np.array([Observed_VEI[0],
                                     Observed_VEI[1],
                                     Observed_VEI[2],
                                     Observed_VEI[3],
                                     Observed_VEI[4],
                                     Observed_VEI[5],
                                     Observed_VEI[6]])

            elif VEI_freq_2["6"].mean() > 0 and VEI_freq_2["7"].mean() == 0:
                print("Using second concentration parameter for Model 2")
                concentration_2 = np.array([VEI_freq_2["1"].mean() * 10000,
                                            VEI_freq_2["2"].mean() * 10000,
                                            VEI_freq_2["3"].mean() * 10000,
                                            VEI_freq_2["4"].mean() * 10000,
                                            VEI_freq_2["5"].mean() * 10000,
                                            VEI_freq_2["6"].mean() * 10000])
                Observed = np.array([Observed_VEI[0],
                                     Observed_VEI[1],
                                     Observed_VEI[2],
                                     Observed_VEI[3],
                                     Observed_VEI[4],
                                     Observed_VEI[5]])

            elif VEI_freq_2["5"].mean() > 0 and VEI_freq_2["6"].mean() == 0:
                print("Using third concentration parameter for Model 2")
                concentration_2 = np.array([VEI_freq_2["1"].mean() * 10000,
                                            VEI_freq_2["2"].mean() * 10000,
                                            VEI_freq_2["3"].mean() * 10000,
                                            VEI_freq_2["4"].mean() * 10000,
                                            VEI_freq_2["5"].mean() * 10000])
                Observed = np.array([Observed_VEI[0],
                                     Observed_VEI[1],
                                     Observed_VEI[2],
                                     Observed_VEI[3],
                                     Observed_VEI[4]])

            elif VEI_freq_2["4"].mean() > 0 and VEI_freq_2["4"].mean() == 0:
                print("Using fourth concentration parameter for Model 2")
                concentration_2 = np.array([VEI_freq_2["1"].mean() * 10000,
                                            VEI_freq_2["2"].mean() * 10000,
                                            VEI_freq_2["3"].mean() * 10000,
                                            VEI_freq_2["4"].mean() * 10000])
                Observed = np.array([Observed_VEI[0],
                                     Observed_VEI[1],
                                     Observed_VEI[2],
                                     Observed_VEI[3]])

            else:
                print("Error with concentration parameter or observed parameter for Model 2")

            # Jenkins_c = Jenkins_1
            shape_2 = concentration_2.size

            with pm.Model() as VEI_model_2:
                parameters = pm.Dirichlet('parameters', a=concentration_2, shape=shape_2)
                # Observed data is from a Multinomial distribution
                observed_data = pm.Multinomial(
                    'observed_data', n=Observed_sum, p=parameters, shape=shape_2, observed=Observed)

            with VEI_model_2:
                # Sample from the posterior
                VEI_trace_2 = pm.sample(draws=10000, chains=2, cores=cores, tune=5000,
                                        discard_tuned_samples=True, target_accept=0.99)

            traces_VEI = [VEI_trace_1, VEI_trace_2]
            model_dict_VEI = dict(zip([VEI_model_1, VEI_model_2], traces_VEI))
            comp_VEI = arviz.compare(model_dict_VEI, method=averaging_method)  # BB-pseudo-BMA
            comp_VEI['index'] = [0, 1]
            comp_VEI.set_index("index", inplace=True)

            ppc_ave_VEI = pm.sample_posterior_predictive_w(traces_VEI, 1000, [VEI_model_1, VEI_model_2],
                                                           weights=comp.weight.sort_index(ascending=True),
                                                           progressbar=False)

            try:
                ppc_ave_VEI_df = pd.DataFrame(ppc_ave_VEI)
                ppc_ave_VEI_df2 = pd.DataFrame(ppc_ave_VEI_df['observed_data'].to_list(),
                                               columns=['VEI1', 'VEI2', 'VEI3', 'VEI4', 'VEI5', 'VEI6', 'VEI7'])
                ppc_ave_VEI_df2 = ppc_ave_VEI_df2.fillna(0)
            except ValueError:
                ppc_ave_VEI_df = pd.DataFrame(ppc_ave_VEI['observed_data'])
                ppc_ave_VEI_df2 = ppc_ave_VEI_df
                if VEI_freq_2["7"].mean() > 0:
                    ppc_ave_VEI_df2.columns = ['VEI1', 'VEI2', 'VEI3', 'VEI4', 'VEI5', 'VEI6', 'VEI7']
                elif VEI_freq_2["5"].mean() > 0 and VEI_freq_2["7"].mean() == 0:
                    if len(ppc_ave_VEI_df2.columns) < 5:
                        print("Number of columns is:", len(ppc_ave_VEI_df2.columns))
                        ppc_ave_VEI_df2.columns = ['VEI1', 'VEI2', 'VEI3', 'VEI4', 'VEI5']
                        ppc_ave_VEI_df2['VEI6'] = 0
                        ppc_ave_VEI_df2['VEI7'] = 0
                    else:
                        ppc_ave_VEI_df2.columns = ['VEI1', 'VEI2', 'VEI3', 'VEI4', 'VEI5', 'VEI6', 'VEI7']

            ppc_ave_VEI_df2['total'] = ppc_ave_VEI_df2['VEI1'] + ppc_ave_VEI_df2['VEI2'] + ppc_ave_VEI_df2['VEI3'] + \
                                       ppc_ave_VEI_df2['VEI4'] + ppc_ave_VEI_df2['VEI5'] + ppc_ave_VEI_df2['VEI6'] + \
                                       ppc_ave_VEI_df2['VEI7']
            ppc_ave_VEI_df2['pVEI1'] = ppc_ave_VEI_df2['VEI1'] / ppc_ave_VEI_df2['total']
            ppc_ave_VEI_df2['pVEI2'] = ppc_ave_VEI_df2['VEI2'] / ppc_ave_VEI_df2['total']
            ppc_ave_VEI_df2['pVEI3'] = ppc_ave_VEI_df2['VEI3'] / ppc_ave_VEI_df2['total']
            ppc_ave_VEI_df2['pVEI4'] = ppc_ave_VEI_df2['VEI4'] / ppc_ave_VEI_df2['total']
            ppc_ave_VEI_df2['pVEI5'] = ppc_ave_VEI_df2['VEI5'] / ppc_ave_VEI_df2['total']
            ppc_ave_VEI_df2['pVEI6'] = ppc_ave_VEI_df2['VEI6'] / ppc_ave_VEI_df2['total']
            ppc_ave_VEI_df2['pVEI7'] = ppc_ave_VEI_df2['VEI7'] / ppc_ave_VEI_df2['total']

            Ave_VEI1 = ppc_ave_VEI_df2['pVEI1'].mean()
            Std_VEI1 = ppc_ave_VEI_df2['pVEI1'].std()

            Ave_VEI2 = ppc_ave_VEI_df2['pVEI2'].mean()
            Std_VEI2 = ppc_ave_VEI_df2['pVEI2'].std()

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

            if Ave_VEI1 > 0:
                a_ave_VEI1 = Ave_VEI1 ** 2 * ((1 - Ave_VEI1) / Std_VEI1 ** 2 - 1 / Ave_VEI1)
                b_ave_VEI1 = a_ave_VEI1 * (1 / Ave_VEI1 - 1)
                data_beta_VEI1 = beta.rvs(a=a_ave_VEI1, b=b_ave_VEI1, size=1000000)
            else:
                a_ave_VEI1 = 0
                b_ave_VEI1 = 0
                data_beta_VEI1 = 0

            if Ave_VEI2 > 0:
                a_ave_VEI2 = Ave_VEI2 ** 2 * ((1 - Ave_VEI2) / Std_VEI2 ** 2 - 1 / Ave_VEI2)
                b_ave_VEI2 = a_ave_VEI2 * (1 / Ave_VEI2 - 1)
                data_beta_VEI2 = beta.rvs(a=a_ave_VEI2, b=b_ave_VEI2, size=1000000)
            else:
                a_ave_VEI2 = 0
                b_ave_VEI2 = 0
                data_beta_VEI2 = 0

            if Ave_VEI3 > 0:
                a_ave_VEI3 = Ave_VEI3 ** 2 * ((1 - Ave_VEI3) / Std_VEI3 ** 2 - 1 / Ave_VEI3)
                b_ave_VEI3 = a_ave_VEI3 * (1 / Ave_VEI3 - 1)
                data_beta_VEI3 = beta.rvs(a=a_ave_VEI3, b=b_ave_VEI3, size=1000000)
            else:
                a_ave_VEI3 = 0
                b_ave_VEI3 = 0
                data_beta_VEI3 = 0

            if Ave_VEI4 > 0:
                a_ave_VEI4 = Ave_VEI4 ** 2 * ((1 - Ave_VEI4) / Std_VEI4 ** 2 - 1 / Ave_VEI4)
                b_ave_VEI4 = a_ave_VEI4 * (1 / Ave_VEI4 - 1)
                data_beta_VEI4 = beta.rvs(a=a_ave_VEI4, b=b_ave_VEI4, size=1000000)
            else:
                a_ave_VEI4 = 0
                b_ave_VEI4 = 0
                data_beta_VEI4 = 0

            if Ave_VEI5 > 0:
                a_ave_VEI5 = Ave_VEI5 ** 2 * ((1 - Ave_VEI5) / Std_VEI5 ** 2 - 1 / Ave_VEI5)
                b_ave_VEI5 = a_ave_VEI5 * (1 / Ave_VEI5 - 1)
                data_beta_VEI5 = beta.rvs(a=a_ave_VEI5, b=b_ave_VEI5, size=1000000)
            else:
                a_ave_VEI5 = 0
                b_ave_VEI5 = 0
                data_beta_VEI5 = 0

            if Ave_VEI6 > 0:
                a_ave_VEI6 = Ave_VEI6 ** 2 * ((1 - Ave_VEI6) / Std_VEI6 ** 2 - 1 / Ave_VEI6)
                b_ave_VEI6 = a_ave_VEI6 * (1 / Ave_VEI6 - 1)
                data_beta_VEI6 = beta.rvs(a=a_ave_VEI6, b=b_ave_VEI6, size=1000000)
            else:
                a_ave_VEI6 = 0
                b_ave_VEI6 = 0
                data_beta_VEI6 = 0

            if Ave_VEI7 > 0:
                a_ave_VEI7 = Ave_VEI7 ** 2 * ((1 - Ave_VEI7) / Std_VEI7 ** 2 - 1 / Ave_VEI7)
                b_ave_VEI7 = a_ave_VEI7 * (1 / Ave_VEI7 - 1)
                data_beta_VEI7 = beta.rvs(a=a_ave_VEI7, b=b_ave_VEI7, size=1000000)
            else:
                a_ave_VEI7 = 0
                b_ave_VEI7 = 0
                data_beta_VEI7 = 0

            VEI_1_mean = Ave_VEI1
            VEI_2_mean = Ave_VEI2
            VEI_3_mean = Ave_VEI3
            VEI_4_mean = Ave_VEI4
            VEI_5_mean = Ave_VEI5
            VEI_6_mean = Ave_VEI6
            VEI_7_mean = Ave_VEI7

            VEI_1_sd = Std_VEI1
            VEI_2_sd = Std_VEI2
            VEI_3_sd = Std_VEI3
            VEI_4_sd = Std_VEI4
            VEI_5_sd = Std_VEI5
            VEI_6_sd = Std_VEI6
            VEI_7_sd = Std_VEI7

            # Combining probability of eruption and relative probability of different VEI
            if VEI_1_mean > 0:
                AP_VEI_1 = []
                Number = 10000
                lower, upper = 0, 1
                for i in range(Number):
                    """Monte Carlo simulation"""
                    AP_VEI1 = (np.random.beta(a_ave, b_ave) * stats.truncnorm.rvs((lower - VEI_1_mean) / VEI_1_sd,
                                                                                  (upper - VEI_1_mean) / VEI_1_sd,
                                                                                  loc=VEI_1_mean,
                                                                                  scale=VEI_1_sd))
                    AP_VEI_1 = AP_VEI_1 + [AP_VEI1]
            else:
                AP_VEI_1 = 0

            if VEI_2_mean > 0:
                AP_VEI_2 = []
                Number = 10000
                lower, upper = 0, 1
                for i in range(Number):
                    """Monte Carlo simulation"""
                    AP_VEI2 = (np.random.beta(a_ave, b_ave) * stats.truncnorm.rvs((lower - VEI_2_mean) / VEI_2_sd,
                                                                                  (upper - VEI_2_mean) / VEI_2_sd,
                                                                                  loc=VEI_2_mean,
                                                                                  scale=VEI_2_sd))
                    AP_VEI_2 = AP_VEI_2 + [AP_VEI2]
            else:
                AP_VEI_2 = 0

            if VEI_3_mean > 0:
                AP_VEI_3 = []
                Number = 10000
                lower, upper = 0, 1
                for i in range(Number):
                    """Monte Carlo simulation"""
                    AP_VEI3 = (np.random.beta(a_ave, b_ave) * stats.truncnorm.rvs((lower - VEI_3_mean) / VEI_3_sd,
                                                                                  (upper - VEI_3_mean) / VEI_3_sd,
                                                                                  loc=VEI_3_mean,
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
                                                                                  (upper - VEI_4_mean) / VEI_4_sd,
                                                                                  loc=VEI_4_mean,
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
                                                                                  (upper - VEI_5_mean) / VEI_5_sd,
                                                                                  loc=VEI_5_mean,
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
                                                                                  (upper - VEI_6_mean) / VEI_6_sd,
                                                                                  loc=VEI_6_mean,
                                                                                  scale=VEI_6_sd))
                    AP_VEI_6 = AP_VEI_6 + [AP_VEI6]
            else:
                AP_VEI_6 = 0

            if VEI_7_mean > 0:
                AP_VEI_7 = []
                Number = 10000
                lower, upper = 0, 1
                for i in range(Number):
                    """Monte Carlo simulation"""
                    AP_VEI7 = (np.random.beta(a_ave, b_ave) * stats.truncnorm.rvs((lower - VEI_7_mean) / VEI_7_sd,
                                                                                  (upper - VEI_7_mean) / VEI_7_sd,
                                                                                  loc=VEI_7_mean,
                                                                                  scale=VEI_7_sd))
                    AP_VEI_7 = AP_VEI_7 + [AP_VEI7]
            else:
                AP_VEI_7 = 0

            #####
            # Visualisation of percentiles and developing csv files.
            #####

            if Uncertainty_schema == 1:
                Percentile_1 = 5
                p1_string = "5th percentile"
                Percentile_2 = 50
                p2_string = "50th percentile"
                Percentile_3 = 95
                p3_string = "95th percentile"
            elif Uncertainty_schema == 2:
                Percentile_1 = 10
                p1_string = "10th percentile"
                Percentile_2 = 50
                p2_string = "50th percentile"
                Percentile_3 = 90
                p3_string = "590th percentile"
            elif Uncertainty_schema == 3:
                Percentile_1 = 25
                p1_string = "25th percentile"
                Percentile_2 = 50
                p2_string = "50th percentile"
                Percentile_3 = 75
                p3_string = "75th percentile"

            VEI_1_Median = stats.scoreatpercentile(AP_VEI_1, Percentile_2)
            VEI_1_Percentile_1 = stats.scoreatpercentile(AP_VEI_1, Percentile_1)
            VEI_1_Percentile_3 = stats.scoreatpercentile(AP_VEI_1, Percentile_3)

            VEI_2_Median = stats.scoreatpercentile(AP_VEI_2, Percentile_2)
            VEI_2_Percentile_1 = stats.scoreatpercentile(AP_VEI_2, Percentile_1)
            VEI_2_Percentile_3 = stats.scoreatpercentile(AP_VEI_2, Percentile_3)

            VEI_3_Median = stats.scoreatpercentile(AP_VEI_3, Percentile_2)
            VEI_3_Percentile_1 = stats.scoreatpercentile(AP_VEI_3, Percentile_1)
            VEI_3_Percentile_3 = stats.scoreatpercentile(AP_VEI_3, Percentile_3)

            VEI_4_Median = stats.scoreatpercentile(AP_VEI_4, Percentile_2)
            VEI_4_Percentile_1 = stats.scoreatpercentile(AP_VEI_4, Percentile_1)
            VEI_4_Percentile_3 = stats.scoreatpercentile(AP_VEI_4, Percentile_3)

            VEI_5_Median = stats.scoreatpercentile(AP_VEI_5, Percentile_2)
            VEI_5_Percentile_1 = stats.scoreatpercentile(AP_VEI_5, Percentile_1)
            VEI_5_Percentile_3 = stats.scoreatpercentile(AP_VEI_5, Percentile_3)

            VEI_6_Median = stats.scoreatpercentile(AP_VEI_6, Percentile_2)
            VEI_6_Percentile_1 = stats.scoreatpercentile(AP_VEI_6, Percentile_1)
            VEI_6_Percentile_3 = stats.scoreatpercentile(AP_VEI_6, Percentile_3)

            VEI_7_Median = stats.scoreatpercentile(AP_VEI_7, Percentile_2)
            VEI_7_Percentile_1 = stats.scoreatpercentile(AP_VEI_7, Percentile_1)
            VEI_7_Percentile_3 = stats.scoreatpercentile(AP_VEI_7, Percentile_3)

            Prob_erupt_Median = VEI_1_Median + VEI_2_Median + VEI_3_Median + VEI_4_Median + VEI_5_Median + VEI_6_Median + VEI_7_Median
            Prob_erupt_Percentile_1 = VEI_1_Percentile_1 + VEI_2_Percentile_1 + VEI_3_Percentile_1 + VEI_4_Percentile_1 + VEI_5_Percentile_1 + VEI_6_Percentile_1 + VEI_7_Percentile_1
            Prob_erupt_Percentile_3 = VEI_1_Percentile_3 + VEI_2_Percentile_3 + VEI_3_Percentile_3 + VEI_4_Percentile_3 + VEI_5_Percentile_3 + VEI_6_Percentile_3 + VEI_7_Percentile_3


            Percentile_1s = [[1, 2, 3, 4, 5, 6, 7],
                              [VEI_1_Percentile_1,
                               VEI_2_Percentile_1,
                               VEI_3_Percentile_1,
                               VEI_4_Percentile_1,
                               VEI_5_Percentile_1,
                               VEI_6_Percentile_1,
                               VEI_7_Percentile_1]]
            Percentile_3s = [[2, 3, 4, 5, 6, 7],
                              [VEI_1_Percentile_3,
                               VEI_2_Percentile_3,
                               VEI_3_Percentile_3,
                               VEI_4_Percentile_3,
                               VEI_5_Percentile_3,
                               VEI_6_Percentile_3,
                               VEI_7_Percentile_3]]

            Probabilities = pd.DataFrame([[1, VEI_1_Median],
                                          [2, VEI_2_Median],
                                          [3, VEI_3_Median],
                                          [4, VEI_4_Median],
                                          [5, VEI_5_Median],
                                          [6, VEI_6_Median],
                                          [7, VEI_7_Median],
                                          [1, VEI_1_Percentile_1],
                                          [2, VEI_2_Percentile_1],
                                          [3, VEI_3_Percentile_1],
                                          [4, VEI_4_Percentile_1],
                                          [5, VEI_5_Percentile_1],
                                          [6, VEI_6_Percentile_1],
                                          [7, VEI_7_Percentile_1],
                                          [1, VEI_1_Percentile_3],
                                          [2, VEI_2_Percentile_3],
                                          [3, VEI_3_Percentile_3],
                                          [4, VEI_4_Percentile_3],
                                          [5, VEI_5_Percentile_3],
                                          [6, VEI_6_Percentile_3],
                                          [7, VEI_7_Percentile_3]],
                                         columns=["VEI", "Annual probability"])

            Volcano_data = GVP_volcanoes[(GVP_volcanoes['Volcano Number'] == volc_num)]
            Volcano_name = Volcano_data.iloc[0]['Volcano Name']
            Volcano_Number = volc_num
            textstr = Volcano_name

            Probabilities_df = pd.DataFrame([[Volcano_Number, Volcano_name, year, eruption_certainty, Analysis_method,
                                              method, averaging_method, Power_law,
                                              Prob_erupt_Percentile_1, Prob_erupt_Median, Prob_erupt_Percentile_3,
                                              VEI_1_Percentile_1, VEI_1_Median, VEI_1_Percentile_3,
                                              VEI_2_Percentile_1, VEI_2_Median, VEI_2_Percentile_3,
                                              VEI_3_Percentile_1, VEI_3_Median, VEI_3_Percentile_3,
                                              VEI_4_Percentile_1, VEI_4_Median, VEI_4_Percentile_3,
                                              VEI_5_Percentile_1, VEI_5_Median, VEI_5_Percentile_3,
                                              VEI_6_Percentile_1, VEI_6_Median, VEI_6_Percentile_3,
                                              VEI_7_Percentile_1, VEI_7_Median, VEI_7_Percentile_3]],
                                            columns=["Volcano Number",
                                                     "Volcano Name",
                                                     "GVP DB Year",
                                                     "Included eruptions",
                                                     "Estimate method",
                                                     "Change point",
                                                     "Average method",
                                                     "Power law",
                                                     "Probability of eruption " + p1_string,
                                                     "Probability of eruption " + p2_string,
                                                     "Probability of eruption " + p3_string,
                                                     "<= VEI 1 " + p1_string,
                                                     "<= VEI 1 " + p2_string,
                                                     "<= VEI 1 " + p3_string,
                                                     "VEI 2 " + p1_string,
                                                     "VEI 2 " + p2_string,
                                                     "VEI 2 " + p3_string,
                                                     "VEI 3 " + p1_string,
                                                     "VEI 3 " + p2_string,
                                                     "VEI 3 " + p3_string,
                                                     "VEI 4 " + p1_string,
                                                     "VEI 4 " + p2_string,
                                                     "VEI 4 " + p3_string,
                                                     "VEI 5 " + p1_string,
                                                     "VEI 5 " + p2_string,
                                                     "VEI 5 " + p3_string,
                                                     "VEI 6 " + p1_string,
                                                     "VEI 6 " + p2_string,
                                                     "VEI 6 " + p3_string,
                                                     "VEI 7 " + p1_string,
                                                     "VEI 7 " + p2_string,
                                                     "VEI 7 " + p3_string])
            ID = str(volc_num)
            if confirmed == "True":
                certainty = "Certain"
            else:
                certainty = "Uncertain"
            if Power_law == "True":
                power = "PowerLaw"
            else:
                power = "NoPowerLaw"

            path_csv = "Probabilities/" + project_folder + "/" + Volcano_name + "_" + certainty + "_" + yearstr + "_" + method + "_" + averaging_method + "_" + power + ".csv"
            Probabilities_df.to_csv(path_csv, index=False)

            if Create_graph == "No":
                print("No graph created due to sampling error")
            else:
                ax = sns.lineplot(x="VEI", y="Annual probability", markers=True, data=Probabilities)
                ax.set(yscale="log")
                ax.set_xticks(range(1, 8))  # <--- set the ticks first
                ax.tick_params(labelsize=10)
                ax.set_xticklabels(['<=1', '2','3', '4', '5', '6', '7'])
                ax.set_xlabel("VEI", fontsize=12)
                ax.set_ylabel("Annual probability", fontsize=12)
                props = dict(boxstyle='round', facecolor='white', alpha=0.5)
                # place a text box in upper left in axes coords
                ax.text(0.80, 0.95, textstr, transform=ax.transAxes, fontsize=12,
                        verticalalignment='top', bbox=props)
                path_fig = "Figures/" + project_folder + "/" + Volcano_name + "_" + certainty + "_" + yearstr + "_" + method + "_" + "_" + averaging_method + "_" + power + ".png"
                plt.savefig(path_fig)
                plt.close()
                # plt.show()

def get_A1_prior_model (GVP_volcanoes, volc_num, Volcano_type_1, Freq_rate_1, Freq_rate_std_1, VEI_freq_1, Power_law, year, confirmed, method, VEI_schema, project_folder, Uncertainty_schema):
    AP_VEI_3 = None
    AP_VEI_4 = None
    AP_VEI_5 = None
    AP_VEI_6 = None
    AP_VEI_7 = None
    yearstr = str(year)
    averaging_method = "None"
    if confirmed == "True":
        eruption_certainty = "Confirmed"
    else:
        eruption_certainty = "Confirmed + Uncertain"
    print("freq_rate_1:", Freq_rate_1[Volcano_type_1])
    print("freq_rate_1_std: ", Freq_rate_std_1[Volcano_type_1])

    alpha_analogue_1 = Freq_rate_1[Volcano_type_1] ** 2 * (
            (1 - Freq_rate_1[Volcano_type_1]) / Freq_rate_std_1[Volcano_type_1] ** 2 - 1 / Freq_rate_1[Volcano_type_1])
    beta_analogue_1 = alpha_analogue_1 * (1 / Freq_rate_1[Volcano_type_1] - 1)
    Analysis_method = "Analogue A1"
    print("Initiating model using: ", Analysis_method)
    a_ave = alpha_analogue_1
    b_ave = beta_analogue_1

    if Power_law == "True":
        power = "PowerLaw"
    else:
        power = "NoPowerLaw"
    if VEI_schema == 1:
        # Combining probability of eruption and relative probability of different VEI
        AP_VEI_3 = []
        Number = 10000
        lower, upper = 0, 1
        for i in range(Number):
            """Monte Carlo simulation"""
            # AP_VEI3 = (np.random.beta(alpha_analogue_1, beta_analogue_1) * stats.truncnorm.rvs((lower - VEI_freq_1[0]) / VEI_freq_1[5],
            #                                                               (upper - VEI_freq_1[0]) / VEI_freq_1[5], loc=VEI_freq_1[0],
            #                                                               scale=VEI_freq_1[5]))
            # AP_VEI_3 = AP_VEI_3 + [AP_VEI3]
            AP_VEI3 = (np.random.beta(alpha_analogue_1, beta_analogue_1) * VEI_freq_1['3'])
            AP_VEI_3 = AP_VEI_3 + [AP_VEI3]

        AP_VEI_4 = []
        Number = 10000
        lower, upper = 0, 1
        for i in range(Number):
            """Monte Carlo simulation"""
            # AP_VEI4 = (np.random.beta(a_ave, b_ave) * stats.truncnorm.rvs((lower - VEI_freq_1[1]) / VEI_freq_1[6],
            #                                                               (upper - VEI_freq_1[1]) / VEI_freq_1[6], loc=VEI_freq_1[1],
            #                                                               scale=VEI_freq_1[6]))
            # AP_VEI_4 = AP_VEI_4 + [AP_VEI4]

            AP_VEI4 = (np.random.beta(alpha_analogue_1, beta_analogue_1) * VEI_freq_1['4'])
            AP_VEI_4 = AP_VEI_4 + [AP_VEI4]

        # if VEI_freq_1['5'] >0:
        AP_VEI_5 = []
        Number = 10000
        lower, upper = 0, 1
        for i in range(Number):
            """Monte Carlo simulation"""
            # AP_VEI5 = (np.random.beta(a_ave, b_ave) * stats.truncnorm.rvs((lower - VEI_freq_1[2]) / VEI_freq_1[7],
            #                                                               (upper - VEI_freq_1[2]) / VEI_freq_1[7], loc=VEI_freq_1[2],
            #                                                               scale=VEI_freq_1[7]))
            # AP_VEI_5 = AP_VEI_5 + [AP_VEI5]

            AP_VEI5 = (np.random.beta(alpha_analogue_1, beta_analogue_1) * VEI_freq_1['5'])
            AP_VEI_5 = AP_VEI_5 + [AP_VEI5]
        # else: AP_VEI_5 = 0

        #if VEI_freq_1['6'] > 0:
        AP_VEI_6 = []
        Number = 10000
        lower, upper = 0, 1
        for i in range(Number):
            """Monte Carlo simulation"""
            # AP_VEI6 = (np.random.beta(a_ave, b_ave) * stats.truncnorm.rvs((lower - VEI_freq_1[3]) / VEI_freq_1[8],
            #                                                               (upper - VEI_freq_1[3]) / VEI_freq_1[8], loc=VEI_freq_1[3],
            #                                                               scale=VEI_freq_1[8]))
            # AP_VEI_6 = AP_VEI_6 + [AP_VEI6]
            AP_VEI6 = (np.random.beta(alpha_analogue_1, beta_analogue_1) * VEI_freq_1['6'])
            AP_VEI_6 = AP_VEI_6 + [AP_VEI6]
        # else:
        #     AP_VEI_6 = 0

        # if VEI_freq_1['7'] > 0:
        AP_VEI_7 = []
        Number = 10000
        lower, upper = 0, 1
        for i in range(Number):
            """Monte Carlo simulation"""
            # AP_VEI7 = (np.random.beta(a_ave, b_ave) * stats.truncnorm.rvs((lower - VEI_freq_1[4]) / VEI_freq_1[9],
            #                                                               (upper - VEI_freq_1[4]) / VEI_freq_1[9], loc=VEI_freq_1[4],
            #                                                               scale=VEI_freq_1[9]))
            # AP_VEI_7 = AP_VEI_7 + [AP_VEI7]
            AP_VEI7 = (np.random.beta(alpha_analogue_1, beta_analogue_1) * VEI_freq_1['7'])
            AP_VEI_7 = AP_VEI_7 + [AP_VEI7]

        # else:
        #     AP_VEI_7 = 0

        if Uncertainty_schema == 1:
            Percentile_1 = 5
            p1_string = "5th percentile"
            Percentile_2 = 50
            p2_string = "50th percentile"
            Percentile_3 = 95
            p3_string = "95th percentile"
        elif Uncertainty_schema == 2:
            Percentile_1 = 10
            p1_string = "10th percentile"
            Percentile_2 = 50
            p2_string = "50th percentile"
            Percentile_3 = 90
            p3_string = "590th percentile"
        elif Uncertainty_schema == 3:
            Percentile_1 = 25
            p1_string = "25th percentile"
            Percentile_2 = 50
            p2_string = "50th percentile"
            Percentile_3 = 75
            p3_string = "75th percentile"


        VEI_3_Median = stats.scoreatpercentile(AP_VEI_3, Percentile_2)
        VEI_3_Percentile_1 = stats.scoreatpercentile(AP_VEI_3, Percentile_1)
        VEI_3_Percentile_3 = stats.scoreatpercentile(AP_VEI_3, Percentile_3)

        VEI_4_Median = stats.scoreatpercentile(AP_VEI_4, Percentile_2)
        VEI_4_Percentile_1 = stats.scoreatpercentile(AP_VEI_4, Percentile_1)
        VEI_4_Percentile_3 = stats.scoreatpercentile(AP_VEI_4, Percentile_3)

        VEI_5_Median = stats.scoreatpercentile(AP_VEI_5, Percentile_2)
        VEI_5_Percentile_1 = stats.scoreatpercentile(AP_VEI_5, Percentile_1)
        VEI_5_Percentile_3 = stats.scoreatpercentile(AP_VEI_5, Percentile_3)

        VEI_6_Median = stats.scoreatpercentile(AP_VEI_6, Percentile_2)
        VEI_6_Percentile_1 = stats.scoreatpercentile(AP_VEI_6, Percentile_1)
        VEI_6_Percentile_3 = stats.scoreatpercentile(AP_VEI_6, Percentile_3)

        VEI_7_Median = stats.scoreatpercentile(AP_VEI_7, Percentile_2)
        VEI_7_Percentile_1 = stats.scoreatpercentile(AP_VEI_7, Percentile_1)
        VEI_7_Percentile_3 = stats.scoreatpercentile(AP_VEI_7, Percentile_3)

        Prob_erupt_Median = VEI_3_Median + VEI_4_Median + VEI_5_Median + VEI_6_Median + VEI_7_Median
        Prob_erupt_Percentile_1 = VEI_3_Percentile_1 + VEI_4_Percentile_1 + VEI_5_Percentile_1 + VEI_6_Percentile_1 + VEI_7_Percentile_1
        Prob_erupt_Percentile_3 = VEI_3_Percentile_3 + VEI_4_Percentile_3 + VEI_5_Percentile_3 + VEI_6_Percentile_3 + VEI_7_Percentile_3

        Median = "50th percentile"
        Percentile_05 = "5th percentile"
        Percentile_95 = "95th percentile"

        Percentile_1s = [[3, 4, 5, 6, 7],
                          [VEI_3_Percentile_1, VEI_4_Percentile_1, VEI_5_Percentile_1, VEI_6_Percentile_1,
                           VEI_7_Percentile_1]]
        Percentile_3s = [[3, 4, 5, 6, 7],
                          [VEI_3_Percentile_3, VEI_4_Percentile_3, VEI_5_Percentile_3, VEI_6_Percentile_3,
                           VEI_7_Percentile_3]]

        Probabilities = pd.DataFrame([[3, VEI_3_Median],
                                      [4, VEI_4_Median],
                                      [5, VEI_5_Median],
                                      [6, VEI_6_Median],
                                      [7, VEI_7_Median],
                                      [3, VEI_3_Percentile_1],
                                      [4, VEI_4_Percentile_1],
                                      [5, VEI_5_Percentile_1],
                                      [6, VEI_6_Percentile_1],
                                      [7, VEI_7_Percentile_1],
                                      [3, VEI_3_Percentile_3],
                                      [4, VEI_4_Percentile_3],
                                      [5, VEI_5_Percentile_3],
                                      [6, VEI_6_Percentile_3],
                                      [7, VEI_7_Percentile_3]],
                                     columns=["VEI", "Annual probability"])

        Volcano_data = GVP_volcanoes[(GVP_volcanoes['Volcano Number'] == volc_num)]
        Volcano_name = Volcano_data.iloc[0]['Volcano Name']
        Volcano_Number = volc_num
        textstr = Volcano_name

        Probabilities_df = pd.DataFrame([[Volcano_Number, Volcano_name, year, eruption_certainty, Analysis_method, method, averaging_method, Power_law,
                                          Prob_erupt_Percentile_1, Prob_erupt_Median, Prob_erupt_Percentile_3,
                                          VEI_3_Percentile_1, VEI_3_Median, VEI_3_Percentile_3,
                                          VEI_4_Percentile_1, VEI_4_Median, VEI_4_Percentile_3,
                                          VEI_5_Percentile_1, VEI_5_Median, VEI_5_Percentile_3,
                                          VEI_6_Percentile_1, VEI_6_Median, VEI_6_Percentile_3,
                                          VEI_7_Percentile_1, VEI_7_Median, VEI_7_Percentile_3]],
                                        columns=["Volcano Number",
                                                 "Volcano Name",
                                                 "GVP DB Year",
                                                 "Included eruptions",
                                                 "Estimate method",
                                                 "Change point",
                                                 "Average method",
                                                 "Power law",
                                                 "Probability of eruption " + p1_string,
                                                 "Probability of eruption " + p2_string,
                                                 "Probability of eruption " + p3_string,
                                                 "<= VEI 3 " + p1_string,
                                                 "<= VEI 3 " + p2_string,
                                                 "<= VEI 3 " + p3_string,
                                                 "VEI 4 " + p1_string,
                                                 "VEI 4 " + p2_string,
                                                 "VEI 4 " + p3_string,
                                                 "VEI 5 " + p1_string,
                                                 "VEI 5 " + p2_string,
                                                 "VEI 5 " + p3_string,
                                                 "VEI 6 " + p1_string,
                                                 "VEI 6 " + p2_string,
                                                 "VEI 6 " + p3_string,
                                                 "VEI 7 " + p1_string,
                                                 "VEI 7 " + p2_string,
                                                 "VEI 7 " + p3_string])
        ID = str(volc_num)
        if confirmed == "True":
            certainty = "Certain"
        else:
            certainty = "Uncertain"
        path_csv = "Probabilities/" + project_folder + "/" + Volcano_name + "_" + certainty + "_" + method + "_" + yearstr + "_" + power + "_A1.csv"
        Probabilities_df.to_csv(path_csv, index=False)

        ax = sns.lineplot(x="VEI", y="Annual probability", markers=True, data=Probabilities)
        ax.set(yscale="log")
        ax.set_xticks(range(3, 8))  # <--- set the ticks first
        ax.tick_params(labelsize=10)
        ax.set_xticklabels(['<=3', '4', '5', '6', '7'])
        ax.set_xlabel("VEI", fontsize=12)
        ax.set_ylabel("Annual probability", fontsize=12)
        props = dict(boxstyle='round', facecolor='white', alpha=0.5)
        # place a text box in upper left in axes coords
        ax.text(0.80, 0.95, textstr, transform=ax.transAxes, fontsize=12,
                verticalalignment='top', bbox=props)
        path_fig = "Figures/" + project_folder + "/" + Volcano_name + "_" + certainty + "_" + method + "_" + yearstr + "_" + power + "_A1.png"
        plt.savefig(path_fig)
        plt.close()
        #plt.show()

    if VEI_schema == 2:
        # VEI <=2, 3, 4, 5, 6, 7
        # Combining probability of eruption and relative probability of different VEI

        AP_VEI_2 = []
        Number = 10000
        lower, upper = 0, 1
        for i in range(Number):
            """Monte Carlo simulation"""
            # AP_VEI3 = (np.random.beta(alpha_analogue_1, beta_analogue_1) * stats.truncnorm.rvs((lower - VEI_freq_1[0]) / VEI_freq_1[5],
            #                                                               (upper - VEI_freq_1[0]) / VEI_freq_1[5], loc=VEI_freq_1[0],
            #                                                               scale=VEI_freq_1[5]))
            # AP_VEI_3 = AP_VEI_3 + [AP_VEI3]
            AP_VEI2 = (np.random.beta(alpha_analogue_1, beta_analogue_1) * VEI_freq_1['2'])
            AP_VEI_2 = AP_VEI_2 + [AP_VEI2]

        AP_VEI_3 = []
        Number = 10000
        lower, upper = 0, 1
        for i in range(Number):
            """Monte Carlo simulation"""
            # AP_VEI3 = (np.random.beta(alpha_analogue_1, beta_analogue_1) * stats.truncnorm.rvs((lower - VEI_freq_1[0]) / VEI_freq_1[5],
            #                                                               (upper - VEI_freq_1[0]) / VEI_freq_1[5], loc=VEI_freq_1[0],
            #                                                               scale=VEI_freq_1[5]))
            # AP_VEI_3 = AP_VEI_3 + [AP_VEI3]
            AP_VEI3 = (np.random.beta(alpha_analogue_1, beta_analogue_1) * VEI_freq_1['3'])
            AP_VEI_3 = AP_VEI_3 + [AP_VEI3]

        AP_VEI_4 = []
        Number = 10000
        lower, upper = 0, 1
        for i in range(Number):
            """Monte Carlo simulation"""
            # AP_VEI4 = (np.random.beta(a_ave, b_ave) * stats.truncnorm.rvs((lower - VEI_freq_1[1]) / VEI_freq_1[6],
            #                                                               (upper - VEI_freq_1[1]) / VEI_freq_1[6], loc=VEI_freq_1[1],
            #                                                               scale=VEI_freq_1[6]))
            # AP_VEI_4 = AP_VEI_4 + [AP_VEI4]

            AP_VEI4 = (np.random.beta(alpha_analogue_1, beta_analogue_1) * VEI_freq_1['4'])
            AP_VEI_4 = AP_VEI_4 + [AP_VEI4]

        # if VEI_freq_1['5'] >0:
        AP_VEI_5 = []
        Number = 10000
        lower, upper = 0, 1
        for i in range(Number):
            """Monte Carlo simulation"""
            # AP_VEI5 = (np.random.beta(a_ave, b_ave) * stats.truncnorm.rvs((lower - VEI_freq_1[2]) / VEI_freq_1[7],
            #                                                               (upper - VEI_freq_1[2]) / VEI_freq_1[7], loc=VEI_freq_1[2],
            #                                                               scale=VEI_freq_1[7]))
            # AP_VEI_5 = AP_VEI_5 + [AP_VEI5]

            AP_VEI5 = (np.random.beta(alpha_analogue_1, beta_analogue_1) * VEI_freq_1['5'])
            AP_VEI_5 = AP_VEI_5 + [AP_VEI5]
        # else: AP_VEI_5 = 0

        # if VEI_freq_1['6'] > 0:
        AP_VEI_6 = []
        Number = 10000
        lower, upper = 0, 1
        for i in range(Number):
            """Monte Carlo simulation"""
            # AP_VEI6 = (np.random.beta(a_ave, b_ave) * stats.truncnorm.rvs((lower - VEI_freq_1[3]) / VEI_freq_1[8],
            #                                                               (upper - VEI_freq_1[3]) / VEI_freq_1[8], loc=VEI_freq_1[3],
            #                                                               scale=VEI_freq_1[8]))
            # AP_VEI_6 = AP_VEI_6 + [AP_VEI6]
            AP_VEI6 = (np.random.beta(alpha_analogue_1, beta_analogue_1) * VEI_freq_1['6'])
            AP_VEI_6 = AP_VEI_6 + [AP_VEI6]
        # else:
        #     AP_VEI_6 = 0

        # if VEI_freq_1['7'] > 0:
        AP_VEI_7 = []
        Number = 10000
        lower, upper = 0, 1
        for i in range(Number):
            """Monte Carlo simulation"""
            # AP_VEI7 = (np.random.beta(a_ave, b_ave) * stats.truncnorm.rvs((lower - VEI_freq_1[4]) / VEI_freq_1[9],
            #                                                               (upper - VEI_freq_1[4]) / VEI_freq_1[9], loc=VEI_freq_1[4],
            #                                                               scale=VEI_freq_1[9]))
            # AP_VEI_7 = AP_VEI_7 + [AP_VEI7]
            AP_VEI7 = (np.random.beta(alpha_analogue_1, beta_analogue_1) * VEI_freq_1['7'])
            AP_VEI_7 = AP_VEI_7 + [AP_VEI7]

        # else:
        #     AP_VEI_7 = 0

        if Uncertainty_schema == 1:
            Percentile_1 = 5
            p1_string = "5th percentile"
            Percentile_2 = 50
            p2_string = "50th percentile"
            Percentile_3 = 95
            p3_string = "95th percentile"
        elif Uncertainty_schema == 2:
            Percentile_1 = 10
            p1_string = "10th percentile"
            Percentile_2 = 50
            p2_string = "50th percentile"
            Percentile_3 = 90
            p3_string = "90th percentile"
        elif Uncertainty_schema == 3:
            Percentile_1 = 25
            p1_string = "25th percentile"
            Percentile_2 = 50
            p2_string = "50th percentile"
            Percentile_3 = 75
            p3_string = "75th percentile"

        VEI_2_Median = stats.scoreatpercentile(AP_VEI_2, Percentile_2)
        VEI_2_Percentile_1 = stats.scoreatpercentile(AP_VEI_2, Percentile_1)
        VEI_2_Percentile_3 = stats.scoreatpercentile(AP_VEI_2, Percentile_3)

        VEI_3_Median = stats.scoreatpercentile(AP_VEI_3, Percentile_2)
        VEI_3_Percentile_1 = stats.scoreatpercentile(AP_VEI_3, Percentile_1)
        VEI_3_Percentile_3 = stats.scoreatpercentile(AP_VEI_3, Percentile_3)

        VEI_4_Median = stats.scoreatpercentile(AP_VEI_4, Percentile_2)
        VEI_4_Percentile_1 = stats.scoreatpercentile(AP_VEI_4, Percentile_1)
        VEI_4_Percentile_3 = stats.scoreatpercentile(AP_VEI_4, Percentile_3)

        VEI_5_Median = stats.scoreatpercentile(AP_VEI_5, Percentile_2)
        VEI_5_Percentile_1 = stats.scoreatpercentile(AP_VEI_5, Percentile_1)
        VEI_5_Percentile_3 = stats.scoreatpercentile(AP_VEI_5, Percentile_3)

        VEI_6_Median = stats.scoreatpercentile(AP_VEI_6, Percentile_2)
        VEI_6_Percentile_1 = stats.scoreatpercentile(AP_VEI_6, Percentile_1)
        VEI_6_Percentile_3 = stats.scoreatpercentile(AP_VEI_6, Percentile_3)

        VEI_7_Median = stats.scoreatpercentile(AP_VEI_7, Percentile_2)
        VEI_7_Percentile_1 = stats.scoreatpercentile(AP_VEI_7, Percentile_1)
        VEI_7_Percentile_3 = stats.scoreatpercentile(AP_VEI_7, Percentile_3)

        Prob_erupt_Median = VEI_2_Median + VEI_3_Median + VEI_4_Median + VEI_5_Median + VEI_6_Median + VEI_7_Median
        Prob_erupt_Percentile_1 = VEI_2_Percentile_1 + VEI_3_Percentile_1 + VEI_4_Percentile_1 + VEI_5_Percentile_1 + VEI_6_Percentile_1 + VEI_7_Percentile_1
        Prob_erupt_Percentile_3 = VEI_2_Percentile_3 + VEI_3_Percentile_3 + VEI_4_Percentile_3 + VEI_5_Percentile_3 + VEI_6_Percentile_3 + VEI_7_Percentile_3

        Median = "50th percentile"
        Percentile_05 = "5th percentile"
        Percentile_95 = "95th percentile"

        Percentile_1s = [[2, 3, 4, 5, 6, 7],
                          [VEI_2_Percentile_1, VEI_3_Percentile_1, VEI_4_Percentile_1, VEI_5_Percentile_1, VEI_6_Percentile_1,
                           VEI_7_Percentile_1]]
        Percentile_3s = [[2, 3, 4, 5, 6, 7],
                          [VEI_2_Percentile_1, VEI_3_Percentile_1, VEI_4_Percentile_1, VEI_5_Percentile_1, VEI_6_Percentile_1,
                           VEI_7_Percentile_1]]

        Probabilities = pd.DataFrame([[2, VEI_2_Median],
                                      [3, VEI_3_Median],
                                      [4, VEI_4_Median],
                                      [5, VEI_5_Median],
                                      [6, VEI_6_Median],
                                      [7, VEI_7_Median],
                                      [2, VEI_2_Percentile_1],
                                      [3, VEI_3_Percentile_1],
                                      [4, VEI_4_Percentile_1],
                                      [5, VEI_5_Percentile_1],
                                      [6, VEI_6_Percentile_1],
                                      [7, VEI_7_Percentile_1],
                                      [2, VEI_2_Percentile_3],
                                      [3, VEI_3_Percentile_3],
                                      [4, VEI_4_Percentile_3],
                                      [5, VEI_5_Percentile_3],
                                      [6, VEI_6_Percentile_3],
                                      [7, VEI_7_Percentile_3]],
                                     columns=["VEI", "Annual probability"])

        Volcano_data = GVP_volcanoes[(GVP_volcanoes['Volcano Number'] == volc_num)]
        Volcano_name = Volcano_data.iloc[0]['Volcano Name']
        Volcano_Number = volc_num
        textstr = Volcano_name

        Probabilities_df = pd.DataFrame([[Volcano_Number, Volcano_name, year, eruption_certainty, Analysis_method,
                                          method, averaging_method, Power_law,
                                          Prob_erupt_Percentile_1, Prob_erupt_Median, Prob_erupt_Percentile_3,
                                          VEI_2_Percentile_1, VEI_2_Median, VEI_2_Percentile_3,
                                          VEI_3_Percentile_1, VEI_3_Median, VEI_3_Percentile_3,
                                          VEI_4_Percentile_1, VEI_4_Median, VEI_4_Percentile_3,
                                          VEI_5_Percentile_1, VEI_5_Median, VEI_5_Percentile_3,
                                          VEI_6_Percentile_1, VEI_6_Median, VEI_6_Percentile_3,
                                          VEI_7_Percentile_1, VEI_7_Median, VEI_7_Percentile_3]],
                                        columns=["Volcano Number",
                                                 "Volcano Name",
                                                 "GVP DB Year",
                                                 "Included eruptions",
                                                 "Estimate method",
                                                 "Change point",
                                                 "Average method",
                                                 "Power law",
                                                 "Probability of eruption " + p1_string,
                                                 "Probability of eruption " + p2_string,
                                                 "Probability of eruption " + p3_string,
                                                 "<= VEI 2 " + p1_string,
                                                 "<= VEI 2 " + p2_string,
                                                 "<= VEI 2 " + p3_string,
                                                 "VEI 3 " + p1_string,
                                                 "VEI 3 " + p2_string,
                                                 "VEI 3 " + p3_string,
                                                 "VEI 4 " + p1_string,
                                                 "VEI 4 " + p2_string,
                                                 "VEI 4 " + p3_string,
                                                 "VEI 5 " + p1_string,
                                                 "VEI 5 " + p2_string,
                                                 "VEI 5 " + p3_string,
                                                 "VEI 6 " + p1_string,
                                                 "VEI 6 " + p2_string,
                                                 "VEI 6 " + p3_string,
                                                 "VEI 7 " + p1_string,
                                                 "VEI 7 " + p2_string,
                                                 "VEI 7 " + p3_string])
        ID = str(volc_num)
        if confirmed == "True":
            certainty = "Certain"
        else:
            certainty = "Uncertain"
        path_csv = "Probabilities/" + project_folder + "/" + Volcano_name + "_" + certainty + "_" + method + "_" + yearstr + "_" + power + "_A1.csv"
        Probabilities_df.to_csv(path_csv, index=False)

        ax = sns.lineplot(x="VEI", y="Annual probability", markers=True, data=Probabilities)
        ax.set(yscale="log")
        ax.set_xticks(range(2, 8))  # <--- set the ticks first
        ax.tick_params(labelsize=10)
        ax.set_xticklabels(['<=2', '3', '4', '5', '6', '7'])
        ax.set_xlabel("VEI", fontsize=12)
        ax.set_ylabel("Annual probability", fontsize=12)
        props = dict(boxstyle='round', facecolor='white', alpha=0.5)
        # place a text box in upper left in axes coords
        ax.text(0.80, 0.95, textstr, transform=ax.transAxes, fontsize=12,
                verticalalignment='top', bbox=props)
        path_fig = "Figures/" + project_folder + "/" + Volcano_name + "_" + certainty + "_" + method + "_" + yearstr + "_" + power + "_A1.png"
        plt.savefig(path_fig)
        plt.close()
        # plt.show()

    if VEI_schema == 3:
        # VEI <=1, 2, 3, 4, 5, 6, 7
        # Combining probability of eruption and relative probability of different VEI

        AP_VEI_1 = []
        Number = 10000
        lower, upper = 0, 1
        for i in range(Number):
            """Monte Carlo simulation"""
            # AP_VEI3 = (np.random.beta(alpha_analogue_1, beta_analogue_1) * stats.truncnorm.rvs((lower - VEI_freq_1[0]) / VEI_freq_1[5],
            #                                                               (upper - VEI_freq_1[0]) / VEI_freq_1[5], loc=VEI_freq_1[0],
            #                                                               scale=VEI_freq_1[5]))
            # AP_VEI_3 = AP_VEI_3 + [AP_VEI3]
            AP_VEI1 = (np.random.beta(alpha_analogue_1, beta_analogue_1) * VEI_freq_1['1'])
            AP_VEI_1 = AP_VEI_1 + [AP_VEI1]

        AP_VEI_2 = []
        Number = 10000
        lower, upper = 0, 1
        for i in range(Number):
            """Monte Carlo simulation"""
            # AP_VEI3 = (np.random.beta(alpha_analogue_1, beta_analogue_1) * stats.truncnorm.rvs((lower - VEI_freq_1[0]) / VEI_freq_1[5],
            #                                                               (upper - VEI_freq_1[0]) / VEI_freq_1[5], loc=VEI_freq_1[0],
            #                                                               scale=VEI_freq_1[5]))
            # AP_VEI_3 = AP_VEI_3 + [AP_VEI3]
            AP_VEI2 = (np.random.beta(alpha_analogue_1, beta_analogue_1) * VEI_freq_1['2'])
            AP_VEI_2 = AP_VEI_2 + [AP_VEI2]

        AP_VEI_3 = []
        Number = 10000
        lower, upper = 0, 1
        for i in range(Number):
            """Monte Carlo simulation"""
            # AP_VEI3 = (np.random.beta(alpha_analogue_1, beta_analogue_1) * stats.truncnorm.rvs((lower - VEI_freq_1[0]) / VEI_freq_1[5],
            #                                                               (upper - VEI_freq_1[0]) / VEI_freq_1[5], loc=VEI_freq_1[0],
            #                                                               scale=VEI_freq_1[5]))
            # AP_VEI_3 = AP_VEI_3 + [AP_VEI3]
            AP_VEI3 = (np.random.beta(alpha_analogue_1, beta_analogue_1) * VEI_freq_1['3'])
            AP_VEI_3 = AP_VEI_3 + [AP_VEI3]

        AP_VEI_4 = []
        Number = 10000
        lower, upper = 0, 1
        for i in range(Number):
            """Monte Carlo simulation"""
            # AP_VEI4 = (np.random.beta(a_ave, b_ave) * stats.truncnorm.rvs((lower - VEI_freq_1[1]) / VEI_freq_1[6],
            #                                                               (upper - VEI_freq_1[1]) / VEI_freq_1[6], loc=VEI_freq_1[1],
            #                                                               scale=VEI_freq_1[6]))
            # AP_VEI_4 = AP_VEI_4 + [AP_VEI4]

            AP_VEI4 = (np.random.beta(alpha_analogue_1, beta_analogue_1) * VEI_freq_1['4'])
            AP_VEI_4 = AP_VEI_4 + [AP_VEI4]

        # if VEI_freq_1['5'] >0:
        AP_VEI_5 = []
        Number = 10000
        lower, upper = 0, 1
        for i in range(Number):
            """Monte Carlo simulation"""
            # AP_VEI5 = (np.random.beta(a_ave, b_ave) * stats.truncnorm.rvs((lower - VEI_freq_1[2]) / VEI_freq_1[7],
            #                                                               (upper - VEI_freq_1[2]) / VEI_freq_1[7], loc=VEI_freq_1[2],
            #                                                               scale=VEI_freq_1[7]))
            # AP_VEI_5 = AP_VEI_5 + [AP_VEI5]

            AP_VEI5 = (np.random.beta(alpha_analogue_1, beta_analogue_1) * VEI_freq_1['5'])
            AP_VEI_5 = AP_VEI_5 + [AP_VEI5]
        # else: AP_VEI_5 = 0

        # if VEI_freq_1['6'] > 0:
        AP_VEI_6 = []
        Number = 10000
        lower, upper = 0, 1
        for i in range(Number):
            """Monte Carlo simulation"""
            # AP_VEI6 = (np.random.beta(a_ave, b_ave) * stats.truncnorm.rvs((lower - VEI_freq_1[3]) / VEI_freq_1[8],
            #                                                               (upper - VEI_freq_1[3]) / VEI_freq_1[8], loc=VEI_freq_1[3],
            #                                                               scale=VEI_freq_1[8]))
            # AP_VEI_6 = AP_VEI_6 + [AP_VEI6]
            AP_VEI6 = (np.random.beta(alpha_analogue_1, beta_analogue_1) * VEI_freq_1['6'])
            AP_VEI_6 = AP_VEI_6 + [AP_VEI6]
        # else:
        #     AP_VEI_6 = 0

        # if VEI_freq_1['7'] > 0:
        AP_VEI_7 = []
        Number = 10000
        lower, upper = 0, 1
        for i in range(Number):
            """Monte Carlo simulation"""
            # AP_VEI7 = (np.random.beta(a_ave, b_ave) * stats.truncnorm.rvs((lower - VEI_freq_1[4]) / VEI_freq_1[9],
            #                                                               (upper - VEI_freq_1[4]) / VEI_freq_1[9], loc=VEI_freq_1[4],
            #                                                               scale=VEI_freq_1[9]))
            # AP_VEI_7 = AP_VEI_7 + [AP_VEI7]
            AP_VEI7 = (np.random.beta(alpha_analogue_1, beta_analogue_1) * VEI_freq_1['7'])
            AP_VEI_7 = AP_VEI_7 + [AP_VEI7]

        # else:
        #     AP_VEI_7 = 0

        if Uncertainty_schema == 1:
            Percentile_1 = 5
            p1_string = "5th percentile"
            Percentile_2 = 50
            p2_string = "50th percentile"
            Percentile_3 = 95
            p3_string = "95th percentile"
        elif Uncertainty_schema == 2:
            Percentile_1 = 10
            p1_string = "10th percentile"
            Percentile_2 = 50
            p2_string = "50th percentile"
            Percentile_3 = 90
            p3_string = "590th percentile"
        elif Uncertainty_schema == 3:
            Percentile_1 = 25
            p1_string = "25th percentile"
            Percentile_2 = 50
            p2_string = "50th percentile"
            Percentile_3 = 75
            p3_string = "75th percentile"

        VEI_1_Median = stats.scoreatpercentile(AP_VEI_1, Percentile_2)
        VEI_1_Percentile_1 = stats.scoreatpercentile(AP_VEI_1, Percentile_1)
        VEI_1_Percentile_3 = stats.scoreatpercentile(AP_VEI_1, Percentile_3)

        VEI_2_Median = stats.scoreatpercentile(AP_VEI_2, Percentile_2)
        VEI_2_Percentile_1 = stats.scoreatpercentile(AP_VEI_2, Percentile_1)
        VEI_2_Percentile_3 = stats.scoreatpercentile(AP_VEI_2, Percentile_3)

        VEI_3_Median = stats.scoreatpercentile(AP_VEI_3, Percentile_2)
        VEI_3_Percentile_1 = stats.scoreatpercentile(AP_VEI_3, Percentile_1)
        VEI_3_Percentile_3 = stats.scoreatpercentile(AP_VEI_3, Percentile_3)

        VEI_4_Median = stats.scoreatpercentile(AP_VEI_4, Percentile_2)
        VEI_4_Percentile_1 = stats.scoreatpercentile(AP_VEI_4, Percentile_1)
        VEI_4_Percentile_3 = stats.scoreatpercentile(AP_VEI_4, Percentile_3)

        VEI_5_Median = stats.scoreatpercentile(AP_VEI_5, Percentile_2)
        VEI_5_Percentile_1 = stats.scoreatpercentile(AP_VEI_5, Percentile_1)
        VEI_5_Percentile_3 = stats.scoreatpercentile(AP_VEI_5, Percentile_3)

        VEI_6_Median = stats.scoreatpercentile(AP_VEI_6, Percentile_2)
        VEI_6_Percentile_1 = stats.scoreatpercentile(AP_VEI_6, Percentile_1)
        VEI_6_Percentile_3 = stats.scoreatpercentile(AP_VEI_6, Percentile_3)

        VEI_7_Median = stats.scoreatpercentile(AP_VEI_7, Percentile_2)
        VEI_7_Percentile_1 = stats.scoreatpercentile(AP_VEI_7, Percentile_1)
        VEI_7_Percentile_3 = stats.scoreatpercentile(AP_VEI_7, Percentile_3)

        Prob_erupt_Median = VEI_1_Median + VEI_2_Median + VEI_3_Median + VEI_4_Median + VEI_5_Median + VEI_6_Median + VEI_7_Median
        Prob_erupt_Percentile_1 = VEI_1_Percentile_1 + VEI_2_Percentile_1 + VEI_3_Percentile_1 + VEI_4_Percentile_1 + VEI_5_Percentile_1 + VEI_6_Percentile_1 + VEI_7_Percentile_1
        Prob_erupt_Percentile_3 = VEI_1_Percentile_3 + VEI_2_Percentile_3 + VEI_3_Percentile_3 + VEI_4_Percentile_3 + VEI_5_Percentile_3 + VEI_6_Percentile_3 + VEI_7_Percentile_3

        Percentile_1s = [[1, 2, 3, 4, 5, 6, 7],
                          [VEI_1_Percentile_1, VEI_2_Percentile_1, VEI_3_Percentile_1, VEI_4_Percentile_1, VEI_5_Percentile_1,
                           VEI_6_Percentile_1,
                           VEI_7_Percentile_1]]
        Percentile_3s = [[1, 2, 3, 4, 5, 6, 7],
                          [VEI_1_Percentile_3, VEI_2_Percentile_3, VEI_3_Percentile_3, VEI_4_Percentile_3, VEI_5_Percentile_3,
                           VEI_6_Percentile_3,
                           VEI_7_Percentile_3]]

        Probabilities = pd.DataFrame([[1, VEI_1_Median],
                                      [2, VEI_2_Median],
                                      [3, VEI_3_Median],
                                      [4, VEI_4_Median],
                                      [5, VEI_5_Median],
                                      [6, VEI_6_Median],
                                      [7, VEI_7_Median],
                                      [1, VEI_1_Percentile_1],
                                      [2, VEI_2_Percentile_1],
                                      [3, VEI_3_Percentile_1],
                                      [4, VEI_4_Percentile_1],
                                      [5, VEI_5_Percentile_1],
                                      [6, VEI_6_Percentile_1],
                                      [7, VEI_7_Percentile_1],
                                      [1, VEI_1_Percentile_3],
                                      [2, VEI_2_Percentile_3],
                                      [3, VEI_3_Percentile_3],
                                      [4, VEI_4_Percentile_3],
                                      [5, VEI_5_Percentile_3],
                                      [6, VEI_6_Percentile_3],
                                      [7, VEI_7_Percentile_3]],
                                     columns=["VEI", "Annual probability"])

        Volcano_data = GVP_volcanoes[(GVP_volcanoes['Volcano Number'] == volc_num)]
        Volcano_name = Volcano_data.iloc[0]['Volcano Name']
        Volcano_Number = volc_num
        textstr = Volcano_name

        Probabilities_df = pd.DataFrame([[Volcano_Number, Volcano_name, year, eruption_certainty, Analysis_method,
                                          method, averaging_method, Power_law,
                                          Prob_erupt_Percentile_1, Prob_erupt_Median, Prob_erupt_Percentile_3,
                                          VEI_1_Percentile_1, VEI_1_Median, VEI_1_Percentile_3,
                                          VEI_2_Percentile_1, VEI_2_Median, VEI_2_Percentile_3,
                                          VEI_3_Percentile_1, VEI_3_Median, VEI_3_Percentile_3,
                                          VEI_4_Percentile_1, VEI_4_Median, VEI_4_Percentile_3,
                                          VEI_5_Percentile_1, VEI_5_Median, VEI_5_Percentile_3,
                                          VEI_6_Percentile_1, VEI_6_Median, VEI_6_Percentile_3,
                                          VEI_7_Percentile_1, VEI_7_Median, VEI_7_Percentile_3]],
                                        columns=["Volcano Number",
                                                 "Volcano Name",
                                                 "GVP DB Year",
                                                 "Included eruptions",
                                                 "Estimate method",
                                                 "Change point",
                                                 "Average method",
                                                 "Power law",
                                                 "Probability of eruption " + p1_string,
                                                 "Probability of eruption " + p2_string,
                                                 "Probability of eruption " + p3_string,
                                                 "<= VEI 1 " + p1_string,
                                                 "<= VEI 1 " + p2_string,
                                                 "<= VEI 1 " + p3_string,
                                                 "VEI 2 " + p1_string,
                                                 "VEI 2 " + p2_string,
                                                 "VEI 2 " + p3_string,
                                                 "VEI 3 " + p1_string,
                                                 "VEI 3 " + p2_string,
                                                 "VEI 3 " + p3_string,
                                                 "VEI 4 " + p1_string,
                                                 "VEI 4 " + p2_string,
                                                 "VEI 4 " + p3_string,
                                                 "VEI 5 " + p1_string,
                                                 "VEI 5 " + p2_string,
                                                 "VEI 5 " + p3_string,
                                                 "VEI 6 " + p1_string,
                                                 "VEI 6 " + p2_string,
                                                 "VEI 6 " + p3_string,
                                                 "VEI 7 " + p1_string,
                                                 "VEI 7 " + p2_string,
                                                 "VEI 7 " + p3_string])
        ID = str(volc_num)
        if confirmed == "True":
            certainty = "Certain"
        else:
            certainty = "Uncertain"
        path_csv = "Probabilities/" + project_folder + "/" + Volcano_name + "_" + certainty + "_" + method + "_" + yearstr + "_" + power + "_A1.csv"
        Probabilities_df.to_csv(path_csv, index=False)

        ax = sns.lineplot(x="VEI", y="Annual probability", markers=True, data=Probabilities)
        ax.set(yscale="log")
        ax.set_xticks(range(1, 8))  # <--- set the ticks first
        ax.tick_params(labelsize=10)
        ax.set_xticklabels(['<=1', '2', '3', '4', '5', '6', '7'])
        ax.set_xlabel("VEI", fontsize=12)
        ax.set_ylabel("Annual probability", fontsize=12)
        props = dict(boxstyle='round', facecolor='white', alpha=0.5)
        # place a text box in upper left in axes coords
        ax.text(0.80, 0.95, textstr, transform=ax.transAxes, fontsize=12,
                verticalalignment='top', bbox=props)
        path_fig = "Figures/" + project_folder + "/" + Volcano_name + "_" + certainty + "_" + method + "_" + yearstr + "_" + power + "_A1.png"
        plt.savefig(path_fig)
        plt.close()
        # plt.show()

def get_A2_prior_model (GVP_volcanoes, volc_num, Volcano_type_2, Freq_rate_2, Freq_rate_std_2, VEI_freq_2, Power_law, year, confirmed, method, VEI_schema, project_folder, Uncertainty_schema):
    AP_VEI_3 = None
    AP_VEI_4 = None
    AP_VEI_5 = None
    AP_VEI_6 = None
    AP_VEI_7 = None
    yearstr = str(year)
    averaging_method = "None"
    if confirmed == "True":
        eruption_certainty = "Confirmed"
    else:
        eruption_certainty = "Confirmed + Uncertain"
    alpha_analogue_2 = Freq_rate_2[Volcano_type_2] ** 2 * (
            (1 - Freq_rate_2[Volcano_type_2]) / Freq_rate_std_2[Volcano_type_2] ** 2 - 1 / Freq_rate_2[Volcano_type_2])
    beta_analogue_2 = alpha_analogue_2 * (1 / Freq_rate_2[Volcano_type_2] - 1)
    Analysis_method = "Analogue A2"
    print("Initiating model using: ", Analysis_method)
    a_ave = alpha_analogue_2
    b_ave = beta_analogue_2

    if Power_law == "True":
        power = "PowerLaw"
    else:
        power = "NoPowerLaw"

    if VEI_schema == 1:
        # Combining probability of eruption and relative probability of different VEI
        AP_VEI_3 = []
        Number = 10000
        lower, upper = 0, 1
        for i in range(Number):
            """Monte Carlo simulation"""
            AP_VEI3 = (np.random.beta(alpha_analogue_2, beta_analogue_2) * VEI_freq_2['3'])
            AP_VEI_3 = AP_VEI_3 + [AP_VEI3]

        AP_VEI_4 = []
        Number = 10000
        lower, upper = 0, 1
        for i in range(Number):
            """Monte Carlo simulation"""
            # if VEI_freq_2[1] > 0:
            AP_VEI4 = (np.random.beta(alpha_analogue_2, beta_analogue_2) * VEI_freq_2['4'])
            AP_VEI_4 = AP_VEI_4 + [AP_VEI4]
            # else:
            #     AP_VEI4 = 0
            # AP_VEI_4 = AP_VEI_4 + [AP_VEI4]

        # if VEI_freq_2[2] > 0:
        AP_VEI_5 = []
        Number = 10000
        lower, upper = 0, 1
        for i in range(Number):
            """Monte Carlo simulation"""
            # if VEI_freq_2[2] > 0:
            AP_VEI5 = (np.random.beta(alpha_analogue_2, beta_analogue_2) * VEI_freq_2['5'])
            AP_VEI_5 = AP_VEI_5 + [AP_VEI5]

        #         else:
        #             AP_VEI5 = 0
        #         AP_VEI_5 = AP_VEI_5 + [AP_VEI5]
        # else:
        #     AP_VEI_5 = 0

        # if VEI_freq_2[3] > 0:
        AP_VEI_6 = []
        Number = 10000
        lower, upper = 0, 1
        for i in range(Number):
            """Monte Carlo simulation"""
            # if VEI_freq_2[3] > 0:
            AP_VEI6 = (np.random.beta(alpha_analogue_2, beta_analogue_2) * VEI_freq_2['6'])
            AP_VEI_6 = AP_VEI_6 + [AP_VEI6]

        #         else:
        #             AP_VEI6 = 0
        #         AP_VEI_6 = AP_VEI_6 + [AP_VEI6]
        # else:
        #     AP_VEI_6 = 0

    #    if VEI_freq_2[4] > 0:
        AP_VEI_7 = []
        Number = 10000
        lower, upper = 0, 1
        for i in range(Number):
            """Monte Carlo simulation"""
            # if VEI_freq_2[4] > 0:
            AP_VEI7 = (np.random.beta(alpha_analogue_2, beta_analogue_2) * VEI_freq_2['7'])
            AP_VEI_7 = AP_VEI_7 + [AP_VEI7]
        #         else:
        #             AP_VEI7 = 0
        #         AP_VEI_7 = AP_VEI_7 + [AP_VEI7]
        # else:
        #     AP_VEI_7 = 0

        if Uncertainty_schema == 1:
            Percentile_1 = 5
            p1_string = "5th percentile"
            Percentile_2 = 50
            p2_string = "50th percentile"
            Percentile_3 = 95
            p3_string = "95th percentile"
        elif Uncertainty_schema == 2:
            Percentile_1 = 10
            p1_string = "10th percentile"
            Percentile_2 = 50
            p2_string = "50th percentile"
            Percentile_3 = 90
            p3_string = "590th percentile"
        elif Uncertainty_schema == 3:
            Percentile_1 = 25
            p1_string = "25th percentile"
            Percentile_2 = 50
            p2_string = "50th percentile"
            Percentile_3 = 75
            p3_string = "75th percentile"

        VEI_3_Median = stats.scoreatpercentile(AP_VEI_3, Percentile_2)
        VEI_3_Percentile_1 = stats.scoreatpercentile(AP_VEI_3, Percentile_1)
        VEI_3_Percentile_3 = stats.scoreatpercentile(AP_VEI_3, Percentile_3)

        VEI_4_Median = stats.scoreatpercentile(AP_VEI_4, Percentile_2)
        VEI_4_Percentile_1 = stats.scoreatpercentile(AP_VEI_4, Percentile_1)
        VEI_4_Percentile_3 = stats.scoreatpercentile(AP_VEI_4, Percentile_3)

        VEI_5_Median = stats.scoreatpercentile(AP_VEI_5, Percentile_2)
        VEI_5_Percentile_1 = stats.scoreatpercentile(AP_VEI_5, Percentile_1)
        VEI_5_Percentile_3 = stats.scoreatpercentile(AP_VEI_5, Percentile_3)

        VEI_6_Median = stats.scoreatpercentile(AP_VEI_6, Percentile_2)
        VEI_6_Percentile_1 = stats.scoreatpercentile(AP_VEI_6, Percentile_1)
        VEI_6_Percentile_3 = stats.scoreatpercentile(AP_VEI_6, Percentile_3)

        VEI_7_Median = stats.scoreatpercentile(AP_VEI_7, Percentile_2)
        VEI_7_Percentile_1 = stats.scoreatpercentile(AP_VEI_7, Percentile_1)
        VEI_7_Percentile_3 = stats.scoreatpercentile(AP_VEI_7, Percentile_3)

        Prob_erupt_Median = VEI_3_Median + VEI_4_Median + VEI_5_Median + VEI_6_Median + VEI_7_Median
        Prob_erupt_Percentile_1 = VEI_3_Percentile_1 + VEI_4_Percentile_1 + VEI_5_Percentile_1 + VEI_6_Percentile_1 + VEI_7_Percentile_1
        Prob_erupt_Percentile_3 = VEI_3_Percentile_3 + VEI_4_Percentile_3 + VEI_5_Percentile_3 + VEI_6_Percentile_3 + VEI_7_Percentile_3

        Percentile_1s = [[3, 4, 5, 6, 7],
                          [VEI_3_Percentile_1, VEI_4_Percentile_1, VEI_5_Percentile_1, VEI_6_Percentile_1,
                           VEI_7_Percentile_1]]
        Percentile_3s = [[3, 4, 5, 6, 7],
                          [VEI_3_Percentile_3, VEI_4_Percentile_3, VEI_5_Percentile_3, VEI_6_Percentile_3,
                           VEI_7_Percentile_3]]

        Probabilities = pd.DataFrame([[3, VEI_3_Median],
                                      [4, VEI_4_Median],
                                      [5, VEI_5_Median],
                                      [6, VEI_6_Median],
                                      [7, VEI_7_Median],
                                      [3, VEI_3_Percentile_1],
                                      [4, VEI_4_Percentile_1],
                                      [5, VEI_5_Percentile_1],
                                      [6, VEI_6_Percentile_1],
                                      [7, VEI_7_Percentile_1],
                                      [3, VEI_3_Percentile_3],
                                      [4, VEI_4_Percentile_3],
                                      [5, VEI_5_Percentile_3],
                                      [6, VEI_6_Percentile_3],
                                      [7, VEI_7_Percentile_3]],
                                     columns=["VEI", "Annual probability"])

        Volcano_data = GVP_volcanoes[(GVP_volcanoes['Volcano Number'] == volc_num)]
        Volcano_name = Volcano_data.iloc[0]['Volcano Name']
        Volcano_Number = volc_num
        print(Volcano_name)
        textstr = Volcano_name

        Probabilities_df = pd.DataFrame([[Volcano_Number, Volcano_name, year, eruption_certainty, Analysis_method, method, averaging_method, Power_law,
                                          Prob_erupt_Percentile_1, Prob_erupt_Median, Prob_erupt_Percentile_3,
                                          VEI_3_Percentile_1, VEI_3_Median, VEI_3_Percentile_3,
                                          VEI_4_Percentile_1, VEI_4_Median, VEI_4_Percentile_3,
                                          VEI_5_Percentile_1, VEI_5_Median, VEI_5_Percentile_3,
                                          VEI_6_Percentile_1, VEI_6_Median, VEI_6_Percentile_3,
                                          VEI_7_Percentile_1, VEI_7_Median, VEI_7_Percentile_3]],
                                        columns=["Volcano Number",
                                                 "Volcano Name",
                                                 "GVP DB Year",
                                                 "Included eruptions",
                                                 "Estimate method",
                                                 "Change point",
                                                 "Average method",
                                                 "Power law",
                                                 "Probability of eruption " + p1_string,
                                                 "Probability of eruption " + p2_string,
                                                 "Probability of eruption " + p3_string,
                                                 "<= VEI 3 " + p1_string,
                                                 "<= VEI 3 " + p2_string,
                                                 "<= VEI 3 " + p3_string,
                                                 "VEI 4 " + p1_string,
                                                 "VEI 4 " + p2_string,
                                                 "VEI 4 " + p3_string,
                                                 "VEI 5 " + p1_string,
                                                 "VEI 5 " + p2_string,
                                                 "VEI 5 " + p3_string,
                                                 "VEI 6 " + p1_string,
                                                 "VEI 6 " + p2_string,
                                                 "VEI 6 " + p3_string,
                                                 "VEI 7 " + p1_string,
                                                 "VEI 7 " + p2_string,
                                                 "VEI 7 " + p3_string])
        ID = str(volc_num)
        if confirmed == "True":
            certainty = "Certain"
        else:
            certainty = "Uncertain"
        path_csv = "Probabilities/" + project_folder + "/" + Volcano_name + "_" + certainty + "_" + method + "_" + yearstr + "_" + power + "_A2.csv"
        Probabilities_df.to_csv(path_csv, index=False)

        ax = sns.lineplot(x="VEI", y="Annual probability", markers=True, data=Probabilities)
        ax.set(yscale="log")
        ax.set_xticks(range(3, 8))  # <--- set the ticks first
        ax.tick_params(labelsize=10)
        ax.set_xticklabels(['<=3', '4', '5', '6', '7'])
        ax.set_xlabel("VEI", fontsize=12)
        ax.set_ylabel("Annual probability", fontsize=12)
        props = dict(boxstyle='round', facecolor='white', alpha=0.5)
        # place a text box in upper left in axes coords
        ax.text(0.80, 0.95, textstr, transform=ax.transAxes, fontsize=12,
                verticalalignment='top', bbox=props)
        path_fig = "Figures/" + project_folder + "/" + Volcano_name + "_" + certainty + "_" + method + "_" + yearstr + "_" + power + "_A2.png"
        plt.savefig(path_fig)
        plt.close()
        #plt.show()

    elif VEI_schema == 2:
        # VEI <=2, 3, 4, 5, 6, 7
        # Combining probability of eruption and relative probability of different VEI

        AP_VEI_2 = []
        Number = 10000
        lower, upper = 0, 1
        for i in range(Number):
            """Monte Carlo simulation"""
            AP_VEI2 = (np.random.beta(alpha_analogue_2, beta_analogue_2) * VEI_freq_2['2'])
            AP_VEI_2 = AP_VEI_2 + [AP_VEI2]

        AP_VEI_3 = []
        Number = 10000
        lower, upper = 0, 1
        for i in range(Number):
            """Monte Carlo simulation"""
            AP_VEI3 = (np.random.beta(alpha_analogue_2, beta_analogue_2) * VEI_freq_2['3'])
            AP_VEI_3 = AP_VEI_3 + [AP_VEI3]

        AP_VEI_4 = []
        Number = 10000
        lower, upper = 0, 1
        for i in range(Number):
            """Monte Carlo simulation"""
            # if VEI_freq_2[1] > 0:
            AP_VEI4 = (np.random.beta(alpha_analogue_2, beta_analogue_2) * VEI_freq_2['4'])
            AP_VEI_4 = AP_VEI_4 + [AP_VEI4]
            # else:
            #     AP_VEI4 = 0
            # AP_VEI_4 = AP_VEI_4 + [AP_VEI4]

        # if VEI_freq_2[2] > 0:
        AP_VEI_5 = []
        Number = 10000
        lower, upper = 0, 1
        for i in range(Number):
            """Monte Carlo simulation"""
            # if VEI_freq_2[2] > 0:
            AP_VEI5 = (np.random.beta(alpha_analogue_2, beta_analogue_2) * VEI_freq_2['5'])
            AP_VEI_5 = AP_VEI_5 + [AP_VEI5]

        #         else:
        #             AP_VEI5 = 0
        #         AP_VEI_5 = AP_VEI_5 + [AP_VEI5]
        # else:
        #     AP_VEI_5 = 0

        # if VEI_freq_2[3] > 0:
        AP_VEI_6 = []
        Number = 10000
        lower, upper = 0, 1
        for i in range(Number):
            """Monte Carlo simulation"""
            # if VEI_freq_2[3] > 0:
            AP_VEI6 = (np.random.beta(alpha_analogue_2, beta_analogue_2) * VEI_freq_2['6'])
            AP_VEI_6 = AP_VEI_6 + [AP_VEI6]

        #         else:
        #             AP_VEI6 = 0
        #         AP_VEI_6 = AP_VEI_6 + [AP_VEI6]
        # else:
        #     AP_VEI_6 = 0

    #    if VEI_freq_2[4] > 0:
        AP_VEI_7 = []
        Number = 10000
        lower, upper = 0, 1
        for i in range(Number):
            """Monte Carlo simulation"""
            # if VEI_freq_2[4] > 0:
            AP_VEI7 = (np.random.beta(alpha_analogue_2, beta_analogue_2) * VEI_freq_2['7'])
            AP_VEI_7 = AP_VEI_7 + [AP_VEI7]
        #         else:
        #             AP_VEI7 = 0
        #         AP_VEI_7 = AP_VEI_7 + [AP_VEI7]
        # else:
        #     AP_VEI_7 = 0

        if Uncertainty_schema == 1:
            Percentile_1 = 5
            p1_string = "5th percentile"
            Percentile_2 = 50
            p2_string = "50th percentile"
            Percentile_3 = 95
            p3_string = "95th percentile"
        elif Uncertainty_schema == 2:
            Percentile_1 = 10
            p1_string = "10th percentile"
            Percentile_2 = 50
            p2_string = "50th percentile"
            Percentile_3 = 90
            p3_string = "590th percentile"
        elif Uncertainty_schema == 3:
            Percentile_1 = 25
            p1_string = "25th percentile"
            Percentile_2 = 50
            p2_string = "50th percentile"
            Percentile_3 = 75
            p3_string = "75th percentile"

        VEI_2_Median = stats.scoreatpercentile(AP_VEI_2, Percentile_2)
        VEI_2_Percentile_1 = stats.scoreatpercentile(AP_VEI_2, Percentile_1)
        VEI_2_Percentile_3 = stats.scoreatpercentile(AP_VEI_2, Percentile_3)

        VEI_3_Median = stats.scoreatpercentile(AP_VEI_3, Percentile_2)
        VEI_3_Percentile_1 = stats.scoreatpercentile(AP_VEI_3, Percentile_1)
        VEI_3_Percentile_3 = stats.scoreatpercentile(AP_VEI_3, Percentile_3)

        VEI_4_Median = stats.scoreatpercentile(AP_VEI_4, Percentile_2)
        VEI_4_Percentile_1 = stats.scoreatpercentile(AP_VEI_4, Percentile_1)
        VEI_4_Percentile_3 = stats.scoreatpercentile(AP_VEI_4, Percentile_3)

        VEI_5_Median = stats.scoreatpercentile(AP_VEI_5, Percentile_2)
        VEI_5_Percentile_1 = stats.scoreatpercentile(AP_VEI_5, Percentile_1)
        VEI_5_Percentile_3 = stats.scoreatpercentile(AP_VEI_5, Percentile_3)

        VEI_6_Median = stats.scoreatpercentile(AP_VEI_6, Percentile_2)
        VEI_6_Percentile_1 = stats.scoreatpercentile(AP_VEI_6, Percentile_1)
        VEI_6_Percentile_3 = stats.scoreatpercentile(AP_VEI_6, Percentile_3)

        VEI_7_Median = stats.scoreatpercentile(AP_VEI_7, Percentile_2)
        VEI_7_Percentile_1 = stats.scoreatpercentile(AP_VEI_7, Percentile_1)
        VEI_7_Percentile_3 = stats.scoreatpercentile(AP_VEI_7, Percentile_3)

        Prob_erupt_Median = VEI_2_Median + VEI_3_Median + VEI_4_Median + VEI_5_Median + VEI_6_Median + VEI_7_Median
        Prob_erupt_Percentile_1 = VEI_2_Percentile_1 + VEI_3_Percentile_1 + VEI_4_Percentile_1 + VEI_5_Percentile_1 + VEI_6_Percentile_1 + VEI_7_Percentile_1
        Prob_erupt_Percentile_3 = VEI_2_Percentile_3 + VEI_3_Percentile_3 + VEI_4_Percentile_3 + VEI_5_Percentile_3 + VEI_6_Percentile_3 + VEI_7_Percentile_3


        Percentile_1s = [[2, 3, 4, 5, 6, 7],
                          [VEI_2_Percentile_1, VEI_3_Percentile_1, VEI_4_Percentile_1, VEI_5_Percentile_1, VEI_6_Percentile_1,
                           VEI_7_Percentile_1]]
        Percentile_3s = [[2, 3, 4, 5, 6, 7],
                          [VEI_2_Percentile_3, VEI_3_Percentile_3, VEI_4_Percentile_3, VEI_5_Percentile_3, VEI_6_Percentile_3,
                           VEI_7_Percentile_3]]

        Probabilities = pd.DataFrame([[2, VEI_2_Median],
                                      [3, VEI_3_Median],
                                      [4, VEI_4_Median],
                                      [5, VEI_5_Median],
                                      [6, VEI_6_Median],
                                      [7, VEI_7_Median],
                                      [2, VEI_2_Percentile_1],
                                      [3, VEI_3_Percentile_1],
                                      [4, VEI_4_Percentile_1],
                                      [5, VEI_5_Percentile_1],
                                      [6, VEI_6_Percentile_1],
                                      [7, VEI_7_Percentile_1],
                                      [2, VEI_2_Percentile_3],
                                      [3, VEI_3_Percentile_3],
                                      [4, VEI_4_Percentile_3],
                                      [5, VEI_5_Percentile_3],
                                      [6, VEI_6_Percentile_3],
                                      [7, VEI_7_Percentile_3]],
                                     columns=["VEI", "Annual probability"])

        Volcano_data = GVP_volcanoes[(GVP_volcanoes['Volcano Number'] == volc_num)]
        Volcano_name = Volcano_data.iloc[0]['Volcano Name']
        Volcano_Number = volc_num
        print(Volcano_name)
        textstr = Volcano_name

        Probabilities_df = pd.DataFrame([[Volcano_Number, Volcano_name, year, eruption_certainty, Analysis_method, method, averaging_method, Power_law,
                                          Prob_erupt_Percentile_1, Prob_erupt_Median, Prob_erupt_Percentile_3,
                                          VEI_2_Percentile_1, VEI_2_Median, VEI_2_Percentile_3,
                                          VEI_3_Percentile_1, VEI_3_Median, VEI_3_Percentile_3,
                                          VEI_4_Percentile_1, VEI_4_Median, VEI_4_Percentile_3,
                                          VEI_5_Percentile_1, VEI_5_Median, VEI_5_Percentile_3,
                                          VEI_6_Percentile_1, VEI_6_Median, VEI_6_Percentile_3,
                                          VEI_7_Percentile_1, VEI_7_Median, VEI_7_Percentile_3]],
                                        columns=["Volcano Number",
                                                 "Volcano Name",
                                                 "GVP DB Year",
                                                 "Included eruptions",
                                                 "Estimate method",
                                                 "Change point",
                                                 "Average method",
                                                 "Power law",
                                                 "Probability of eruption " + p1_string,
                                                 "Probability of eruption " + p2_string,
                                                 "Probability of eruption " + p3_string,
                                                 "<= VEI 2 " + p1_string,
                                                 "<= VEI 2 " + p2_string,
                                                 "<= VEI 2 " + p3_string,
                                                 "VEI 3 " + p1_string,
                                                 "VEI 3 " + p2_string,
                                                 "VEI 3 " + p3_string,
                                                 "VEI 4 " + p1_string,
                                                 "VEI 4 " + p2_string,
                                                 "VEI 4 " + p3_string,
                                                 "VEI 5 " + p1_string,
                                                 "VEI 5 " + p2_string,
                                                 "VEI 5 " + p3_string,
                                                 "VEI 6 " + p1_string,
                                                 "VEI 6 " + p2_string,
                                                 "VEI 6 " + p3_string,
                                                 "VEI 7 " + p1_string,
                                                 "VEI 7 " + p2_string,
                                                 "VEI 7 " + p3_string])
        ID = str(volc_num)
        if confirmed == "True":
            certainty = "Certain"
        else:
            certainty = "Uncertain"
        path_csv = "Probabilities/" + project_folder + "/" + Volcano_name + "_" + certainty + "_" + method + "_" + yearstr + "_" + power + "_A2.csv"
        Probabilities_df.to_csv(path_csv, index=False)

        ax = sns.lineplot(x="VEI", y="Annual probability", markers=True, data=Probabilities)
        ax.set(yscale="log")
        ax.set_xticks(range(2, 8))  # <--- set the ticks first
        ax.tick_params(labelsize=10)
        ax.set_xticklabels(['<=2', '3', '4', '5', '6', '7'])
        ax.set_xlabel("VEI", fontsize=12)
        ax.set_ylabel("Annual probability", fontsize=12)
        props = dict(boxstyle='round', facecolor='white', alpha=0.5)
        # place a text box in upper left in axes coords
        ax.text(0.80, 0.95, textstr, transform=ax.transAxes, fontsize=12,
                verticalalignment='top', bbox=props)
        path_fig = "Figures/" + project_folder + "/" + Volcano_name + "_" + certainty + "_" + method + "_" + yearstr + "_" + power + "_A2.png"
        plt.savefig(path_fig)
        plt.close()
        #plt.show()

    elif VEI_schema == 3:
        # VEI <=1, 2, 3, 4, 5, 6, 7
        # Combining probability of eruption and relative probability of different VEI
        AP_VEI_1 = []
        Number = 10000
        lower, upper = 0, 1
        for i in range(Number):
            """Monte Carlo simulation"""
            AP_VEI1 = (np.random.beta(alpha_analogue_2, beta_analogue_2) * VEI_freq_2['1'])
            AP_VEI_1 = AP_VEI_1 + [AP_VEI1]

        AP_VEI_2 = []
        Number = 10000
        lower, upper = 0, 1
        for i in range(Number):
            """Monte Carlo simulation"""
            AP_VEI2 = (np.random.beta(alpha_analogue_2, beta_analogue_2) * VEI_freq_2['2'])
            AP_VEI_2 = AP_VEI_2 + [AP_VEI2]

        AP_VEI_3 = []
        Number = 10000
        lower, upper = 0, 1
        for i in range(Number):
            """Monte Carlo simulation"""
            AP_VEI3 = (np.random.beta(alpha_analogue_2, beta_analogue_2) * VEI_freq_2['3'])
            AP_VEI_3 = AP_VEI_3 + [AP_VEI3]

        AP_VEI_4 = []
        Number = 10000
        lower, upper = 0, 1
        for i in range(Number):
            """Monte Carlo simulation"""
            # if VEI_freq_2[1] > 0:
            AP_VEI4 = (np.random.beta(alpha_analogue_2, beta_analogue_2) * VEI_freq_2['4'])
            AP_VEI_4 = AP_VEI_4 + [AP_VEI4]
            # else:
            #     AP_VEI4 = 0
            # AP_VEI_4 = AP_VEI_4 + [AP_VEI4]

        # if VEI_freq_2[2] > 0:
        AP_VEI_5 = []
        Number = 10000
        lower, upper = 0, 1
        for i in range(Number):
            """Monte Carlo simulation"""
            # if VEI_freq_2[2] > 0:
            AP_VEI5 = (np.random.beta(alpha_analogue_2, beta_analogue_2) * VEI_freq_2['5'])
            AP_VEI_5 = AP_VEI_5 + [AP_VEI5]

        #         else:
        #             AP_VEI5 = 0
        #         AP_VEI_5 = AP_VEI_5 + [AP_VEI5]
        # else:
        #     AP_VEI_5 = 0

        # if VEI_freq_2[3] > 0:
        AP_VEI_6 = []
        Number = 10000
        lower, upper = 0, 1
        for i in range(Number):
            """Monte Carlo simulation"""
            # if VEI_freq_2[3] > 0:
            AP_VEI6 = (np.random.beta(alpha_analogue_2, beta_analogue_2) * VEI_freq_2['6'])
            AP_VEI_6 = AP_VEI_6 + [AP_VEI6]

        #         else:
        #             AP_VEI6 = 0
        #         AP_VEI_6 = AP_VEI_6 + [AP_VEI6]
        # else:
        #     AP_VEI_6 = 0

        #    if VEI_freq_2[4] > 0:
        AP_VEI_7 = []
        Number = 10000
        lower, upper = 0, 1
        for i in range(Number):
            """Monte Carlo simulation"""
            # if VEI_freq_2[4] > 0:
            AP_VEI7 = (np.random.beta(alpha_analogue_2, beta_analogue_2) * VEI_freq_2['7'])
            AP_VEI_7 = AP_VEI_7 + [AP_VEI7]
        #         else:
        #             AP_VEI7 = 0
        #         AP_VEI_7 = AP_VEI_7 + [AP_VEI7]
        # else:
        #     AP_VEI_7 = 0

        if Uncertainty_schema == 1:
            Percentile_1 = 5
            p1_string = "5th percentile"
            Percentile_2 = 50
            p2_string = "50th percentile"
            Percentile_3 = 95
            p3_string = "95th percentile"
        elif Uncertainty_schema == 2:
            Percentile_1 = 10
            p1_string = "10th percentile"
            Percentile_2 = 50
            p2_string = "50th percentile"
            Percentile_3 = 90
            p3_string = "590th percentile"
        elif Uncertainty_schema == 3:
            Percentile_1 = 25
            p1_string = "25th percentile"
            Percentile_2 = 50
            p1_string = "50th percentile"
            Percentile_3 = 75
            p1_string = "75th percentile"

        VEI_1_Median = stats.scoreatpercentile(AP_VEI_1, Percentile_2)
        VEI_1_Percentile_1 = stats.scoreatpercentile(AP_VEI_1, Percentile_1)
        VEI_1_Percentile_3 = stats.scoreatpercentile(AP_VEI_1, Percentile_3)

        VEI_2_Median = stats.scoreatpercentile(AP_VEI_2, Percentile_2)
        VEI_2_Percentile_1 = stats.scoreatpercentile(AP_VEI_2, Percentile_1)
        VEI_2_Percentile_3 = stats.scoreatpercentile(AP_VEI_2, Percentile_3)

        VEI_3_Median = stats.scoreatpercentile(AP_VEI_3, Percentile_2)
        VEI_3_Percentile_1 = stats.scoreatpercentile(AP_VEI_3, Percentile_1)
        VEI_3_Percentile_3 = stats.scoreatpercentile(AP_VEI_3, Percentile_3)

        VEI_4_Median = stats.scoreatpercentile(AP_VEI_4, Percentile_2)
        VEI_4_Percentile_1 = stats.scoreatpercentile(AP_VEI_4, Percentile_1)
        VEI_4_Percentile_3 = stats.scoreatpercentile(AP_VEI_4, Percentile_3)

        VEI_5_Median = stats.scoreatpercentile(AP_VEI_5, Percentile_2)
        VEI_5_Percentile_1 = stats.scoreatpercentile(AP_VEI_5, Percentile_1)
        VEI_5_Percentile_3 = stats.scoreatpercentile(AP_VEI_5, Percentile_3)

        VEI_6_Median = stats.scoreatpercentile(AP_VEI_6, Percentile_2)
        VEI_6_Percentile_1 = stats.scoreatpercentile(AP_VEI_6, Percentile_1)
        VEI_6_Percentile_3 = stats.scoreatpercentile(AP_VEI_6, Percentile_3)

        VEI_7_Median = stats.scoreatpercentile(AP_VEI_7, Percentile_2)
        VEI_7_Percentile_1 = stats.scoreatpercentile(AP_VEI_7, Percentile_1)
        VEI_7_Percentile_3 = stats.scoreatpercentile(AP_VEI_7, Percentile_3)

        Prob_erupt_Median = VEI_1_Median + VEI_2_Median + VEI_3_Median + VEI_4_Median + VEI_5_Median + VEI_6_Median + VEI_7_Median
        Prob_erupt_Percentile_1 = VEI_1_Percentile_1 + VEI_2_Percentile_1 + VEI_3_Percentile_1 + VEI_4_Percentile_1 + VEI_5_Percentile_1 + VEI_6_Percentile_1 + VEI_7_Percentile_1
        Prob_erupt_Percentile_3 = VEI_1_Percentile_3 + VEI_2_Percentile_3 + VEI_3_Percentile_3 + VEI_4_Percentile_3 + VEI_5_Percentile_3 + VEI_6_Percentile_3 + VEI_7_Percentile_3

        Percentile_1s = [[1, 2, 3, 4, 5, 6, 7],
                          [VEI_1_Percentile_1, VEI_2_Percentile_1, VEI_3_Percentile_1, VEI_4_Percentile_1, VEI_5_Percentile_1,
                           VEI_6_Percentile_1,
                           VEI_7_Percentile_1]]
        Percentile_3s = [[1, 2, 3, 4, 5, 6, 7],
                          [VEI_1_Percentile_3, VEI_2_Percentile_3, VEI_3_Percentile_3, VEI_4_Percentile_3, VEI_5_Percentile_3,
                           VEI_6_Percentile_3,
                           VEI_7_Percentile_3]]

        Probabilities = pd.DataFrame([[1, VEI_1_Median],
                                      [2, VEI_2_Median],
                                      [3, VEI_3_Median],
                                      [4, VEI_4_Median],
                                      [5, VEI_5_Median],
                                      [6, VEI_6_Median],
                                      [7, VEI_7_Median],
                                      [1, VEI_1_Percentile_1],
                                      [2, VEI_2_Percentile_1],
                                      [3, VEI_3_Percentile_1],
                                      [4, VEI_4_Percentile_1],
                                      [5, VEI_5_Percentile_1],
                                      [6, VEI_6_Percentile_1],
                                      [7, VEI_7_Percentile_1],
                                      [1, VEI_1_Percentile_3],
                                      [2, VEI_2_Percentile_3],
                                      [3, VEI_3_Percentile_3],
                                      [4, VEI_4_Percentile_3],
                                      [5, VEI_5_Percentile_3],
                                      [6, VEI_6_Percentile_3],
                                      [7, VEI_7_Percentile_3]],
                                     columns=["VEI", "Annual probability"])

        Volcano_data = GVP_volcanoes[(GVP_volcanoes['Volcano Number'] == volc_num)]
        Volcano_name = Volcano_data.iloc[0]['Volcano Name']
        Volcano_Number = volc_num
        print(Volcano_name)
        textstr = Volcano_name

        Probabilities_df = pd.DataFrame([[Volcano_Number, Volcano_name, year, eruption_certainty, Analysis_method,
                                          method, averaging_method, Power_law,
                                          Prob_erupt_Percentile_1, Prob_erupt_Median, Prob_erupt_Percentile_3,
                                          VEI_1_Percentile_1, VEI_1_Median, VEI_1_Percentile_3,
                                          VEI_2_Percentile_1, VEI_2_Median, VEI_2_Percentile_3,
                                          VEI_3_Percentile_1, VEI_3_Median, VEI_3_Percentile_3,
                                          VEI_4_Percentile_1, VEI_4_Median, VEI_4_Percentile_3,
                                          VEI_5_Percentile_1, VEI_5_Median, VEI_5_Percentile_3,
                                          VEI_6_Percentile_1, VEI_6_Median, VEI_6_Percentile_3,
                                          VEI_7_Percentile_1, VEI_7_Median, VEI_7_Percentile_3]],
                                        columns=["Volcano Number",
                                                 "Volcano Name",
                                                 "GVP DB Year",
                                                 "Included eruptions",
                                                 "Estimate method",
                                                 "Change point",
                                                 "Average method",
                                                 "Power law",
                                                 "Probability of eruption " + p1_string,
                                                 "Probability of eruption " + p2_string,
                                                 "Probability of eruption " + p3_string,
                                                 "<= VEI 1 " + p1_string,
                                                 "<= VEI 1 " + p2_string,
                                                 "<= VEI 1 " + p3_string,
                                                 "VEI 2 " + p1_string,
                                                 "VEI 2 " + p2_string,
                                                 "VEI 2 " + p3_string,
                                                 "VEI 3 " + p1_string,
                                                 "VEI 3 " + p2_string,
                                                 "VEI 3 " + p3_string,
                                                 "VEI 4 " + p1_string,
                                                 "VEI 4 " + p2_string,
                                                 "VEI 4 " + p3_string,
                                                 "VEI 5 " + p1_string,
                                                 "VEI 5 " + p2_string,
                                                 "VEI 5 " + p3_string,
                                                 "VEI 6 " + p1_string,
                                                 "VEI 6 " + p2_string,
                                                 "VEI 6 " + p3_string,
                                                 "VEI 7 " + p1_string,
                                                 "VEI 7 " + p2_string,
                                                 "VEI 7 " + p3_string])
        ID = str(volc_num)
        if confirmed == "True":
            certainty = "Certain"
        else:
            certainty = "Uncertain"
        path_csv = "Probabilities/" + project_folder + "/" + Volcano_name + "_" + certainty + "_" + method + "_" + yearstr + "_" + power + "_A2.csv"
        Probabilities_df.to_csv(path_csv, index=False)

        ax = sns.lineplot(x="VEI", y="Annual probability", markers=True, data=Probabilities)
        ax.set(yscale="log")
        ax.set_xticks(range(1, 8))  # <--- set the ticks first
        ax.tick_params(labelsize=10)
        ax.set_xticklabels(['<=1', '2', '3', '4', '5', '6', '7'])
        ax.set_xlabel("VEI", fontsize=12)
        ax.set_ylabel("Annual probability", fontsize=12)
        props = dict(boxstyle='round', facecolor='white', alpha=0.5)
        # place a text box in upper left in axes coords
        ax.text(0.80, 0.95, textstr, transform=ax.transAxes, fontsize=12,
                verticalalignment='top', bbox=props)
        path_fig = "Figures/" + project_folder + "/" + Volcano_name + "_" + certainty + "_" + method + "_" + yearstr + "_" + power + "_A2.png"
        plt.savefig(path_fig)
        plt.close()
        # plt.show()

# fix this function with the new additions
def get_Bayes_update (GVP_volcanoes, volc_num, Volcano_type_1, Observed_rate, Observed_VEI, confirmed_VEI_eruptions, Freq_rate_1, Freq_rate_std_1, VEI_freq_1, Power_law, averaging_method, year, confirmed, method, analogue_model):
    AP_VEI_3 = None
    AP_VEI_4 = None
    AP_VEI_5 = None
    AP_VEI_6 = None
    AP_VEI_7 = None
    yearstr = str(year)
    if confirmed == "True":
        eruption_certainty = "Confirmed"
    else:
        eruption_certainty = "Confirmed + Uncertain"
    alpha_analogue_1 = Freq_rate_1[Volcano_type_1] ** 2 * (
            (1 - Freq_rate_1[Volcano_type_1]) / Freq_rate_std_1[Volcano_type_1] ** 2 - 1 / Freq_rate_1[Volcano_type_1])
    beta_analogue_1 = alpha_analogue_1 * (1 / Freq_rate_1[Volcano_type_1] - 1)

    if confirmed_VEI_eruptions == 0:
        print("No confirmed VEI eruptions")

    else:
        Analysis_method = "Bayesian update " + analogue_model
        print("Initiating model using: ", Analysis_method)

        ###
        # Model averaging for frequency of eruption of any magnitude
        # Probability of an eruption of any VEI
        cores = 1
        y = np.array([Observed_rate * 10000])
        n = np.array([10000])
        N = len(n)
        # if __name__ == '__main__':
        # Probability of eruption model
        with pm.Model() as AP_model_1:
            theta = pm.Beta('theta', alpha=alpha_analogue_1, beta=beta_analogue_1)
            p = pm.Binomial('y', p=theta, observed=y, n=n)
            AP_trace_1 = pm.sample(20000, tune=10000, chains=2, cores=cores, target_accept=0.99)
        ###

        Ave = AP_trace_1['theta'].mean()
        Std = AP_trace_1['theta'].std()
        a_ave = Ave**2 * ((1-Ave)/Std**2 -1/Ave)
        b_ave = a_ave * (1/Ave-1)
        data_beta = beta.rvs(a=a_ave, b=b_ave, size=1000000)

        # # Comparing conditional probability of each VEI eruption

        # # # Model 1

        ###
        if VEI_freq_1["7"].mean() > 0:
            print("Using first concentration parameter for Model 1")
            concentration_1 = np.array([VEI_freq_1["3"].mean() * 10000,
                                        VEI_freq_1["4"].mean() * 10000,
                                        VEI_freq_1["5"].mean() * 10000,
                                        VEI_freq_1["6"].mean() * 10000,
                                        VEI_freq_1["7"].mean() * 10000])
            Observed = np.array([Observed_VEI[0],
                                 Observed_VEI[1],
                                 Observed_VEI[2],
                                 Observed_VEI[3],
                                 Observed_VEI[4]])

        elif VEI_freq_1["6"].mean() > 0 and VEI_freq_1["7"].mean() == 0:
            print("Using second concentration parameter for Model 1")
            concentration_1 = np.array([VEI_freq_1["3"] * 10000,
                                        VEI_freq_1["4"] * 10000,
                                        VEI_freq_1["5"] * 10000,
                                        VEI_freq_1["6"] * 10000])
            Observed = np.array([Observed_VEI[0],
                                 Observed_VEI[1],
                                 Observed_VEI[2],
                                 Observed_VEI[3]])

        elif VEI_freq_1["5"].mean() > 0 and VEI_freq_1["6"].mean() == 0:
            print("Using third concentration parameter for Model 1")
            concentration_1 = np.array([VEI_freq_1["3"].mean() * 10000,
                                        VEI_freq_1["4"].mean() * 10000,
                                        VEI_freq_1["5"].mean() * 10000])
            Observed = np.array([Observed_VEI[0],
                                 Observed_VEI[1],
                                 Observed_VEI[2]])

        elif VEI_freq_1["4"].mean() > 0 and VEI_freq_1["5"].mean() == 0:
            print("Using fourth concentration parameter for Model 1")
            concentration_1 = np.array([VEI_freq_1["3"].mean() * 10000,
                                        VEI_freq_1["4"].mean() * 10000])
            Observed = np.array([Observed_VEI[0],
                                 Observed_VEI[1]])

        else:
            print("Error with concentration parameter or observed parameter for Model 1")
            # Jenkins_c = Jenkins_1
        shape_1 = concentration_1.size

        Observed_sum = np.sum(Observed)

        ###

        # if VEI_freq_1[4] > 0:
        #     print("Using first concentration parameter for Model 1")
        #     concentration_1 = np.array([VEI_freq_1[0] * 10000,
        #                                 VEI_freq_1[1] * 10000,
        #                                 VEI_freq_1[2] * 10000,
        #                                 VEI_freq_1[3] * 10000,
        #                                 VEI_freq_1[4] * 10000])
        #     Observed = np.array([Observed_VEI[0],
        #                          Observed_VEI[1],
        #                          Observed_VEI[2],
        #                          Observed_VEI[3],
        #                          Observed_VEI[4]])
        #
        # elif VEI_freq_1[3] > 0 and VEI_freq_1[4] == 0:
        #     print("Using second concentration parameter for Model 1")
        #     concentration_1 = np.array([VEI_freq_1[0] * 10000,
        #                                 VEI_freq_1[1] * 10000,
        #                                 VEI_freq_1[2] * 10000,
        #                                 VEI_freq_1[3] * 10000])
        #     Observed = np.array([Observed_VEI[0],
        #                          Observed_VEI[1],
        #                          Observed_VEI[2],
        #                          Observed_VEI[3]])
        #
        # elif VEI_freq_1[2] > 0 and VEI_freq_1[3] == 0:
        #     print("Using third concentration parameter for Model 1")
        #     concentration_1 = np.array([VEI_freq_1[0] * 10000,
        #                                 VEI_freq_1[1] * 10000,
        #                                 VEI_freq_1[2] * 10000])
        #     Observed = np.array([Observed_VEI[0],
        #                          Observed_VEI[1],
        #                          Observed_VEI[2]])
        #
        # elif VEI_freq_1[1] > 0 and VEI_freq_1[2] == 0:
        #     print("Using fourth concentration parameter for Model 1")
        #     concentration_1 = np.array([VEI_freq_1[0] * 10000,
        #                                 VEI_freq_1[1] * 10000])
        #     Observed = np.array([Observed_VEI[0],
        #                          Observed_VEI[1]])
        #
        # else:
        #     print("Error with concentration parameter or observed parameter for Model 1")
        # # Jenkins_c = Jenkins_1
        # shape_1 = concentration_1.size
        #
        #
        # Observed_sum = np.sum(Observed)

        # Create model 1

        Create_graph = "Yes"
        with pm.Model() as VEI_model_1:
            parameters = pm.Dirichlet('parameters', a=concentration_1, shape=shape_1)
            # Observed data is from a Multinomial distribution
            observed_data = pm.Multinomial(
                'observed_data', n=Observed_sum, p=parameters, shape=shape_1, observed=Observed)

        with VEI_model_1:
            # Sample from the posterior
            VEI_trace_1 = pm.sample(draws=10000, chains=2, cores=cores, tune=5000,
                                      discard_tuned_samples=True, target_accept=0.99)
            VEI_trace_data = VEI_trace_1.get_values(parameters, burn=5000)
            VEI_trace_data_df = pd.DataFrame(VEI_trace_data)
            ppc_VEI = pm.sample_posterior_predictive(VEI_trace_1, var_names=['observed_data'],
                                                     random_seed=42)

        if shape_1 == 5:
            try:
                ppc_ave_VEI_df = pd.DataFrame(ppc_VEI)
                ppc_ave_VEI_df2 = pd.DataFrame(ppc_ave_VEI_df['observed_data'].to_list(),
                                               columns=['VEI3', 'VEI4', 'VEI5', 'VEI6', 'VEI7'])
                ppc_ave_VEI_df2 = ppc_ave_VEI_df2.fillna(0)
            except ValueError:
                ppc_ave_VEI_df = pd.DataFrame(ppc_VEI['observed_data'])
                ppc_ave_VEI_df2 = ppc_ave_VEI_df
                ppc_ave_VEI_df2.columns = ['VEI3', 'VEI4', 'VEI5', 'VEI6', 'VEI7']

            # arviz.plot_density(trace_VEI_GVP)
            # plt.show()

            Ave_VEI3 = VEI_trace_data_df.iloc[:, 0].mean()
            Std_VEI3 = VEI_trace_data_df.iloc[:, 0].std()

            Ave_VEI4 = VEI_trace_data_df.iloc[:, 1].mean()
            Std_VEI4 = VEI_trace_data_df.iloc[:, 1].std()

            Ave_VEI5 = VEI_trace_data_df.iloc[:, 2].mean()
            Std_VEI5 = VEI_trace_data_df.iloc[:, 2].std()

            Ave_VEI6 = VEI_trace_data_df.iloc[:, 3].mean()
            Std_VEI6 = VEI_trace_data_df.iloc[:, 3].std()

            Ave_VEI7 = VEI_trace_data_df.iloc[:, 4].mean()
            Std_VEI7 = VEI_trace_data_df.iloc[:, 4].std()

        elif shape_1 == 3:
            try:
                ppc_ave_VEI_df = pd.DataFrame(ppc_VEI)
                ppc_ave_VEI_df2 = pd.DataFrame(ppc_ave_VEI_df['observed_data'].to_list(),
                                               columns=['VEI3', 'VEI4', 'VEI5'])
                ppc_ave_VEI_df2 = ppc_ave_VEI_df2.fillna(0)
            except ValueError:
                ppc_ave_VEI_df = pd.DataFrame(ppc_VEI['observed_data'])
                ppc_ave_VEI_df2 = ppc_ave_VEI_df
                ppc_ave_VEI_df2.columns = ['VEI3', 'VEI4', 'VEI5']

            # arviz.plot_density(trace_VEI_GVP)
            # plt.show()

            Ave_VEI3 = VEI_trace_data_df.iloc[:, 0].mean()
            Std_VEI3 = VEI_trace_data_df.iloc[:, 0].std()

            Ave_VEI4 = VEI_trace_data_df.iloc[:, 1].mean()
            Std_VEI4 = VEI_trace_data_df.iloc[:, 1].std()

            Ave_VEI5 = VEI_trace_data_df.iloc[:, 2].mean()
            Std_VEI5 = VEI_trace_data_df.iloc[:, 2].std()

            Ave_VEI6 = 0
            Std_VEI6 = 0

            Ave_VEI7 = 0
            Std_VEI7 = 0

        elif shape_1 == 2:
            try:
                ppc_ave_VEI_df = pd.DataFrame(ppc_VEI)
                ppc_ave_VEI_df2 = pd.DataFrame(ppc_ave_VEI_df['observed_data'].to_list(),
                                               columns=['VEI3', 'VEI4'])
                ppc_ave_VEI_df2 = ppc_ave_VEI_df2.fillna(0)
            except ValueError:
                ppc_ave_VEI_df = pd.DataFrame(ppc_VEI['observed_data'])
                ppc_ave_VEI_df2 = ppc_ave_VEI_df
                ppc_ave_VEI_df2.columns = ['VEI3', 'VEI4']

            # arviz.plot_density(trace_VEI_GVP)
            # plt.show()

            Ave_VEI3 = VEI_trace_data_df.iloc[:, 0].mean()
            Std_VEI3 = VEI_trace_data_df.iloc[:, 0].std()

            Ave_VEI4 = VEI_trace_data_df.iloc[:, 1].mean()
            Std_VEI4 = VEI_trace_data_df.iloc[:, 1].std()

            Ave_VEI5 = 0
            Std_VEI5 = 0

            Ave_VEI6 = 0
            Std_VEI6 = 0

            Ave_VEI7 = 0
            Std_VEI7 = 0

        if Ave_VEI3 > 0:
            a_ave_VEI3 = Ave_VEI3 ** 2 * ((1 - Ave_VEI3) / Std_VEI3 ** 2 - 1 / Ave_VEI3)
            b_ave_VEI3 = a_ave_VEI3 * (1 / Ave_VEI3 - 1)
            data_beta_VEI3 = beta.rvs(a=a_ave_VEI3, b=b_ave_VEI3, size=1000000)
        else:
            a_ave_VEI3 = 0
            b_ave_VEI3 = 0
            data_beta_VEI3 = 0

        if Ave_VEI4 > 0:
            a_ave_VEI4 = Ave_VEI4 ** 2 * ((1 - Ave_VEI4) / Std_VEI4 ** 2 - 1 / Ave_VEI4)
            b_ave_VEI4 = a_ave_VEI4 * (1 / Ave_VEI4 - 1)
            data_beta_VEI4 = beta.rvs(a=a_ave_VEI4, b=b_ave_VEI4, size=1000000)
        else:
            a_ave_VEI4 = 0
            b_ave_VEI4 = 0
            data_beta_VEI4 = 0

        if Ave_VEI5 > 0:
            a_ave_VEI5 = Ave_VEI5 ** 2 * ((1 - Ave_VEI5) / Std_VEI5 ** 2 - 1 / Ave_VEI5)
            b_ave_VEI5 = a_ave_VEI5 * (1 / Ave_VEI5 - 1)
            data_beta_VEI5 = beta.rvs(a=a_ave_VEI5, b=b_ave_VEI5, size=1000000)
        else:
            a_ave_VEI5 = 0
            b_ave_VEI5 = 0
            data_beta_VEI5 = 0

        if Ave_VEI6 > 0:
            a_ave_VEI6 = Ave_VEI6 ** 2 * ((1 - Ave_VEI6) / Std_VEI6 ** 2 - 1 / Ave_VEI6)
            b_ave_VEI6 = a_ave_VEI6 * (1 / Ave_VEI6 - 1)
            data_beta_VEI6 = beta.rvs(a=a_ave_VEI6, b=b_ave_VEI6, size=1000000)
        else:
            a_ave_VEI6 = 0
            b_ave_VEI6 = 0
            data_beta_VEI6 = 0

        if Ave_VEI7 > 0:
            a_ave_VEI7 = Ave_VEI7 ** 2 * ((1 - Ave_VEI7) / Std_VEI7 ** 2 - 1 / Ave_VEI7)
            b_ave_VEI7 = a_ave_VEI7 * (1 / Ave_VEI7 - 1)
            data_beta_VEI7 = beta.rvs(a=a_ave_VEI7, b=b_ave_VEI7, size=1000000)
        else:
            a_ave_VEI7 = 0
            b_ave_VEI7 = 0
            data_beta_VEI7 = 0

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

        # Combining probability of eruption and relative probability of different VEI


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
            AP_VEI_6 = 0

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

        Volcano_data = GVP_volcanoes[(GVP_volcanoes['Volcano Number'] == volc_num)]
        Volcano_name = Volcano_data.iloc[0]['Volcano Name']
        Volcano_Number = volc_num
        textstr = Volcano_name

        Probabilities_df = pd.DataFrame([[Volcano_Number, Volcano_name, year, eruption_certainty, Analysis_method, method, averaging_method, Power_law,
                                          Prob_erupt_Percentile_05, Prob_erupt_Median, Prob_erupt_Percentile_95,
                                          VEI_3_Percentile_05, VEI_3_Median, VEI_3_Percentile_95,
                                          VEI_4_Percentile_05, VEI_4_Median, VEI_4_Percentile_95,
                                          VEI_5_Percentile_05, VEI_5_Median, VEI_5_Percentile_95,
                                          VEI_6_Percentile_05, VEI_6_Median, VEI_6_Percentile_95,
                                          VEI_7_Percentile_05, VEI_7_Median, VEI_7_Percentile_95]],
                                        columns=["Volcano Number",
                                                 "Volcano Name",
                                                 "GVP DB Year",
                                                 "Included eruptions",
                                                 "Estimate method",
                                                 "Change point",
                                                 "Average method",
                                                 "Power law",
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
        ID = str(volc_num)
        if confirmed == "True":
            certainty = "Certain"
        else:
            certainty = "Uncertain"
        if Power_law == "True":
            power = "PowerLaw"
        else:
            power = "NoPowerLaw"

        path_csv = "Probabilities/Full/" + Volcano_name + "_" + certainty + "_" + yearstr + "_" + method + "_" + Analysis_method + "_" + averaging_method + "_" + power + ".csv"
        Probabilities_df.to_csv(path_csv,  index=False)

        if Create_graph == "No":
            print("No graph created due to sampling error")
        else:
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
            path_fig = "Figures/Full/" + Volcano_name + "_" + certainty + "_" + yearstr + "_" + method + "_" + Analysis_method + "_" + averaging_method + "_" + power + ".png"
            plt.savefig(path_fig)
            plt.close()
            #plt.show()