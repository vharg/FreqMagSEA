import pymc3 as pm
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from scipy.stats import beta
import arviz
# import kneebow
# from kneebow.rotor import Rotor
import csv
from datetime import datetime
from datetime import date

from Functions import assign_A1_classification
from Functions import get_complete_record
from Functions import calc_mean_annual_frequency_analogue
from Functions import calc_std_annual_frequency_analogue
from Functions import calc_relative_probability
from Functions import get_VEI_eruptions_count
from Functions import get_observed_eruption_rate
from Functions import get_observed_relative_frequency
from Functions import get_freq_mag_dirichlet
from Functions import set_analogue_eruption_rate
from Functions import get_observed_eruption_rate_full_record
from Functions import get_observed_relative_frequency_full_record
from Functions import get_master_csv
from Functions import get_model_average
from Functions import get_Bayes_update
from Functions import get_A1_prior_model
from Functions import get_A2_prior_model

startTime = datetime.now()

######
GVP_DB_Year = [2019]
######

selected_volcanoes = pd.read_csv("Selected_volcanoes.csv")
volcanoes = list(selected_volcanoes['volc_num'])
Total_volcanoes = len(volcanoes)
Runstart = date.today()
save_file_name = f"FM_SEA_probabilities_{Runstart}.csv"
project_folder = "Test_run"



#### For system 1
system_1 = "A1"
GVP_volcanoes = pd.read_csv("Volcano_list.csv")
A1_volcanoes = assign_A1_classification(GVP_volcanoes)
A2_volcanoes = pd.read_csv("Whelley_SEA.csv")


# Full_record = pd.read_csv("GVPDB2019.csv")
# GVP_df = Full_record
change_points = pd.read_csv("Change_points_all.csv")

Change_point_method = ["MeadMagill50", "MeadMagill5", "MeadMagill95"]
Eruption_certainty = ["True"]
Power_law_used = ["True", "False"]
Uncertainty_schema = 1
VEI_schema = 1

x = 0

for y in GVP_DB_Year:
    year = y
    print("Year: ", year)
    if year == 2019:
        Full_record = pd.read_csv("GVPDB2019.csv")
    else:
        print("Year not valid")
    GVP_df = Full_record
    for i in volcanoes:
        volc_num = i
        x += 1
        percentage_complete = (x / Total_volcanoes)
        Estimated_time_to_completion_min = (Total_volcanoes - x) * 3
        Estimated_time_to_completion_max = (Total_volcanoes - x) * 4
        print("Percengtage of analysis complete = ", percentage_complete)
        print("Estimation duration until completion: ", Estimated_time_to_completion_min, "-",
              Estimated_time_to_completion_max, " hours")
        print("Starting analysis for: ", volc_num)
        # confirmed eruptions
        for j in Eruption_certainty:
            confirmed = j
            print("Confirmed eruptions:", confirmed)
            for k in Change_point_method:
                method = k
                print("Using changepoint:", method)
                complete_record = get_complete_record(GVP_df, GVP_volcanoes, change_points, confirmed, method)
                #complete_record = pd.read_csv("complete_record_MeadMagill_50_confirmed.csv") # Add in the function
                Volcano_data = GVP_volcanoes[(GVP_volcanoes['Volcano Number'] == volc_num)]
                Volcano_name = Volcano_data.iloc[0]['Volcano Name']
                # # Analogue system 1
                Analogue_volcanoes = A1_volcanoes
                analogue_annual_frequency_1 = calc_mean_annual_frequency_analogue(Analogue_volcanoes, complete_record, confirmed, method, year)
                std_annual_frequency_analogue_1 = calc_std_annual_frequency_analogue(Analogue_volcanoes, complete_record, confirmed, method, year)
                analogue_freq_mag_1 = calc_relative_probability(Analogue_volcanoes, Full_record, confirmed, VEI_schema)
                average_eruptions_per_year_1, std_eruptions_per_year_1 = set_analogue_eruption_rate(volc_num, Analogue_volcanoes, Freq_rate=analogue_annual_frequency_1, Freq_rate_std=std_annual_frequency_analogue_1)
                Volcano_type_1 = Analogue_volcanoes.loc[Analogue_volcanoes['Volcano Number'] == volc_num, 'Volcano type'].iloc[0]  # A1_volcanoes.loc[A1_volcanoes['Volcano Number'] == volc_num, 'Volcano type'].iloc[0]
                Freq_rate_1 = analogue_annual_frequency_1
                Freq_rate_std_1 = std_annual_frequency_analogue_1

                # # Analogue system 2
                Analogue_volcanoes = A2_volcanoes
                analogue_annual_frequency_2 = calc_mean_annual_frequency_analogue(Analogue_volcanoes, complete_record, confirmed, method, year)
                std_annual_frequency_analogue_2 = calc_std_annual_frequency_analogue(Analogue_volcanoes, complete_record, confirmed, method, year)
                analogue_freq_mag_2 = calc_relative_probability(Analogue_volcanoes, Full_record, confirmed, VEI_schema)
                average_eruptions_per_year_2, std_eruptions_per_year_2 = set_analogue_eruption_rate(volc_num, Analogue_volcanoes, Freq_rate=analogue_annual_frequency_2, Freq_rate_std=std_annual_frequency_analogue_2)
                Volcano_type_2 = Analogue_volcanoes.loc[Analogue_volcanoes['Volcano Number'] == volc_num, 'Volcano type'].iloc[0]  # A1_volcanoes.loc[A1_volcanoes['Volcano Number'] == volc_num, 'Volcano type'].iloc[0]
                Freq_rate_2 = analogue_annual_frequency_2
                Freq_rate_std_2 = std_annual_frequency_analogue_2

                # Reported eruptions for volcano
                VEI_count = get_VEI_eruptions_count(volc_num, complete_record)
                Observed_rate = get_observed_eruption_rate(volc_num, complete_record, confirmed, year, Full_record, method)
                Observed_VEI = get_observed_relative_frequency(volc_num, complete_record, confirmed, Full_record, VEI_schema)
                # Observed_rate_full_record = get_observed_eruption_rate_full_record(volc_num, year, confirmed, Full_record, VEI_schema)
                # Observed_VEI_full_record = get_observed_relative_frequency_full_record(volc_num, Full_record, confirmed, VEI_schema)
                confirmed_VEI_eruptions = VEI_count
                Observed_rate = Observed_rate
                Observed_VEI = Observed_VEI


                for l in Power_law_used:
                    Power_law = l
                    print("Power law: ", Power_law)
                    Analogue_volcanoes = A1_volcanoes
                    freq_mag_dirichlet_1 = get_freq_mag_dirichlet(analogue_freq_mag=analogue_freq_mag_1, Analogue_volcanoes=Analogue_volcanoes, volc_num=volc_num, Power_law=Power_law, VEI_schema=VEI_schema)
                    VEI_freq_1 = freq_mag_dirichlet_1
                    Analogue_volcanoes = A2_volcanoes
                    freq_mag_dirichlet_2 = get_freq_mag_dirichlet(analogue_freq_mag=analogue_freq_mag_2, Analogue_volcanoes=Analogue_volcanoes, volc_num=volc_num, Power_law=Power_law, VEI_schema=VEI_schema)
                    VEI_freq_2 = freq_mag_dirichlet_2
                    if confirmed_VEI_eruptions > 1:
                        print("Using model averaging")
                        averaging_method = "Stacking"
                        print("averaging method: ", averaging_method)
                        Model_average_stacking = get_model_average(GVP_volcanoes,
                                                                   volc_num,
                                                                   Volcano_type_1,
                                                                   Volcano_type_2,
                                                                   Observed_rate,
                                                                   Observed_VEI,
                                                                   confirmed_VEI_eruptions,
                                                                   Freq_rate_1,
                                                                   Freq_rate_2,
                                                                   Freq_rate_std_1,
                                                                   Freq_rate_std_2,
                                                                   VEI_freq_1,
                                                                   VEI_freq_2,
                                                                   Power_law,
                                                                   averaging_method,
                                                                   year,
                                                                   confirmed,
                                                                   method,
                                                                   project_folder,
                                                                   VEI_schema,
                                                                   Uncertainty_schema)
                        averaging_method = "pseudo-BMA"
                        print("averaging method: ", averaging_method)
                        Model_average_pseudo_BMA = get_model_average(GVP_volcanoes,
                                                                     volc_num,
                                                                     Volcano_type_1,
                                                                     Volcano_type_2,
                                                                     Observed_rate,
                                                                     Observed_VEI,
                                                                     confirmed_VEI_eruptions,
                                                                     Freq_rate_1,
                                                                     Freq_rate_2,
                                                                     Freq_rate_std_1,
                                                                     Freq_rate_std_2,
                                                                     VEI_freq_1,
                                                                     VEI_freq_2,
                                                                     Power_law,
                                                                     averaging_method,
                                                                     year,
                                                                     confirmed,
                                                                     method,
                                                                     project_folder,
                                                                     VEI_schema,
                                                                     Uncertainty_schema)
                        averaging_method = "BB-pseudo-BMA"
                        print("averaging method: ", averaging_method)
                        Model_average_BB_pseudo_BMA = get_model_average(GVP_volcanoes,
                                                                        volc_num,
                                                                        Volcano_type_1,
                                                                        Volcano_type_2,
                                                                        Observed_rate,
                                                                        Observed_VEI,
                                                                        confirmed_VEI_eruptions,
                                                                        Freq_rate_1,
                                                                        Freq_rate_2,
                                                                        Freq_rate_std_1,
                                                                        Freq_rate_std_2,
                                                                        VEI_freq_1,
                                                                        VEI_freq_2,
                                                                        Power_law,
                                                                        averaging_method,
                                                                        year,
                                                                        confirmed,
                                                                        method,
                                                                        project_folder,
                                                                        VEI_schema,
                                                                        Uncertainty_schema)
                        averaging_method = "None"
                        analogue_model = "A1"
                        Bayesian_update_A1 = get_Bayes_update(GVP_volcanoes, volc_num, Volcano_type_1, Observed_rate,
                                                           Observed_VEI, confirmed_VEI_eruptions, Freq_rate_1,
                                                           Freq_rate_std_1, VEI_freq_1, Power_law, averaging_method,
                                                           year, confirmed, method, analogue_model)

                        analogue_model = "A2"
                        Bayesian_update_A2 = get_Bayes_update(GVP_volcanoes, volc_num, Volcano_type_2, Observed_rate,
                                                           Observed_VEI, confirmed_VEI_eruptions, Freq_rate_2,
                                                           Freq_rate_std_2, VEI_freq_2, Power_law, averaging_method,
                                                           year, confirmed, method, analogue_model)
                    else:
                        print("Not updated")
                        print("Using prior analogues")
                        averaging_method = "None"
                        print("averaging method: ", averaging_method)
                        A1_FM = get_A1_prior_model(GVP_volcanoes,
                                                   volc_num,
                                                   Volcano_type_1,
                                                   Freq_rate_1,
                                                   Freq_rate_std_1,
                                                   VEI_freq_1,
                                                   Power_law,
                                                   year,
                                                   confirmed,
                                                   method,
                                                   VEI_schema,
                                                   project_folder,
                                                   Uncertainty_schema)
                        A2_FM = get_A2_prior_model(GVP_volcanoes,
                                                   volc_num,
                                                   Volcano_type_2,
                                                   Freq_rate_2,
                                                   Freq_rate_std_2,
                                                   VEI_freq_2,
                                                   Power_law,
                                                   year,
                                                   confirmed,
                                                   method,
                                                   VEI_schema,
                                                   project_folder,
                                                   Uncertainty_schema)


path_csv = "Probabilities/" + project_folder + "/"
master_csv = get_master_csv(path_csv, save_file_name)
print("script finished in: ", datetime.now() - startTime)



