#!/usr/bin/env python
# -*- coding: utf-8 -*-

# SPDX-FileCopyrightText: : 2017-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""

"""


import pandas as pd
import logging
import yaml


from _helpers import configure_logging, set_scenario_config, get_snapshots



logger = logging.getLogger(__name__)





if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("build_electricity_demand")

    configure_logging(snakemake)
    set_scenario_config(snakemake)



    #################### Handle inputs and parameters

    ##### Unwrap annual electricity demand
    annual_electricity_demand = snakemake.params.annual_electricity_demand

    ##### Load profiles
    df_profiles = pd.read_csv(snakemake.params.profiles, index_col=0)#.fillna(0, inplace=True)


    ##### Load percentages
    df_percentages = pd.read_csv(snakemake.params.percentages, index_col=0)


    ##### Load dic_datadis, that contains relevant info
    with open(snakemake.input.dic_datadis, 'r') as file:
        dic_datadis = yaml.safe_load(file)





    #################### Initialise output
    
    df_output = pd.DataFrame()





##### REMOVE porque los códigos NUTS ya van en el nombre de las columnas
    #################### Load dic with communities, provinces or country, according to the number of columns in df_percentages
    #if df_percentages.shape[1] == 16:
    #    dic_region = dic_datadis['dic_community_NUTSid']
    #if df_percentages.shape[1] == 48:
    #    dic_region = dic_datadis['dic_province_NUTSid']
    #if df_percentages.shape[1] == 1:
    #    dic_region = dic_datadis['dic_province_NUTSid']
   



    ############### Loop to combine profiles, percentages and annual demand



    ##### Loop over rr and ff, multiply each profile by corresponding factor

    for rr in df_percentages.columns:  # rr are the NUTS_ID


        for ff in df_percentages.index:


            ### *1e6 because annaul_electricity_demand is provided in TWh, but the time series is in MWh
            ### /8760 because hourly profiles have been obtained with the mean load, not the sum
            ### /100 because percentages are over 100
            factor = (annual_electricity_demand * 1e6 / 8760 )* df_percentages.at[ff,rr] / 100

            df_profiles.loc[:,f'{rr}-{ff}'] = df_profiles[f'{rr}-{ff}'] * factor



        ### Aggregate df_profiles according to rr
        df_output[rr] = df_profiles.filter(like=rr).sum(axis=1)






    snapshots = get_snapshots(
        snakemake.params.snapshots, snakemake.params.drop_leap_day
    )


    df_output.index=snapshots



    #assert not load.isna().any().any(), (
    #    "Load data contains nans. Adjust the parameters "
    #    "`time_shift_for_large_gaps` or modify the `manual_adjustment` function "
    #    "for implementing the needed load data modifications."
    #)



    df_output.round(4).to_csv(snakemake.output[0], float_format='%.4f')



