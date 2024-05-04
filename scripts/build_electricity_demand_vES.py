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
    df_profiles = pd.read_csv(snakemake.params.profiles, index_col=0)
    df_profiles.fillna(0, inplace=True)

    ##### Load percentages
    df_percentages = pd.read_csv(snakemake.params.percentages, index_col=0)




    #################### Initialise output
    
    df_output = pd.DataFrame()





    #################### Load dic with communities or provinces if required

    ##### If communities
    if 'Galicia' in df_percentages.columns:
        with open(snakemake.input.dic_communities, 'r') as file:
            dic_regions = yaml.safe_load(file)




    ############### Loop to combine profiles, percentages and annual demand

    ##### Check if dic_regions exists (use the dic name as string, otherwise it does not work)
    if 'dic_regions' in locals():



        ##### Loop over rr and ff, multiply each profile by corresponding factor

        for rr in dic_regions.values():



            for ff_values in df_percentages.index:


                ### *1e6 because annaul_electricity_demand is provided in TWh, but the time series is in MWh
                ### /8760 because hourly profiles have been obtained with the mean load, not the sum
                ### /100 because percentages are over 100
                factor = (annual_electricity_demand * 1e6 / 8760 )* df_percentages.at[ff_values,rr] / 100

                df_profiles.loc[:,f'{rr}-{ff_values}'] = df_profiles[f'{rr}-{ff_values}'] * factor



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
