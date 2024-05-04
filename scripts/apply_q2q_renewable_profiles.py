#!/usr/bin/env python
# -*- coding: utf-8 -*-

# SPDX-FileCopyrightText: : 2017-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""

"""


import logging
import pickle
import xarray as xr 


from _helpers import configure_logging, set_scenario_config



##### Just during testing
# import pandas as pd
# import io                ########## Esto es solo para exportar fácilmente un string vía df


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("apply_q2q_renewable_profiles", technology="onwind")
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    

    #################### Allocate relevant variables
    ### parameters
    q2q_transforms = snakemake.params.q2q_transform
    ### inputs
    profile_initial = snakemake.input.profile_initial
    ### outputs
    profile = snakemake.output.profile
    ### wildcards
    technology = snakemake.wildcards.technology



    #################### Load initial profile
    profile_data = xr.open_dataset(profile_initial)



    #################### Retrieve the name of the file with the q2q transform from config.yaml
    name_of_function=q2q_transforms[technology]



    ##### If name_of_function is empty, do nothing, just generate the output
    if not name_of_function:
        print(f'########## pypsa-es INFO: No q2q transform available for {technology}..')
        profile_data.to_netcdf(profile)


    ##### If name_of_function exists, apply q2q trasform and save output
    else:
        print(f'########## pypsa-es INFO: Applying q2q transform {name_of_function} to {technology}..')
        with open(name_of_function, 'rb') as f:
            interp_func = pickle.load(f)

            profile_data["profile"][:, :, :] = interp_func(profile_data.variables["profile"])

            profile_data.to_netcdf(profile)
    

