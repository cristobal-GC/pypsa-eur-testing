#!/usr/bin/env python
# -*- coding: utf-8 -*-

# SPDX-FileCopyrightText: : 2017-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""

"""

import pandas as pd
import geopandas as gpd
import pypsa
import yaml


from _helpers import configure_logging, set_scenario_config



##### Just during testing
# import pandas as pd
# import io                ########## Esto es solo para exportar fácilmente un string vía df


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("fit_elec_capacities")
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    

    ############################## Handle inputs and parameters

    ##### params
    update_with_REE = snakemake.params.update_with_REE
    dic_carriers_to_update = snakemake.params.carriers_to_update
    ##### inputs
    elec_initial = snakemake.input.elec_initial
    nuts2_ES_file = snakemake.input.nuts2_ES
    dic_nuts_file = snakemake.input.dic_nuts



    ############################## Load and define relevant variables:

    ##### initial net
    n = pypsa.Network(elec_initial)
    # related dfs
    df_buses = n.buses
    df_generators = n.generators


    ##### nuts2_ES
    nuts2_ES = gpd.read_file(nuts2_ES_file)


    ##### dic_nuts < is employed to loop over the regions considered in pypsa-es
    with open(dic_nuts_file, "r") as archivo:
        dic_nuts = yaml.safe_load(archivo)

    dic_nuts2 = dic_nuts['dic_NUTS2']








    ############################################################ Update not requested
    if not update_with_REE:

        print('########## Pypsa-ES [rule: update_elec_capacities]: No update requested')








    ############################################################ Update requested
    if update_with_REE:

        print('########## Pypsa-ES [rule: update_elec_capacities]: Updating generation capacities according to REE...')


        

        #################### Loop over carriers to update
        for cc, ff in dic_carriers_to_update.items():


            ##### Load file with REE data
            df_installed_capacity = pd.read_csv(ff, index_col='datetime')






            #################### Loop over NUTS2 regions included in df_installed_capacity columns
            for rr in df_installed_capacity.columns:



                ##### check that rr is in dic_NUTS2 (to avoid 'total', 'Melilla',... that ar in rr, and 'Canarias' that is not in dic_nuts)
                if rr in dic_nuts2.values():


                    rr_name = nuts2_ES.loc[ nuts2_ES["NUTS_ID"]==rr , "NUTS_NAME"].values[0]

                    
                    ########## get a list with local buses located in region rr
                    geometry_buses = gpd.points_from_xy(df_buses["x"], df_buses["y"])
                    gdf_buses = gpd.GeoDataFrame(df_buses,geometry=geometry_buses, crs=4326)

                    gdf_region = nuts2_ES[nuts2_ES['NUTS_ID']==rr]


                    ### Get intersection
                    gdf_buses_local = gpd.sjoin( gdf_buses , gdf_region , how="inner", predicate="within")

                    ### Get the list of local buses
                    list_buses_local = gdf_buses_local.index.to_list()



                    ########## Get the generators in local buses and carrier cc, with capacity>0.01 (to avoid everywhere generators). Compute total installed capacity
                    df_generators_local_cc = df_generators.loc[ (df_generators['bus'].isin(list_buses_local)) & (df_generators['carrier']==cc) & (df_generators['p_nom']>0.01)]

                    initial_installed_capacity = df_generators_local_cc['p_nom'].sum()



                    ########## Get the real capacity reported by ESIOS in that region and carrier for the desired year        

                    ### Check if region rr has installed capacity of carrier cc, otherwise assign 0
                    required_installed_capacity = df_installed_capacity[rr].mean()



                    ############### Modify the network:

                    ### If capacity needs to be increased: add same capacity to each bus
                    if int(required_installed_capacity) > int(initial_installed_capacity):

                        capacity_added_at_each_bus = (required_installed_capacity-initial_installed_capacity)/df_generators_local_cc.shape[0]


                        print(f'########## Pypsa-ES [rule "update_elec_capacities"]: {cc} capacity in {rr_name} was increased from {initial_installed_capacity:.2f} MW to {required_installed_capacity:.2f} MW ({capacity_added_at_each_bus:.2f} MW were added to each of the {df_generators_local_cc.shape[0]}/{len(list_buses_local)} buses).')

                        df_generators_local_cc.loc[:, 'p_nom'] += capacity_added_at_each_bus


                    ### If capacity needs to be reduced: reduce capacity proportionally at each bus
                    if int(required_installed_capacity) < int(initial_installed_capacity):

                        factor = required_installed_capacity / initial_installed_capacity

                        print(f'########## Pypsa-ES [rule "update_elec_capacities"]: {cc} capacity in {rr_name} was reduced from {initial_installed_capacity:.2f} MW to {required_installed_capacity:.2f} MW (a factor of {factor:.2f} is applied to each of the {df_generators_local_cc.shape[0]}/{len(list_buses_local)} buses).')

                        df_generators_local_cc.loc[:, 'p_nom'] *= factor


                    ### If capacity matches, do nothing
                    if int(required_installed_capacity) == int(initial_installed_capacity):
                        print(f'########## Pypsa-ES [rule "update_elec_capacities"]: {cc} capacity perfectly matches in {rr_name}! ({df_generators_local_cc.shape[0]}/{len(list_buses_local)} buses).')




                    ### Integrate the modifications
                    n.generators.update(df_generators_local_cc)





    ############################################################ Save output net

    n.export_to_netcdf(snakemake.output.elec)
    

