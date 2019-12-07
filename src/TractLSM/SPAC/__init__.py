try:
    from hydraulics import f, Weibull_params, k_regulate, hydraulics
    from canatm import vpsat, slope_vpsat, LH_water_vapour, psychometric
    from canatm import emissivity, net_radiation, absorbed_radiation_2_leaves
    from leaf import conductances, leaf_temperature, leaf_energy_balance
    from leaf import foliar_resp, calc_photosynthesis, rubisco_limit
    from soil import fwsoil, wetness, water_potential

except (ImportError, ModuleNotFoundError):
    from TractLSM.SPAC.hydraulics import f, Weibull_params, k_regulate
    from TractLSM.SPAC.hydraulics import hydraulics
    from TractLSM.SPAC.canatm import vpsat, slope_vpsat, LH_water_vapour
    from TractLSM.SPAC.canatm import psychometric, emissivity, net_radiation
    from TractLSM.SPAC.canatm import absorbed_radiation_2_leaves
    from TractLSM.SPAC.leaf import conductances, leaf_temperature
    from TractLSM.SPAC.leaf import leaf_energy_balance, foliar_resp
    from TractLSM.SPAC.leaf import calc_photosynthesis, rubisco_limit
    from TractLSM.SPAC.soil import fwsoil, wetness, water_potential
