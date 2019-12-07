try:
    from USO import solve_std, set_trans_std
    from ProfitMax import profit_psi

except (ImportError, ModuleNotFoundError):
    from TractLSM.CH2OCoupler.USO import solve_std, set_trans_std
    from TractLSM.CH2OCoupler.ProfitMax import profit_psi
