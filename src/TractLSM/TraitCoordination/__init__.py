try:
    from optimise_kmax import optimal_kmax

except (ImportError, ModuleNotFoundError):
    from TractLSM.TraitCoordination.optimise_kmax import optimal_kmax
