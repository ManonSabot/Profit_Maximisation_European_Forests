try:
    from general_utils import get_script_dir, get_main_dir
    from general_utils import retrieve_class, read_csv, read_netcdf
    from calculate_solar_geometry import cos_zenith
    from cru_climate_lat_lon import main as cru_climate

except (ImportError, ModuleNotFoundError):
    from TractLSM.Utils.general_utils import get_script_dir, get_main_dir
    from TractLSM.Utils.general_utils import retrieve_class
    from TractLSM.Utils.general_utils import read_csv, read_netcdf
    from TractLSM.Utils.calculate_solar_geometry import cos_zenith
    from TractLSM.Utils.cru_climate_lat_lon import main as cru_climate
