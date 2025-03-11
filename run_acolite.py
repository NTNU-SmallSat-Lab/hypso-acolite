import numpy as np
from importlib.resources import files
from urllib.parse import urlparse, unquote

from pathlib import Path
import netCDF4 as nc

import sys
import os

from hypso import Hypso
from hypso.write import write_l1c_nc_file


def main(l1a_nc_file_path: str, acolite_path: str) -> np.ndarray:
    """
    Run the ACOLITE correction model. Adjustments of the original files are made to ensure they work for HYPSO

    :param hypso_info: Dictionary containing the hypso capture information
    :param atmos_dict: Dictionary containing the information required for the atmospheric correction
    :param nc_file_acoliteready: Absolute path for the .nc file for ACOLITE (L1B)

    :return: Returns surface reflectanc corrected spectral image
    """
    # https://odnature.naturalsciences.be/remsem/acolite-forum/viewtopic.php?t=238
    
    # add acolite clone to Python path and import acolite

    acolite_path = Path(acolite_path)
    acolite_path = os.path.join(acolite_path, "acolite")
    acolite_path = Path(acolite_path).absolute()

    #print(sys.path)
    sys.path.append(str(acolite_path))
    #print(sys.path)

    import acolite as ac


    # Check if the first file exists
    if not os.path.isfile(l1a_nc_file_path):
        print(f"Error: The file '{l1a_nc_file_path}' does not exist.")
        return

    # Process the first file
    print(f"Processing file: {l1a_nc_file_path}")

    l1a_nc_file_path = Path(l1a_nc_file_path)

    satobj = Hypso(path=l1a_nc_file_path, verbose=True)

    satobj.generate_l1b_cube()
    satobj.generate_l1c_cube()

    write_l1c_nc_file(satobj, overwrite=True, datacube=True)


    l1c_nc_file_path = satobj.l1c_nc_file

    output_path = satobj.parent_dir



    # optional file with processing settings
    # if set to None defaults will be used
    settings_file = None

    # import settings
    settings = ac.acolite.settings.load(settings_file)

    # set settings provided above
    settings['inputfile'] = str(l1c_nc_file_path)
    settings['output'] = str(output_path)

    settings['polygon'] = None
    settings['l2w_parameters'] = None
    settings['rgb_rhot'] = True
    settings['rgb_rhos'] = True
    settings['map_l2w'] = False

    # user and password from https://urs.earthdata.nasa.gov/profile
    # optional but good
    #settings['EARTHDATA_u'] = "alvarof"
    #settings['EARTHDATA_p'] = "nwz7xmu8dak.UDG9kqz"

    # other settings can also be provided here, e.g.
    # settings['s2_target_res'] = 60
    # settings['dsf_path_reflectance'] = 'fixed'
    # settings['l2w_parameters'] = ['t_nechad', 't_dogliotti']

    # process the current bundle
    processed = ac.acolite.acolite_run(settings=settings)

    acolite_l2_file = processed[0]['l2r'][0]

    # Maintainer comment:
    # Source: https://odnature.naturalsciences.be/remsem/acolite-forum/viewtopic.php?t=311
    # - L1R, containing the top-of-atmosphere reflectance (rhot_*) as converted to the ACOLITE format from the sensor specific L1 files
    # - L2R, containing the top-of-atmosphere reflectance (rhot_*) and the surface-level reflectance after atmospheric correction (rhos_*)
    # - L2W, containing user requested parameters, e.g. water-leaving radiance reflectance (rhow_*), Remote sensing reflectance (Rrs_*), or outputs from any of the included parameter retrieval algorithms

    # Read .nc
    '''
    final_acolite_l2 = None

    with nc.Dataset(acolite_l2_file, format="NETCDF4") as f:
        group = f
        keys = [i for i in f.variables.keys()]

        toa_keys = [k for k in keys if 'rhos' not in k]
        surface_keys = [kk for kk in keys if 'rhot' not in kk]

        # Add Cube

        for i, k in enumerate(surface_keys):
            current_channel = np.array(group.variables[k][:])
            if final_acolite_l2 is None:
                final_acolite_l2 = np.empty(
                    (current_channel.shape[0], current_channel.shape[1], len(surface_keys)))

            final_acolite_l2[:, :, i] = current_channel

        # TODO: Confirm if zeros should be appended at the beginning or end
        # ACOLITE returns 118 bands
        # If number of bands less that 120, append zeros to the end
        delta = int(120 - final_acolite_l2.shape[2])
        if delta > 0:
            for _ in range(delta):
                zeros_arr = np.zeros((final_acolite_l2.shape[0], final_acolite_l2.shape[1]), dtype=float)
                final_acolite_l2 = np.dstack((final_acolite_l2, zeros_arr))

    return final_acolite_l2
    '''

if __name__ == "__main__":
    if len(sys.argv) < 3 or len(sys.argv) > 3:
        print("Usage: python script.py <path_to_nc_file> <path_to_acolite>")
        sys.exit(1)

    l1a_nc_file_path = sys.argv[1]
    acolite_path = sys.argv[2]

    main(l1a_nc_file_path, acolite_path)
