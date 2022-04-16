import MDSplus as mds
import numpy as np


def get_outergap_efit(shot_number, tree="ANALYSIS"):
    """
    Extracts the data on the outergap determined by EFIT.
    The outergap is the distance of limiter surface relative to the separatrix.
    More details about EFIT are here: https://omfit.io/modules/mod_EFIT.html

    Args:
        shot_number: Shot number(s) of interest.
        tree: Which EFIT data you want.
              By default, this is set to "ANALYSIS" which is a lower resolution EFIT.
              The other option is "EFIT19" which is a higher resolution EFIT.
              Not all shots have EFIT19 data!
              This is case sensisitve, so use capitals.

    Returns:
        time_outergap_EFIT: Time data of the EFIT data in seconds.
        outergap_EFIT: Outergap determined by EFIT in millimetres.

    """

    c = mds.Connection("alcdata")

    if tree == "EFIT19":
        c.openTree("efit19", shot_number)
        tree_path = "\EFIT19::TOP.RESULTS.A_EQDSK:ORIGHT"
    else:
        c.openTree("analysis", shot_number)
        tree_path = "\ANALYSIS::EFIT_AEQDSK:ORIGHT"

    outergap_efit = c.get(tree_path)
    time_outergap_efit = c.get(f"dim_of({tree_path})").data()

    return time_outergap_efit, outergap_efit
