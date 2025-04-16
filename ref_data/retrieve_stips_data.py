### Instructions ###
# 1. In your ~/.bash_profile file, copy the following environment variables into
#    the same file:
#       export stips_data="<absolute_path_to_this_folder>/ref_data/stips_data"
#       export STPSF_PATH="<absolute_path_to_this_folder>/ref_data/stpsf-data"
#       export PYSYN_CDBS="<absolute_path_to_this_folder>/ref_data/grp/redcat/trds"
#       export pandeia_refdata="<absolute_path_to_this_folder>/ref_data/pandeia_data-
#       2024.12-roman"
#
# 2. In the ref_data folder, run the script: python retrieve_stips_data.py

import stips

stips.DownloadReferenceData()
