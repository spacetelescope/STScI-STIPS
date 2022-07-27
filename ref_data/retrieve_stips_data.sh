### Instructions ###
# 1. In your ~/.bash_profile file, copy the following environment variables into
#    the same file (making sure to update 1.x.x version of Pandeia you're using
#    –– for STIPS this should be 1.7):
#       export stips_data="<absolute_path_to_this_folder>/ref_data/stips_data"
#       export WEBBPSF_PATH="<absolute_path_to_this_folder>/ref_data/webbpsf-data"
#       export PYSYN_CDBS="<absolute_path_to_this_folder>/ref_data/grp/redcat/trds"
#       export pandeia_refdata="<absolute_path_to_this_folder>/ref_data/pandeia_data-
#       1.7_roman"
# 2. In the ref_data folder, run the script: bash retrieve_stips_data.sh

### STIPS Supporting Data Files ###
curl -OL https://stsci.box.com/shared/static/4nebx2ndxr7c77lgocfbvxo7c2hyd3in.tgz
tar -xzvf 4nebx2ndxr7c77lgocfbvxo7c2hyd3in.tgz
rm 4nebx2ndxr7c77lgocfbvxo7c2hyd3in.tgz

### Synphot and STSynphot Reference Data ###
curl -o synphot1.tar.gz https://ssb.stsci.edu/trds/tarfiles/synphot1.tar.gz
curl -o synphot2.tar.gz https://ssb.stsci.edu/trds/tarfiles/synphot2.tar.gz
curl -o synphot3.tar.gz https://ssb.stsci.edu/trds/tarfiles/synphot3.tar.gz
curl -o synphot4.tar.gz https://ssb.stsci.edu/trds/tarfiles/synphot4.tar.gz
curl -o synphot5.tar.gz https://ssb.stsci.edu/trds/tarfiles/synphot5.tar.gz
curl -o synphot6.tar.gz https://ssb.stsci.edu/trds/tarfiles/synphot6.tar.gz
curl -o synphot7.tar.gz https://ssb.stsci.edu/trds/tarfiles/synphot7.tar.gz

tar -xzvf synphot1.tar.gz
tar -xzvf synphot2.tar.gz
tar -xzvf synphot3.tar.gz
tar -xzvf synphot4.tar.gz
tar -xzvf synphot5.tar.gz
tar -xzvf synphot6.tar.gz
tar -xzvf synphot7.tar.gz

rm synphot1.tar.gz
rm synphot2.tar.gz
rm synphot3.tar.gz
rm synphot4.tar.gz
rm synphot5.tar.gz
rm synphot6.tar.gz
rm synphot7.tar.gz

### Pandeia Reference Data ###
curl -OL https://stsci.box.com/shared/static/ycbm34uxhzafgb7te74vyl2emnr1mdty.gz
tar -xvzf ycbm34uxhzafgb7te74vyl2emnr1mdty.gz
rm ycbm34uxhzafgb7te74vyl2emnr1mdty.gz

### WebbPSF Supporting Data Files ###
curl -OL https://stsci.box.com/shared/static/34o0keicz2iujyilg4uz617va46ks6u9.gz
tar -xvzf 34o0keicz2iujyilg4uz617va46ks6u9.gz
rm 34o0keicz2iujyilg4uz617va46ks6u9.gz
