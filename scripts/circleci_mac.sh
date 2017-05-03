export UVCDAT_ANONYMOUS_LOG=False
export PATH=${HOME}/miniconda/bin:${PATH}
python run_tests.py -v2 
rc=$?
if [ ${rc} -ne 0 ]; then exit ${rc} ; fi
if [ ${rc} -eq 0 -a $CIRCLE_BRANCH == "master" ]; then bash ./scripts/conda_upload.sh ; fi
