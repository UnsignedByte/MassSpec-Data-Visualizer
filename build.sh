# @Author: UnsignedByte
# @Date:   15:02:23, 09-Jun-2020
# @Last Modified by:   UnsignedByte
# @Last Modified time: 17:41:00, 22-Oct-2020

# create venv if not exist
if [ ! -d "$.venv" ]; then
  virtualenv --python=python3.7 .venv
fi


# activate venv
source .venv/bin/activate
pip install -r requirements.txt

# add matlab alias
alias matlab='/Applications/MATLAB_R2019b.app/bin/matlab -nodesktop -nosplash $*'

# declare -a langs=("matlab" "r")
# for i in "${langs[@]}"
# do
# mkdir -p "${PWD}/utils/cpptemplates/${i}"
# echo "#define ${i}" > "${PWD}/utils/cpptemplates/${i}/t.h";
# done