# @Author: UnsignedByte
# @Date:   15:02:23, 09-Jun-2020
# @Last Modified by:   UnsignedByte
# @Last Modified time: 13:36:06, 10-Jan-2021

# create venv if not exist
if [ ! -d "./.venv" ]; then
  virtualenv --python=python3 --system-site-packages .venv
fi


# activate venv
source .venv/bin/activate
pip install -r requirements.txt

# declare -a langs=("matlab" "r")
# for i in "${langs[@]}"
# do
# mkdir -p "${PWD}/utils/cpptemplates/${i}"
# echo "#define ${i}" > "${PWD}/utils/cpptemplates/${i}/t.h";
# done