no_conda_msg1="Anaconda is not installed. Install it and try again!"
no_conda_msg2="To install Anaconda, head to: https://docs.anaconda.com/anaconda/install/"

all: check_conda create_env

check_conda:
	@conda info --envs || (echo ${no_conda_msg1}; echo ${no_conda_msg2}; exit)

create_env:
	@conda env create -f src/extra/envs/TractLSM.yml; conda env create -f src/extra/envs/TractLSM37.yml; conda env create -f src/extra/envs/TractLSM27.yml || echo "please, activate either TractLSM or TractLSM37 or TractLSM27"
