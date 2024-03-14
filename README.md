# QCML
Quantum Computing and Machine Learning


# IBM Account
visit https://quantum.ibm.com to open an account in IBM.
You'll need the token from IBM homepage to run the examples in IBM.
For this repository the toten is assuemd to be stored in "~/Documents/keys/ibm-token.txt" file.


# Install conda 

# Qiskit 1.0 Configuration
$ conda create --name quantum 

$ conda activate quantum 

$ conda install jupyter python=3.8

$ pip install qiskit  qiskit-ibm-runtime qiskit-aer

### In MacOS
$ pip install 'qiskit[visualization]'

### In Windows/Linux
$ pip install qiskit[visualization]

# Nvidia Tensor Core to simulate
$ conda install conda-forge::cuquantum-python


# Removing conda environment
If you want to remove the quantum environment 

$ conda remove --name quantum --all





