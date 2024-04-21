###################################

# This quick script will run the toy example provided by Qiskit using the quantum protein folding
# algorithm proposed in the article 'Resource-efficient quantum algorithm for protein folding' [1]. 
# This algorithm from my understanding is for determining/simulating the tertiary structure, as in
# a single chain protein (within the tertirary structure monomers are amino acids; within the quaternary
# structure the monomers would be the subunits, aka the invidual proteins that form the complex). 
#   ~ Andrew Jordan Siciliano

###################################

# git clone https://github.com/qiskit-research/qiskit-research.git
# cd qiskit-research
# pip install .

# git clone https://github.com/qiskit-community/quantum-protein-folding

###################################

# Note this code is fetched from the qiskit-research repository on protein folding!

# References: 
#   [1] https://www.nature.com/articles/s41534-021-00368-4
#   [2] https://www.sciencedirect.com/science/article/abs/pii/S002228369690114X?via%3Dihub
#   [3] https://doi.org/10.1016/j.sbi.2006.06.013
#   [4] https://www.nature.com/articles/s41598-020-73300-z
#   [5] https://www.sciencedirect.com/science/article/abs/pii/B9780080951676001099
#   [6] https://arxiv.org/pdf/2312.00875.pdf
#   [Code] https://github.com/qiskit-community/quantum-protein-folding/blob/main/docs/protein_folding.ipynb
#   [Source] https://qiskit-community.github.io/qiskit-research/protein_folding/protein_folding.html

# [1] Futher Resources:
#   > Hamiltonian Operator Creation See: 
#       > https://github.com/qiskit-community/quantum-protein-folding/blob/main/src/protein_folding/qubit_op_builder.py
#           > For backbone interaction: _create_h_bbbb
#   > https://www.mathworks.com/help/matlab/math/ground-state-protein-folding-using-variational-quantum-eigensolver-vqe.html

# Resources on Tetrahedral (Diamond) Lattice Structures: 
#   > https://edu.itp.phys.ethz.ch/15FS/sst/solution2.pdf
#   > https://unlcms.unl.edu/cas/physics/tsymbal/teaching/SSP-927/Section%2001_Crystal%20Structure.pdf

# Resources for using real quantum device w/IBM: 
#   > https://medium.com/@FelA350/how-to-run-a-quantum-computer-with-qiskit-d539d1f2b8e5

# Hamiltonians and interactions: https://compphysics.github.io/ComputationalPhysics2/doc/LectureNotes/_build/html/basicmanybody.html

# HUBO (Higher-Order Unconstrained Binary Optimization) & Qiskit Pauli Ops:
#   > https://quantumcomputing.stackexchange.com/questions/32288/does-qaoa-require-that-the-problem-hamiltonian-be-an-ising-hamiltonian-as-a-quad
#   > https://quantumcomputing.stackexchange.com/questions/16020/solving-higher-order-unconstrained-binary-optimization-problems-with-qaoa-with
#   > https://arxiv.org/pdf/2009.07309.pdf
#       - this article provides ways of constructing hamiltonians for HUBO!
#   > D-Wave requires QUBO format b/c of quantum hardware limitatations
#   > VQE & QOAO can I think work with higher order polynomials, so as long as you can properly 
#     convert your HUBO into a hamiltonian (not trivial...). Another option is using quadratization but 
#     requires additional qubits bits...
#   > https://docs.quantum.ibm.com/api/qiskit/qiskit.quantum_info.Pauli

# Peptide Database:
#   > https://www.iedb.org/result_v3.php?cookie_id=ee5fc5
#       > https://www.iedb.org/epitope/155

###################################

import sys

sys.path.append("./quantum-protein-folding/src/")

from protein_folding.interactions.miyazawa_jernigan_interaction import MiyazawaJerniganInteraction

from protein_folding.peptide.peptide import Peptide
from protein_folding.protein_folding_problem import ProteinFoldingProblem

from protein_folding.penalty_parameters import PenaltyParameters

from qiskit.utils import algorithm_globals, QuantumInstance

algorithm_globals.random_seed = 23

###################################

import argparse

parser = argparse.ArgumentParser()

parser.add_argument("-seq","--protein_sequence",action="store",type=str,default="APRLRFY", help="protein sequence")
parser.add_argument("-pb","--penalty_back",action="store",type=int,default=10)
parser.add_argument("-po","--penalty_overlap",action="store",type=int,default=10)
parser.add_argument("-pc","--penalty_chiral",action="store",type=int,default=10)
parser.add_argument("-i","--iterations",action="store",type=int,default=100)
parser.add_argument("-o","--out",action="store",type=str,default="./results/")

args = parser.parse_args()

###################################

N = len(args.protein_sequence) # protein sequence length
side_chains = [""] * N

assert N > 3, "Too small of a protein!"

C = 2 * (len(args.protein_sequence) - 3)  # Configuration Qubits

###################################

print("*"*15)
print("'Resource-efficient quantum algorithm for protein folding' Qiskit Implementation")
print("*"*15)
print("sequence:",args.protein_sequence,"->",N,"residues")
print("Chirality Penalty:",args.penalty_chiral)
print("Back (turns on same axis) Penalty :",args.penalty_back)
print("Overlap Penalty:",args.penalty_overlap)
print("Iterations:",args.iterations)
print("Outdir:",args.out)

###################################

import os

save_dir = args.out + "/" + args.protein_sequence + "/"

if not os.path.exists(save_dir): os.makedirs(save_dir)

print(vars(args))    

###################################    

mj_interaction = MiyazawaJerniganInteraction() # [2]

# "Knowledge-based potentials are statistical parameters derived from databases of known protein properties that empirically
# capture aspects of the physical chemistry of protein structure and function." [3]

###################################

penalty_terms = PenaltyParameters(args.penalty_chiral, args.penalty_back, args.penalty_overlap)

# "Except for glycine, the amino acids exist as chiral molecules in two forms, the l- and d-enantiomers." [4]
# "Chirality plays an important role in the recognition phenomenon between the biologically active 
#  molecule and its target;" [5]
# further of chirality -> https://pl.khanacademy.org/test-prep/mcat/chemical-processes/stereochemistry/a/chiral-drugs

###################################

peptide = Peptide(args.protein_sequence, side_chains)

# a peptide is a chain of amino acids that is smaller in length
# a protein (not a complex... single chain) is a medium-long length chain of amino acids

###################################

protein_folding_problem = ProteinFoldingProblem(peptide, mj_interaction, penalty_terms)
qubit_op = protein_folding_problem.qubit_op() # hamiltonian given by a 'PauliSumOp'

# note that this operator only uses identity and pauli Z operators, so it is equivilant to a binary function by linear transformation!

print("*"*15)
print("Total Qubits:",qubit_op.num_qubits)
print("Configuration Qubits:",C) # this I am unsure of at the moment
print("Interaction Qubits:",qubit_op.num_qubits - C) # this I am unsure of at the moment

###################################

from qiskit.circuit.library import RealAmplitudes
from qiskit.algorithms.optimizers import COBYLA
from qiskit.algorithms import NumPyMinimumEigensolver
from qiskit.algorithms.minimum_eigensolvers import SamplingVQE
from qiskit import execute, Aer
from qiskit.primitives import Sampler

optimizer = COBYLA(maxiter=args.iterations)
ansatz = RealAmplitudes(reps=1)

###################################

import json
import numpy as np

if not os.path.exists(save_dir + "iterations/"): os.makedirs(save_dir + "iterations/")

def store_intermediate_result(eval_count, parameters, mean, std):
    global storage
    
    # note this will overwrite existing iterations!
    
    with open(save_dir + "iterations/" + str(eval_count) + ".json", 'w') as outfile: 
        json.dump({
            "iter": int(eval_count),
            "mean-eng": np.real(mean).item()
        }, outfile)
        
vqe = SamplingVQE(
    Sampler(),
    ansatz=ansatz,
    optimizer=optimizer,
    aggregation=0.1, # CVaR: alpha = 0.1
    callback=store_intermediate_result,
)

###################################

raw_result = vqe.compute_minimum_eigenvalue(qubit_op)
result = protein_folding_problem.interpret(raw_result=raw_result)

print("*"*15)
print("The bitstring representing the shape of the protein during optimization is: ",result.turn_sequence)
print("The expanded expression is:", result.get_result_binary_vector())
print(f"The folded protein's main sequence of turns is: {result.protein_shape_decoder.main_turns}")
print("*"*15)

final_result = {
    
}

for i, entry in enumerate(result.protein_shape_file_gen.get_xyz_data()):
    assert entry[0] == args.protein_sequence[i], "ERROR: output prot sequence at position " +\
        str(i) + "not matching input prot sequence > " + args.sequence[i] + "," + entry[0]
    print("Residue:",i+1,"| Amino Acid:",entry[0],"| (X,Y,Z):",tuple(entry[1:]))
    final_result[i+1] = {
        "AA": entry[0],
        "coordinate": tuple(np.real(entry[1:]))
    }
    
with open(save_dir + "protein.json", 'w') as outfile: json.dump(final_result, outfile)
print("Saved Protein Coordinates @:",save_dir + "protein.json")

###################################

energies = {}

for fl in os.listdir(save_dir + "iterations/"): 
    metrics_iter = json.load(open(save_dir + "iterations/" + fl,'r'))
    energies[metrics_iter["iter"]] = metrics_iter["mean-eng"]
    
iterations = list(sorted(energies.keys()))
energies = [energies[i] for i in iterations]

###################################

import matplotlib.pyplot as plt

fig = plt.figure()

plt.plot(iterations, energies)
plt.ylabel("Conformation Energy")
plt.xlabel("VQE Iterations")

plt.savefig(save_dir + "curve.png")

plt.clf()

fig = result.get_figure(title="Protein Structure", ticks=False, grid=True)
fig.get_axes()[0].view_init(10, 70)

plt.savefig(save_dir + "protein.png")

print("*"*15)
print("Saved Curve Image @:",save_dir + "curve.png")
print("Saved Protein Image @:",save_dir + "protein.png")

###################################

print("*"*15)
print("Done!")
print("*"*15)

###################################










