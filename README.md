# Z3 Installation

```
git clone https://github.com/Z3Prover/z3.git
cd pk-sat
python3 -m venv venv
source venv/bin/activate
cd ../z3
python scripts/mk_make.py --python
cd build
make
make install
# You will find Z3 and the Python bindings installed in the virtual environment
venv/bin/z3 -h
```
## LinearPartition Installation

From the /SAT directory do:

```
cd ..
git clone https://github.com/LinearFold/LinearPartition.git
cd LinearPartition/
make
cat ../SAT/CRW_00774.fa | ./linearpartition --prefix ./prob_matrix
```

