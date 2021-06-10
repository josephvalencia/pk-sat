# Installing Z3

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
## Installing LinearPartition

From the /pk-sat directory do:

```
cd ..
git clone https://github.com/LinearFold/LinearPartition.git
cd LinearPartition/
make
```
To verify installation run:

```
cat ../pk-sat/CRW_00774.fa | ./linearpartition --prefix ./prob_matrix
```
You should see output in the text file prob_matrix_1

## Testing Baseline
To test the baseline model (Ganesh et al 2012), run
```
python pksat.py
```
## Testing Advanced
To test the advanced model (Sato et al 2011), run:
```
python ipknot.py
```
