#! /bin/bash

cd trioasm/whatshap
$PYTHON -m pip install . --no-deps --ignore-installed --no-cache-dir -vvv
#python setup.py build_ext -i
cd ../whatshap_trioasm
$PYTHON -m pip install . --no-deps --ignore-installed --no-cache-dir -vvv
#python setup.py build_ext -i
cd ../../

mkdir -p $PREFIX/bin/
cp $RECIPE_DIR/Aligner $PREFIX/bin/Aligner
cp $RECIPE_DIR/vg $PREFIX/bin/vg
cp $RECIPE_DIR/art_illumina $PREFIX/bin/art_illumina
if [ ! -d spades ]; then
git clone https://github.com/ablab/spades.git
fi
cd spades
git checkout spades_3.13.0
git apply ../patches/spade_muteBulgeRemover.patch
cd assembler/
./spades_compile.sh 
#ln -n spades.py $PREFIX/bin/spades.py
cd ../../

#python setup.py install --single-version-externally-managed --record=record.txt build_ext -i
#python setup.py install build_ext -i
$PYTHON -m pip install . --no-deps --ignore-installed --no-cache-dir -vvv
