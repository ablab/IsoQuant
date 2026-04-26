
############################################################################
# Copyright (c) 2022-2024 University of Helsinki
# All Rights Reserved
# See file LICENSE for details.
############################################################################

VERSION="$(cat ./VERSION)"
TARGET_DIR=IsoQuant-$VERSION
rm -rf $TARGET_DIR
mkdir $TARGET_DIR

# cleaning .pyc and .pyo
rm -f */*.pyc
rm -f */*.pyo
rm -rf */__pycache__/
rm -f */*/*.pyc
rm -f */*/*.pyo
rm -rf */*/__pycache__/
rm -rf ./isoquant_tests/out*/
rm -rf ./isoquant_tests/.config/

cp -r ./isoquant_lib $TARGET_DIR/
cp -r ./isoquant_tests $TARGET_DIR/
cp -r ./docs $TARGET_DIR/


rm -fr $TARGET_DIR/isoquant_tests/short_reads_toy/
rm -fr $TARGET_DIR/isoquant_tests/simple_data/
rm -fr $TARGET_DIR/isoquant_tests/universal_data/

cp isoquant.py $TARGET_DIR/
cp isoquant_detect_barcodes.py $TARGET_DIR/
cp isoquant_visualize.py $TARGET_DIR/
cp README.md $TARGET_DIR/
cp VERSION $TARGET_DIR/
cp LICENSE $TARGET_DIR/
cp changelog.md $TARGET_DIR/
cp requirements.txt $TARGET_DIR/
cp CODE_OF_CONDUCT.md $TARGET_DIR/
cp GPL2.txt $TARGET_DIR/
cp pyproject.toml $TARGET_DIR/
cp MANIFEST.in $TARGET_DIR/

tar -pczf $TARGET_DIR.tar.gz $TARGET_DIR
rm -r $TARGET_DIR
