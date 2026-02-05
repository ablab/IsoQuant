
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
rm -rf ./tests/out*/
rm -rf ./tests/.config/

cp -r ./src $TARGET_DIR/
cp -r ./tests $TARGET_DIR/
cp -r ./docs $TARGET_DIR/


rm -fr $TARGET_DIR/tests/short_reads_toy/

cp isoquant.py $TARGET_DIR/
cp isoquant_detect_barcodes.py $TARGET_DIR/
cp visualize.py $TARGET_DIR/
cp README.md $TARGET_DIR/
cp VERSION $TARGET_DIR/
cp LICENSE $TARGET_DIR/
cp changelog.md $TARGET_DIR/
cp requirements.txt $TARGET_DIR/
cp CODE_OF_CONDUCT.md $TARGET_DIR/

tar -pczf $TARGET_DIR.tar.gz $TARGET_DIR
rm -r $TARGET_DIR
