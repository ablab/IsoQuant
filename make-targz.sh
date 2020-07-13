VERSION="$(cat ./VERSION)"
TARGET_DIR=IsoQuant-$VERSION
rm -rf $TARGET_DIR
SRC_DIR=$TARGET_DIR/src
mkdir -p $SRC_DIR

cp -r ./src $SRC_DIR/

# cleaning .pyc and .pyo
rm -f $SRC_DIR/*.pyc
rm -f $SRC_DIR/*.pyo
rm -rf $SRC_DIR/__pycache__/

cp isoquant.py $TARGET_DIR/
cp README.md $TARGET_DIR/
cp VERSION $TARGET_DIR/
cp LICENSE $TARGET_DIR/
cp changelog.html $TARGET_DIR/
cp requirements.txt $TARGET_DIR/

tar -pczf $TARGET_DIR.tar.gz $TARGET_DIR
rm -r $TARGET_DIR
