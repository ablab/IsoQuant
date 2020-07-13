VERSION="$(cat ./VERSION)"
TARGET_DIR=IsoQuant-$VERSION
rm -rf $TARGET_DIR
mkdir $TARGET_DIR

cp -r ./src $TARGET_DIR/
cp -r ./tests $TARGET_DIR/

# cleaning .pyc and .pyo
rm -f */*.pyc
rm -f */*.pyo
rm -rf */__pycache__/
rm -f */*/*.pyc
rm -f */*/*.pyo
rm -rf */*/__pycache__/

cp isoquant.py $TARGET_DIR/
cp README.md $TARGET_DIR/
cp VERSION $TARGET_DIR/
cp LICENSE $TARGET_DIR/
cp changelog.html $TARGET_DIR/
cp requirements.txt $TARGET_DIR/

tar -pczf $TARGET_DIR.tar.gz $TARGET_DIR
rm -r $TARGET_DIR
