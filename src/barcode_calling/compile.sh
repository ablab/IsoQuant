g++ -c -fPIC cseqlib.cpp -o cseqlib.o; g++ -shared -Wl,-soname,cseqlib.so -o cseqlib.so  cseqlib.o
