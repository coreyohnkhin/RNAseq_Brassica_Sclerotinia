cd software
wget https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2
tar -xf samtools-1.21.tar.bz2
cd samtools-1.21
./configure --prefix=$(pwd)/../samtools
make
make install
cd ..
rm samtools-1.21.tar.bz2
cd ..
