SUBDIRS = SeqLib/htslib SeqLib/bwa SeqLib/fermi-lite SeqLib/src src

apophenia:
	cd apophenia-pkg && ./configure --prefix=`pwd` && make && make install

install:
	mkdir -p bin && mv src/ginseng bin
