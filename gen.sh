#!/bin/zsh
scriptpath=`dirname $(readlink -f "$0")` && cd "$scriptpath"
if [[ ! -r "$1" ]]; then
	echo "File $1 not readable"
	exit 1
fi

# escapes symbols that cannot be treated by pSAscan
if [[ ! -r "$1.0" ]]; then
	echo ./standardize "${1}"
	time ./standardize "${1}"
fi
# generates SA from text
if [[ ! -r "$1.0.sa5" ]]; then
	echo ./psascan "${1}.0"
	time ./psascan "${1}.0"
fi
# generates BWT and ISA
if [[ ! -r "$1.0.isa5" ]]; then
	echo ./isaandbwt "${1}.0"
	time ./isaandbwt "${1}.0"
fi
# generates PLCP from BWT and SA
if [[ ! -r "$1.0.plcp" ]]; then
	echo ./construct_plcp_sequential "${1}.0"
	time ./construct_plcp_sequential "${1}.0"
fi
# generates LCP from text and SA
if [[ ! -r "$1.0.lcp5" ]]; then
	echo ./construct_lcp "${1}.0"
	time ./construct_lcp "${1}.0"
fi

#runs plcp, has to be generated with 'make plcp' in tudocomp
echo "./plcp ${1}.0 ${1}.0.tdc"
./plcp "${1}.0" "${1}.0.tdc"

# echo "./tdc -a 'lcpcomp(coder=bit)' -d --raw ${1}.0.tdc -f -o ${1}.0.orig"
# ./tdc -a 'lcpcomp(coder=bit)' -d --raw "${1}.0.tdc" -f -o "${1}.0.orig"

#decompress plcp-compressed file
echo "./tdc -a 'lcpcomp(coder=huff)' -d --raw ${1}.0.tdc -f -o ${1}.0.orig"
./tdc -a 'lcpcomp(coder=huff,dec=ext)' -d --raw "${1}.0.tdc" -f -o "${1}.0.orig"

# un-escape decompressed file
echo "./unstandardize ${1}.0.orig"
./unstandardize "${1}.0.orig"

# compare decompressed file with original file
echo "diff ${1} ${1}.0.orig.1"
diff "${1}" "${1}.0.orig.1"
[[ $? -eq 0 ]] && echo 'PERFECT'

# echo ./em_lpf "${1}.0" /local2/human.lz
# time ./em_lpf 800 "${1}.0" /local2/human.lz
# echo ./em_lzend "${1}.0"
# time ./em_lzend "${1}.0"
# echo ./em_lzscan "${1}.0"
# time ./em_lzscan "${1}.0" 800 /local2/human.lz

#/bighome/workspace/eSAIS/src/test/a.out "${1}.0.sa5"
#/bighome/workspace/eSAIS/src/test/a.out "$1.isa5"
#/bighome/workspace/stxxltools/build/readplcp "${1}.0"

