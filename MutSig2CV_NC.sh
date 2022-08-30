maf=$1
out=$2

MCRROOT=/share/public/software/Onc_Soft/MATLAB_Compiler_Runtime/v901
export LD_LIBRARY_PATH=.:${MCRROOT}/runtime/glnxa64
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/bin/glnxa64
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/sys/os/glnxa64
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/sys/opengl/lib/glnxa64


bin/MutSig2CV_NC $maf $out run/params/Liver-HCC_alt_prom.params.txt
