VPATH=src include
NVCC=nvcc  
CC=g++
obj_dir=obj
bin_dir=bin
INCLUDE_PATH=include
GCCFLAGS= ${addprefix -I,${INCLUDE_PATH}} 

ROOTFLAGS=`root-config --cflags`
ROOTLIBS=-L${ROOTSYS}/lib  -lCore -lRIO -lNet -lHist \
			-lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript \
			-lMatrix -lPhysics -lMathCore -lThread  -lm -ldl -lMinuit -m64  -I${ROOTSYS}/include \
					-lRooFit -lRooFitCore -w \
					-L /home/liul/software/myMinuit -I /home/liul/software/myMinuit \
					-lmyMinuit  -lMinuit2

obj1=${obj_dir}/rootfile.o
obj2=${obj_dir}/main.o
obj3=${obj_dir}/Amplitude.o
obj4=${obj_dir}/RooArgusPoly.o
obj5=${obj_dir}/RooArgusGauss.o
obj6=${obj_dir}/floatreduce.o
OBJGECTS=${obj1} ${obj2} ${obj3}  ${obj4} ${obj5} ${obj6}
TARGET=${bin_dir}/MLL

all : prepare ${TARGET}
${TARGET}: ${OBJGECTS}
	${NVCC} -o $@ $? ${GCCFLAGS} ${ROOTLIBS}
${obj2}:main.cpp
	${CC} -c $< -o $@ ${GCCFLAGS} ${ROOTLIBS}
${obj1}: rootfile.cu
	 ${NVCC} -c $< -o $@ ${GCCFLAGS} ${ROOTLIBS}
${obj3}: Amplitude.cu
	 ${NVCC} -c $< -o $@ ${GCCFLAGS} ${ROOTLIBS}
${obj4}: RooArgusPoly.cxx
	 ${CC} -c $< -o $@ ${GCCFLAGS} ${ROOTLIBS}
${obj5}: RooArgusGauss.cxx
	 ${CC} -c $< -o $@ ${GCCFLAGS} ${ROOTLIBS}
${obj6}: floatreduce.cu
	 ${NVCC} -c $< -o $@ ${GCCFLAGS} ${ROOTLIBS}


prepare:
	@if [ ! -d ${bin_dir} ]; then \
			mkdir -p ${bin_dir} ${obj_dir};  \
	fi
	

.PHONY:all clean

clean:
	-rm -rf  ${bin_dir} ${obj_dir}

