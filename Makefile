VPATH=src include
CC=nvcc  
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
					-lmyMinuit

obj1=${obj_dir}/rootfile.o
obj2=${obj_dir}/main.o
obj3=${obj_dir}/Amplitude.o
OBJGECTS=${obj1} ${obj2} ${obj3}
TARGET=${bin_dir}/MLL

all : prepare ${TARGET}
${TARGET}: ${OBJGECTS}
	${CC} -o $@ $? ${GCCFLAGS} ${ROOTLIBS}
${obj2}:main.cu  
	${CC} -c $< -o $@ ${GCCFLAGS} ${ROOTLIBS}
${obj1}: rootfile.cu
	 ${CC} -c $< -o $@ ${GCCFLAGS} ${ROOTLIBS}
${obj3}: Amplitude.cu
	 ${CC} -c $< -o $@ ${GCCFLAGS} ${ROOTLIBS}


prepare:
	@if [ ! -d ${bin_dir} ]; then \
			mkdir -p ${bin_dir} ${obj_dir};  \
	fi
	

.PHONY:all clean

clean:
	-rm -rf  ${bin_dir} ${obj_dir}

