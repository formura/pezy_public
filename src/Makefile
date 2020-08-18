CXX = mpicxx
CXXFLAG = -std=c++11 -O2 -Wall -Iinclude
LDFLAG = -lm -lgsl -lgslcblas -lfftw3 -lOpenCL

PZ_OLEVEL = -O2

all: main

run: main
	mkdir -p data
	./main

clean:
	rm -rf main *.tar.bz2 *.o kernel.pz pzcc-temps

dist: main
	tar -cjf main.tar.bz2 *.cpp *.h kernel.pz include Makefile

main: main.o filter.o update_sc.o kernel.o
	$(CXX) $(CXXFLAG) $^ $(LDFLAG) -o $@

%.s: %.cpp
	$(CXX) $(CXXFLAG) -S $< -o $@

%.o: %.cpp cache_profiler.h config.h update.h
	$(CXX) $(CXXFLAG) -c $< -o $@

# pzcc-temps/kernel.s: kernel_gen.pzc deriv.txt config.h Makefile
# 	$(PZSDK_HOME)/bin/clang -target pz64 -mcpu=sc2 -xc++ -std=c++11 -I$(PZSDK_HOME)/inc -Iinclude -E $< > kernel.pzc
# 	rm -rf pzcc-tmps

# 	pzcc kernel.pzc --release --save-temps --verbose-pz-hexobj $(PZ_OLEVEL) -o annot.pz && (rg 'i64\.addi\tsp sp [^-]' pzcc-temps/annot.s | sort -V | tail -3)

# 	$(PZSDK_HOME)/bin/llvm-dis pzcc-temps/annot.bc
# 	mv pzcc-temps/annot.ll .

# 	pzcc kernel.pzc --release --save-temps $(PZ_OLEVEL)
# 	rm kernel.pz

pzcc-temps/kernel.s: kernel.pzc config.h update.h Makefile
	rm -rf pzcc-tmps

	pzcc kernel.pzc --release --save-temps --verbose-pz-hexobj $(PZ_OLEVEL) -o annot.pz && (rg 'i64\.addi\tsp sp [^-]' pzcc-temps/annot.s | sort -V | tail -3)

	$(PZSDK_HOME)/bin/llvm-dis pzcc-temps/annot.bc
	mv pzcc-temps/annot.ll .

	pzcc kernel.pzc --release --save-temps $(PZ_OLEVEL)
	rm kernel.pz

pzcc-temps/kernel.opt.s: pzcc-temps/kernel.s
	# cd ../scheduler/ && cargo run --release ../src/$<
	cp $< $@

kernel.pz: pzcc-temps/kernel.opt.s
	$(PZSDK_HOME)/bin/clang "-target" "pz64" "-mcpu=sc2" "-I$(PZSDK_HOME)/inc" "-E" "-DPZ_LLVMMC_ASM" "-DPZC_RELEASE_BUILD" "-xc" "-o" pzcc-temps/kernel.opt.expanded.s $< "-v"

	$(PZSDK_HOME)/bin/llvm-mc "-filetype=obj" "-arch=pz64" "-mcpu=sc2" "-verbose-pz-hexobj" "-o" pzcc-temps//kernel.hex pzcc-temps/kernel.opt.expanded.s

	$(PZSDK_HOME)/bin/pzlink "-o" "kernel.pz" pzcc-temps/kernel.hex "-v"

%.o: %.pz
	objcopy -I binary -O elf64-x86-64 -B i386 $< $@
