TARGETS=O0 O2 sse-aligned sse-unaligned avx-aligned avx-unaligned

all:$(TARGETS)
ARGS=-DSIZE=64 -DLOOPTIME=1 -DSIMD_OPT
COMPILER_OPT=-DSIZE=64 -DLOOPTIME=1

O0:
	g++ -o $@ main.cpp -mavx --std=c++17 -O0 $(COMPILER_OPT)
O2:
	g++ -o $@ main.cpp -mavx --std=c++17 -O2 $(COMPILER_OPT)
sse-aligned:
	g++ -o $@ $@.cpp -msse3 --std=c++17 -O0 $(ARGS) 
sse-unaligned:
	g++ -o $@ $@.cpp -msse3 --std=c++17 -O0 $(ARGS) 
avx-aligned:
	g++ -o $@ $@.cpp -mavx --std=c++17 -O0 $(ARGS) 
avx-unaligned:
	g++ -o $@ $@.cpp -mavx --std=c++17 -O0 $(ARGS) 

disasm:
	objdump -dSl main > disasm

clean: 
	rm -f disasm $(TARGETS)

.PHONY: $(TARGETS)

