simd-gc:
	g++ -o main main.cpp -msse3 --std=c++17 -O3 $(ARGS) 

disasm:
	objdump -dSl main > disasm

clean: 
	rm main disasm

.PHONY: simd-gc disasm
