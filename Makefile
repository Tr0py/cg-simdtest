simd-gc:
	g++ -o main main.cpp --std=c++17 -mavx -O3

disasm:
	objdump -dSl main > disasm

clean: 
	rm main disasm

.PHONY: simd-gc disasm
