simd-gc:
	g++ -o main main.cpp -g --std=c++17

disasm:
	objdump -dSl main > disasm

clean: 
	rm main disasm

.PHONY: simd-gc disasm
