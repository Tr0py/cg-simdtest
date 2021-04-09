simd-gc:
	g++ -o main main.cpp -mavx --std=c++17 -O0 $(ARGS)  -g

disasm:
	objdump -dSl main > disasm

clean: 
	rm main disasm

.PHONY: simd-gc disasm
