simd-gc:
	g++ -o main main.cpp -g -O0 --std=c++17

disasm:
	objdump -dSl main > disasm
clean: 
	rm main disasm
