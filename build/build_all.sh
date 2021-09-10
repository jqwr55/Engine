nasm ../src/core/asm/funcs.asm -g -o ../bin/asm_func.o -felf64
./core_build.sh &
./renderer_build.sh &
rm ../src/core/cpp/main.S &
wait
objdump -drwC -Mintel ../bin/core >> ../src/core/cpp/main.S
