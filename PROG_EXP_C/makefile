
prog.exe : mainc.o subf.o VTSW.o 
	gcc -g -o prog.exe mainc.o subf.o VTSW.o -lm 

mainc.o : mainc.c
	gcc -g -c mainc.c 	

VTSW.o : VTSW.c
	gcc -g -c VTSW.c
subf.o : subf.c
	gcc -g -c subf.c


