all: myprogram.cpp
		g++ –o strassen strassen.cpp
clean:
		$(RM) strassen