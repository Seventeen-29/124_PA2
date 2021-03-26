all: strassen
strassen:
	g++ -o strassen strassen.cpp
clean:
	$(RM) strassen && $(RM) opt