all: strassen opt
strassen:
	g++ -o strassen strassen.cpp
opt:
	g++ -o opt opt_strassen.cpp
clean:
	$(RM) strassen && $(RM) opt