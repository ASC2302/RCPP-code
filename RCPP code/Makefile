CXX = C:/msys64/ucrt64/bin/g++.exe
CXXFLAGS = -fdiagnostics-color=always -g -I"C:\Users\chiti\AppData\Local\R\win-library\4.4\Rcpp\include" -I"C:\Program Files\R\R-4.4.2\include"
SOURCES = test2.cpp
EXECUTABLE = test2.exe

all: $(EXECUTABLE)

$(EXECUTABLE): $(SOURCES)
	$(CXX) $(CXXFLAGS) -o $@ $^

clean:
	rm -f $(EXECUTABLE)


	"C:\Users\chiti\Desktop"