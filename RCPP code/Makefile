CXX = C:/msys64/ucrt64/bin/g++.exe
CXXFLAGS = -fdiagnostics-color=always -g -I"\path\to\Rcpp\include" -I"\path\to\R\R-4.4.2\include"
SOURCES = test2.cpp
EXECUTABLE = test2.exe

all: $(EXECUTABLE)

$(EXECUTABLE): $(SOURCES)
	$(CXX) $(CXXFLAGS) -o $@ $^

clean:
	rm -f $(EXECUTABLE)


	"C:\Users\chiti\Desktop"
