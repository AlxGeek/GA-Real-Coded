CXX = g++
LINKER = g++

LFLAGS = -Wall 
CXXFLAGS = -Wall -c -std=c++14

OBJDIR := obj
SRCDIR := src
BINDIR := bin

TEXDIR := latex
DOCDIR := docs

TARGET = ga-rc.out

INCLUDES := $(wildcard $(SRCDIR)/*.hpp*)
SOURCES  := $(wildcard $(SRCDIR)/*.cpp)

OBJECTS := $(SOURCES:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)

rm = rm -f

$(BINDIR)/$(TARGET): $(OBJECTS)
	@mkdir -p $(BINDIR)
	@$(LINKER) $(LFLAGS) -o $@ $(OBJECTS) -lm
	@echo "Linking complete!"

$(OBJECTS): $(OBJDIR)/%.o : $(SRCDIR)/%.cpp 
	@mkdir -p $(OBJDIR)
	@$(CXX) $(CXXFLAGS) $< -o $@
	@echo "Compiled "$<" successfully!"

optimized:
	make CXXFLAGS="-O2 $(CXXFLAGS) "; LFLAGS="-O2 -Wall"

debug:
	make CXXFLAGS="-g $(CXXFLAGS) "; LFLAGS="-g -Wall"

.PHONY: latex
latex:
	@mkdir $(DOCDIR)
	pdflatex -output-directory $(DOCDIR) $(TEXDIR)/reporte.tex
	@echo "Reporte generado"

.PHONY: help
help:
	@echo "Se usan los siguientes headers: "$(INCLUDES)
	@echo "Se compilan los siguientes archivos: "$(SOURCES)
	@echo "Se genera el programa: "$(TARGET)" en "$(BINDIR)
	@echo "Se puede usar la opcion de compilacion debug para compilar con la opcion de debug."
	@echo "Se puede usar la opcion de compilacion optimized para compilar una version optimizada."
	@echo "La opcion clean borra los .o."
	@echo "La opcion cleanall borra las carpetas obj, bin y docs."
	@echo "La opcion latex genera el reporte."

.PHONY: clean
clean:
	@$(rm) $(OBJECTS) 
	@echo "Clean complete!"

.PHONY: cleanall
cleanall: clean
	@$(rm) -r $(BINDIR)	
	@$(rm) -r $(DOCDIR)
	@$(rm) -r $(OBJDIR)
	@echo "Cleanall complete!"





