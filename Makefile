# Compile the project
graph:
	g++ -std=c++11 -fopenmp -I lib/Eigen src/graph.cpp src/utils.cpp src/main.cpp -o bin/graph

# Make the pdf slides from the markdown files
SLIDES := $(wildcard slides/*.md)
slides: $(SLIDES:md=pdf)

slides/%.pdf: slides/%.md
	pandoc $< --from markdown --to beamer -o $@

# Download and munge the data sets
data: init
	./setup/data.sh
	python src/make-datasets.py

init:
	./setup/init.sh
