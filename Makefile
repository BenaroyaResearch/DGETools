all: doc

doc:
	RScript -e 'devtools::document()'

install:
	Rscript -e "devtools::install()"

