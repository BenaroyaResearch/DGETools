all: doc

doc:
	RScript -e 'devtools::document()'

