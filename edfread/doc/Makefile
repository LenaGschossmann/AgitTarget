all: dev.pdf user.pdf

install:
	cp dev.pdf ../current
	cp user.pdf ../current

user.pdf: user.sgml
	docbook2pdf user.sgml

dev.pdf: latex/refman.pdf
	mv latex/refman.pdf dev.pdf

latex/refman.pdf: Doxyfile
	doxygen
	(cd latex && make)

clean:
	rm -f dev.pdf
	rm -f user.pdf
	rm -rf latex
