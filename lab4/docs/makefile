pan: Lab3.md
	pandoc --pdf-engine=xelatex -V mainfont="Liberation Serif" Lab3.md -o lab3.pdf

title: title.pdf
	pdftk title.pdf cat 1 output title_o.pdf
	rm title.pdf
	mv title_o.pdf title.pdf

merge: title pan
	pdftk title.pdf Lab4.pdf cat output merged.pdf


# very dirty
# merge: title.pdf Lab1.pdf Lab1_no_code.pdf
#	convert title.pdf Lab1.pdf lab1_Соколов_Б20_505.pdf
#	convert title.pdf Lab1_no_code.pdf lab1_1_Соколов_Б20_505.pdf
