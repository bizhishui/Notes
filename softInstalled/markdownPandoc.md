#Convert chinese markdown file to pdf using pancoc#
Some sites useful

[pandoc中文pdf转换攻略](http://afoo.me/posts/2013-07-10-how-to-transform-chinese-pdf-with-pandoc.html)

[Convert Markdown files into PDF via LaTeX using Pandoc](http://felixfan.github.io/Pandoc-LaTeX-Chinese/)

[How-To Generate PDF From Markdown](http://kevin.deldycke.com/2012/01/how-to-generate-pdf-markdown/)

[Pandoc Demos](http://pandoc.org/demos.html)

[让pandoc输出pdf时支持中文](http://www.verydemo.com/demo_c173_i98333.html)

在正确安装Texlive (with chinese support)和pandoc后，可以使用pandoc将含中文的markdown文档转为pdf文档。具体包含如下几点

- [生成合适的template.tex文件](http://www.verydemo.com/demo_c173_i98333.html)，并指定中文字体。
- 使用命令pandoc -o test.pdf test.md  --latex-engine=xelatex --template=template.tex生成pdf。需要注意template.tex的路径信息。
- 为方便可以将template.tex放置在\$HOME/Templates文件夹下，并在\$HOME/.bashrc中添加
  alias pandocCN="pandoc --template=\$HOME/Templates/template.tex --latex-engine=xelatex"。
  这样就可以使用pandocCN -s test.md -o test.pdf
