#TexLive 2015 安装#
##方法一：##
sudo apt-get purge texlive*

    rm -rf /usr/local/texlive/2012 and rm -rf ~/.texlive2012

    rm -rf /usr/local/share/texmf

    rm -rf /var/lib/texmf

    rm -rf /etc/texmf

    sudo apt-get remove tex-common --purge

    rm -rf ~/.texlive


##方法二：##
[Ref. Site](http://seisman.info/install-texlive-under-linux.html)

1. 安装TexLive+Latex+CJK：

如果硬盘充裕的话，直接完整安装也可以：

sudo apt-get install texlive-full latex-beamer

**也可下载完整镜像，挂载，sudo ./instal-tl -gui, 并将其中NO->YES**

安装完后，就可以安装CJK的相关软件包了，如果只需要获得中文支持，那么执行：

sudo apt-get install latex-cjk-chinese ttf-arphic-* hbf-*

否则，建议安装latex-cjk-all以获取完整支持。

2. 环境变量

\# TeX Live 2015

export MANPATH=\${MANPATH}:/usr/local/texlive/2015/texmf-dist/doc/man

export INFOPATH=\${INFOPATH}:/usr/local/texlive/2015/texmf-dist/doc/info

export PATH=\${PATH}:/usr/local/texlive/2015/bin/x86_64-linux

3. 中文

[Ref. Site2](http://www.cnblogs.com/JackieMe/p/4677555.html)

xeCJK将默认使用TeXLive自带的Fandole字体。要能够编译通过，需要将TeXLive自带的中文字体安装到系统中，最简单的办法是在~/.fonts目录下建一个软链接:

ln -s /opt/texlive/2014/texmf-dist/fonts/opentype/public/fandol ~/.fonts/

使用时在文档前面加入

\\documentclass{article}

\\usepackage{xeCJK}

即可使用xelatex编译
