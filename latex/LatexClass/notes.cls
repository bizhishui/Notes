% notes.cls
% LaTeX document class for class notes
%
% Author: May Lynn Forssen
% mforssen@hmc.edu

\ProvidesClass{notes}[11/5/2013]
\NeedsTeXFormat{LaTeX2e}

% OPTIONS

\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}
\ProcessOptions \relax

% GEOMETRY

\RequirePackage[paper=a4paper]{geometry}

% LOAD PACKAGES

\RequirePackage{ifpdf}
\RequirePackage{amsmath}
\RequirePackage{amssymb}
\RequirePackage{fancyhdr}
\RequirePackage{framed}
\RequirePackage[usenames,dvipsnames]{color}
\RequirePackage{xcolor}
\RequirePackage{etoolbox}
\RequirePackage{paralist}


\LoadClass{article}

% COLORS

\definecolor{thmcolor}{RGB}{219,255,217}
\definecolor{proofcolor}{RGB}{186,235,183}
\definecolor{importantcolor}{RGB}{250,250,170}
\definecolor{keywordcolor}{RGB}{245,127,81}
\definecolor{headercolor}{RGB}{150,150,150}

% HEADER AND FOOTER

\fancypagestyle{firstpage}{%
  \fancyhead[LE,LO]{\color{headercolor}\@name}
  \fancyhead[CE,CO]{\color{headercolor}\@subject}
  \fancyhead[RE,RO]{\color{headercolor}\@date}
  \fancyfoot[CE,CO]{\color{headercolor}\thepage}
  \patchcmd{\headrule}{\hrule}{\color{headercolor}\hrule}{}{}
  \patchcmd{\footrule}{\hrule}{\color{headercolor}\hrule}{}{}
  \renewcommand{\headrulewidth}{0.4pt}
  \renewcommand{\footrulewidth}{0.4pt}
}



\pagestyle{firstpage}

\renewcommand{\footnoterule}{\relax}

% VARIABLES

\global\let\@subject\@empty
\newcommand{\subject}[1]{\def\@subject{#1}}

\global\let\@date\@empty
\renewcommand{\date}[1]{\def\@date{#1}}

\global\let\@name\@empty
\newcommand{\name}[1]{\def\@name{#1}}

% COMMANDS

\newcommand{\topic}[1]{\begin{center} \LARGE \bfseries #1 \end{center} \vspace{10pt}}

\newcommand{\subtopic}[1]{{\large\bfseries #1 \newline}}

\newcommand{\key}[1]{{\bfseries\color{keywordcolor} #1}}

\newcommand{\bul}[1]{\begin{compactitem} \item #1 \end{compactitem}}

\newcommand{\RedBold}[1]{{\bfseries\color{red} #1}}

\newcommand{\VioletBold}[1]{{\bfseries\color{Violet} #1}}

\newcommand{\Bold}[1]{{\bfseries #1}}


% ENVIRONMENTS

\newenvironment{important}{%
    \colorlet{shadecolor}{importantcolor}
    \begin{shaded}
        \par\noindent
    } {%
    \end{shaded}
}

\newenvironment{thm}[1]{%
    \colorlet{shadecolor}{thmcolor}
    \begin{shaded} 
        \par\noindent {\bfseries #1:}
    } {%
    \end{shaded}
}

\newenvironment{proof}{%
    \colorlet{shadecolor}{proofcolor}
    \begin{shaded}
        \par\noindent {\bfseries Proof:}
    }   {%
    \hfill $\square$
    \end{shaded}
}

\newenvironment{definition}[1]{%
    \par\noindent\key{#1:}
    } {%
}
