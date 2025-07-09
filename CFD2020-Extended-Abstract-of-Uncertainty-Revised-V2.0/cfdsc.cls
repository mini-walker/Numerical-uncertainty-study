\def\update{03/01/23 v1.0}

\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{cfdsc}[\update\ Societe Canadienne de CFD Society of Canada]

%=================================================
% Classe de base - Basic class
%=================================================
\LoadClass[twocolumn,letterpaper,10pt]{article}
\RequirePackage{pslatex} % Gives cleaner translation to PDF
\RequirePackage{graphicx}
\RequirePackage{graphics}

%======================================================================
% Title and author stuff styled after the A&A Document class, designed
% for the Astronomy and Astrophysics Journal
%======================================================================
\newcounter{inst}
\newcounter{auth}
\let\@date\undefined
\def\institute#1{\gdef\@institute{#1}}
\def\email#1{\gdef\@email{#1}}

\def\institutename{\par
 \begingroup
 \newdimen\instindent
 \instindent=0pt
 \parindent=0pt
 \parskip=0pt
 \setcounter{inst}{1}%
 \def\and{\par\stepcounter{inst}%
 \hangindent\instindent
 \noindent
 \hbox to\instindent
 {{}$^{\theinst}$\enspace}\ignorespaces}%
 \setbox0=\vbox{\def\thanks##1{}\@institute}
 \ifnum\value{inst}>9\relax\setbox0=\hbox{$^{88}$\enspace}%
                \else\setbox0=\hbox{$^{8}$\enspace}\fi
 \instindent=\wd0\relax
 \ifnum\value{inst}=1\relax
 \else
   \setcounter{inst}{1}%
   \hangindent\instindent
   \noindent
   \hbox to\instindent
   {{}$^{\theinst}$\enspace}\ignorespaces
 \fi
 \ignorespaces
 \@institute\par
 \endgroup}

\def\inst#1{\unskip$^{#1}$}
\def\fnmsep{\unskip$^,$}

\def\maketitle{\par
\begingroup
%
    \renewcommand\thefootnote{\@fnsymbol\c@footnote}%
    \def\@makefnmark{$^{\@thefnmark}$}%
%
 \def\lastand{\ifnum\value{inst}=2\relax
                 \unskip{} and\
              \else
                 \unskip, and\
              \fi}%
 \def\and{\stepcounter{auth}\relax
          \ifnum\value{auth}=\value{inst}%
             \lastand
          \else
             \unskip,
          \fi}%
 \ifnum \col@number=\@ne
   \@maketitle
 \else
   \twocolumn[\@maketitle]%
 \fi
 \thispagestyle{empty}
 \global\@topnum\z@
 \@thanks
 \endgroup
 \setcounter{footnote}{0}
 \setcounter{footnote}{0}%
 \global\let\thanks\relax
 \global\let\maketitle\relax
 \global\let\@maketitle\relax
 \global\let\@thanks\@empty
 \global\let\@author\@empty
 \global\let\@date\@empty
 \global\let\@title\@empty
 \global\let\@subtitle\@empty
 \global\let\@abstract\@empty
 \global\let\title\relax
 \global\let\author\relax
 \global\let\aabstract\relax
 \global\let\date\relax
 \global\let\and\relax
}

\def\@maketitle{
 \begin{center}
  \vskip24pt
 {\LARGE \bfseries\boldmath
  \@title \par}\vskip 14pt
 \setbox0=\vbox{\setcounter{auth}{1}\def\and{\stepcounter{auth} }%
                \hfuzz=2\textwidth\def\thanks##1{}\@author}%
 \value{inst}=\value{auth}%
 \setcounter{auth}{1}%
 {\large \@author}
 \end{center}
 \vskip 12pt
 {\centering \normalsize\itshape\institutename}
 \vskip 12pt
 \begin{center}
 Email: \normalsize\itshape\@email
 \end{center}
 \vskip 12pt
}

%=================================================
% Redéfinition des sections et sous-sections
% Redefinition of the sections and subsections
%=================================================
\renewcommand\section{\@startsection
  {section}%                    % section name
  {1}%                          % level
  {\z@}%                        % indentation of heading
  {3ex}%                        % before skip (neg, no parindent)
  {0.1ex}%                      % after skip (neg, run-on heading space)
  {\normalfont\large\scshape\bfseries\centering}}
%=================================================
% Mot clé pour identifier la section références:
% en majuscules
% Keyword to identify the references section:
% in cap. letters
%=================================================
\renewcommand\refname{\sc{References}}

%=================================================
% Distances pour la mise en page - For the layout
%=================================================
\setlength{\parindent}{0mm}
\setlength{\parskip}{6pt}
\usepackage{geometry}
\geometry{tmargin=0.75in,bmargin=1in,lmargin=0.75in,rmargin=0.75in,headheight=0in,headsep=0in,footskip=0in}
\setlength{\columnsep}{0.25in}

%=================================================
% Macros
%=================================================
\newcommand{\bm}[1]{\boldsymbol{#1}}

%=================================================
% Pas de pagination: \pagestyle{empty}
%=================================================
\pagestyle{empty}
\flushbottom

\bibliographystyle{abbrv}
