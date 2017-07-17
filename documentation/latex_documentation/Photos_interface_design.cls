\ProvidesClass{my_report}
\LoadClass[11pt, titlepage]{article}
\usepackage[paperwidth=21cm,paperheight=29.7cm,includehead,headheight=2.0cm,pdftex,hmargin={3cm,2.5cm},vmargin={0cm,2cm},]{geometry} 
\setlength{\parindent}{0cm} 
\renewcommand{\baselinestretch}{1.2}
\parskip=1em

\renewcommand{\contentsname}{Table of Contents}

\renewcommand{\normalsize}{\fontsize{11pt}{11pt}\selectfont}
\renewcommand{\title}[1]{\textsf{\Huge{#1}}}
\renewcommand{\author}[1]{\textsf{\Huge{#1}}}
\newcommand{\subtitle}[1]{\textsf{\textbf{\Large{#1}}}}
\newcommand{\abstractheading}[1]{\textsf{\textbf{\LARGE{#1}}}}
\newcommand{\code}[1]{\texttt{\footnotesize{#1}}}

\renewcommand{\section}{\@startsection
{section}%                    % the name
{1}%                          % the level
{0mm}%                        % the indent
{10.6mm}%             % the beforeskip
{4.2mm}%           % the afterskip
{\LARGE\bfseries\sffamily}}  % the style

\renewcommand{\subsection}{\@startsection
{subsection}%                    % the name
{2}%                          % the level
{0mm}%                        % the indent
{6.4mm}%             % the beforeskip
{1.1mm}%           % the afterskip
{\Large\bfseries\sffamily}}  % the style

\renewcommand{\subsubsection}{\@startsection
{subsubsection}%                    % the name
{3}%                          % the level
{0mm}%                        % the indent
{4.2mm}%             % the beforeskip
{1.1mm}%           % the afterskip
{\normalsize\bfseries\sffamily}}  % the style

%% Chapter headings should be centered, uppercase, and at the top of the page.
\def\@makechapterhead#1{%
  { \parindent \z@ \raggedright \normalfont
    %\centering
    \ifnum \c@secnumdepth >\m@ne
        \huge\textsf{\@chapapp\space \thechapter:}
        % \par\nobreak
        %\vskip 20\p@
    \fi
    \interlinepenalty\@M
    \huge \bfseries \textsf{#1}\par\nobreak
    \rule{5cm}{0.5pt}
    \vskip 20\p@
  } }
  
\def\@makeschapterhead#1{%
  %\vspace*{50\p@}%
  { \parindent \z@ \raggedright
    %\centering
    \normalfont
    \interlinepenalty\@M
    \huge \bfseries  \textsf{#1}\par\nobreak
    \rule{5cm}{0.5pt}
    \vskip 20\p@

  }}
  
 \renewenvironment{abstract}{%
      \chapter*{\abstractname}%
      \addcontentsline{toc}{chapter}{\abstractname}
 }
     
     
\makeatletter
\renewcommand{\l@section}{\bfseries\@dottedtocline{1}{0em}{2.3em}}
\renewcommand{\l@subsection}{\normalfont\@dottedtocline{2}{2em}{2.3em}}
\renewcommand{\l@subsubsection}{\normalfont\@dottedtocline{3}{2em}{2.3em}}
%\renewcommand{\l@subsubsection}{\normalfont\@dottedtocline{4}{2em}{2.3em}}
\makeatother

\def\maketitle{\begin{titlepage}
\thispagestyle{empty}
\let\footnotesize\small \let\footnoterule\relax \setcounter{page}{0}
\null
\vfil
\begin{flushright}
{\bf  CERN-PH-TH/2010-261 \\
  IFJPAN-IV-2010-6  \\
 Updated version; for the program as of  April 2015}
\end{flushright}

\begin{center}
\title{PHOTOS Interface in C++\\[0.5cm]{} 
\textbf{\small Technical and Physics Documentation  \\[0.1cm]\rule{4cm}{1pt}\\[0.3cm]}}

\author{\LARGE N. Davidson$^{a}$, T. Przedzinski$^{b,}$%
\footnote{Corresponding author; e-mail: {\tt tomasz.przedzinski@cern.ch}; phone no.: $+48 503651129$ },
Z. Was$^{c,d}$}\\[0.3cm]
\textsf{\Large{\em {\em $^a$ Murdoch Childrens Research Institute, Royal
Children’s Hospital, \\ Flemington Road, Parkville 3052, Melbourne, VIC,
Australia
 }\\
{\em $^b $ The Faculty of Physics, Astronomy and Applied Computer
Science,\\
Jagellonian University, Stanislawa
Lojasiewicza 11,
30-348 Cracow, Poland}
 \\
{\em $^c$ Institute of Nuclear Physics, 
Polish Academy of Sciences,\\
         ul. Radzikowskiego 152, 31-342 Cracow, Poland}\\
{\em $^d $ Theory Group, Physics Department, CERN, CH-1211, Geneva 23, Switzerland}}}
\\[0.3cm]

%\includegraphics[height=3cm]{logo.eps}\\[1cm]
{\bf \Large Abstract} \vskip 2 mm 

\end{center}

{ For five years now, {\tt PHOTOS} Monte Carlo for
bremsstrahlung in the decay of particles and resonances
has been available with an interface to the {\tt C++ HepMC} event record. 
The main purpose of the present
paper is to document the technical aspects of the {\tt PHOTOS} Monte Carlo 
installation and present version use. A multitude of test results and examples are  
distributed together 
with the program code.

The
{\tt PHOTOS} {\tt C++} physics precision is   better than  its {\tt FORTRAN}
predecessor and more convenient steering options are also available. 
An algorithm for the event record interface
necessary for process dependent 
photon emission kernel is implemented. 
It is used  in $Z$ and $W$ decays for kernels
of complete first order matrix elements of the decays.
Additional  emission of  final state lepton  pairs  is also available.

Physics assumptions used in the program and properties of the
solution are reviewed. 
In particular, it is explained how
the  second order matrix elements were used in design and 
validation of the program iteration procedure. Also, it is explained
that the phase space parameterization used in the program is exact. 
 }

\centerline{\normalsize Submitted to Comp. Phys. Commun.}

\vfil

\footnoterule
\noindent
{\footnotesize
\begin{itemize}
\item[${\dag}$]

This project is financed in part from funds of Polish National Science
Centre under decisions  DEC-2011/03/B/ST2/00107 (TP), 
%DEC-2011/03/B/ST2/02632 (ZW), 
DEC-2011/03/B/ST2/00220 (ZW).
This research was supported in part by
the Research Executive Agency (REA) of the European Union under
the Grant Agreement PITN–GA–2012–316704 (HiggsTools).

% \\
%{\em $^d$ Until September 1 2013 Grup de F\'{\i}sica Te\`orica, Institut de F\'{\i}sica d'Altes Energies, 
%Universitat Aut\`onoma de Barcelona, E-08193 Bellaterra, Barcelona, Spain} \\
%MS DEC-2012/04/M/ST2/00240"
%ER DEC-2011/03/B/ST2/00220
%ZW DEC-2011/03/B/ST2/00107"
%SJ DEC-2011/03/B/ST2/02632"
\end{itemize}
}
\null
\end{titlepage}}

\pagestyle{plain}

