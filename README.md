# DiphotonAnalysis
  2
  3 ## Installation
  4
  5 To install this package you must compile a new library into ROOT for it to work.
  6 Begin by installing some python packages. I recommend setting up in a virtualenv:
  7
  8 ```bash
  9 virtualenv testGP/
 10 cd testGP
 11 source bin/activate
 12 ```
 13
 14 Then clone the git repo and install the needed python packages:
 15
 16 ```bash
 17 git clone git@github.com:Rob-Fletcher/DiphotonAnalysis.git
 18 pip install numpy
 19 pip install scipy
 20 pip install scikit-learn
 21 ```
 22
 23 Next we will compile some libraries into ROOT. This is not the best way to do this
 24 but it was the easiest way for me to get this working on my local installation.
 25 To start we will copy some files into your local ROOT installation.
 26
 27 ```bash
 28 cp DiphotonAnalysis/src/*.cxx <path-to-ROOT>/bindings/pyroot/src/
 29 cp DiphotonAnalysis/inc/*.h <path-to-ROOT>/bindings/pyroot/inc/
 30 ```
 31
 32 Then recompile root with whatever method you normally use. For cmake this looks something like:
 33 ```bash
 34 cd rootbuild
 35 cmake ../root
 36 cmake --build .
 37 ```
