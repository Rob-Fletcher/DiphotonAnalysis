# DiphotonAnalysis

## Installation

To install this package you must compile a new library into ROOT for it to work.
Begin by installing some python packages. I recommend setting up in a virtualenv:

```bash
virtualenv testGP/
cd testGP
source bin/activate
```

Then clone the git repo and install the needed python packages:

```bash
git clone git@github.com:Rob-Fletcher/DiphotonAnalysis.git
pip install numpy
pip install scipy
pip install scikit-learn
```

Next we will compile some libraries into ROOT. This is not the best way to do this
but it was the easiest way for me to get this working on my local installation.
To start we will copy some files into your local ROOT installation.

```bash
cp DiphotonAnalysis/src/*.cxx <path-to-ROOT>/bindings/pyroot/src/
cp DiphotonAnalysis/inc/*.h <path-to-ROOT>/bindings/pyroot/inc/
```

Then recompile root with whatever method you normally use. For cmake this looks something like:
```bash
cd rootbuild
cmake ../root
cmake --build .
```
