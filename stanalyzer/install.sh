#!/bin/bash
if hash pip 2>/dev/null; then
   pip install Django 
   pip install MDAnalysis 
   pip install Pyhull 
else
   echo "pip not does not exist or permission denied!"
fi

if hash easy_install 2>/dev/null; then
   easy_install install Django 
   easy_install install MDAnalysis 
   easy_install install Pyhull 
else
   echo "easy_install does not exsit or permission denied!" 
fi

cd "$1"

if hash git 2>/dev/null; then
   git clone git@github.com:stanalyzer/ST-analyzer.git
else
   echo "git does not exist or permission denied!"
fi
