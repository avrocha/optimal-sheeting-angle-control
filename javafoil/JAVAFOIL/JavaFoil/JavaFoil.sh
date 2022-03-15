#!/bin/sh

cd ~/MH-AeroTools/JavaFoil

ARCHS=javafoil.jar:mhclasses.jar
CLASS=MH.JavaFoil.JavaFoil

java -cp $ARCHS $CLASS
