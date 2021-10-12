# ProBoundTools

This repository contains the source code for ProBoundTools (which manipulates and applies ProBound models) and buildCountTable.py (which builds count tables). 

To compile ProBoundTools, first install Maven and then execute:

cd ProBoundTools/
mvn package
cd ..

To run ProBoundTools, run:

java -cp ProBoundTools/target/ProBound-jar-with-dependencies.jar proBoundTools/App -c 'loadMotifCentralModel(15413).addNScoring().inputTXT(demo/initialPool.30mer1.R0.txt).bindingModeScores(/dev/stdout)'

