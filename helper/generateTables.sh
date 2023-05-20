#!/bin/sh


# Download CRPropa3-data repository.
# This will provide additional tools to compute the LIV tables.
if [ -d "CRPropa3-data" ] 
then
	echo "CRPropa3-data already exists. Continuing..." 
else
	echo "CRPropa3-data already exists. Downloading it..."
	git clone https://github.com/CRPropa/CRPropa3-data.git
fi



# Compute the tables
$PYTHON_EXECUTABLE computeElectromagneticInteractions.py