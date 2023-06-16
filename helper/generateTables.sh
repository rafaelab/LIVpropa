#!/bin/sh


# download CRPropa3-data repository
# this will provide additional tools to compute the LIV tables
if [ -d "CRPropa3-data" ] 
then
	echo "CRPropa3-data already exists. Continuing..." 
else
	echo "CRPropa3-data already exists. Downloading it..."
	git clone https://github.com/CRPropa/CRPropa3-data.git
fi



# compute the tables
echo "Generating the interaction tables..."
$PYTHON_EXECUTABLE computeElectromagneticInteractions.py

# copy to the installation directory (assumed to be build/data)
if [ -d "../build/data" ]
then
	cp -r ../data/PairProductionLIV/* ../build/data/PairProductionLIV/.
	cp -r ../data/PairProductionLIV/* ../build/share/livpropa/PairProductionLIV/.
	echo "Files copied to the appropriate folder." 
else
	echo "Installation path is not standard." 
	echo "You should mannually copy the data files from ../data/ to the folder where LIVPropa was installed."
fi