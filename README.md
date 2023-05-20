# LIVpropa

This is a plugin for the CRPropa code.
It enables simulations including Lorentz invariance violation (LIV).
For now, only Breit-Wheeler pair production is supported, but we will soon extend the code to include other processes such as inverse Compton scattering.


## Science


## Installation procedure

To install *livpropa*, you will need to have CRPropa 3 installed. 
Go to https://github.com/CRPropa/CRPropa3/ and follow the instructions.

Now proceed with the installation of this module.

1. Download the latest version of the code.
```
git clone https://github.com/rafaelab/livpropa.git
```

2. Navigate into the downloaded folder and create a folder called "build/".

3. Install the code with CMake and then make:

```
cmake ..
make
make install
```

4. If the code was properly compiled, you are ready to go!
Make sure to add the path where livpropa.py is created to your PYTHONPATH.
Alternatively, this can be hard-coded in your python script.

For consistency with CRPropa, the flag `-DCMAKE_CXX_FLAGS="-std=c++11"` may be required to force compilers to adopt C++11.


## Disclaimer

This program is provided 'as is', without warranties of any kind. 
Please use your discernement to interpret the results obtained with it.


## Acknowledgements