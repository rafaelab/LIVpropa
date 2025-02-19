# LIVpropa

This is a plugin for the CRPropa code.
It enables simulations including Lorentz invariance violation (LIV).

For now, only the following (QED) processes are considered:
- Breit-Wheeler pair production 
- inverse Compton scattering
- photon decay
- vacuum Cherenkov
Other processes will be added later.



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

To install the code locally, add the flag `-DCMAKE_INSTALL_PREFIX=$PWD` when running `cmake` from the `build/` folder (important!), as follows:
```
cmake .. -DCMAKE_INSTALL_PREFIX=$PWD
```
Then the python file `livpropa.py` will be generated. You can add this directory to your `PYTHONPATH` or do:
```
import sys
sys.path.append('/path/to/livpropa/build/')
from livpropa import *
```

IMPORTANT NOTICE: For some reason, gcc complains about the templates, whereas clang does not. Therefore, use clang. (I will probably not fix this.)


## Citation

If you use this code in your work, please use the following citation:
```
@Article{saveliev2024,
    author = "{Saveliev}, Andrey and {Alves Batista}, Rafael",
    title = "{Simulating electromagnetic cascades with Lorentz invariance violation}",
    doi = "10.1088/1361-6382/ad40f1",
    eid = "115011",
    eprint = "2312.10803",
    number = "11",
    pages = "115011",
    volume = "41",
    adsurl = "https://ui.adsabs.harvard.edu/abs/2024CQGra..41k5011S",
    archiveprefix = "arXiv",
    journal = "Classical and Quantum Gravity",
    month = "June",
    primaryclass = "astro-ph.HE",
    year = "2024"
}
```

## Disclaimer

This program is provided 'as is', without warranties of any kind. 
Please use your discernement to interpret the results obtained with it.


## Acknowledgements