# ior-curved
 A library to trace curved trajectories of light according to variations of the index of refraction.
 
## Installation

To download this library, you can use the command line:
```
git clone https://github.com/adolfomunoz/ior-curved.git
```

To compile all of the executables you can use [cmake](https://cmake.org/), as follows:
```
mkdir build
cd build
cmake -G "Unix Makefiles" .. -DCMAKE_BUILD_TYPE=Release
make
```

This downloads (or updates) all the required dependencies into the `external` subfolder and compiles all the provided examples and tests and moves the executable files to the `bin` subfolder. Other cmake generators (such as Visual Studio) have not been tested but might work as well. 

Once the external dependencies are downloaded, they won't be updated automatically, but they can be updated manually as follows
```
make update
```