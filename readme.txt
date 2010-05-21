Getting Started with the 1dflameV2 Code

[Last updated 5/21/2009]

0. Dependinces: These directions assume you are using a Linux system
with the following installed:
   * Matlab
   * Boost.Build
   * Cantera
   * Sundials
   * Libconfig

1. The 1dflameV2 code is stored in a Git repository on Pharos. To
check out a copy of the code, run:

    $ git clone pharos:/var/cache/git/1dflameV2.git

2. Edit the file "1dflameV2/Jamroot" to point to the correct include
and library directories for your system.

3. Compiling the Code:

    $ cd 1dflameV2
    $ bjam

This produces the binary "bin/1dflameV2".

4. Prepare your kinetics mechanism. If you're starting with a
mechanism that's in the Chemkin format, it will need to be converted:

    $ source /usr/local/bin/setup_cantera
    $ ck2cti -i <mech> -t <thermo> -tr <transport>

This will produce mech.cti, which needs to be further converted:

    $ cti2ctml mech.cti

to produce mech.xml.

5. Prepare the input file, based on the one specified in

    1dflameV2/input/example-input.txt

Specify mech.xml (which needs to be in whatever directory is specified
as the "input" directory in input file) as the mechanismFile and "gas"
as the phaseID. The other parameters you may want to change are in the
InitialCondition and StrainParameters section of the input file.

6. Run the code: The code takes the path to the input file as an
argument. The code relies on the Matlab libraries for saving data, so
the path to those libraries must be specified.

    $ export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/matlab_r2007a/bin/glnxa64
    $ ./1dflameV2 input/myinputfile.txt

7. Examine the output files:
   * integral.mat contains integral flame properties (flame speed, e.g.) as a function of time
   * profNNNNNN.mat contain the temperature & species profiles output periodically.
