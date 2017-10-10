# QC-FPT
Quantum Chemical Fragment Parent Tests
This set of Perl scripts implements the QC-FPT algorithm described in "Quantum Chemical Fragment Parent Tests", B. G. Janesko, L. Li, R. Mensing, Anal. Chim. Acta, in press 

The scripts depend on the PerlMol Chemistry module. 

Use of the scripts is illustrated with the m/z 46 fragment obtained from the [M+H](+) parent ion of diethylamine. Begin with diethylamine.com and diethylamine.log, input and output files from a Gaussian 09 geometry optimization of the parent ion.  

First run 
perl DriverRoutine.prl diethylamine Write 46 

This creates a new subdirectory ./diethylamine_46/ and fills it with Gaussian input files for all the fragments

Next run Gaussian 09 on all of the input files in ./diethylamine_46/

Next run 
perl DriverRoutine.prl diethylamine Read 46 

This reads all of the fragment output files and writes the predicted fragment filenames and fragmentation energies in diethylamine_46.res 

Finally run 
cat diethylamine_46.res|grep "^ *46" |sort -g -k2 | head -n4 

This prints the four most stable fragmentations giving fragment ions of the desired m/z. I find the following output: 

# Mass DE Efragment Etotal Filename
    46.092    28.1  -135.554132  -214.153777 001.log,007.log
    46.092    28.1  -135.554132  -214.153777 001.log,007.log

