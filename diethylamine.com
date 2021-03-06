%chk=diethylamine.chk
%mem=4GB
%nprocshared=6
#p opt=maxcycles=100 b3lyp/6-31+g(d,p) nosymm geom=connectivity
scf=(fermi,xqc)

Title Card Required

1 1
 C                  1.13024300    0.14301600   -0.01390700
 H                  1.45310000   -0.90279300   -0.03726300
 H                  1.45295100    0.64654200   -0.93083000
 H                  0.03646500    0.14501700   -0.01012100
 C                  1.61795600    0.86164500    1.23159400
 H                  1.29602400    1.90332200    1.26334700
 H                  1.29550000    0.36837500    2.14945100
 N                  3.16679100    0.89017300    1.28134700
 H                  3.53368700    1.84742900    1.28885300
 H                  3.57011200    0.41739600    0.46467000
 C                  3.69297606    0.21609042    2.47702336
 H                  3.28143071    0.67157446    3.35341637
 H                  3.42136596   -0.81855407    2.45175032
 C                  5.22745039    0.34346695    2.50465973
 H                  5.49919088    1.37823906    2.52208916
 H                  5.60982670   -0.14048537    3.37900435
 H                  5.63949387   -0.11885103    1.63208752

 1 2 1.0 3 1.0 4 1.0 5 1.0
 2
 3
 4
 5 6 1.0 7 1.0 8 1.0
 6
 7
 8 9 1.0 10 1.0 11 1.0
 9
 10
 11 12 1.0 13 1.0 14 1.0
 12
 13
 14 15 1.0 16 1.0 17 1.0
 15
 16
 17

