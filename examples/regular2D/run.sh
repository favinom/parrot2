nx=80      

ref=9

np=1


time mpirun -n ${np} ../../parrot2-opt -i regular2D_error.i  nxMoose=${nx}  nxStringMoose=${nx} ref=${ref}
# time mpirun -n ${np} ../../parrot2-opt -i smallFeatures2D_GMG.i    nxMoose=${nx} nyMoose=${ny} nxStringMoose=${nxString}
# time mpirun -n ${np} ../../parrot2-opt -i smallFeatures2D_SSFBBG.i nxMoose=${nx} nyMoose=${ny} nxStringMoose=${nxString}

#/Applications/ParaView-5.9.0-RC4.app/Contents/bin/pvpython postprocessor.py output_${nxString}.e.${np}.0
# /Applications/ParaView-5.9.0-RC4.app/Contents/bin/pvpython postprocessor.py output_GMG_${nxString}.e.${np}.0    GMG_
# /Applications/ParaView-5.9.0-RC4.app/Contents/bin/pvpython postprocessor.py output_SSFBBG_${nxString}.e.${np}.0 SSFBBG_

rm *_0.csv
