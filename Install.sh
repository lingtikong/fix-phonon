# Install/unInstall package classes in LAMMPS

if (test $1 = 1) then

  cp -p fix_phonon.h ..
  cp -p fix_phonon.cpp ..

elif (test $1 = 0) then

  rm ../fix_phonon.h
  rm ../fix_phonon.cpp

fi
