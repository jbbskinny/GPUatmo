# GPUatmo

## What to do..
To build the code on your system you will need SCons:
```
$ brew install scons
```

Now, all you have to do is:
```
$ python pywrap/main.py
```
If you would like to adjust the initial conditions see the beginning of pywrap/main.py.

## Description

This program was built to model radiative transfer in neutron star atmospheres. Neutron stars are very small, dense stars.  The atmospheres are on the order of a centimeter thick, as compared to the Earth's atmosphere being on the order of miles.  The high temperatures ionize the atmosphere leaving a dense plasma of ions and electrons.  This plasma gives rise to some very interesting quantum mechanics and electromagnetic spectrums.

We have assumed local thermodynamic, hydrostatic and radiative equilibrium.  Another simplification was the assumption of plane-parallel geometry.  This allowed us to exclude time-dependency and exclude solid-angle when calculating the extinction coefficients.

To see further details see slide.pdf.

