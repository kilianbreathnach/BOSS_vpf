Finding the Void Probability Function of galaxies in BOSS data (CMASS and LOWZ)
===============================================================================

Kilian Walsh and Jeremy Tinker
(arxiv article: ...)


Introduction
------------

> And I, infinitesimal being,
> drunk with the great starry
> void,
> likeness, image of
> mystery,
> felt myself a pure part
> of the abyss,
> I wheeled with the stars,
> my heart broke loose on the wind.
>
> - Pablo Neruda

The Void Probability Function (VPF) is a probability density function over the
radius of randomly placed spheres in space, showing the probability that a
sphere of a given radius will be empty - in this case, empty of galaxies. We
calculate the VPF for galaxies in the CMASS and LOWZ surveys of BOSS and
compare it with a VPF calculated from a HOD-populated dark matter simulation.

The Halo Occupation Distribution is a model for nonlinear galaxy bias that
assigns galaxies to dark matter halos based solely on the mass of the halo.
By fitting the HOD parameters to have a mock galaxy catalogue that agrees with
the observed galaxy number density and 2-point projected correlation function,
we show the consistency of the HOD with the BOSS data via comparison of VPFs.


Code
----

The code in this repository is designed to work as a standalone automated
reproduction of the exact methodology of our experiment. Simply by running the
script "run_expt.py" you should reproduce all the results from our paper.

Before this is possible, a few configurations may be necessary. Firstly, all
the necessary python packages will need to be available. An easy way to ensure
this is through the command

```
pip install -r requirements.txt
```

which should install the nonstandard python package this experiment uses. If
you need to make a local install, consult the pip documentation here:

http://www.pip-installer.org/en/latest/cookbook.html

As for having all of the necessary data in the right place, you will need a lot
of free space to store the necessary files, in particular the simulation output
and mock galaxy catalogues. Most of the observational data is downloaded
automatically by the main script but for custom directory organisation in the
case of server space restrictions, you will have to symlink accordingly. These
details are treated in the DATA.doc readme

TODO

For access to as of yet proprietary data, you will have to contact the code
authors.
have to symlink some directories

