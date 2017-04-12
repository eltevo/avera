# N-body code for Average Expansion Rate Approximation (AvERA) cosmology

Authors: Gábor Rácz:1, László Dobos:1, István Csabai:1

1:Department of Physics of Complex Systems, Eötvös Loránd University, Pf. 32, H-1518 Budapest, Hungary

This is the source code used in the following paper:

# Concordance cosmology without dark energy

Gábor Rácz:1, László Dobos:1, Róbert Beck:1, István Szapudi:2 and István Csabai:1

1: Department of Physics of Complex Systems, Eötvös Loránd University, Pf. 32, H-1518 Budapest, Hungary
2: Institute for Astronomy, University of Hawaii, 2680 Woodlawn Drive, Honolulu, HI 96822, USA

Abstract

According to the separate universe conjecture, spherically symmetric sub-regions in an
isotropic universe behave like mini-universes with their own cosmological parameters. This is
an excellent approximation in both Newtonian and general relativistic theories. We estimate
local expansion rates for a large number of such regions, and use a scale parameter calculated
from the volume-averaged increments of local scale parameters at each time step in an otherwise
standard cosmological N-body simulation. The particle mass, corresponding to a coarse
graining scale, is an adjustable parameter. This mean field approximation neglects tidal forces
and boundary effects, but it is the first step towards a non-perturbative statistical estimation
of the effect of non-linear evolution of structure on the expansion rate. Using our algorithm,
a simulation with an initial Omega_m = 1 Einstein–de Sitter setting closely tracks the expansion
and structure growth history of the Lambda cold dark matter (LCDM) cosmology. Due to small but
characteristic differences, our model can be distinguished from the LCDM model by future
precision observations. Moreover, our model can resolve the emerging tension between local
Hubble constant measurements and the Planck best-fitting cosmology. Further improvements
to the simulation are necessary to investigate light propagation and confirm full consistency
with cosmic microwave background observations.

Mon Not R Astron Soc Lett slx026.
DOI: 10.1093/mnrasl/slx026
https://academic.oup.com/mnrasl/article/2982870/Concordance

Paper web site: http://www.vo.elte.hu/papers/2016/avera

This work was supported by NKFI NN 114560. IS acknowledges
NASA grants NNX12AF83G and NNX10AD53G for support. RB
was supported through the New National Excellence Program of
the Ministry of Human Capacities, Hungary.

---

Gábor Rácz, 2016
Department of Physics of Complex Systems, Eötvös Loránd University
ragraat@caesar.elte.hu

Cosmological simulation code with backreaction.
- written in C++
- parallelized with openmp and CUDA
- read Gadget2 and ascii IC formats
- output in ascii format
- Can handle standard homogenious and statistical backreaction cosmological simulations
- in this early version the code does not make difference between baryonic and dark matter (there is no baryonic effects now)

---

To compile the code, first install the following libraries:
-voro++ (http://math.lbl.gov/voro++/download/)
-DTFE (http://www.astro.rug.nl/~voronoi/DTFE/download.html)

After installing these, you should modify the Makefile, so the compiler can find the necessary libraries. After the modification, simply type:

	make

---

Once you compiled the code, you can simply run it by typing:

	./CCLEA <parameterfile>

where the parameterfile specifies the parameters of the simulation (see the EXAMPLE file for more details)

If you want to run the code, for example in 8 processors, type:

	export OMP_NUM_THREADS=8
	./CCLEA <parameterfile>

---

Output format for the Logfile:
Logfile.dat:
	t[Gy] error h[Gy] a z H[km/s/Mpc] q Omega_m_eff

Output format for the particle data files:
t*.dat:
	x[Mpc]	y[Mpc]	z[Mpc]	v_x[0.0482190732645461km/s] v_y[0.0482190732645461km/s] v_z[0.0482190732645461km/s]

(The code uses internal units. The distance unit is in Mpc, and the velocity unit is 0.0482190732645461km/s)

---

The source code provided here is for reference purposes.
