# ICVL_Fortran
FCC QD heterostructures.
The Fortran code designs the position of atoms in an un-relaxed quantum dot (QD) standing on a wetting layer (WL) . The QD+WL ensemble is embedded in a matrix.  Both QD+WL and matrix have FCC symmetry. Two QD shapes are considered, cone and hemi-torus.

Parameters:

	a1, a2 are the lattice constants.
	ratio is used in layer=a2*ratio=layer_thicknes.
	a1, a2, ratio may be changed it Parameters.inp
	Nw is #mono-layers. 
	nxmin, nxmax, nymin, nymax, nzmin, nzmax are #atoms on x,y,z axes.
	nxmin is  #atoms on negative part of x axis.
	nxmax is  #atoms on the positive part of x axis; similar for nymin, nymax, nzmin, nzmax.

For cone:
	Rc is the base radius of the regular cone and h is the cone height.
For HT:
	Rt is the distance between the center of the tube and the center of the torus, Rq is the radius of the tube.

The parameters Rc, h, Rt, Rq, Nw, nxmin, nxmax, nymin, nymax, nzmin, nzmax may be changed by typing their values in QD.f90 in the sector:
	if (QD==1)then …  end if
