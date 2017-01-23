#!/bin/csh
#------------------------------------------------------------------------------
# Version 9-November-2000
#------------------------------------------------------------------------------
# To install xnirspec properly, you must follow these steps:
# make include
# make xnirspec
# make clean
#------------------------------------------------------------------------------
# define here maximum image dimensions
# NXMAX must be the image size in the wavelength direction
# NYMAX must be the image size in the spatial direction
# Note1: xnirspec can rotate the image after reading it (see SLEEFITS)
# Note2: NXYMAX is the largest of NXMAX and NYMAX
# Note3: NXMAXB9.LE.NXMAX and NYMAXB9.LE.NYMAX
# --------------
NXMAX   = 4300
NYMAX   = 4300
NXYMAX  = 4300
NXMAXB9 =  900
NYMAXB9 =  900
NOVERSAMPMAX = 10
NZMAX   = 10000
#------------------------------------------------------------------------------
# include here the whole path to SExtractor
#SEXPATH = /home/ncl/SExtractor/sextractor2.1.6/source
#SEXPATH = /home/scisoft/sextractor/sextractor2.1.6/source
SEXPATH= /Users/cardiel/s/sextractor/sextractor-2.5.0/src
#------------------------------------------------------------------------------
# verify that the required libraries are properly set
#PGPDIR  = /usr/local/pgplot
PGPDIR  = /opt/local/lib
#------------------------------------------------------------------------------
#X11DIR  = /usr/openwin/lib
#X11DIR  = /usr/X11R6/lib
X11DIR  = /opt/local/lib
#------------------------------------------------------------------------------
#FIODIR  = /usr/local/cfitsio
FIODIR  = /opt/local/lib
#FIODIR  = /home/cardiel/w/WinterSchool2005/CDROM/software/galfit/cfitsio
# Set the appropiate compilers (note: in hana.keck.hawaii.edu it is necessary
# to use g77 instead of cc; otherwise it is not possible to mix fortran and C
# code).
#------------------------------------------------------------------------------
#FCOMPIL = f77 -g -ftrap=common,underflow
#FCOMPIL = gfortran -O3 -g -Wall
#FCOMPIL = g95 -O2 -g -Wall
#FCOMPIL = gfortran-mp-4.8 -O2 -g -Wall
FCOMPIL = gfortran-mp-5 -O2 -g -Wall
#------------------------------------------------------------------------------
# Nothing SHOULD be modified below this comment line
#------------------------------------------------------------------------------
NMAXBUFF = 12
# macro definitions
FSOURCE = button.f buttqbr.f buttqcf.f buttqch.f buttqex.f \
          buttqit.f buttqpr.f buttqxb.f buttqyb.f buttqytext.f buttsbr.f \
          buttscf.f buttsch.f buttsex.f buttsit.f buttspr.f buttsxb.f \
          buttsyb.f buttsytext.f ifbutton.f rpgband.f rpgbegin.f rpgbegok.f \
          rpgenv.f rpgeras.f rpgerasb.f rpgerasw.f \
          activebut.f \
          arclength.f \
          binsearch.f \
          boundary.f \
          centroid.f \
          combine.f \
          combpf.f \
          creatergb.f \
          cubspl.f \
          cubsplx.f \
          distpp.f \
          distpr.f \
          downhill.f \
          drawpolxy.f \
          drawpolyx.f \
          expline.f \
          factorialpf.f \
          fderi.f \
          findmax.f \
          findmml.f \
          fitpol.f \
          fmap.f \
          fmean0.f \
          fmean0e.f \
          fmean0w.f \
          fmean1.f \
          fmean2.f \
          fmean2e.f \
          fmedian1.f \
          fmedian1e.f \
          fpoly.f \
          fpolyinv.f \
          fstatistic.f \
          gauscfit.f \
          gauss2dfit.f \
          gaussfit.f \
          gaussfita.f \
          guessef.f \
          histogram.f \
          hplot.f \
          intersec.f \
          iofunctions.f \
          leecolumn.f \
          leeonefile.f \
          lrombo.f \
          ludcmp.f \
          ludcmpd.f \
          lusolv.f \
          lusolvd.f \
          mapping.f \
          maxdimb9.f \
          mideoff.f \
          mideoffset.f \
          operbox9.f \
          ordena1f.f \
          ordena1f1i.f \
          ordena3f3i.f \
          palette.f \
          pmapp.f \
          pbound.f \
          plot3dbars.f \
          polfit.f \
          polfitsig.f \
          poly2.f \
          pcut.f \
          pwave.f \
          pxmap.f \
          randomnumber.f \
          readpwave.f \
          rectify.f \
          refineba.f \
          refinebl.f \
          sescrfits.f \
          sescrredu.f \
          sextractor.f \
          skysub.f \
          srcpat.f \
          sfitpol.f \
          shiftb9.f \
          show_fitpol.f \
          showperc.f \
          shiftimage.f \
          sleefits.f \
          slit_aperture.f \
          sortlines.f \
          splfit.f \
          statistic.f \
          statisticb9.f \
          subimath.f \
          sublook.f \
          subplot.f \
          subplotbis.f \
          subzoom.f \
          systemfunction.f \
          try_reduceme.f \
          unrectify.f \
          whereami.f \
          x2arc.f \
          xcut.f \
          ycut.f \
          zoomhistog.f \
          xextract.f \
          xnirspec.f \
          zscale.f
FOBJECT = $(FSOURCE:.f=.o)
# Default rule to create program
xnirspec:  $(FOBJECT)
#	$(FCOMPIL) -o $@ $(FOBJECT) -L$(PGPDIR) -L$(FIODIR) -L$(X11DIR) -lpgplot -lcfitsio -lX11 -lnsl -lsocket
	$(FCOMPIL) -o $@ $(FOBJECT) -L$(PGPDIR) -L$(FIODIR) -L$(X11DIR) -lpgplot -lcfitsio -lX11
# Target to clean object modules
clean:    $(FOBJECT)
	rm -f $(FOBJECT)
	rm -f dimensions.inc
	rm -f largest.inc
	rm -f sexpath.inc
# Target to touch source modules
touch:
	touch $(FSOURCE)
# Target to create the file dimensions.inc
include:
	rm -f dimensions.inc
	echo "        INTEGER NXMAX" > dimensions.inc
	echo "        PARAMETER (NXMAX=$(NXMAX))" >> dimensions.inc
	echo "        INTEGER NYMAX" >> dimensions.inc
	echo "        PARAMETER (NYMAX=$(NYMAX))" >> dimensions.inc
	echo "        INTEGER NMAXBUFF" >> dimensions.inc
	echo "        PARAMETER (NMAXBUFF=$(NMAXBUFF))" >> dimensions.inc
	echo "        INTEGER NXMAXB9" >> dimensions.inc
	echo "        PARAMETER (NXMAXB9=$(NXMAXB9))" >> dimensions.inc
	echo "        INTEGER NYMAXB9" >> dimensions.inc
	echo "        PARAMETER (NYMAXB9=$(NYMAXB9))" >> dimensions.inc
	echo "        INTEGER NZMAX" >> dimensions.inc
	echo "        PARAMETER (NZMAX=$(NZMAX))" >> dimensions.inc
	rm -f largest.inc
	echo "        INTEGER NXYMAX" > largest.inc
	echo "        PARAMETER (NXYMAX=$(NXYMAX))" >> largest.inc
	echo "        INTEGER NOVERSAMPMAX" >> largest.inc
	echo "        PARAMETER (NOVERSAMPMAX=$(NOVERSAMPMAX))" >> largest.inc
	rm -f sexpath.inc
	echo "        CHARACTER*255 SEXPATH" > sexpath.inc
	echo "        PARAMETER (SEXPATH=" >> sexpath.inc
	echo "     +   '$(SEXPATH)')" >> sexpath.inc
	touch $(FSOURCE)
# second level dependencies
.f.o: $(FSOURCE)
	$(FCOMPIL) -c $?
# definitions
.PRECIOUS: xnirspec
