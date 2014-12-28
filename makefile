# define here your compiler:
#general use:

compile = gcc -O4 -Wall  


#-----------

all: 	pips tools 

pips: 	begin circ_ap rect_ap file_int file_pha forvard tilt gauss \
b_split b_mix lens absorber normal cros_out lens_forvard \
interpol interp1 convert circ_screen rect_screen Strehl \
file_ps   Zernike forward    Fresnel l_amplif \
lens_fresn pip_fft gauss_screen fil_ter random \
 tor_lens steps  

tools: 	file_pgm unf3 unf4 a_phase negate mask c_matlab c_xplot 4_int gamma scale

small: 	begin rect_ap file_ps gauss 

clean: 
	rm begin circ_ap rect_ap file_int file_pha \
forvard tilt gauss b_split b_mix lens absorber normal \
cros_out lens_forvard interpol interp1 convert circ_screen \
rect_screen Strehl file_ps   Zernike forward   \
Fresnel l_amplif lens_fresn pip_fft gauss_screen fil_ter random \
file_pgm unf3 unf4 a_phase negate tor_lens steps mask  c_matlab c_xplot 4_int\
scale gamma

begin:  begin.c pipes.h 
	$(compile) -o begin  begin.c -lm

circ_ap: circ_ap.c pipes.h
	$(compile) -o circ_ap circ_ap.c -lm

rect_ap: rect_ap.c pipes.h
	$(compile) -o rect_ap rect_ap.c -lm

file_int: file_int.c pipes.h
	$(compile) -o file_int  file_int.c -lm

file_pha: file_pha.c pipes.h
	$(compile) -o file_pha  file_pha.c -lm

forvard: forvard.c fftn.c pipes.h fft_prop.c
	$(compile) -o forvard  forvard.c fftn.c -lm

tilt: tilt.c pipes.h
	$(compile) -o tilt  tilt.c -lm

gauss: gauss.c pipes.h
	$(compile) -o gauss gauss.c -lm

b_split: b_split.c pipes.h
	$(compile) -o b_split b_split.c -lm

b_mix: b_mix.c pipes.h
	$(compile) -o b_mix b_mix.c -lm

lens: lens.c pipes.h
	$(compile) -o lens lens.c -lm

absorber: absorber.c pipes.h
	$(compile) -o absorber absorber.c -lm

normal: normal.c pipes.h
	$(compile) -o normal normal.c -lm

cros_out: cros_out.c pipes.h
	$(compile) -o cros_out cros_out.c -lm

lens_forvard: lens_forvard.c fftn.c pipes.h fft_prop.c
	$(compile) -o lens_forvard lens_forvard.c fftn.c -lm

interpol: interpol.c pipes.h
	$(compile) -o interpol interpol.c -lm

interp1: interp1.c pipes.h
	$(compile) -o interp1 interp1.c -lm

convert: convert.c pipes.h
	$(compile) -o convert convert.c -lm

circ_screen: circ_screen.c pipes.h
	$(compile) -o circ_screen circ_screen.c -lm

rect_screen: rect_screen.c pipes.h
	$(compile) -o rect_screen rect_screen.c -lm

Strehl: strehl.c pipes.h
	$(compile) -o Strehl strehl.c -lm

file_ps: file_ps.c pipes.h
	$(compile) -o file_ps file_ps.c -lm


Zernike: zernike.c pipes.h
	$(compile) -o Zernike zernike.c -lm

forward: forward.c pipes.h
	$(compile) -o forward forward.c -lm


Fresnel: fresnel.c frsn_prop.c pipes.h fftn.c
	$(compile) -o Fresnel fresnel.c  fftn.c -lm

l_amplif: l_amplif.c pipes.h
	$(compile) -o l_amplif l_amplif.c -lm

lens_fresn: lens_fresn.c frsn_prop.c pipes.h fftn.c
	$(compile) -o lens_fresn lens_fresn.c  fftn.c -lm


pip_fft: pip_fft.c fftn.c pipes.h
	$(compile) -o pip_fft  pip_fft.c  fftn.c -lm

gauss_screen: gauss_screen.c pipes.h
	$(compile) -o gauss_screen  gauss_screen.c  -lm

fil_ter:	fil_ter.c pipes.h
	$(compile) -o fil_ter fil_ter.c -lm

random:		random.c pipes.h
	$(compile) -o random random.c -lm

file_pgm:	file_pgm.c pipes.h
	$(compile) -o file_pgm file_pgm.c -lm

unf3:	unf3.c
	$(compile) -o unf3 unf3.c -lm

unf4:	unf4.c
	$(compile) -o unf4 unf4.c -lm

a_phase:	a_phase.c
	$(compile) -o a_phase a_phase.c -lm	

negate:		negate.c
	$(compile) -o negate negate.c -lm

tor_lens:	tor_lens.c pipes.h
	$(compile) -o tor_lens tor_lens.c -lm

steps:		steps.c pipes.h
	$(compile) -o steps steps.c -lm


mask:	mask.c 
	$(compile)   -o mask mask.c -lm

c_matlab:	c_matlab.c 
	$(compile)   -o c_matlab c_matlab.c -lm

c_xplot:	c_xplot.c 
	$(compile)   -o c_xplot c_xplot.c -lm

4_int:		4_int.c
	$(compile)   -o 4_int 4_int.c -lm

gamma:		gamma.c
	$(compile) -o gamma gamma.c -lm

scale:  	scale.c
	$(compile) -o scale scale.c -lm
