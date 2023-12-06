**What the code is doing**
	preliminary: I compute the probabilities for N particles in an event (with simplification that Nmax=5)

For 1 event: 
	-I generate N particles, in random position (X, Y) according to the beam profile (Umberto's results)
	-Particles only have Pz componant (simplification)
	-Using the position (X,Y), I take the efficiency of each 2s Modules at this position (Riccardo's results).
	-check if the particle give a track in each of the stations (2X, 2Y and at least 1U/V stubs)
	-try match tracks in stations 1 and tracks in stations 2


**Run the code:**
to build the library:
	$ make
to run macros:
	cd macros/
	root efficiency.cpp
	(You have to be in macros folder to use the macro, it won't use rootlogon.C otherwise)