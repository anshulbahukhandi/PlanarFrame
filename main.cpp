/*********************************************
Planar Frame Analysis Program
Copyright(c) 2000-13, S. D. Rajan
All rights reserved

Introduction to Structural Analysis and Design
Object-Oriented Numerical Analysis
*********************************************/
#include "frame.h"
#include "F:\ASU Courses\cee532\Library\Library\clock.h"

int main (int argc, char *argv[])
{
    CFrame TheFrame; // the one and only frame!

	// show program banner
    TheFrame.Banner (std::cout);

	// Prepare for I/O
	TheFrame.PrepareIO (argc, argv);
	
    // start the timer --------------------------------------------------------
    CClock Timer;
    std::string szDateTime;
    Timer.GetDateTime (szDateTime);
    std::cout << "\nStarting out at : " << szDateTime << "\n";

	// read the problem size
	TheFrame.ReadProblemSize ();

	// set problem size
	TheFrame.SetSize ();						//******STEP 3*******

	// read nodal and element data
	TheFrame.ReadFrameModel ();

	// construct system equations
	TheFrame.ConstructK ();
	TheFrame.ConstructF ();

	// impose boundary conditions
	TheFrame.ImposeBC ();

	// solve for the nodal displacements
	TheFrame.Solve ();

	// compute secondary unknowns
	TheFrame.Response ();
	//computes errors							changed this
	TheFrame.ResidualError();

	// create output file
	TheFrame.CreateOutput ();						//***STEP 28***		

	// Close input and output files
	TheFrame.TerminateProgram ();

    // get the current date and time
    Timer.GetDateTime (szDateTime);
    std::cout << "\n      Ending at : " << szDateTime << "\n";
    // compute the elapsed time -----------------------------------------------
    std::cout << "Elapsed wall clock time: " << Timer.DiffTime ()
              << " seconds\n";

	return 0;
}