/*********************************************
Program Planar Frame
Copyright(c) 2000-08, S. D. Rajan
All rights reserved

Introduction to Structural Analysis and Design
Object-Oriented Numerical Analysis
*********************************************/
#ifndef __RAJAN_NODALLOADS_H__
#define __RAJAN_NODALLOADS_H__

#include "F:\ASU Courses\cee532\Library\Library\vectortemplate.h"
#include "constants.h"

class CNodalLoads
{
    public:
        CNodalLoads ();                   // default ctor 
        CNodalLoads (const CNodalLoads&); //copy ctor
        ~CNodalLoads ();

        // accessor functions
        void GetValues (CVector<float>&);
		void GetTemp (float&);
        // modifier functions
        void SetValues (const CVector<float>&);
		void SetTemp (const float&);

    private:
        CVector<float> m_fVLoads;	// nodal loads
		float m_ftemperature;   // Nodal temerature						******changed this for regrading ******
};

#endif