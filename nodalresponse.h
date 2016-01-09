/*********************************************
Program Planar Frame
Copyright(c) 2000-08, S. D. Rajan
All rights reserved

Introduction to Structural Analysis and Design
Object-Oriented Numerical Analysis
*********************************************/
#ifndef __RAJAN_NODALRESPONSE_H__
#define __RAJAN_NODALRESPONSE_H__

#include "F:\ASU Courses\cee532\Library\Library\vectortemplate.h"
#include "F:\ASU Courses\cee532\Library\Library\matrixtemplate.h"
#include "constants.h"

class CNodalResponse
{
    public:
        CNodalResponse ();   // ctor
        ~CNodalResponse ();  // dtor

        // accessor functions
        void GetValues (CVector<float>& fVDisp) const;
		
        // modifier functions
        void SetValues (const CVector<float>& fVDisp);
		void SetValues (const CVector<double>& fVDisp);
            
private:
        CVector<float> m_fVDisp; // nodal displacements	
		CMatrix<double>NSupportReac; // support reactions at each node
};

#endif