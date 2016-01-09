/*********************************************
Program Planar Frame
Copyright(c) 2000-08, S. D. Rajan
All rights reserved

Introduction to Structural Analysis and Design
Object-Oriented Numerical Analysis
*********************************************/
#ifndef __RAJAN_RECTSOLID_H__
#define __RAJAN_RECTSOLID_H__

#include "F:\ASU Courses\cee532\Library\Library\vectortemplate.h"
#include "xstype.h"
const int numRectDimensions = 2;

class CRectSolid: public CXSType
{
    public:
        CRectSolid (const CVector<float>& fV);
        CRectSolid (const CRectSolid&);
        ~CRectSolid ();

        // helper functions
        virtual void ComputeProperties ();

    private:
};

#endif