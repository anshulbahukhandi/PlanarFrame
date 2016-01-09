#include <cmath>
#include <iostream>
#include "CCircsolid.h"

CCircSolid::CCircSolid (const CVector<float>& fV) 
                     : CXSType (numCircDimensions)
// ---------------------------------------------------------------------------
// Function: overloaded ctor
// Input:    vector with rectangular solid dimensions
// Output:   none
// ---------------------------------------------------------------------------
{
    assert (fV.GetSize() >= numCircDimensions);
    for (int i=1; i <= numCircDimensions; i++)
        m_fVDimensions(i) = fV(i);
    ComputeProperties ();
	m_szID="circs";
}

CCircSolid::~CCircSolid ()
// ---------------------------------------------------------------------------
// Function: dtor
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
}

void CCircSolid::ComputeProperties ()
// ---------------------------------------------------------------------------
// Function: computes the Circular solid properties
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    // height
    float fR = m_fVDimensions(1); 
     // cross-sectional area
    m_fArea = static_cast<float>(3.1415*fR*fR);
    // MOI y-axis
    m_fIyy = static_cast<float>(pow(fR, 4.0f)*3.1415/4.0f);
    // MOI z-axis
    m_fIyy = static_cast<float>(pow(fR, 4.0f)*3.1415/4.0f);
}