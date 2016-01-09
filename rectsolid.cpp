/*********************************************
Program Planar Frame
Copyright(c) 2000-08, S. D. Rajan
All rights reserved

Introduction to Structural Analysis and Design
Object-Oriented Numerical Analysis

Derived class from CXSType base class
to store rectangular solid data and properties
*********************************************/
#include <cmath>
#include <iostream>
#include "rectsolid.h"

CRectSolid::CRectSolid (const CVector<float>& fV) 
                     : CXSType (numRectDimensions)
// ---------------------------------------------------------------------------
// Function: overloaded ctor
// Input:    vector with rectangular solid dimensions
// Output:   none
// ---------------------------------------------------------------------------
{
    assert (fV.GetSize() >= numRectDimensions);
    for (int i=1; i <= numRectDimensions; i++)
        m_fVDimensions(i) = fV(i);
    ComputeProperties ();
	m_szID="rects";
}

CRectSolid::~CRectSolid ()
// ---------------------------------------------------------------------------
// Function: dtor
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
}

void CRectSolid::ComputeProperties ()
// ---------------------------------------------------------------------------
// Function: computes the rectangular solid properties
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    // height
    float fH = m_fVDimensions(1); 
    // width
    float fW = m_fVDimensions(2); 

    // cross-sectional area
    m_fArea = fH * fW;
    // MOI y-axis
    m_fIzz = pow(fH, 3.0f)*fW/12.0f;
    // MOI z-axis
    m_fIyy = pow(fW, 3.0f)*fH/12.0f;
}