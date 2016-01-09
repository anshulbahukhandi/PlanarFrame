/*********************************************
Program Planar Frame
Copyright(c) 2000-08, S. D. Rajan
All rights reserved

Introduction to Structural Analysis and Design
Object-Oriented Numerical Analysis

Derived class from CXSType base class
to store I-section data and properties
*********************************************/
#include <cmath>
#include <iostream>
#include "isection.h"

CISection::CISection (const CVector<float>& fV) : CXSType (numISDimensions)
// ---------------------------------------------------------------------------
// Function: overloaded ctor
// Input:    vector with I-section dimensions
// Output:   none
// ---------------------------------------------------------------------------
{
    assert (fV.GetSize() >= numISDimensions);
    for (int i=1; i <= numISDimensions; i++)
        m_fVDimensions(i) = fV(i);
    ComputeProperties ();
	m_szID="isection";
}

CISection::~CISection ()
// ---------------------------------------------------------------------------
// Function: dtor
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
}

void CISection::ComputeProperties ()
// ---------------------------------------------------------------------------
// Function: computes the I-section properties
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    // web height
    float wH = m_fVDimensions(1); 
    // web thickness
    float wT = m_fVDimensions(2); 
    // flange width
    float fW = m_fVDimensions(3); 
    // flange thickness
    float fT = m_fVDimensions(4); 

    // cross-sectional area
    m_fArea = wH*wT + 2.0f*fW*fT;
    // MOI y-axis
    float fOI, fII;
    fOI = (pow(wH+2.0f*fT, 3.0f) * fW)/12.0f;
    fII = 2.0f*(pow(wH, 3.0f) * (0.5f*(fW-wT)))/12.0f;
    m_fIyy = fOI - fII;
    // MOI z-axis
    float fIzz = 2.0f*(pow(fW,3.0f)*fT)/12.0f + 
                 pow(wT,3)*wH/12.0f;
    m_fIzz = fIzz;
}