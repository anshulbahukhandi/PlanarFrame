/*********************************************
Program Planar Frame
Copyright(c) 2000-08, S. D. Rajan
All rights reserved

Introduction to Structural Analysis and Design
Object-Oriented Numerical Analysis

Base class for element cross-sectional properties
*********************************************/
#include <cassert>
#include "xstype.h"

CXSType::CXSType ()
// ---------------------------------------------------------------------------
// Function: default ctor
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    Initialize ();
}

CXSType::CXSType (int numDimensions)
// ---------------------------------------------------------------------------
// Function: overloaded ctor
// Input:    # of cross-sectional dimensions
// Output:   none
// ---------------------------------------------------------------------------
{
    Initialize ();
    m_numDimensions = numDimensions;
    m_fVDimensions.SetSize (m_numDimensions);
	
}

void CXSType::Initialize ()
// ---------------------------------------------------------------------------
// Function: initializes all the member variables with default values
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    m_numDimensions = 0;
    m_fArea = 0.0f;
    m_fIyy = 0.0f;
    m_fIzz = 0.0f;
}

CXSType::~CXSType ()
// ---------------------------------------------------------------------------
// Function: dtor
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
}
//To set the shape of the elements
//void CXSType::SetShape(const std::string shape)
//{
//	m_szID=shape;
//}
//


//To get the hape of the cross section 
//Returns the cross section shape of the element
void CXSType::GetShape(std::string& sz)
{
	sz=m_szID;
}

void CXSType::GetProperties (float& fArea, float& fIyy, float& fIzz)
// ---------------------------------------------------------------------------
// Function: gets the computed cross-sectional properties
// Input:    variables to hold area, Iyy, Izz
// Output:   area, Iyy, Izz values
// ---------------------------------------------------------------------------
{
    ComputeProperties ();
    fArea = m_fArea;
    fIyy = m_fIyy;
    fIzz = m_fIzz;
}

void CXSType::GetDimensions (CVector<double>& fV) const
// ---------------------------------------------------------------------------
// Function: gets the cross-sectional dimensions
// Input:    vector to hold x/s dimensions
// Output:   x/s dimensions
// ---------------------------------------------------------------------------
{
    assert (fV.GetSize() >= m_numDimensions);
    for (int i=1; i <= m_numDimensions; i++)
        fV(i) = m_fVDimensions(i);
}
