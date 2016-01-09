/*********************************************
Program Planar Frame
Copyright(c) 2000-08, S. D. Rajan
All rights reserved

Introduction to Structural Analysis and Design
Object-Oriented Numerical Analysis

Implementation of the CElementLoads class.

*********************************************/
#include "ElementLoads.h"

CElementLoads::CElementLoads ()
// ---------------------------------------------------------------------------
// Function: default constructor
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    m_nElement = 0;
    m_fValue1 = m_fValue2 = 0.0f;
	m_Type= LinearlyDistributed;
}

CElementLoads::CElementLoads (const CElementLoads& ELO)				//NOT GETTING..why do we need this here ?
// ---------------------------------------------------------------------------
// Function: copy constructor
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    m_nElement = ELO.m_nElement;
    m_Type = ELO.m_Type;
    m_fValue1 = ELO.m_fValue1;
    m_fValue2 = ELO.m_fValue2;
}

CElementLoads::~CElementLoads ()
// ---------------------------------------------------------------------------
// Function: destructor
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
}

void CElementLoads::GetValues (int& nE, ELType& type,
                               float& f1, float& f2) const
// ---------------------------------------------------------------------------
// Function: gets the element load information
// Input:    element #, load type, 2 load associated values
// Output:   values for all these variables
// ---------------------------------------------------------------------------
{
    nE = m_nElement;
    type = m_Type;
    f1 = m_fValue1;
    f2 = m_fValue2;
}

void CElementLoads::SetValues (const int nE, const ELType type,
                               const float f1, const float f2)
// ---------------------------------------------------------------------------
// Function: sets the element load information
// Input:    element #, load type, 2 load associated values
// Output:   none
// ---------------------------------------------------------------------------
{
    m_nElement = nE;
	m_Type=type;
    m_fValue1 = f1;
    m_fValue2 = f2;
}

void CElementLoads::GetENF (CVector<float>& fV1, CVector<float>& fV2,
                            const float fLength) const
// ---------------------------------------------------------------------------
// Function: gets the equivalent nodal forces
// Input:    vectors to store ENF at start node and end node, element length
// Output:   the two vectors suitably populated
// ---------------------------------------------------------------------------
{
}

std::ostream &operator<< (std::ostream& ofs, const CElementLoads::ELType& ELT)
// ---------------------------------------------------------------------------
// Function: non-class function to display the Load Type as a text string
// Input:    output stream, Load Type
// Output:   none
// ---------------------------------------------------------------------------
{
	if (ELT ==CElementLoads::LinearlyDistributed)
        ofs << "Linearly Distributed";
	if (ELT == CElementLoads::ConcentratedX)
        ofs << "Concenterated X";
	if (ELT == CElementLoads::ConcentratedY)
        ofs << "Concenterated Y";
	if (ELT == CElementLoads::Moment)
        ofs << "Concenterated Moment";
return ofs;
}