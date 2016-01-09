/*********************************************
Program Planar Frame
Copyright(c) 2000-08, S. D. Rajan
All rights reserved

Introduction to Structural Analysis and Design
Object-Oriented Numerical Analysis

Implementation of the CNodalResponse class.

*********************************************/
#include "NodalResponse.h"

CNodalResponse::CNodalResponse ()
// ---------------------------------------------------------------------------
// Function: default constructor
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    m_fVDisp.SetSize (DOFPN);
}

CNodalResponse::~CNodalResponse ()
// ---------------------------------------------------------------------------
// Function: destructor
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
}

void CNodalResponse::GetValues (CVector<float>& fV) const
// ---------------------------------------------------------------------------
// Function: gets the nodal displacements at the current node
// Input:    vector to hold the nodal displacement values
// Output:   nodal displacement values
// ---------------------------------------------------------------------------
{
    for (int i=1; i <= DOFPN; i++)
        fV(i) = m_fVDisp(i);
}

void CNodalResponse::SetValues (const CVector<float>& fV)
// ---------------------------------------------------------------------------
// Function: sets the nodal displacements at the current node
// Input:    vector holding the nodal displacement values
// Output:   none
// ---------------------------------------------------------------------------
{
    for (int i=1; i <= DOFPN; i++)
        m_fVDisp(i) = fV(i);
}

void CNodalResponse::SetValues (const CVector<double>& fV)
// ---------------------------------------------------------------------------
// Function: sets the nodal displacements at the current node
// Input:    vector holding the nodal displacement values
// Output:   none
// ---------------------------------------------------------------------------
{
    for (int i=1; i <= DOFPN; i++)
        m_fVDisp(i) = static_cast<float>(fV(i));
}

