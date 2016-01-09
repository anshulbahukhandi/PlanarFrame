/*********************************************
Program Planar Frame
Copyright(c) 2000-08, S. D. Rajan
All rights reserved

Introduction to Structural Analysis and Design
Object-Oriented Numerical Analysis

Implementation of the CElementResponse class.

*********************************************/
#include "ElementResponse.h"

CElementResponse::CElementResponse ()
// ---------------------------------------------------------------------------
// Function: default constructor
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    m_fVFElement.SetSize (DOFPE,1);  
    m_fVFElement.Set(0.0f);
	Esupportreaction.SetSize(DOFPE,1);
	Esupportreaction.Set(0);
}

CElementResponse::~CElementResponse ()
// ---------------------------------------------------------------------------
// Function: destructor
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
}

void CElementResponse::GetValues (CMatrix<double>& fVFE) const
// ---------------------------------------------------------------------------
// Function: gets the element force response
// Input:    vectors to hold element forces at start and end nodes
// Output:   the element force values
// ---------------------------------------------------------------------------
{
    for (int i=1; i <= DOFPE; i++)
    {
        fVFE(i,1) = m_fVFElement(i,1);
        
    }
}

void CElementResponse::SetValues (const CMatrix<double>& fVFE)
// ---------------------------------------------------------------------------
// Function: gets the element force response
// Input:    vectors to hold element forces at start and end nodes
// Output:   none
// ---------------------------------------------------------------------------
{
    for (int i=1; i <= DOFPE; i++)
    {  m_fVFElement(i,1) = fVFE(i,1);
    }
}

//Gets The Maximum tensile,compressive and shear stress in the element
void CElementResponse:: GetStress(double& maxt ,double& maxc,double& maxs)
{
	maxt=MaxTensile;
	maxc=MaxComp;
	maxs=MaxShear;
}
// Sets the maximum tensile , compressive and shear stress
void CElementResponse:: SetStress(const double& maxt ,const double& maxc,const double& maxs)
{
	MaxComp=maxc;
	MaxTensile=maxt;
	MaxShear=maxs;
}
void CElementResponse::SetESupportReaction(const CMatrix<double>& SR)
{
	for(int i=1;i<=DOFPE;i++)
	{
		Esupportreaction(i,1)=SR(i,1);
	}

}
//get support reactions of the elements
void CElementResponse::GetESupportReaction(CMatrix<double>& SR)
{for(int i=1;i<=DOFPE;i++)
{
	SR(i,1)=Esupportreaction(i,1);
}
}