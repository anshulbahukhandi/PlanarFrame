/*********************************************
Program Planar Frame
Copyright(c) 2000-08, S. D. Rajan
All rights reserved

Introduction to Structural Analysis and Design
Object-Oriented Numerical Analysis
*********************************************/
#ifndef __RAJAN_ELEMENTRESPONSE_H__
#define __RAJAN_ELEMENTRESPONSE_H__

#include "F:\ASU Courses\cee532\Library\Library\matrixtemplate.h"
#include "constants.h"

class CElementResponse
{
    public:
        CElementResponse ();    // default ctor
        ~CElementResponse ();   // dtor

        // accessor functions
        void GetValues (CMatrix<double>&) const;
		void GetStress(double& ,double&,double&);
        // modifier functions
        void SetValues (const CMatrix<double>&);
		void SetStress(const double&, const double&,const double&);
		
void SetESupportReaction(const CMatrix<double>&);
void GetESupportReaction(CMatrix<double>&);
private:
        CMatrix<double> m_fVFElement;					 // Forces in an Element
        double MaxTensile;
		double MaxComp;
		double MaxShear;
		CMatrix<double> Esupportreaction;			//Element support Reactions
};

#endif