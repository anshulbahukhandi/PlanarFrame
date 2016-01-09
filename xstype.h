/*********************************************
Program Planar Frame
Copyright(c) 2000-08, S. D. Rajan
All rights reserved

Introduction to Structural Analysis and Design
Object-Oriented Numerical Analysis
*********************************************/
#ifndef __RAJAN__CXSTYPE__H__
#define __RAJAN__CXSTYPE__H__

#include <string>
#include "F:\ASU Courses\cee532\Library\Library\vectortemplate.h"

class CXSType
{
    public:
        CXSType ();
        CXSType (int);
        virtual ~CXSType ();

        // helper function
        void DisplayProperties ();

        // accessor functions
        void GetProperties (float&, float&, float&);
        void GetDimensions (CVector<double>&) const;
        virtual void ComputeProperties () = 0;
		void GetShape(std::string&);	//changed this
		//void SetShape(const std::string);		//changed this
    
private:
        void Initialize ();

    protected:
        std::string m_szID;            // identification tag
        float m_fArea;                 // x/s area
        float m_fIyy;                  // MOI y-axis
        float m_fIzz;                  // MOI z-axis
        int   m_numDimensions;         // number of dimensions
		CVector<float> m_fVDimensions; // the dimensions
};

#endif