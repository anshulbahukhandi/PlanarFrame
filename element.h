/*********************************************
Program Planar Frame
Copyright(c) 2000-08, S. D. Rajan
All rights reserved

Introduction to Structural Analysis and Design
Object-Oriented Numerical Analysis
*********************************************/
#ifndef __RAJAN_ELEMENT_H__
#define __RAJAN_ELEMENT_H__

#include "xstype.h"
class CElement
{
    public:
        CElement ();                  // default ctor
        CElement(const CElement& EO); // copy ctor
        ~CElement ();                 // dtor

        // accessor functions
        void GetENodes (int& nSN, int& nEN) const;
        int  GetMatPropertyGroup () const;
        CXSType* GetEPropertyGroup () const;
        
        // modifier functions
        void SetENodes (const int nSN, const int nEN);
        void SetEPropertyGroup (CXSType*);
        void SetMatPropertyGroup (const int);

    private:
        int m_nSN;		     // start node
        int m_nEN;		     // end node
        int m_nMPGroup;      // material property group
        CXSType* m_pEPGroup; // element property group
		
};

#endif