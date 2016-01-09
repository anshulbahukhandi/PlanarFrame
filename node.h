/*********************************************
Program Planar Frame
Copyright(c) 2000-08, S. D. Rajan
All rights reserved

Introduction to Structural Analysis and Design
Object-Oriented Numerical Analysis
*********************************************/
#ifndef __RAJAN_NODE_H__
#define __RAJAN_NODE_H__

#include <iostream>
#include "F:\ASU Courses\cee532\Library\Library\vectortemplate.h"
#include "constants.h"

class CNode
{
    public:
        enum Fixity {FREE, SPECIFIED};
        friend std::ostream &operator<< (std::ostream&, 
                                         const CNode::Fixity&);			//NOT GETTING
        CNode ();             // ctor
        CNode (const CNode&); // copy ctor
        ~CNode ();            // dtor

        // accessor functions
        void GetCoords (CVector<float>& fVC) const;
        void GetFixity (CVector<Fixity>& VFC) const;
		void GetDisp (CVector<float>& VD)const;								//changed this						
        // modifier functions
        void SetCoords (const CVector<float>& fVC);
        void SetFixity (const CVector<Fixity>& VFC);
		void SetDisp (const CVector<float>& VD);

    private:
        CVector<float>	m_fVCoor;	// coordinates
        CVector<Fixity>	m_VFC;		// fixity conditions
        CVector<float>	m_fVDisp;	// known displacements
		
};

#endif	