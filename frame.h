/*********************************************
Planar Frame Analysis Program
Copyright(c) 2000-08, S. D. Rajan
All rights reserved

Introduction to Structural Analysis and Design
Object-Oriented Numerical Analysis
*********************************************/
#ifndef __RAJAN_FRAME_H__
#define __RAJAN_FRAME_H__

#include <fstream>
#include <iostream>
#include "constants.h"
#include "F:\ASU Courses\cee532\Library\Library\vectortemplate.h"
#include "F:\ASU Courses\cee532\Library\Library\matrixtemplate.h"
#include "node.h"
#include "element.h"
#include "material.h"
#include "xstype.h"
#include "rectsolid.h"
#include "isection.h"
#include "nodalloads.h"
#include "elementloads.h"
#include "nodalresponse.h"
#include "elementresponse.h"
#include "CCircSolid.h"							//CHANGED THIS

class CFrame
{
    public:
        CFrame ();    // ctor
        ~CFrame ();   // dtor
        enum ERRORCODE {NUMNODES, NUMELEMENTS, DEBUGCODE,
                        NODENUMBER, ELEMENTNUMBER, XSAREA,
                        YOUNGSMODULUS, UNSTABLE, INVALIDINPUT,
                        INVALIDLASTLINE, NODALFIXITY, NUMMATGROUPS,
                        NUMEPROPGROUPS, EPROPNUMBER, XSDIMENSION,
                        MATGROUPNUMBER, MATPROPERTY, ELOADTYPE,
                        XSTYPE, INVALIDCOMMANDLINE,
		                CANNOTOPENIFILE, CANNOTOPENOFILE};

        // helper functions
        void Banner (std::ostream& OF);
        void PrepareIO (int argc, char *argv[]);
        void ReadProblemSize ();
        void ReadFrameModel ();
        void ConstructK ();
        void ConstructF ();
        void ImposeBC ();
        void Solve ();
        void Response ();
        void CreateOutput ();
        void TerminateProgram ();
        void SuppressDOF (const int);
		void ResidualError();			//Residual error  changed this

        // modifier functions
        void SetSize ();

    private:
        
		int m_nNodes;		   // number of nodes
        int m_nElements;       // number of elements
        int m_nEPGroups;       // number of x/s properties
        int m_nMatGroups;      // number of material groups
        int m_nElementLoads;   // number of element loads

        int m_nDOF;			   // total degrees-of-freedom
        int m_nDebugLevel;	   // debugging level
        int m_nLineNumber;	   // current line number in input file
		double RError;				//changed this
		double AError;				//changed this

		CMatrix<double> supportreactions;        //matrix to hold support reactions at each node of the frame 
        CVector<CNode>				m_NodalData;
        CVector<CElement>			m_ElementData;
        CVector<CMaterial>			m_MaterialData;
        CVector<CXSType*>			m_EPData;				//NOT GETTING.. how is pointer useful here??
        CVector<CNodalLoads>		m_NodalLoadData;
        CVector<CElementLoads>		m_ElementLoadData;
        CVector<CNodalResponse>		m_NodalResponseData;
        CVector<CElementResponse>	m_ElementResponseData;						//All six forces in an element     changed this
        CMatrix<double> m_SSM;	 // structural stiffness matrix
        CMatrix<double> m_SND;	 // structural nodal displacements
        CMatrix<double> m_SNF;	 // structural nodal forces
		CMatrix<double> m_ELL;		//Element load local coordinates							
		CMatrix<double> m_ELG;		//Element load local coordinates
        std::ifstream m_FileInput;	// File Input
        std::ofstream m_FileOutput;	// File Output

        // error handlers
        void ErrorHandler (ERRORCODE nCode) const; 
        void IOErrorHandler (ERRORCODE nCode) const;
};

#endif