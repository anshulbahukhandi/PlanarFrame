/*********************************************
Planar Frame Analysis Program
Copyright(c) 2000-13, S. D. Rajan
All rights reserved

Introduction to Structural Analysis and Design
Object-Oriented Numerical Analysis
*********************************************/
#include "frame.h"
#include "F:\ASU Courses\cee532\ANSHUL_BAHUKHANDI_PLANAR FRAME\PlanarFrame\MatToolBox.h"

/* ==================================================================
   ======================= CFrame class =============================
   ================================================================== */

CFrame::CFrame ()
// ---------------------------------------------------------------------------
// Function: default constructor
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    m_nNodes = m_nElements = m_nEPGroups = m_nMatGroups
             = m_nDOF = m_nDebugLevel = m_nLineNumber
             = m_nElementLoads = 0;
	supportreactions.Set(0);
}

CFrame::~CFrame ()
// ---------------------------------------------------------------------------
// Function: destructor
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    // deallocate ... deallocate ... deallocate
    for (int i=1; i <= m_nEPGroups; i++)
    {
        delete m_EPData(i); 
    }
}

void CFrame::SetSize ()									//***STEP 3*****
// ---------------------------------------------------------------------------
// Function: memory allocation for all major arrays in the program
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
	//allocate soace for element load data
	m_ElementLoadData.SetSize(m_nElements);					//**CHANGED THIS FOR REGRADING******
    // allocate space for nodal data
    m_NodalData.SetSize (m_nNodes);
    // allocate space for nodal loads data
    m_NodalLoadData.SetSize (m_nNodes);
    // allocate space for material data
    m_MaterialData.SetSize (m_nMatGroups);
    // allocate space for x/s data
    m_EPData.SetSize (m_nEPGroups);
    // allocate space for nodal response data
    m_NodalResponseData.SetSize (m_nNodes);
    // allocate space for element data
    m_ElementData.SetSize (m_nElements);
    // allocate space for element response data
    m_ElementResponseData.SetSize (m_nElements);
    // allocate space for element loads, if required
    if (m_nElementLoads >= 0)					//*****CHANGED THIS FOR REGRADING******				
    {
        m_ELL.SetSize (DOFPE, m_nElements); m_ELL.Set(0.0f);
        m_ELG.SetSize (DOFPE, m_nElements); m_ELG.Set(0.0f);
	
	}

    // allocate and initialize major matrices
    m_nDOF = DOFPN*m_nNodes;
    m_SSM.SetSize (m_nDOF, m_nDOF);  m_SSM.Set (0.0);
    m_SND.SetSize (m_nDOF, 1);       m_SND.Set (0.0);
    m_SNF.SetSize (m_nDOF, 1);       m_SNF.Set (0.0);
//allocate size for support reactions
	supportreactions.SetSize(m_nDOF,1);												//changed this
}

void CFrame::ConstructF ()
// ---------------------------------------------------------------------------
// Function: constructs the system load vector
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{	CVector<int>nE(DOFPN);					//adding nodal loads into structural load vector
	CVector<float>fNLD(DOFPN);
	for(int i=1;i<=m_nNodes;i++)
	{																				//***STEP 8****

		m_NodalLoadData(i).GetValues(fNLD);
		nE(1)=3*i-2;nE(2)=3*i-1;nE(3)=3*i;
		m_SNF(nE(1),1)+=fNLD(1);
		m_SNF(nE(2),1)+=fNLD(2);
		m_SNF(nE(3),1)+=fNLD(3);
	
	}
	// converting ELL into ELG  AND storing  ELG into SNF					*****STEP 9*******
	CVector<float> fsn(NDIM);												//changed this
	CVector<float>fen(NDIM);
	
	for(int i=1;i<=m_nElements;i++)
	{   CMatToolBox<double> MTB;
		int sn,en;
		m_ElementData(i).GetENodes(sn,en);
		m_NodalData(sn).GetCoords(fsn);
		m_NodalData(en).GetCoords(fen);
    double length=sqrt((fen(1)-fsn(1))*(fen(1)-fsn(1)) + (fen(2)-fsn(2))*(fen(2)-fsn(2)));
	CMatrix<double>dT(6,6);		
	CMatrix<double>dTT(6,6);									//changed this
	dT.Set(0.0);
	double dl= dT(1,1)=(fen(1)-fsn(1))/length;
	double dm=dT(1,2)=(fen(2)-fsn(2))/length;
	dT(2,1)=-(fen(2)-fsn(2))/length;
	dT(2,2)=(fen(1)-fsn(1))/length;
	dT(3,3)=1;
	dT(6,6)=1;
	dT(4,4)=dl;
	dT(4,5)=dm;
	dT(5,4)=-dm;
	dT(5,5)=dl;	
	
	CVector<int>nE(6);
		m_ElementData(i).GetENodes(sn,en);
		nE(1)=3*sn-2;nE(2)=3*sn-1;nE(3)=3*sn;
		nE(4)=3*en-2;nE(5)=3*en-1;nE(6)=3*en;
	
	if(MTB.Transpose(dT,dTT))
	{
		if(MTB.Multiply(dTT,m_ELL,m_ELG))
		{
			for(int j=1;j<=6;j++)
			{
			m_SNF(nE(j),1)+=m_ELG(j,i);						//Structural force matrix updated with element loads
			}
			}
		}
	}
}

void CFrame::ConstructK ()
// ---------------------------------------------------------------------------
// Function: constructs the system stiffness matrix
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
	CMatToolBox<double> MTB;
    int nSN, nEN, nMatGrp;
    CXSType* pXSGrp;
    float fArea, fIyy, fIzz;
    CVector<float> fVC1(NDIM), fVC2(NDIM);
    CVector<CNode::Fixity> VFC1(DOFPN), VFC2(DOFPN);
	CMatrix<double>dKEL(6,6);
	dKEL.Set(0);
		CMatrix<double>dT(6,6);		
	CMatrix<double>dTT(6,6);									//changed this
	dT.Set(0.0);
		CMatrix<double>dTemp(6,6);
		CMatrix<double>dKEG(6,6);
		dTemp.Set(0.0);
		dKEG.Set(0.0);

	for (int i=1; i <= m_nElements; i++)
    {
        m_ElementData(i).GetENodes (nSN, nEN);
        m_NodalData(nSN).GetCoords (fVC1);
        m_NodalData(nEN).GetCoords (fVC2);
        m_NodalData(nSN).GetFixity (VFC1);
        m_NodalData(nEN).GetFixity (VFC2);
        pXSGrp = m_ElementData(i).GetEPropertyGroup ();
        nMatGrp = m_ElementData(i).GetMatPropertyGroup ();
        pXSGrp->GetProperties(fArea, fIyy, fIzz);
		
        float fYM = m_MaterialData(nMatGrp).GetYM ();
        float fCTE = m_MaterialData(nMatGrp).GetCTE ();
	double length=sqrt((fVC2(1)-fVC1(1))*(fVC2(1)-fVC1(1)) + (fVC2(2)-fVC1(2))*(fVC2(2)-fVC1(2)));
	
	double dl= dT(1,1)=(fVC2(1)-fVC1(1))/length;			//***STEP 10***
	double dm=dT(1,2) =(fVC2(2)-fVC1(2))/length;
	dT(2,1)=-dm;
	dT(2,2)=dl;
	dT(3,3)=1;
	dT(6,6)=1;
	dT(4,4)=dl;
	dT(4,5)=dm;
	dT(5,4)=-dm;
	dT(5,5)=dl;	
	
	 double dAEOL = double((fArea*fYM)/length);
	 double d1= (12*fYM*fIzz)/(length*length*length);
	 double d2= (6*fYM*fIzz)/(length*length);
	 double d3=(4*fYM*fIzz)/(length);	
	 double d4 = (2*fYM*fIzz)/length;

	    dKEL(1,1) = dKEL(4,4) = dAEOL;							//***STEP 11***
		dKEL(1,4) = dKEL(4,1) = -dAEOL;
		dKEL(2,2) = dKEL(5,5) = d1;
		dKEL(2,5) = dKEL(5,2) = -d1;
		dKEL(2,3) = dKEL(2,6) =  dKEL(3,2) = dKEL(6,2) = d2;
		dKEL(3,5) = dKEL(5,3) = dKEL(5,6) = dKEL(6,5) = -d2;
		dKEL(3,6) = dKEL(6,3) = d4;
		dKEL(3,3) = dKEL(6,6) = d3;
		
		

		CVector<int>nE(6);
		nE(1) = 3*nSN-2; nE(2) = 3*nSN-1; nE(3) = 3*nSN;
		nE(4) = 3*nEN-2; nE(5) = 3*nEN-1; nE(6) = 3*nEN;

		if(MTB.Transpose(dT,dTT))
		{
			if(MTB.Multiply(dKEL,dT,dTemp))					//***STEP 12***
			{
			if(MTB.Multiply(dTT,dTemp,dKEG))
			{		
				for(int o=1 ;o<=6;o++)
				{
				
				int row=nE(o);								//***STEP 13***
				for(int j=1;j<=6;j++)
			{
				int column=nE(j);
				m_SSM(row,column)+=dKEG(o,j);
			}
				}
			}
		}
	}
	}
	
			
			
}

void CFrame::ImposeBC ()
// ---------------------------------------------------------------------------
// Function: imposes the homogenous and non-homogenous EBC
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
	
	CVector<CNode::Fixity> VFC(DOFPN);
	CVector<float>VD(DOFPN);
	// loop thro' all nodes
	for (int i=1; i <= m_nNodes; i++)
	{

		//Imposing boundary conditions on nodal displacements

		m_NodalData(i).GetFixity (VFC);
		m_NodalData(i).GetDisp(VD);
		if (VFC(1) == 1)
		{
			int nGDOF = 3*i-2;
			SuppressDOF (nGDOF);

			for (int j=nGDOF ;j<=m_nDOF;j++)				//Changing the structural load matrix according to B.C
			{ 
				if (j==nGDOF)
				{
					m_SNF(j,1)=VD(1);
				}

				else
				{
					m_SNF(j,1)=m_SNF(j,1)-m_SSM(j,nGDOF)*VD(1);					//*****STEP 15****
				}
			}
		}
		if (VFC(2) == 1)
		{
			int nGDOF = 3*i-1;
			SuppressDOF (nGDOF);
			for (int k=nGDOF ;k<=m_nDOF;k++)			//changing the structural load matrix according to B.C
			{ if (k==nGDOF)
			{
				m_SNF(k,1)=VD(2);
			}
			else
			{
				m_SNF(k,1)=m_SNF(k,1)-m_SSM(k,nGDOF)*VD(2);
			}
			}

		}
		if (VFC(3) == 1)
		{
			int nGDOF = 3*i;
			SuppressDOF (nGDOF);
			for (int n=nGDOF ;n<=m_nDOF;n++)            //changing the structural load matrix according to B.C
			{ if (n==nGDOF)
			m_SNF(n,1)=VD(3);
			else
				m_SNF(n,1)=m_SNF(n,1)-m_SSM(n,nGDOF)*VD(3);
			}
		}
	}

}

void CFrame::SuppressDOF (const int nEqn)
// ---------------------------------------------------------------------------
// Function: works in conjunction with ImposeBC
// Input:    global equation # to impose EBC
// Output:   none
// ---------------------------------------------------------------------------
{
		for (int j=1; j <= m_nDOF; j++)
	{
		// zero out the row
		m_SSM(nEqn, j) = 0.0;								//***STEP 14*****
		// zero out the column
		m_SSM(j, nEqn) = 0.0;
	}
	// set diagonal to 1  
	m_SSM(nEqn, nEqn) = 1.0;

		
}

void CFrame::Solve ()
// ---------------------------------------------------------------------------
// Function: solves for the nodal displacements
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
	CMatToolBox<double> MTB;
	double TOL = 1.0e-6;
	CVector<double>fVD(3);
	if(MTB.LDLTFactorization (m_SSM, TOL) )       //Factorizing the stiffness matrix      //***STEP 16***
	{
		
		if(MTB.LDLTSolve(m_SSM,m_SND,m_SNF))     //solving for the structural nodal displacements  //***STEP 17***
		{	
						for ( int i=1; i <= m_nNodes; i++)
			{
				fVD(1) = static_cast<double>(m_SND(3*i-2,1));
				fVD(2) = static_cast<double>(m_SND(3*i-1,1));
				fVD(3) = static_cast<double>(m_SND(3*i,1));
				m_NodalResponseData(i).SetValues(fVD);          //Storing values from SND to the respective nodes
			}			
		}
		else ErrorHandler(UNSTABLE);
	}
	else ErrorHandler(UNSTABLE);					//changed this
}


void CFrame::Response ()
// ---------------------------------------------------------------------------
// Function: computes the element nodal forces and support reactions
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
	//calculating support reactions
	CMatToolBox<double> MTB;
	if(MTB.Multiply(m_SSM,m_SND,m_SNF))
	{

	}																//PROBLEM PROBLEM PROBLEM

	int nSN,nEN;	
	CVector<float>fCS(NDIM);
	CVector<float>fCE(NDIM);
	double length;
	CMatrix<double>dT(6,6);
	CMatrix<double>dTT(6,6);
	CMatrix<double>dEDG(6,1);
	CMatrix<double>dEDL(6,1);
	CMatrix<double>dELsingle(6,1);
	dELsingle.Set(0.0);
	double maxT=0;
	double maxC=0;
	double maxS=0;
	CMatrix<double>dTempa(6,1);
	CMatrix<double>dTempb(6,1);
	CMatrix<double>dTempc(6,1);
	CMatrix<double>dKEL(6,6);
	dTempb.Set(0);
	dT.Set(0);
	dTT.Set(0);
	dKEL.Set(0);
	dTempa.Set(0);
	dTempc.Set(0);
	for (int i=1; i <= m_nElements; i++)
	{

		m_ElementData(i).GetENodes(nSN,nEN);
		m_NodalData(nSN).GetCoords (fCS);
		m_NodalData(nEN).GetCoords (fCE);
		length = sqrt((fCE(1)-fCS(1))*(fCE(1)-fCS(1)) +
			(fCE(2)-fCS(2))*(fCE(2)-fCS(2)));

		// local-to-global transformation matrix

		double dl=dT(1,1) = double((fCE(1)-fCS(1))/length);
		double dm=dT(1,2) = double((fCE(2)-fCS(2))/length);
		dT(2,1)=-dm;
		dT(2,2)=dl;
		dT(3,3)=1;
		dT(6,6)=1;
		dT(4,4)=dl;
		dT(4,5)=dm;
		dT(5,4)=-dm;
		dT(5,5)=dl;	
		CVector<int> nE(6);

		// get element nodal displacements in global coordinates                      
		nE(1) = 3*nSN-2; nE(2) = 3*nSN-1;nE(3) = 3*nSN;
		nE(4) = 3*nEN-2; nE(5) = 3*nEN-1;nE(6) = 3*nEN;
		for (int j=1; j <= 6; j++)
		{
			dEDG(j,1) = m_SND(nE(j),1);				//***STEP 18***
			dELsingle(j,1)=m_ELL(j,i);				//***STEP 19***
		}
		//constructing K local for each element
		std::string shape;							//changed this
		CVector<double> Dims(MAXEPDIM);					//changed this

		CXSType* pXSGrp;
		float fArea, fIyy, fIzz;
		pXSGrp = m_ElementData(i).GetEPropertyGroup ();
		int nMatGrp = m_ElementData(i).GetMatPropertyGroup ();
		pXSGrp->GetProperties(fArea, fIyy, fIzz);

		pXSGrp->GetShape(shape);							//changed this
		pXSGrp->GetDimensions(Dims);						//changed this
		float fYM = m_MaterialData(nMatGrp).GetYM ();
		float fCTE = m_MaterialData(nMatGrp).GetCTE ();

		double dAEOL = double((fArea*fYM)/length);
		double d1= (12*fYM*fIzz)/(length*length*length);
		double d2= (6*fYM*fIzz)/(length*length);
		double d3=(4*fYM*fIzz)/(length);	
		double d4 = (2*fYM*fIzz)/length;
		dKEL(1,1) = dKEL(4,4) = dAEOL;
		dKEL(1,4) = dKEL(4,1) = -dAEOL;
		dKEL(2,2) = dKEL(5,5) = d1;
		dKEL(2,5) = dKEL(5,2) = -d1;
		dKEL(2,3) = dKEL(2,6) =  dKEL(3,2) = dKEL(6,2) = d2;
		dKEL(3,5) = dKEL(5,3) = dKEL(5,6) = dKEL(6,5) = -d2;
		dKEL(3,6) = dKEL(6,3) = d4;
		dKEL(3,3) = dKEL(6,6) = d3;

		//Forces

		if(MTB. Multiply ( dT, dEDG,dEDL))	                                 	//***STEP 20***					
		{
			if(MTB.Multiply(dKEL,dEDL, dTempa))						//***STEP 21***
			{

				//subtracting Element loads to get initial loading conditions
				if(MTB.Subtract(dTempa,dELsingle,dTempb))			//***STEP 22***
				{

					m_ElementResponseData(i).SetValues(dTempb);			//Storing Element Nodal Forces in Local coordinates		
				}
				else ErrorHandler(CANNOTOPENOFILE);					//changed this
				if(MTB.Transpose(dT,dTT))						//***STEP 23***
				{
					if(MTB.Multiply(dTT,dTempb,dTempc))		//storing element  Nodal Forces in global cordinates in dTempc
					{
					
						m_ElementResponseData(i).SetESupportReaction(dTempc);
					}
					else ErrorHandler(CANNOTOPENOFILE);  //changed this
				}

			}
		}
		//Calculating Maximum Tensile Compressive and Shear stress

		CMatrix<double> Forces(DOFPE,1);
		m_ElementResponseData(i).GetValues(Forces);
		
		double NMax;
		double MMax;
		double VMax;
		if (abs(Forces(1,1))>abs(Forces(4,1)))  NMax=Forces(1,1);  else NMax=Forces(4,1);		//***STEP 24***
		if(abs(Forces(3,1))>abs(Forces(6,1)))   MMax=Forces(3,1); else MMax=Forces(6,1);
		if(abs(Forces(2,1))>abs(Forces(5,1)))   VMax=Forces(2,1); else VMax=Forces(5,1);

		if(shape=="rects")					//***STEP 25***
		{
			maxT=(NMax/Dims(1)*Dims(2))+(6*MMax/(Dims(1)*Dims(1)*Dims(2)));				//***STEP 26***
			maxC= (NMax/Dims(1)*Dims(2))-(6*MMax/(Dims(1)*Dims(1)*Dims(2)));
			maxS=3*VMax/(2*Dims(1)*Dims(2));

		}
		else if(shape=="isection")
		{ 
			float wH=Dims(1);
			float wT=Dims(2);
			float wF=Dims(3);
			float tF=Dims(4);
			double A=2*wF*tF +wT*wH;
			double I=(wF*tF*tF*tF/6)+(wF*tF*wH*wH/4)+(wT*wH*wH*wH/12);
			double S=(2*I)/(wH+2*tF);
			double SF=8*I*wT/((4*wF*tF*(wH+tF))+(wT*wH*wH));
			maxT=abs(NMax/A+ MMax/S);
			maxC= NMax/A- MMax/S;
			maxS=abs(VMax/SF);

		}
		else if (shape=="circs")
		{
			maxT=abs((NMax/PIE*Dims(1)*Dims(1))+(4*MMax/(PIE*Dims(1)*Dims(1)*Dims(1))));
			maxC= (NMax/PIE*Dims(1)*Dims(1))-(4*MMax/(PIE*Dims(1)*Dims(1)*Dims(1)));
			maxS=abs(4*VMax/3*PIE*Dims(1)*Dims(1));
		}
		else ErrorHandler(XSTYPE);
		m_ElementResponseData(i).SetStress(maxT,maxC,maxS);
	}
	

	//calculating support reactions at each node.
	supportreactions.Set(0);
	for(int i=1;i<=m_nElements;i++)
		
	{
		int SN,EN;
		CVector<int> nE(6);
		m_ElementData(i).GetENodes(SN,EN);
		nE(1) = 3*SN-2; nE(2) = 3*SN-1;nE(3) = 3*SN;				//Respective positions 
		nE(4) = 3*EN-2; nE(5) = 3*EN-1;nE(6) = 3*EN;
		CMatrix<double> ESR(DOFPE,1);

		m_ElementResponseData(i).GetESupportReaction(ESR);
				for(int j=1;j<=DOFPE;j++)
		{
		
			int nrow=nE(j);
			supportreactions(nrow,1)+=ESR(j,1);							//***STEP 27***
			
		}
	}
	
}
//Computes residual and absolute error
void CFrame::ResidualError()
{
	CMatToolBox<double> MTB;
	CMatrix<double> dTempd;
	CMatrix<double> dTempe;
	dTempd.Set(0);
	dTempe.Set(0);
	dTempd.SetSize(m_nDOF,1);
	dTempe.SetSize(m_nDOF,1);
	double sum=0;
	if(MTB.Multiply(m_SSM,m_SND,dTempd))
	{
		
	if(	MTB.Subtract(dTempd,m_SNF,dTempe))
	{
		for(int i=1;i<=m_nDOF;i++)
		{
			sum+=dTempe(i,1);
		}
		RError=sum;
		AError=sum/m_nDOF;
	}

	}

}

void CFrame::ErrorHandler (ERRORCODE nCode) const
// ---------------------------------------------------------------------------
// Function: displays error messages related to frame analysis
// Input:    error code
// Output:   none
// ---------------------------------------------------------------------------
{
    std::cerr << '\n';
	if(nCode==UNSTABLE)
		std::cerr<<"Unstable Trusss";
	else if(nCode==INVALIDCOMMANDLINE)
		std::cerr<<"Invalid command Line!!!!";
	else if(nCode==CANNOTOPENIFILE)
		std::cerr<<"Cannot open Input File.File Being used or does not exists!!";
	else if(nCode== CANNOTOPENOFILE)
		std::cerr<<"Cannot open Input File.File Being used or does not exists!!";
	else std::cerr<<"UNKNOWN ERROR......!!!!!";
	
	

    std::cerr << std::endl;
    exit (1);
}