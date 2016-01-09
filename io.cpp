/*********************************************
Planar Frame Analysis Program
Copyright(c) 2000-13, S. D. Rajan
All rights reserved

Introduction to Structural Analysis and Design
Object-Oriented Numerical Analysis
*********************************************/
#include <iomanip>
#include <sstream>
#include"elementloads.h"
#include "frame.h"
#include "F:\ASU Courses\cee532\Library\Library\parser.h"
#include "F:\ASU Courses\cee532\Library\Library\fileio.h"
#include "F:\ASU Courses\cee532\Programs Chapter 2 through 18\Example10_3_1\MatToolBox.h"

const int MAXCHARS = 80;
std::string szInputString;
std::string szComment = "**";

void CFrame::Banner (std::ostream& OF)
// ---------------------------------------------------------------------------
// Function: Prints program banner
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    OF << '\n';
    OF << "\t\t--------------------------------------------" << '\n';
    OF << "\t\t         Planar Frame Analysis Program      " << '\n';
    OF << "\t\tIntroduction to Structural Analysis & Design" << '\n';
    OF << "\t\t           (c) 2000-14, S. D. Rajan         " << '\n';
    OF << "\t\t         Enhanced By: Anshul Bahukhandi        " << '\n';
    OF << "\t\t--------------------------------------------" << '\n';
}

void CFrame::PrepareIO (int argc, char *argv[])
// ---------------------------------------------------------------------------
// Function: Obtains file names and opens input/output files
// Input:    command line arguments
// Output:   none
// ---------------------------------------------------------------------------
{
    if (argc == 1)
    {
        // open the input file
        OpenInputFileByName ("Complete input file name: ", m_FileInput,						//*******STEP 1*********
                             std::ios::in);

        // open the output file
        OpenOutputFileByName ("Complete output file name: ", m_FileOutput,
                              std::ios::out);
    }
    else if (argc == 3) // spacetruss input_file output_file
    {
        m_FileInput.open (argv[1], std::ios::in);
        if (!m_FileInput)
            ErrorHandler (CANNOTOPENIFILE);
        m_FileOutput.open (argv[2], std::ios::out);
        if (!m_FileOutput)
            ErrorHandler (CANNOTOPENOFILE);
        std::cout << "\n";
        std::cout << argv[1] << " opened as input file.\n";
        std::cout << argv[2] << " opened as output file.\n";
    }
	else
    {
        ErrorHandler (INVALIDCOMMANDLINE);
    }

    // print banner
    Banner (m_FileOutput);
}

void CFrame::ReadProblemSize ()
// ---------------------------------------------------------------------------
// Function: reads the size of the problem from input file but does not
//           store data
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    CParser Parse;																	//*****STEP 2******

    // read the problem description
    if (!Parse.ReadNextLine (m_FileInput, m_nLineNumber, szInputString, 
                             MAXCHARS, szComment))
        IOErrorHandler (INVALIDINPUT);
    if (!Parse.ReadNextLine (m_FileInput, m_nLineNumber, szInputString, 
                             MAXCHARS, szComment))
        IOErrorHandler (INVALIDINPUT);

    // nodal coordinates
    if (!Parse.ReadNextLine (m_FileInput, m_nLineNumber, szInputString, 
                             MAXCHARS, szComment))
        IOErrorHandler (INVALIDINPUT);
    for (;;)
    {
        if (!Parse.ReadNextLine (m_FileInput, m_nLineNumber, szInputString, 
                                 MAXCHARS, szComment))
            IOErrorHandler (INVALIDINPUT);
        if (szInputString.substr(0,13) == "*nodal fixity")
            break;
        ++m_nNodes;
    }

    // nodal fixity
    for (;;)
    {
        if (!Parse.ReadNextLine (m_FileInput, m_nLineNumber, szInputString, 
                                 MAXCHARS, szComment))
            IOErrorHandler (INVALIDINPUT);
        if (szInputString.substr(0,12) == "*nodal loads")
            break;
    }

    // nodal loads
    for (;;)
    {
        if (!Parse.ReadNextLine (m_FileInput, m_nLineNumber, szInputString, 
                                 MAXCHARS, szComment))
            IOErrorHandler (INVALIDINPUT);
        if (szInputString.substr(0,14) == "*material data")
            break;
    }

    // material data
    for (;;)
    {
        if (!Parse.ReadNextLine (m_FileInput, m_nLineNumber, szInputString, 
                                 MAXCHARS, szComment))
            IOErrorHandler (INVALIDINPUT);
        if (szInputString.substr(0,21) == "*cross-sectional data")
            break;
        ++m_nMatGroups;
    }

    // cross-sectional data
    for (;;)
    {
        if (!Parse.ReadNextLine (m_FileInput, m_nLineNumber, szInputString, 
                                 MAXCHARS, szComment))
            IOErrorHandler (INVALIDINPUT);
        if (szInputString.substr(0,13) == "*element data")
            break;
        ++m_nEPGroups;
    }

    // element data
    for (;;)
    {
        if (!Parse.ReadNextLine (m_FileInput, m_nLineNumber, szInputString, 
                                 MAXCHARS, szComment))
            IOErrorHandler (INVALIDINPUT);
        if (szInputString.substr(0,14) == "*element loads")
            break;
        ++m_nElements;
    }

    // element loads
		
    for (;;)
    {
        if (!Parse.ReadNextLine (m_FileInput, m_nLineNumber, szInputString, 
                                 MAXCHARS, szComment))
            IOErrorHandler (INVALIDINPUT);
            if (szInputString.substr(0,4) == "*end")
                break;
        m_nElementLoads++;
    }
	
	m_ElementLoadData.SetSize(m_nElements);				//set the vector size..used to store the number of elements
	
    // check data for validity
    if (m_nNodes <= 1) 
        IOErrorHandler (NUMNODES);
    if (m_nElements <= 0) 
        IOErrorHandler (NUMELEMENTS);
    if (m_nMatGroups <= 0) 
        IOErrorHandler (NUMMATGROUPS); 
    if (m_nEPGroups <= 0) 
        IOErrorHandler (NUMEPROPGROUPS); 
    if (m_nDebugLevel < 0 || m_nDebugLevel > 1) 
        IOErrorHandler (DEBUGCODE);
}

void CFrame::ReadFrameModel ()
// ---------------------------------------------------------------------------
// Function: reads the rest of the frame data from input file
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    int i, nTag;
    CParser Parse;

    // rewind the file to read the input file again
    Rewind (m_FileInput);														//*****STEP 4*******
    
    // header line
    if (!Parse.ReadNextLine (m_FileInput, m_nLineNumber, szInputString, 
                             MAXCHARS, szComment))
        IOErrorHandler (INVALIDINPUT);
    if (!Parse.ReadNextLine (m_FileInput, m_nLineNumber, szInputString, 
                             MAXCHARS, szComment))
        IOErrorHandler (INVALIDINPUT);

    // read nodal coordinates														//*****STEP 5********
    CVector<float> fVC(NDIM);
    if (!Parse.ReadNextLine (m_FileInput, m_nLineNumber, szInputString, 
                             MAXCHARS, szComment))
        IOErrorHandler (INVALIDINPUT);
    for (i=1; i <= m_nNodes; i++)
    {
        if (!Parse.ReadNextLine (m_FileInput, m_nLineNumber, szInputString, 
                                 MAXCHARS, szComment))
            IOErrorHandler (INVALIDINPUT);
        else
        {
            std::istringstream szFormatString (szInputString);
            szFormatString >> nTag >> fVC(1) >> fVC(2);
            if (szFormatString.fail() || szFormatString.bad())
                IOErrorHandler (INVALIDINPUT);
        }
        if (nTag <= 0 || nTag > m_nNodes) 
            IOErrorHandler (NODENUMBER);
        m_NodalData(nTag).SetCoords (fVC);
    }

    // read nodal fixity conditions
    if (!Parse.ReadNextLine (m_FileInput, m_nLineNumber, szInputString, 
                             MAXCHARS, szComment))
        IOErrorHandler (INVALIDINPUT);
    std::string szXFC, szYFC, szZFC;
	CVector<float>VD(DOFPN);
    CVector<CNode::Fixity> VFC(DOFPN);
    for (;;)
    {
        if (!Parse.ReadNextLine (m_FileInput, m_nLineNumber, szInputString, 
                                 MAXCHARS, szComment))
            IOErrorHandler (INVALIDINPUT);
        if (szInputString.substr(0,12) == "*nodal loads")
            break;

        std::istringstream szFormatString (szInputString);
        szFormatString >> nTag;
        if (szFormatString.fail() || szFormatString.bad())
            IOErrorHandler (INVALIDINPUT);
        if (nTag <= 0)
            break;
        szFormatString >> szXFC  >> szYFC  >> szZFC
                       >> VD(1) >> VD(2) >> VD(3);
        if (szFormatString.fail() || szFormatString.bad())
            IOErrorHandler (INVALIDINPUT);

        if (nTag > m_nNodes) 
            IOErrorHandler (NODENUMBER);
        if (szXFC != "free" && szXFC != "specified")
            IOErrorHandler (NODALFIXITY);
        if (szYFC != "free" && szYFC != "specified")
            IOErrorHandler (NODALFIXITY);
        if (szZFC != "free" && szZFC != "specified")
            IOErrorHandler (NODALFIXITY);

        VFC(1) = (szXFC == "specified" ? CNode::SPECIFIED : CNode::FREE);
        VFC(2) = (szYFC == "specified" ? CNode::SPECIFIED : CNode::FREE);
        VFC(3) = (szZFC == "specified" ? CNode::SPECIFIED : CNode::FREE);
        m_NodalData(nTag).SetFixity (VFC);
		m_NodalData(nTag).SetDisp(VD);
}
	
    // read nodal loads
	CVector<float> fVForce(DOFPN);
	float temperature;
	for (;;)
    {
        if (!Parse.ReadNextLine (m_FileInput, m_nLineNumber, szInputString, 
                                 MAXCHARS, szComment))
            IOErrorHandler (INVALIDINPUT);
        if (szInputString.substr(0,14) == "*material data")
            break;

        std::istringstream szFormatString (szInputString);
        szFormatString >> nTag;
        if (szFormatString.fail() || szFormatString.bad())
            IOErrorHandler (INVALIDINPUT);
        if (nTag <= 0)
            break;
		szFormatString >> fVForce(1) >> fVForce(2) >> fVForce(3)>>temperature;
        if (szFormatString.fail() || szFormatString.bad())
            IOErrorHandler (INVALIDINPUT);
        if (nTag > m_nNodes) 
            IOErrorHandler (NODENUMBER);																//Removed the THE WAY  OF INDEXING
	m_NodalLoadData(nTag).SetValues(fVForce);									
	m_NodalLoadData(nTag).SetTemp(temperature);				//*****CHANGED THIS FOR REGRADING********
	}											

    // read material data
    for (i=1; i <= m_nMatGroups; i++)
    {
        if (!Parse.ReadNextLine (m_FileInput, m_nLineNumber, szInputString, 
                                 MAXCHARS, szComment))
         
								 IOErrorHandler (INVALIDINPUT);
        else
        {
            std::istringstream szFormatString (szInputString);
            szFormatString >> nTag >> fVC(1);
            if (szFormatString.fail() || szFormatString.bad())
                IOErrorHandler (INVALIDINPUT);
            if (fVC(1) <= 0.0f)
                IOErrorHandler (MATPROPERTY);
        }
        if (nTag <= 0 || nTag > m_nMatGroups)
            IOErrorHandler (MATGROUPNUMBER);
        m_MaterialData(nTag).SetYM (fVC(1));
    }

    // read cross-sectional data
    std::string szTag;
    CVector<float> fVXSDims(MAXEPDIM);
    if (!Parse.ReadNextLine (m_FileInput, m_nLineNumber, szInputString, 
                             MAXCHARS, szComment))
        IOErrorHandler (INVALIDINPUT);
    for (i=1; i <= m_nEPGroups; i++)
    {
        if (!Parse.ReadNextLine (m_FileInput, m_nLineNumber, szInputString, 
                                 MAXCHARS, szComment))
            IOErrorHandler (INVALIDINPUT);
        else
        {
            std::istringstream szFormatString (szInputString);
            szFormatString >> nTag >> szTag;
            if (nTag > m_nNodes) 
                IOErrorHandler (NODENUMBER);
            if (szTag == "rects")
            {
                szFormatString >> fVXSDims(1) >> fVXSDims(2);
                if (szFormatString.fail() || szFormatString.bad())
                    IOErrorHandler (INVALIDINPUT);
                if (fVXSDims(1) <= 0.0f || fVXSDims(2) <= 0.0f)
                    IOErrorHandler (XSDIMENSION);
                m_EPData(nTag) = new CRectSolid (fVXSDims);						//Not getting
				
			}
            else if (szTag == "isection")
            {
                szFormatString >> fVXSDims(1) >> fVXSDims(2) 
                               >> fVXSDims(3) >> fVXSDims(4);
                if (szFormatString.fail() || szFormatString.bad())
                    IOErrorHandler (INVALIDINPUT);
                if (fVXSDims(1) <= 0.0f || fVXSDims(2) <= 0.0f ||
                    fVXSDims(3) <= 0.0f || fVXSDims(4) <= 0.0f)
                    IOErrorHandler (XSDIMENSION);
                m_EPData(nTag) = new CISection (fVXSDims);
            }

			else if (szTag == "circs")								//CHANGED THIS
            {
                szFormatString >> fVXSDims(1); 
                if (szFormatString.fail() || szFormatString.bad())
                    IOErrorHandler (INVALIDINPUT);
                if (fVXSDims(1) <= 0.0f)
                    IOErrorHandler (XSDIMENSION);
                m_EPData(nTag) = new CCircSolid (fVXSDims);
            }
            else
                IOErrorHandler (XSTYPE);
        }
    }

    // read element data
    if (!Parse.ReadNextLine (m_FileInput, m_nLineNumber, szInputString, 
                             MAXCHARS, szComment))
        IOErrorHandler (INVALIDINPUT);
    int nSN, nEN, nMatGrp, nEPGrp;
    for (i=1; i <= m_nElements; i++)
    {
        if (!Parse.ReadNextLine (m_FileInput, m_nLineNumber, szInputString, 
                                 MAXCHARS, szComment))
            IOErrorHandler (INVALIDINPUT);
        else
        {
            std::istringstream szFormatString (szInputString);
            szFormatString >> nTag >> nSN  >> nEN
                           >> nMatGrp >> nEPGrp;
            if (szFormatString.fail() || szFormatString.bad())
                IOErrorHandler (INVALIDINPUT);
        }
        if (nTag <= 0 || nTag > m_nElements)
            IOErrorHandler (ELEMENTNUMBER);
        if (nMatGrp <= 0 || nMatGrp > m_nMatGroups)
            IOErrorHandler (MATGROUPNUMBER);
        if (nEPGrp <= 0 || nEPGrp > m_nEPGroups)
            IOErrorHandler (EPROPNUMBER);
        if (nSN <= 0 || nSN > m_nNodes)
            IOErrorHandler (NODENUMBER);
        if (nEN <= 0 || nEN > m_nNodes)
            IOErrorHandler (NODENUMBER);
        m_ElementData(nTag).SetMatPropertyGroup (nMatGrp);
        m_ElementData(nTag).SetEPropertyGroup (m_EPData(nEPGrp));
        m_ElementData(nTag).SetENodes (nSN, nEN);
}
    // read element loads
    if (!Parse.ReadNextLine (m_FileInput, m_nLineNumber, szInputString, 
                             MAXCHARS, szComment))
        IOErrorHandler (INVALIDINPUT);
	std::string szEType;
    CVector<float> fVELoads(2);
	CElementLoads::ELType ELT;					// made this!!
	

    for (;;)
    {
		
        if (!Parse.ReadNextLine (m_FileInput, m_nLineNumber, szInputString, 
                                 MAXCHARS, szComment))
            IOErrorHandler (INVALIDINPUT);
        if (szInputString.substr(0,4) == "*end")
            break;
        std::istringstream szFormatString (szInputString);
        szFormatString >> nTag;
        if (szFormatString.fail() || szFormatString.bad())
            IOErrorHandler (INVALIDINPUT);
        if (nTag <= 0)
            break;
        if (nTag > m_nElements) 
            IOErrorHandler (ELEMENTNUMBER);
        szFormatString >> szEType >> fVELoads(1) >> fVELoads(2);
        if (szFormatString.fail() || szFormatString.bad())
            IOErrorHandler (INVALIDINPUT);
        if (szEType != "dly'" && szEType != "ploadx'" &&
			szEType != "ploady'"&& szEType !="cmoment")				//changed this
            IOErrorHandler (ELOADTYPE);
      
		
		// TO DO: construct equivalent nodal loads and update the matrices
		
		
		
		
	CVector<float> fsn(NDIM);												
	CVector<float>fen(NDIM);
	int sn,en;
	float length;
	m_ElementData(nTag).GetENodes(sn , en);
	m_NodalData(sn).GetCoords(fsn);
	m_NodalData(en).GetCoords(fen);
	length=sqrt((fen(1)-fsn(1))*(fen(1)-fsn(1)) + (fen(2)-fsn(2))*(fen(2)-fsn(2)));
	if(szEType=="dly'")														//*******STEP 6*******
	{
			ELT=CElementLoads::LinearlyDistributed;	
		m_ELL(1,nTag)+=0;														//***STEP 7******
		m_ELL(2,nTag)+=length*(7*fVELoads(1)+3*fVELoads(2))/20;
		m_ELL(3,nTag)+=length*length*(3*fVELoads(1)+2*fVELoads(2))/60;
		m_ELL(4,nTag)+=0;
		m_ELL(5,nTag)+=length*(7*fVELoads(2)+3*fVELoads(1))/20;
		m_ELL(6,nTag)+=-length*length*(3*fVELoads(2)+2*fVELoads(1))/60;
	m_ElementLoadData(nTag).SetValues(nTag,ELT,fVELoads(1),fVELoads(2));			
	
	}
	
	if(szEType==" ploadx' ")
	{
		
		ELT=CElementLoads::ConcentratedX;
		m_ELL(1,nTag)+=fVELoads(2)*(length-fVELoads(1))/length;
		
		m_ELL(2,nTag)+=0;
		m_ELL(3,nTag)+=0;
		m_ELL(4,nTag)+=fVELoads(2)*fVELoads(1)/length;
		m_ELL(5,nTag)+=0;
		m_ELL(6,nTag)+=0;
	
	 	m_ElementLoadData(nTag).SetValues(nTag,ELT,fVELoads(1),fVELoads(2));
	}
	if(szEType=="cmoment")
	{
		
		ELT=CElementLoads::Moment;
			m_ELL(1,nTag)+=0;
			m_ELL(2,nTag)+=(-6*fVELoads(1)*(length-fVELoads(1)))/(fVELoads(2)*fVELoads(2)*fVELoads(2));
			m_ELL(3,nTag)+=(-fVELoads(2)*(length-fVELoads(1))*(3*fVELoads(1)-length))/(length*length);
			m_ELL(4,nTag)+=0;
			m_ELL(5,nTag)+=6*fVELoads(1)*fVELoads(2)*(length-fVELoads(1))/(length*length*length);
			m_ELL(6,nTag)+=fVELoads(1)*fVELoads(2)*(3*fVELoads(1)-2*length)/(length * length);
	m_ElementLoadData(nTag).SetValues(nTag,ELT,fVELoads(1),fVELoads(2));	
	}
	
	if(szEType=="ploady'")
	{
		ELT=CElementLoads::ConcentratedY;
		m_ELL(1,nTag)+=0;
		m_ELL(2,nTag)+=(fVELoads(2)*(length-fVELoads(1))*(length-fVELoads(1))*(length+2*fVELoads(1)));
		m_ELL(3,nTag)+=(fVELoads(2)*fVELoads(1)*(length-fVELoads(1))*(length-fVELoads(1))/(length*length));
		m_ELL(4,nTag)+=0;
		m_ELL(5,nTag)+=(fVELoads(2)*fVELoads(1)*fVELoads(1)*(3*length-2*fVELoads(1))/(length*length*length));
		m_ELL(6,nTag)+=(-fVELoads(2)*fVELoads(1)*fVELoads(1)*(length-fVELoads(1))/(length*length));
	m_ElementLoadData(nTag).SetValues(nTag,ELT,fVELoads(1),fVELoads(2));
	}
	
	}

		}
	



void CFrame::CreateOutput ()
// ---------------------------------------------------------------------------
// Function: creates the output file
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{	
	m_FileOutput<<"\n \n"<<"Number of Nodes :"<<m_nNodes<<std::endl;
	m_FileOutput<<"Number of Elements:"<<m_nElements<<"\n\n";
	m_FileOutput<<"____________________________________________________________________________________________\n";
	
	//Nodal cordinates

	m_FileOutput<<"\n\nNodal cordinates\n";
	m_FileOutput<<std::setw(10)<<"# Node"<<std::setw(10)<<"X Cordinate"<<std::setw(10)<<"Y coordinate"<<'\n';	
			
	for(int i=1;i<=m_nNodes;i++)
		{
			CVector<float>a(2);
			m_NodalData(i).GetCoords(a);
			m_FileOutput<<std::setw(15)<<i<<std::setw(15)<<a(1)<<std::setw(15)<<a(2)<<'\n';
		}
	m_FileOutput<<"____________________________________________________________________________________________\n";
	//Nodal fixities
	m_FileOutput<<"\n\n# Node"<<"    "<<"X-FIXITY "<<"    "<<"Y-FIXITY"<<"    "<<"MZ FIXITY"<<'\n';	
	for (int i=1; i <= m_nNodes; i++)
    {
		int nSN,nEN;
		CVector<CNode::Fixity> VFC(DOFPN);
		m_NodalData(i).GetFixity (VFC);
		m_FileOutput<<"   "<<i<<"       "<<VFC(1)<<"       "<<VFC(2)<<"       "<<VFC(3)<<'\n';	
		}
	m_FileOutput<<"____________________________________________________________________________________________\n";
	
	//NODAL  LOADS
	CVector<float>fNL(DOFPN);
	float ftemp;
	m_FileOutput<<"\n\n"<<"# Node"<<"    "<<"X-Force"<<"    "<<"Y-Force"<<"    "<<"Z Moment"<<"         "<<"Temperature"<<'\n';			
	for(int i=1;i<=m_nNodes;i++)
	{m_NodalLoadData(i).GetValues(fNL);
	m_NodalLoadData(i).GetTemp(ftemp);
	m_FileOutput<<"   "<<i<<"        "<<fNL(1)<<"        "<<fNL(2)<<"        "<<fNL(3)<<"         "<<ftemp; //CHANGED THIS
	m_FileOutput<<'\n';
	}
	m_FileOutput<<"____________________________________________________________________________________________\n";	
	//Element loads

	CElementLoads::ELType  ELT;
	int n;
	float f1,f2;
	m_FileOutput<<"\nElement"<<"    "<<"Load Type"<<std::setw(25)<<"    "<<"Value1"<<"    "<<"value2"<<'\n';		
	for(int i=1;i<=m_nElementLoads;i++)
	{
	m_ElementLoadData(i).GetValues(n,ELT,f1,f2);
	m_FileOutput<<n<<"           "<<ELT<<"            "<<f1<<"           "<<f2<<"\n";		
	}
	m_FileOutput<<"____________________________________________________________________________________________\n";
	//Element  material data
	int nSN=0;
	int nEN=0;
	int nMatGrp=0;
	m_FileOutput <<"\n\n"<<std::setw(5)<< "#Element"<<std::setw(5)<<"  START NODE"<<std::setw(5)<<"  END NODE"<<std::setw(5)<<"  ELASTICITY"<<std::setw(5)<<"  THERMAL COEFFICIENT"<<'\n';
	for (int i=1; i <= m_nElements; i++)
	{
		m_ElementData(i).GetENodes(nSN,nEN);
		nMatGrp = m_ElementData(i).GetMatPropertyGroup ();
        float fYM = m_MaterialData(nMatGrp).GetYM ();
        float fCTE = m_MaterialData(nMatGrp).GetCTE ();
		m_FileOutput<<std::setw(10)<< i<<std::setw(10)<<nSN<<std::setw(10)
			<<nEN<<std::setw(15)<<fYM<<std::setw(10)<<fCTE<<'\n';         
                                            
                  
	}
	m_FileOutput<<"____________________________________________________________________________________________\n";
	//Element property data             
	m_FileOutput << "\n\n"<<"#ELEMENT"<<"   X/S AREA"<<"   MOI about y axis:"<< "  MOI about z axis:"<< '\n';
	for (int i=1; i <= m_nElements; i++)
	{
		float fArea,fIyy,fIzz;
		CXSType *pXSGrp;
        pXSGrp = m_ElementData(i).GetEPropertyGroup ();
        pXSGrp->GetProperties(fArea, fIyy, fIzz);					//not getting why not similar to above line??
		m_FileOutput<<i<<std::setw(15)<<fArea<<std::setw(20)<<fIyy<<std::setw(20)<<fIzz<<'\n';
		}
	m_FileOutput<<"____________________________________________________________________________________________\n";
	//support reactions
	m_FileOutput<<"\n\n#NODE"<<"   "<<"X Reaction"<<"   "<<"Y Reaction"<<"   "<<"Moment"<<'\n';

	for(int i=1;i<=m_nNodes;i++)
	{
		int a=3*i-2;
		int b=3*i-1;
		int c=3*i;
		m_FileOutput<<i<<std::setw(10)<<supportreactions(a,1)<<"  "<<std::setw(15)<<supportreactions(b,1)<<"  "<<std::setw(10)<<supportreactions(c,1)<<'\n'; 

	}
	m_FileOutput<<"____________________________________________________________________________________________\n";
//Nodal DIsplacements
	m_FileOutput<<"\n\n"<<"#NODE"<<"   "<<"X Displacement"<<"   "<<"Y Displacement"<<"   "<<"Rotation about Z"<<'\n';
	for(int i=1;i<=m_nNodes;i++)
	{
		m_FileOutput<<i<<std::setw(18)<<m_SND(3*i-2,1)<<std::setw(18)<<m_SND(3*i-1,1)<<std::setw(18)<<m_SND(3*i,1)<<'\n'; 

	}	
	m_FileOutput<<"____________________________________________________________________________________________\n";
	//ELLEMENT NODAL FORCES										!!!!!!!!!!!!!!!!!!CHANGED THIS FOR REGRADING!!!!!!!!!!!!!!!!!!!
	CMatrix<double> ENF(DOFPE,1);
	m_FileOutput<<std::setw(20)<<"\nELEMENT NODAL FORCES IN GLOBAL CORDINATES\n";
	m_FileOutput<<std::setw(20)<<"\n-------------------------------------------\n";
	m_FileOutput<<"Fx start node"<<"  "<<"Fy start node"<<"  "<<"Mz start node"
				<<"    "<<"Fx end node"<<"    "<<"Fy end node"<<"    "<<"Mz end node\n";
	for (int i=1;i<=m_nElements;i++)
	{
		m_FileOutput<<"\nELEMENT : "<<i<<std::endl;
		m_ElementResponseData(i).GetESupportReaction(ENF);
		for(int j=1;j<=DOFPE;j++)
		{
			m_FileOutput<<ENF(j,1)<<std::setw(17);
		}
		m_FileOutput<<std::endl;
	}
	m_FileOutput<<"____________________________________________________________________________________________\n";
		// Element STRESSES
	for (int i=1;i<=m_nElements;i++)
	{
		double MT,MC,MS;
		m_ElementResponseData(i).GetStress(MT,MC,MS);
			m_FileOutput<<"\n\nMaximum Tensile Stress in Element"<<" "<<i<<" is :  "<<MT<<std::endl;
				m_FileOutput<<"Maximum Compressive Stress in Element"<<" "<<i<<" is :  "<<MC<<std::endl;
					m_FileOutput<<"Maximum Shear Stress in Element"<<" "<<i<<" is :  "<<MS<<"\n\n\n";
					
	}
	m_FileOutput<<"____________________________________________________________________________________________\n";
	//ERRORS
	m_FileOutput<<"\n\n Residual Error: "<<RError<<"\n";
	m_FileOutput<<"Absolute Error: "<<AError<<"\n";
m_FileOutput<<"____________________________________________________________________________________________\n";
}

void CFrame::IOErrorHandler (ERRORCODE ECode) const
// ---------------------------------------------------------------------------
// Function: displays error messages related to input data
// Input:    error code
// Output:   none
// ---------------------------------------------------------------------------
{
    std::cerr << '\n';

    if (ECode == NUMNODES) // invalid number of nodes
        std::cerr << "Number of nodes must be >= 2.";
    else if (ECode == NUMELEMENTS) // invalid number of elements
        std::cerr << "Number of elements must be >= 1.";
    else if (ECode == DEBUGCODE) // invalid debug level
        std::cerr << "Debug level must be 0 or 1.";
    else if (ECode == NODENUMBER) // invalid node number
        std::cerr << "Invalid node number.";
    else if (ECode == ELEMENTNUMBER) // invalid element number
        std::cerr << "Invalid element number.";
    else if (ECode == XSAREA) // invalid x/s area
        std::cerr << "Area must be positive.";
    else if (ECode == YOUNGSMODULUS) // invalid E
        std::cerr << "Modulus of elasticity must be positive.";
    else if (ECode == INVALIDINPUT) // invalid input
        std::cerr << "Invalid input.";
    else if (ECode == INVALIDLASTLINE) // invalid input
        std::cerr << "Input file needs *end as last line.";
    else if (ECode == NODALFIXITY) // invalid fixity code
        std::cerr << "Nodal fixity code must be 'free' or 'specified'.";
    else if (ECode == NUMMATGROUPS) // invalid number of material groups
        std::cerr << "Number of material groups must be >= 1.";
    else if (ECode == NUMEPROPGROUPS) // invalid number of property groups
        std::cerr << "Number of element property groups must be >= 1.";
    else if (ECode == EPROPNUMBER) // invalid element property group
        std::cerr << "Invalid element property group number.";
    else if (ECode == XSDIMENSION) // invalid x/s dimension
        std::cerr << "Invalid cross-section dimension.";
    else if (ECode == MATGROUPNUMBER) // invalid material group
        std::cerr << "Invalid material property group number.";
    else if (ECode == MATPROPERTY) // invalid material property
        std::cerr << "Invalid material property.";
    else if (ECode == ELOADTYPE) // invalid element load type
        std::cerr << "Invalid element load type.";
    else if (ECode == XSTYPE) // invalid x/s type
        std::cerr << "Invalid cross-section type.";
    else
        std::cerr << "Unknown error ...?";

    std::cerr << '\n' << "Error in input file line : " << m_nLineNumber;
    std::cerr << std::endl;

    exit (1);
}

void CFrame::TerminateProgram ()
// ---------------------------------------------------------------------------
// Function: terminates the program steps by closing the input/output files
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    // close the input and output files
    m_FileInput.close ();
    m_FileOutput.close ();
}
