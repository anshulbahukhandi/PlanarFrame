#ifndef __ANSHUL_CIRCSOLID_H__
#define __ANSHUL_CIRCSOLID_H__

#include "F:\ASU Courses\cee532\Library\Library\vectortemplate.h"
#include "xstype.h"
const int numCircDimensions = 1;

class CCircSolid: public CXSType
{
    public:
        CCircSolid (const CVector<float>& fV);
        CCircSolid (const CCircSolid&);
        ~CCircSolid ();

        // helper functions
        virtual void ComputeProperties ();

    private:
};

#endif