#ifndef MEGAISTATE_H
#define MEGAISTATE_H

#include "ReadWriteBaseState.h"

class MegaInterleavedState : public ReadWriteBaseState
{
public:
    
    MegaInterleavedState(ReadWriteMS* MachineState) { Machine = MachineState; name="mega_interleaved"; extension="mega"; canLoad=true; };
    
    virtual int CheckAlignment(istream* origin);
    virtual newAlignment* LoadAlignment(string filename);
    virtual bool SaveAlignment(newAlignment* alignment, ostream* output, std::string* FileName);
    virtual bool RecognizeOutputFormat(std::string FormatName);
     
};

#endif // MEGAISTATE_H