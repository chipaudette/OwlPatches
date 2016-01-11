////////////////////////////////////////////////////////////////////////////////////////////////////

/*
Title: EmptyPatch
Author: Chip Audette
Created: Jan 2016
License: The MIT License
*/


////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef __EmptyPatch_h__
#define __EmptyPatch_h__

#include "StompBox.h"

class EmptyPatch : public Patch {
public:
  EmptyPatch(){
    registerParameter(PARAMETER_A, "");    
    registerParameter(PARAMETER_B, "");    
    registerParameter(PARAMETER_C, "");    
    registerParameter(PARAMETER_D, "");    
  }
  void processAudio(AudioBuffer &buffer){
    //float gain = getParameterValue(PARAMETER_A)*2;
    int size = buffer.getSize();
    for(int ch=0; ch<buffer.getChannels(); ++ch){
      float* buf = buffer.getSamples(ch);
      for(int i=0; i<size; ++i)
	//buf[i] = gain*buf[i];
    }
  }
};

#endif // __EmptyPatch_h__


////////////////////////////////////////////////////////////////////////////////////////////////////


