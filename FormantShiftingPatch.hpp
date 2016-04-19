#ifndef __FastFourierTestPatch_hpp__
#define __FastFourierTestPatch_hpp__

#include "StompBox.h"
#include "FastFourierTransform.h"

class FastFourierTestPatch : public Patch {
private:
  FastFourierTransform transform;
  ComplexFloatArray ca;
  float getShiftScaleFac(PatchParameterId id) {
    float f = getParameterValue(id);  //0.0 to 1.0
    return (f*0.67)+0.33;  //0.33 to 1.0
  }
  
public:
  FormantShiftingPatch(){
    int fftSize = getBlockSize();
    transform.init(fftSize);
    ca = ComplexFloatArray::create(fftSize);
    registerParameter(PARAMETER_A, "Formant_Shift", "Formant_Shift");
  }
  ~FastFourierTestPatch(){
    ComplexFloatArray::destroy(ca);
  }
  void processAudio(AudioBuffer &buffer){
    float scale_fac = getShiftScaleFac(PARAMETER_A);
    
    FloatArray buf = buffer.getSamples(0);
    transform.fft(buf, ca);
    transform.ifft(ca, buf);
  }
};

#endif // __FastFourierTestPatch_hpp__
