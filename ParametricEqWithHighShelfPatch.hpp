////////////////////////////////////////////////////////////////////////////////////////////////////

/*
 
 
 LICENSE:
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 
 */
 
 /*
     registerParameter(PARAMETER_A, "PE_Freq", "PE_Freq");
    registerParameter(PARAMETER_B, "PE_Q", "PE_Q");
    registerParameter(PARAMETER_C, "PE_Gain","PE_Gain);
    registerParameter(PARAMETER_D, "Treb_Gain", "Treb_Gain");
    registerParameter(PARAMETER_E, "FreqPedal", "FreqPedal");
    */

/* created by the OWL team 2013 */

////////////////////////////////////////////////////////////////////////////////////////////////////


#ifndef __ParametricEqWithHighShelfPatch_hpp__
#define __ParametricEqWithHighShelfPatch_hpp__


enum filterType {
    PEQ, // Parametric EQ
    HSH, // High Shelf
    LSH  // Low SHelf
};
#define Q_BUTTERWORTH   0.707

class BiquadDF1 {
public:
    BiquadDF1() {}
    ~BiquadDF1() {}
    
  void initStateVariables(){
    x1=0.f;
    x2=0.f;
    y1=0.f;
    y2=0.f;
  }
    
  // function used for PEQ, HSH, LSH
  void setCoeffs(float normalizedFrequency, float Q, float dbGain){        
    float alpha, c, omega, d, e, gamma, beta ;        
    omega = 2*M_PI*normalizedFrequency ;
    c = cosf(omega) ;
    alpha = sinf(omega)/(2*Q);
    d = powf(10,dbGain/40.f);
    gamma = alpha*powf(10,fabsf(dbGain)/40.f);
    e = powf(10,fabsf(dbGain)/20.f);
    beta = 2*alpha*powf(e,0.5f);        
    switch (fType)
      {
      case PEQ: // Parametric EQ
	a[0]=1+gamma/d;
	a[1]=-2*c/a[0];
	a[2]=(1-gamma/d)/a[0];
	b[0]=(1+gamma*d)/a[0];
	b[1]=a[1];
	b[2]=(1-gamma*d)/a[0];
	a[0]=1;
	break;                
      case HSH: // High Shelf
	if (dbGain >0){
	  a[0]=2*(1+alpha);
	  a[1]=-4*c/a[0];
	  a[2]=2*(1-alpha)/a[0];
	  b[0]=((1+e)-(1-e)*c+beta)/a[0];
	  b[1]=2*((1-e)-(1+e)*c)/a[0];
	  b[2]=((1+e)-(1-e)*c-beta)/a[0];
	  a[0]=1;
	}
	else {
	  a[0]=(1+e)-(1-e)*c+beta;
	  a[1]=2*((1-e)-(1+e)*c)/a[0];
	  a[2]=((1+e)-(1-e)*c-beta)/a[0];
	  b[0]=2*(1+alpha)/a[0];
	  b[1]=-4*c/a[0];
	  b[2]=2*(1-alpha)/a[0];
	  a[0]=1;
	}
	break;                
      case LSH: // Low Shelf
	if (dbGain >0){
	  a[0]=2*(1+alpha);
	  a[1]=-4*c/a[0];
	  a[2]=2*(1-alpha)/a[0];
	  b[0]=((1+e)+(1-e)*c+beta)/a[0];
	  b[1]=-(2*((1-e)+(1+e)*c))/a[0];
	  b[2]=((1+e)+(1-e)*c-beta)/a[0];
	  a[0]=1;                    
	}
	else {
	  a[0]=(1+e)+(1-e)*c+beta;
	  a[1]=-2*((1-e)+(1+e)*c)/a[0];
	  a[2]=((1+e)+(1-e)*c-beta)/a[0];
	  b[0]=(2*(1+alpha))/a[0];
	  b[1]=-4*c/a[0];
	  b[2]=2*(1-alpha)/a[0];
	  a[0]=1;
	}
	break;
      }        
  }

  void copyCoeffs(BiquadDF1& other){
    memcpy(a, other.a, sizeof(float)*3);
    memcpy(b, other.b, sizeof(float)*3);
  }

  void process (int numSamples, float* buf){
    float out;
    for (int i=0;i<numSamples;i++){
      out = b[0]*buf[i]+b[1]*x1+b[2]*x2-a[1]*y1-a[2]*y2 ;
      y2 = y1;
      y1 = out;
      x2 = x1;
      x1 = buf[i];
      buf[i]=out;
    }
  }
    
  void setType (filterType typ){
    fType = typ;
  }

private:
    float a[3] ; // ai coefficients
    float b[3] ; // bi coefficients
    float x1, x2, y1, y2 ; // state variables to compute samples
    filterType fType;
};


/**
 * Biquad Parametric EQ filter class
 */
class Biquad1 {
public:
  Biquad1() {}
  ~Biquad1() {}
    
  void initStateVariables(){
        x1=0.f;
        x2=0.f;
        y1=0.f;
        y2=0.f;
    }
    
  void setCoeffsPEQ(float normalizedFrequency, float Q, float dbGain) {
    // Compute the filters coefficients a[i] and b[i];
    float omega, c, alpha, d, gamma;
    omega = 2*M_PI*normalizedFrequency ;
    c = cosf(omega) ;
    alpha = sinf(omega)/(2*Q);
    d = powf(10,dbGain/40);
    gamma = alpha*powf(10,fabsf(dbGain)/40);
      
    a[0] = 1+gamma/d;
    a[1] = -2*c/a[0];
    a[2] = (1-gamma/d)/a[0];
    b[0] = (1+gamma*d)/a[0];
    b[1] = a[1];
    b[2] = (1-gamma*d)/a[0];
    a[0] = 1.0;
  }
    
  void process(int numSamples, float* input, float* out){
    // process a block of more than 2 samples. Basic implementation without coeffs interpolation.
    out[0] = b[0]*input[0]+b[1]*x1+b[2]*x2-a[1]*y1-a[2]*y2 ;
    out[1] = b[0]*input[1]+b[1]*input[0]+b[2]*x1-a[1]*out[0]-a[2]*y1 ;
    for(int i=2; i<numSamples; i++){
      out[i] = b[0]*input[i]+b[1]*input[i-1]+b[2]*input[i-2]-a[1]*out[i-1]-a[2]*out[i-2] ;
    }
      
    // store values for next block
    x1 = input[numSamples-1];
    x2 = input[numSamples-2];
    y1 = out[numSamples-1];
    y2 = out[numSamples-2];
  }
    
  void process (int numSamples, float* buf){
    float out;
    for (int i=0;i<numSamples;i++){
        out = b[0]*buf[i]+b[1]*x1+b[2]*x2-a[1]*y1-a[2]*y2 ;
        y2 = y1;
        y1 = out;
        x2 = x1;
        x1 = buf[i];
        buf[i]=out;
    }
  }
    
private:
  float a[3] ; // ai coefficients
  float b[3] ; // bi coefficients
  float x1, x2, y1, y2 ; // state variables to compute samples
};

class FourBandsEq {
private:
    BiquadDF1 band1; //, band2, band3, band4; // filters
    float fn1; //, fn2, fn3, fn4; // cutoffs frequencies, normalized
public:
  FourBandsEq(double samplerate) {
    band1.initStateVariables();
    band1.setType(HSH);
    fn1=3000/samplerate;
      
  /*  band2.initStateVariables();
    band2.setType(PEQ);
    fn2=250/samplerate;
      
    band3.initStateVariables();
    band3.setType(PEQ);
    fn3=1500/samplerate;
      
    band4.initStateVariables();
    band4.setType(PEQ);
    fn4=4000/samplerate;*/
  }

  void setCoeffs(float a){
    // update filter coefficients
    band1.setCoeffs(fn1, Q_BUTTERWORTH, a);
    //band2.setCoeffs(fn2, Q_BUTTERWORTH, b);
    //band3.setCoeffs(fn3, Q_BUTTERWORTH, c);
    //band4.setCoeffs(fn4, Q_BUTTERWORTH, d);
  }

  void copyCoeffs(FourBandsEq& other){
    band1.copyCoeffs(other.band1);
    //band2.copyCoeffs(other.band2);
    //band3.copyCoeffs(other.band3);
    //band4.copyCoeffs(other.band4);
  }

  void process(int numSamples, float* buf){      
    // process
    band1.process(numSamples, buf);
    //band2.process(numSamples, buf);
    //band3.process(numSamples, buf);
    //band4.process(numSamples, buf);
  }
};

/**
 * Parametric EQ OWL Patch
 */
class ParametricEqWithHighShelfPatch : public Patch {
public:
  ParametricEqWithHighShelfPatch() : eqL(getSampleRate()), eqR(getSampleRate()) {
    registerParameter(PARAMETER_A, "PE_Freq", "PE_Freq");
    registerParameter(PARAMETER_B, "PE_Q", "PE_Q");
    registerParameter(PARAMETER_C, "PE_Gain","PE_Gain");
    registerParameter(PARAMETER_D, "Treb_Gain", "Treb_Gain"); //aka "HIGH"s
    registerParameter(PARAMETER_E, "FreqPedal", "FreqPedal");
    peqL.initStateVariables();
    peqR.initStateVariables();
    
  }    
     

  void processAudio(AudioBuffer &buffer){
    // update filter coefficients
    float fn = getFrequency()/getSampleRate();
    float Q = getQ();
    float g = getDbGain();
    float a= getDbGain2(PARAMETER_D);
    
    peqL.setCoeffsPEQ(fn, Q, g) ;
    peqR.setCoeffsPEQ(fn, Q, g) ;
    
    eqL.setCoeffs(a);
    eqR.copyCoeffs(eqL);
      
    // process
    int size = buffer.getSize();
    float* left = buffer.getSamples(0);
    peqL.process(size, left);
    float* right = buffer.getSamples(1);
    peqR.process(size, right);
    
    eqL.process(size, left);
    eqR.process(size, right);
  }
    
private:
  Biquad1 peqL ; // PEQ filter
  Biquad1 peqR ; // PEQ filter
  
  FourBandsEq eqL;
  FourBandsEq eqR;

  float getFrequency() {
    //float f = getParameterValue(PARAMETER_A)+getParameterValue(PARAMETER_E)/2;
    // param_A = 0    <-> f=50;
    // param_A = 1    <-> f=10050;
//      return powf(10,3*f+1)+40;
//      return (f*8000)+50;
    float f = getParameterValue(PARAMETER_A);
    return (f*4000)+100;  //100 to 4000
  }
        
  float getQ(){
    float q = getParameterValue(PARAMETER_B);
    // param_B = 0    <-> Q=0.5
    // param_B = 1    <-> Q=10
    return q*9.5+0.5;
  }
    
  float getDbGain(){
    float linGain = getParameterValue(PARAMETER_C);
    // linGain = 0    <-> -15 dB
    // linGain = 0.5  <-> 0dB
    // linGain = 1    <-> 15dB
    //return (linGain-0.5)*30;
    return (linGain-1.0)*30;  //-30 to 30 dB
  }
  
  float getDbGain2(PatchParameterId id){
    float linGain = getParameterValue(id);
    // linGain = 0    <-> -15 dB
    // linGain = 0.5  <-> 0dB
    // linGain = 1    <-> 15dB
    return (linGain-0.5)*30;
  }
};

#endif // __ParametricEqPatch_hpp__
