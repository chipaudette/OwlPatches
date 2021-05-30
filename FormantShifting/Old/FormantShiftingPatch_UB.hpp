#ifndef __FastFourierTestPatch_hpp__
#define __FastFourierTestPatch_hpp__

#include "StompBox.h"
#include "FastFourierTransform.h"

#ifndef PI
#define PI 3.141592653589793238
#endif

#define N_BLOCK_BUFF (4)
class FormantShiftingPatch : public Patch {
private:
  int blockSize,fftSize;
  int blockBuff_counter;
  FloatArray blockBuff[N_BLOCK_BUFF], outBlockBuff[N_BLOCK_BUFF];
  FloatArray window, full_data, orig_mag;
  FastFourierTransform transform;
  ComplexFloatArray ca;
  float getShiftScaleFac(PatchParameterId id) {
    float f = getParameterValue(id);  //0.0 to 1.0
    return (f*0.67)+0.33;  //0.33 to 1.0
  }
  float getGainScale(PatchParameterId id){
    float gain = getParameterValue(id);
    gain = (gain*3.0)+1.0;
    return gain*gain;  //1.0 to 16.0
  }
  
  float getTrebleCutoff(PatchParameterId id) {
  	  float treb_cutoff = getParameterValue(id); //0.0 to 1.0
  	  return (treb_cutoff*4000.0)+0000.0; //0000 to 40000
  }
  float getTrebGain(PatchParameterId id) {
     float treb_gain = getParameterValue(id);
     treb_gain = (treb_gain*3.0)+1.0;
     return treb_gain*treb_gain;
  }
  
public:
  FormantShiftingPatch(){
  	  
    blockSize = getBlockSize();
    fftSize = 4*getBlockSize();
    blockBuff_counter = 0;
    for (int Ibuff = 0; Ibuff < N_BLOCK_BUFF; Ibuff++) {
        blockBuff[Ibuff] = FloatArray::create(blockSize);
        outBlockBuff[Ibuff] = FloatArray::create(blockSize);
    }
    transform.init(fftSize);
    ca = ComplexFloatArray::create(fftSize);
    window = FloatArray::create(fftSize);
    full_data = FloatArray::create(fftSize);
    orig_mag = FloatArray::create(fftSize/2+1);
    registerParameter(PARAMETER_A, "Formant_Shift", "Formant_Shift");
    registerParameter(PARAMETER_B, "Treb_Freq","Treb_Freq");
    registerParameter(PARAMETER_C, "Treb_Gain","Treb_Gain");
    registerParameter(PARAMETER_D, "Gain"); 
    
    //fill window with values
    for (int i=0; i<fftSize; i++) {
      window[i] = -0.5*cos(2.0*PI*i/fftSize) + 0.5;
    }
  }
  ~FormantShiftingPatch(){
    ComplexFloatArray::destroy(ca);
    FloatArray::destroy(window);
    FloatArray::destroy(full_data);
    FloatArray::destroy(orig_mag);
    for (int Ibuff = 0; Ibuff < N_BLOCK_BUFF; Ibuff++) {
        FloatArray::destroy(blockBuff[Ibuff]);
        FloatArray::destroy(outBlockBuff[Ibuff]);
    }
  }
  void processAudio(AudioBuffer &buffer){
    //get values from knobs
    float shift_scale_fac = getShiftScaleFac(PARAMETER_A);
    float treb_cutoff_Hz = getTrebleCutoff(PARAMETER_B);
    int treb_cutoff_ind = (int)(((float)fftSize)*treb_cutoff_Hz/48000.0);
    float treb_gain = getTrebGain(PARAMETER_C);
    float gain_scale = getGainScale(PARAMETER_D);   
    
    
    //get the latest data
    FloatArray buf = buffer.getSamples(0);
    
    //store this buffer
    blockBuff_counter++;  //increment
    if (blockBuff_counter >= N_BLOCK_BUFF) blockBuff_counter = 0;  //wrap around
    buf.copyTo(blockBuff[blockBuff_counter]); //copy
    
    //form full record for sending to FFT
    int block_ind;
    for (int i = 0; i < N_BLOCK_BUFF; i++) {
    	  block_ind = (blockBuff_counter+1+i) % N_BLOCK_BUFF;  //oldest first
    	  for (int j = 0; j < blockSize; j++) {
    	      full_data[i*blockSize+j] = blockBuff[block_ind][j];
    	  }
    }
    
    //apply window
    full_data.multiply(window);
    
    //compute FFT
    transform.fft(full_data, ca);
    
    //get the original magnitudes
    int N_2 = fftSize/2+1;
    for (int i=0; i<N_2; i++) orig_mag[i] = ca.mag(i);
    
    //compute the new magnitude and apply
    int source_ind, neg_dest_ind;
    float new_mag, scale;
    //int max_source_ind = N_2;
    int max_source_ind = (int)(((float)N_2)*(10000.0 / (48000.0/2.0))); //highest frequency bin to grab from (Assuming 48kHz sample rate)
    for (int dest_ind=1; dest_ind < N_2; dest_ind++) {
      //what is the source bin for the new magnitude
      source_ind = (int)((((float)dest_ind)/shift_scale_fac) + 0.5);  //round
      
      //what is the new magnitude
      new_mag = 0.0;scale = 0.0;
      if (source_ind < max_source_ind) {
      
        //new_mag=orig_mag[source_ind];  //the magnitude that we desire
        //scale = new_mag / orig_mag[dest_ind];//compute the scale factor
        scale = orig_mag[source_ind]/orig_mag[dest_ind];
        
        //also gain some gain to the higher frequencies
        if (dest_ind > treb_cutoff_ind) {
            scale = scale * treb_gain;
        }
        
        
        //apply scale factor
        ca[dest_ind].re *= scale;
        ca[dest_ind].im *= scale;
      } else {
        ca[dest_ind].re = 0.0;
        ca[dest_ind].im = 0.0;
      }
      
      //set the complex conjugate for the negative freuqency side
      neg_dest_ind = fftSize - dest_ind - 1;
      if (neg_dest_ind >= N_2) {
          ca[neg_dest_ind].re = ca[dest_ind].re;
          ca[neg_dest_ind].im = ca[dest_ind].im;
      }
    }    
    
    //compute IFFT
    transform.ifft(ca, full_data); //"buf" is the output
    
    //apply window again
    full_data.multiply(window);
    
    //overlap and add
    int count=0;
    for (int i=0; i < blockSize; i++) {
         buf[i] = gain_scale* (outBlockBuff[0][i]+full_data[count]); //this goes out
         count++;
    }
    
    //store output blocks for next round
    for (int Iblock = 0; Iblock < (N_BLOCK_BUFF-1); Iblock++) {
        for (int i=0; i < blockSize; i++) {
            outBlockBuff[Iblock][i] = outBlockBuff[Iblock+1][i]+full_data[count]; //overlap and add
            count++;
        }
    }
    
    
  }
};

#endif // __FastFourierTestPatch_hpp__
