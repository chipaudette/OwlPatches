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
    //return (f*0.67)+0.33;  //0.33 to 1.0
    //return (f*1.67)+0.33;  //0.33 to 1.0
    float scale_fac = 1.0;
    if (f < 0.5) {
    	return (f*2.0)*0.67 + 0.33; //0.33 to 1.0
    } else {
    	return ((f-0.5)*2.0)*1.5+1.0;  //1.0 to 2.5
    }
  }
  float getGainScale(PatchParameterId id){
    return getTrebGain(id)/2.0; //should be 1/4.5 to 1.0 to 4.5
  }
  
  float getTrebleCutoff(PatchParameterId id) {
  	  float treb_cutoff = getParameterValue(id); //0.0 to 1.0
  	  return (treb_cutoff*4000.0)+0000.0; //0000 to 40000
  }
  float getTrebGain(PatchParameterId id) {
     float val = getParameterValue(id); //0 to +1.0
     bool is_boost = true;
     if (val < 0.5) is_boost = false;
     float treb_gain = abs((val - 0.5)*2.0); //0.0 to +1.0
     treb_gain *= 2.0;  //0.0 to 2.0
     treb_gain += 1.0;  //1.0 to 3.0;
     treb_gain *= treb_gain; //1 to 9
     if (is_boost) {
     	 return max(1.0,min(9.0,treb_gain));
     } else {
     	 return min(1./treb_gain,1.0);
     }
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
        
    //define some variables
    int N_2 = fftSize/2+1;
    int source_ind, neg_dest_ind;
    float source_ind_float, interp_fac;
    float new_mag, scale;
    int max_source_ind = (int)(((float)N_2)*(10000.0 / (48000.0/2.0))); //highest frequency bin to grab from (Assuming 48kHz sample rate)
    
    //check to see if bypass button is pressed
    if (isButtonPressed(PUSHBUTTON)) {
    	//bypass all freq-domain prcessing
    	
    } else {
    	//do the freq-domain processing
    	
		//get the original magnitudes
		for (int i=0; i<N_2; i++) orig_mag[i] = ca.mag(i);
    
		//compute the new magnitude and apply
		for (int dest_ind=1; dest_ind < N_2; dest_ind++) {
		  //what is the source bin for the new magnitude
		  source_ind_float = (((float)dest_ind)/shift_scale_fac) + 0.5;
		  //source_ind = (int)(source_ind_float+0.5);  //no interpolation but round to the neariest index
		  //source_ind = min(max(source_ind,1),N_2-1);
		  source_ind = min(max(1,(int)source_ind_float),N_2-2);
		  interp_fac = source_ind_float - (float)((int)source_ind);
		  interp_fac = max(0.0,interp_fac);
	
		  
		  //what is the new magnitude
		  new_mag = 0.0;scale = 0.0;
		  if (source_ind < max_source_ind) {
		  
			//new_mag=orig_mag[source_ind];  //the magnitude that we desire
			//scale = new_mag / orig_mag[dest_ind];//compute the scale factor
			new_mag = orig_mag[source_ind];
			new_mag += interp_fac*(orig_mag[source_ind]-orig_mag[source_ind+1]);
			scale = new_mag/orig_mag[dest_ind];
			
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
