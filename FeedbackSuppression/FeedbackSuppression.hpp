#ifndef __FastFourierTestPatch_hpp__
#define __FastFourierTestPatch_hpp__

#include "StompBox.h"
#include "FastFourierTransform.h"

#ifndef PI
#define PI 3.141592653589793238
#endif

#define N_BLOCK_BUFF (4)
class FeedbackSuppression : public Patch {
private:
  int blockSize,fftSize;
  int blockBuff_counter;
  int lower_freq_limit_bin=0, upper_freq_limit_bin = 0;
  int bin_to_cut = 0;
  float sampleRate_Hz = 48000.0;  //is reset to correct value in constructor
  FloatArray blockBuff[N_BLOCK_BUFF], outBlockBuff[N_BLOCK_BUFF];
  FloatArray window, full_data, orig_mag, ave_mag;
  FastFourierTransform transform;
  ComplexFloatArray ca;
  float getSNRThresh(PatchParameterId id) {
    float val = getParameterValue(id);
    val = (val*5.0)+1.0; //1.0 to 6.0
    return val;  //1.0 to 36.0
  }
  float getGainScale(PatchParameterId id){
    float gain = getParameterValue(id);
    gain = (gain*3.0)+1.0; //1.0 to 4.0
    return gain*gain;  //1.0 to 16.0
  }
  
  float getTrebleCutoff(PatchParameterId id) {
  	  float treb_cutoff = getParameterValue(id); //0.0 to 1.0
  	  return (treb_cutoff*4000.0)+0000.0; //0000 to 40000
  }
  float getTrebGain(PatchParameterId id) {
     float treb_gain = getParameterValue(id);
     treb_gain = (treb_gain*3.0)+1.0;   //1.0 to 4.0
     return treb_gain*treb_gain; //1.0 to 16.0
  }
  
public:
  FeedbackSuppression(){
    //constructor
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
    ave_mag = FloatArray::create(fftSize/2+1);
    registerParameter(PARAMETER_A, "SNR Thresh", "SNR Thresh");
    registerParameter(PARAMETER_B, "Treb_Freq","Treb_Freq");
    registerParameter(PARAMETER_C, "Treb_Gain","Treb_Gain");
    registerParameter(PARAMETER_D, "Gain"); 
    
    sampleRate_Hz = getSampleRate();
    lower_freq_limit_bin = (int)(fftSize * (500.0 / sampleRate_Hz)); //500 Hz
    upper_freq_limit_bin = (int)(fftSize * (5000.0 / sampleRate_Hz)); //5000 Hz
    
    //fill window with values
    for (int i=0; i<fftSize; i++) {
      window[i] = -0.5*cos(2.0*PI*i/fftSize) + 0.5;
    }
  }
  ~FeedbackSuppression(){
    ComplexFloatArray::destroy(ca);
    FloatArray::destroy(window);
    FloatArray::destroy(full_data);
    FloatArray::destroy(ave_mag);
    FloatArray::destroy(orig_mag);
    for (int Ibuff = 0; Ibuff < N_BLOCK_BUFF; Ibuff++) {
			FloatArray::destroy(blockBuff[Ibuff]);
			FloatArray::destroy(outBlockBuff[Ibuff]);
    }
  }
  void processAudio(AudioBuffer &buffer){
    //get values from knobs
    float SNR_thresh = getSNRThresh(PARAMETER_A);
    float treb_cutoff_Hz = getTrebleCutoff(PARAMETER_B);
    int treb_cutoff_ind = (int)(((float)fftSize)*treb_cutoff_Hz/sampleRate_Hz);
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
    
    //loop over the bins
    int N_2 = fftSize/2+1;
    float learn_fac = 0.05;
    float unlearn_fac = 1.0 - 0.05 - 1e-6;
    float max_SNR = 0.0;
    int max_SNR_ind = 0;
    float SNR;
    for (int i=0; i<N_2; i++) {
    	//get the current magnitude
    	orig_mag[i] = ca.mag(i);
    	
    	//if we're within the allowed freq bounds, compute what's needed to see if we cut this bin
    	if ((i >= lower_freq_limit_bin) && (i <= upper_freq_limit_bin)) { 
				//compute SNR-like value relative to the running ave
				SNR = orig_mag[i]/ave_mag[i];
				
				//is this the biggest SNR?
				if (SNR > max_SNR) {
					max_SNR = SNR;
					max_SNR_ind = i;
				}
			}
    
    	//update the running average (mean) magnitude
    	ave_mag[i] = unlearn_fac * ave_mag[i] + learn_fac * orig_mag[i];
    }
    
    //decide whether to update which bin is being notched
    if (max_SNR > SNR_thresh) bin_to_cut = max_SNR_ind;
      
    //compute the new magnitude and apply
    int source_ind, neg_dest_ind;
    float new_mag, scale;
    int max_source_ind = (int)(((float)N_2)*(10000.0 / (sampleRate_Hz/2.0))); //highest frequency bin to grab from (Assuming 48kHz sample rate)
    for (int dest_ind=1; dest_ind < N_2; dest_ind++) { //skip the zero bin
      source_ind = dest_ind;
      
      //set the new magnitude
      if (source_ind < max_source_ind) {
        //by default no boost or cut, yet
        scale = 1.0;
        
        //is this the bin to cut?
        if (source_ind == (bin_to_cut - 1)) scale = 0.3*scale;  //cut it a bit (10 dB)
        if (source_ind == bin_to_cut) scale = 0.03*scale;       //cut it big time  (33 dB)
        if (source_ind == (bin_to_cut + 1)) scale = 0.3*scale;  //cut it a bit (10 dB)
        
        //also gain some gain to the higher frequencies
        if (dest_ind > treb_cutoff_ind) scale = scale * treb_gain;  //apply boost to treble
        
        //apply the scale factor
        ca[dest_ind].re *= scale;
        ca[dest_ind].im *= scale;
      } else {
      	//zero everything else
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

#endif 
