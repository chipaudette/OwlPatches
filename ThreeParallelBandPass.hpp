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


/* Original "StateVariableFilterPatch" created by the OWL team 2013 */
/* Adapted by Chip Audette to "ThreeParallelBandPass" 2021 */


////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef __ThreeParallelBandPass_hpp__
#define __ThreeParallelBandPass_hpp__

//include "SampleBasedPatch.hpp"

/**
State variable Filter
http://musicdsp.org/showArchiveComment.php?ArchiveID=23

Type : 12db resonant low, high or bandpass
References : Effect Deisgn Part 1, Jon Dattorro, J. Audio Eng. Soc., Vol 45, No. 9, 1997 September

Notes : 
Digital approximation of Chamberlin two-pole low pass. Easy to calculate coefficients, easy to process algorithm.

Code : 
cutoff = cutoff freq in Hz
fs = sampling frequency //(e.g. 44100Hz)
f = 2 sin (pi * cutoff / fs) //[approximately]
q = resonance/bandwidth [0 < q <= 1]  most res: q=1, less: q=0
low = lowpass output
high = highpass output
band = bandpass output
notch = notch output
scale = q
low=high=band=0;
//--beginloop
low = low + f * band;
high = scale * input - low - q*band;
band = f * high + band;
notch = high + low;
//--endloop
*/

class SampleBasedPatch : public Patch {
public:
  virtual void prepare() = 0;
  virtual float processSample(float sample) = 0;
  void processAudio(AudioBuffer &buffer){
    prepare();
    int size = buffer.getSize();
    float* samples = buffer.getSamples(0); // This Class is Mono (1in, 1out)
      for(int i=0; i<size; ++i){
          samples[i] = processSample(samples[i]);
      }
  }
};

class StateVariableFilterPatch : public SampleBasedPatch {
private:
  float low[3], band[3];
  float f[3], q;
  float gain;
public:
  StateVariableFilterPatch() {
    registerParameter(PARAMETER_A, "Fc1");
    registerParameter(PARAMETER_B, "Fc2");
    registerParameter(PARAMETER_C, "Fc3");
    registerParameter(PARAMETER_D, "Q");
    //registerParameter(PARAMETER_B, "Q");
    //registerParameter(PARAMETER_C, "");
    //registerParameter(PARAMETER_D, "Gain");
	
	//initialize states
	for (int i=0; i<3; i++) {
		low[i] = 0.0;
		band[i]=0.0;
		f[i]=1000.0;
	}
  }
  void prepare(){
    float fc[3];
    fc[0] = getParameterValue(PARAMETER_A);
    fc[1] = getParameterValue(PARAMETER_B);
    fc[2] = getParameterValue(PARAMETER_C);
    q = getParameterValue(PARAMETER_D);
    //gain = getParameterValue(PARAMETER_D); // get gain value
    gain = 1.0;

    // fc = cutoff freq in Hz
    // fs = sampling frequency //(e.g. 44100Hz)
    // q = resonance/bandwidth [0 < q <= 1]  most res: q=1, less: q=0

	for (int i=0; i<3; i++) {
		fc[i] /= 2;
		f[i] = sin(M_PI * fc[i]);
	}

    q = 1 - q;
  }
  float bandpass(float sample, int ind) {
    low[ind] = low[ind] + f[ind] * band[ind];
    float high[ind] = q * sample - low[ind] - q*band[ind];
    band = f[ind] * high[ind] + band[ind];
    //return gain*low[ind];
	return gain*band[ind];
  }	
  float processSample(float sample){
	float out_val = 0.0;
	for (int i=0; i<3; i++) {
		out_val += bandpass(sample, i);
	}
    return out_val;
  }
};

#endif /* __StateVariableFilterPatch_hpp__ */
