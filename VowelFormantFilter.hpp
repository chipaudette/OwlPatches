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
/* Adapted by Chip Audette to "VowelFormantFilter" 2021 */


////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef __VowelFormantFilter_hpp__
#define __VowelFormantFilter_hpp__


/**

Built from parallel bandpass filters

Bandpass filters are built from State variable Filter:
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

class VowelFormantFilter : public SampleBasedPatch {
	public:
		VowelFormantFilter(void) {
			registerParameter(PARAMETER_A, "Vowel"); //will be 0.0 to 1.0
			registerParameter(PARAMETER_B, "Q"); //will be 0.0 to 1.0
			registerParameter(PARAMETER_C, "Model"); //will be 0.0 to 1.0
			registerParameter(PARAMETER_D, "gain");

			//initialize states
			overall_gain = 1.0;
			model = -1;
			for (int i=0; i<3; i++) {
				low[i] = 0.0;
				band[i]=0.0;
				f[i]=1000.0;
				gain[i]=0.0;
			}
		}
		void prepare(void){
			float fc[3];
			float vowel = getParameterValue(PARAMETER_A); //a value of 1.0 means fc = sample rate
			q = getParameterValue(PARAMETER_B); //a value of 1.0 means fc = sample rate
			float new_model = getParameterValue(PARAMETER_C);  
			overall_gain = getParameterValue(PARAMETER_D);
			
			//update the formant model
			chooseModel(new_model);
			

			// fc = cutoff freq in Hz 
			// fs = sampling frequency //(e.g. 44100Hz)
			// q = resonance/bandwidth [0 < q <= 1]  most res: q=1, less: q=0

			//map vowel knob to formant frequencies
			float frac = vowel * (N_table-1);	
			int ind_low = (int)(frac);
			int ind_high = (int)ceil(frac);
			frac = frac - ind_low;
			fc[0] = frac*(table_F1[ind_high]-table_F1[ind_low]) + table_F1[ind_low];
			fc[1] = frac*(table_F2[ind_high]-table_F2[ind_low]) + table_F2[ind_low];
			fc[2] = frac*(table_F3[ind_high]-table_F3[ind_low]) + table_F3[ind_low];
			gain[0] = frac*(table_gain_F1[ind_high]-table_gain_F1[ind_low]) + table_gain_F1[ind_low];
			gain[1] = frac*(table_gain_F2[ind_high]-table_gain_F2[ind_low]) + table_gain_F2[ind_low];
			gain[2] = frac*(table_gain_F3[ind_high]-table_gain_F3[ind_low]) + table_gain_F3[ind_low];

			for (int i=0; i<3; i++) { //only do two formants (two bandpass filters
				fc[i] = fc[i] / 44100.f;  //normalize by the sample rate
				f[i] = sin(M_PI * fc[i]);
			}

			//convert q into the format that the algorithm needs
			q = 1 - q;
			
			//convert overall gain into logarithmic
			overall_gain = overall_gain * 2.0;  //make the center of the dial be zero gain.  max will be G=2 => 6dB
			overall_gain = overall_gain * overall_gain;  //max gain will be 4 => 12 dB
		}
		
		//choose which formant model to use
		int chooseModel(float new_model_f) {
			
			int new_model = max(1,model); //model is negative when starting up 
			if (model < 1) {
				if (new_model_f > 0.5) new_model = 2;
			} else if (model == 1) {
				if (new_model_f > 0.52) new_model = 2;
			} else {
				if (new_model_f < 0.48) new_model = 1;
			}
			
			//decide if we reset the model or not
			if (model != new_model) {
					 
				model = new_model;
				
				switch (model) {
					case 1:
						N_bandpass = N_bandpass_1;
						N_table = N_table_1;
						table_F1 = table_F1_1; 
						table_F2 = table_F2_1; 
						table_F3 = table_F3_1;
						table_gain_F1 = table_gain_F1_1; 
						table_gain_F2 = table_gain_F2_1; 
						table_gain_F3 = table_gain_F3_1;
						break;
						
					case 2:
						N_bandpass = N_bandpass_2;
						N_table = N_table_2;
						table_F1 = table_F1_2; 
						table_F2 = table_F2_2; 
						table_F3 = table_F3_2;
						table_gain_F1 = table_gain_F1_2; 
						table_gain_F2 = table_gain_F2_2; 
						table_gain_F3 = table_gain_F3_2;
						break;
				}
			}
			return model;
		}
		
		float bandpass(float sample, int ind) {
		
			low[ind] = low[ind] + f[ind] * band[ind];
			float high = q * sample - low[ind] - q*band[ind];
			band[ind] = f[ind] * high + band[ind];
			return gain[ind]*band[ind];
		}	
		float processSample(float sample){
			float out_val = 0.0;
			for (int i=0; i<N_bandpass; i++) {  //only do the 2 bandpass filters
				out_val += bandpass(sample, i);
			}
			return max(-1.0f, min(1.0f, overall_gain * out_val));
		}
  
  private:
		float low[3], band[3];
		float f[3], q;
		float gain[3];
		float overall_gain;
		int model;

		  
		#define MAX_TABLE 16
		int N_bandpass, N_table;
		float *table_F1, *table_F2, *table_F3;
		float *table_gain_F1, *table_gain_F2, *table_gain_F3;

		//Canadian english vowels https://home.cc.umanitoba.ca/~krussll/phonetics/acoustic/formants.html
		const int N_bandpass_1 = 2; //how many bandpass filters to use
		const int N_table_1 = 11; //how long are the arrays below
		float table_F1_1[MAX_TABLE] = {280.,	370.,	405.,	600.,	860.,	830.,	560.,	430.,	400.,	330.,	680. };
		float table_F2_1[MAX_TABLE] = {2230.,	2090.,	2080.,	1930.,	1550.,	1170.,	820.,	980.,	1100.,	1260.,	1310. };
		float table_F3_1[MAX_TABLE] = {0., 	0.,		0.,		0.,		0.,		0.,		0.,		0.,		0.,		0.,		0.,};
		float table_gain_F1_1[MAX_TABLE] = {1.0,	1.0,	1.0,	1.0,	1.0,	1.0,	1.0,	1.0,	1.0,	1.0,	1.0};	  
		float table_gain_F2_1[MAX_TABLE] = {1.0,	1.0,	1.0,	1.0,	1.0,	1.0,	1.0,	1.0,	1.0,	1.0,	1.0};
		float table_gain_F3_1[MAX_TABLE] = {0.0,	0.0,	0.0,	0.0,	0.0,	0.0,	0.0,	0.0,	0.0,	0.0,	0.0};

		//wikipedia: https://en.wikipedia.org/wiki/Formant
		const int N_bandpass_2 = 2; //how many bandpass filters to use
		const int N_table_2 = 16-3-3; //how long are the arrays below
		float table_F1_2[MAX_TABLE] = {240., 390., 610., 850., 820., 750., 700., 500., 360., 250.};
		float table_F2_2[MAX_TABLE] = {2400., 2300., 1900., 1610., 1530., 940., 760., 700., 640., 595.};
		float table_F3_2[MAX_TABLE] = {0., 	0.,		0.,		0.,		0.,		0.,		0.,		0.,		0.,		0.,		0.,};
		float table_gain_F1_2[MAX_TABLE] = {1.0,	1.0,	1.0,	1.0,	1.0,	1.0,	1.0,	1.0,	1.0,	1.0,	1.0};	  
		float table_gain_F2_2[MAX_TABLE] = {1.0,	1.0,	1.0,	1.0,	1.0,	1.0,	1.0,	1.0,	1.0,	1.0,	1.0};
		float table_gain_F3_2[MAX_TABLE] = {0.0,	0.0,	0.0,	0.0,	0.0,	0.0,	0.0,	0.0,	0.0,	0.0,	0.0};		 
	
};

#endif /* __StateVariableFilterPatch_hpp__ */


/* Formant Tables

https://en.wikipedia.org/wiki/Formant

Vowel
(IPA)	Formant F1
(Hz)	Formant F2
(Hz)	Difference
F1 – F2
(Hz)
IPA	F1	F2		Diff
i	240	2400	2160
y	235	2100	1865 //don't use
e	390	2300	1910 
ø	370	1900	1530 //don't use
ɛ	610	1900	1290
œ	585	1710	1125  //don't use
a	850	1610	760
ɶ	820	1530	710
ɑ	750	940		190
ɒ	700	760		60
ʌ	600	1170	570   //don't use
ɔ	500	700		200
ɤ	460	1310	850  //don't use
o	360	640		280
ɯ	300	1390	1090  //don't use
u	250	595		345

https://home.cc.umanitoba.ca/~krussll/phonetics/acoustic/formants.html
Canadian english vowels
Vowel	[i]	[ɪ]	[e]	[ɛ]	[æ]	[ɑ]	[ɔ]	[o]	[ʊ]	[u]	[ʌ]
F1	280	370	405	600	860	830	560	430	400	330	680
F2	2230	2090	2080	1930	1550	1170	820	980	1100	1260	1310

*/
