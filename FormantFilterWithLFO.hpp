
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
/* Adapted by Chip Audette to "FormatFilterWithLFO" 2021 */


////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef __FormatFilterWithLFO_hpp__
#define __FormatFilterWithLFO_hpp__



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

class FormantFilterWithLFO : public SampleBasedPatch {
	public:
		FormantFilterWithLFO(void) {
			registerParameter(PARAMETER_A, "Vowel1"); //will be 0.0 to 1.0
			registerParameter(PARAMETER_B, "Vowel2"); //will be 0.0 to 1.0
			registerParameter(PARAMETER_C, "Speed"); //will be 0.0 to 1.0
			registerParameter(PARAMETER_D, "Gain");

			//initialize states
			overall_gain = 1.0;
			model = -1;
			for (int i=0; i<3; i++) {
				low[i] = 0.0;
				band[i]=0.0;
				
				f[i]=1000.0;
				gain[i]=0.0;
			}
			
			//choose baseline formant model
			chooseModel(3);  //this code has four models to choose from?

		}
		
		void updateFilters(float vowel, float *_f, float *_gain) {
			float fc[3];
			
			// fc = cutoff freq in Hz 
			// fs = sampling frequency //(e.g. 44100Hz)
			// q = resonance/bandwidth [0 < q <= 1]  most res: q=1, less: q=0
			
			float frac = vowel * (N_table-1);	
			int ind_low = (int)(frac);
			int ind_high = (int)ceil(frac);
			frac = frac - ind_low;
			fc[0] = frac*(table_F1[ind_high]-table_F1[ind_low]) + table_F1[ind_low];
			fc[1] = frac*(table_F2[ind_high]-table_F2[ind_low]) + table_F2[ind_low];
			fc[2] = frac*(table_F3[ind_high]-table_F3[ind_low]) + table_F3[ind_low];
			
			_gain[0] = frac*(table_gain_F1[ind_high]-table_gain_F1[ind_low]) + table_gain_F1[ind_low];
			_gain[1] = frac*(table_gain_F2[ind_high]-table_gain_F2[ind_low]) + table_gain_F2[ind_low];
			_gain[2] = frac*(table_gain_F3[ind_high]-table_gain_F3[ind_low]) + table_gain_F3[ind_low];

			for (int i=0; i<3; i++) { //only do two formants (two bandpass filters
				fc[i] = fc[i] / (44100.f / 2.0f);  //normalize by the nyquist
				_f[i] = sinf(M_PI * fc[i]);
			}
		}
		
		void prepare(void){
			vowel1 = getParameterValue(PARAMETER_A); //a value of 1.0 means fc = sample rate
			vowel2 = getParameterValue(PARAMETER_B); //a value of 1.0 means fc = sample rate
			//q = getParameterValue(PARAMETER_C); //a value of 1.0 means fc = sample rate
			float speed_frac = getParameterValue(PARAMETER_C);
			overall_gain = getParameterValue(PARAMETER_D);
			
			//choose the formant model to use
			//chooseModel(3);  //this code has four models to choose from?

			//convert q into the format that the algorithm needs
			q = 0.75;
			q = 1 - q;
			
			//convert the speed into an lfo increment
			if (speed_frac < 0.025) {
				//turn off the lfo
				lfo_increment = 0.0;
				lfo_val = 0.0;
			} else {
				lfo_increment = lfo_speed_scale * (speed_frac*speed_frac);  //squaring the speed_frac gives better access to smaller values
			}
			
			//convert overall gain into logarithmic
			float desired_mid_point = 0.7f;  //without scaling, neutral volume appears to be about 75% of the knob
			if (overall_gain < 0.5f) {
				overall_gain = overall_gain / 0.5 * desired_mid_point;
			} else {
				overall_gain = ((overall_gain - 0.5) / 0.5f) * (1.0-desired_mid_point) + desired_mid_point;
			}
			overall_gain = overall_gain * 3.0;  //make the center of the dial be zero gain.  max will be G=2 => 6dB
			overall_gain = overall_gain * overall_gain;  //max gain will be 4 => 12 dB
		}
			
		float processSample(float sample){
			
			//update the lfo
			lfo_val += (lfo_sign*lfo_increment);
			if ((lfo_sign > 0.1)  && (lfo_val > (1.0f - max(0.01f,1.5f*lfo_increment)))) lfo_sign = -lfo_sign; //flip the direction
			if ((lfo_sign < -0.1) && (lfo_val < (0.0f + max(0.01f,1.5f*lfo_increment)))) lfo_sign = -lfo_sign; //flip the direction
			lfo_val = max(0.0,min(1.0,lfo_val)); //limit the value
			
			//update the filter parameters
			float frac = lfo_val;
			float vowel = frac*(vowel2-vowel1)+vowel1;
			updateFilters(vowel, f, gain);
			
			
			//apply the bandpass filters
			float out_val = 0.0; //initialize our output variable
			for (int i=0; i<N_bandpass; i++) {
				out_val += bandpass(sample, i); //compute each bandpass filter in parallel and sum
			}
			out_val *= overall_gain; //apply overall gain
			out_val = max(-1.0f, min(1.0f, out_val));  //saturate whenever the amplitude is too large
			return out_val;
		}
		
		float bandpass(float sample, int ind) {
			low[ind] = low[ind] + f[ind] * band[ind];
			float high = q * sample- low[ind]- q*band[ind];
			band[ind] = f[ind] * high + band[ind];
			return (gain[ind]*band[ind]);
		}
		
		//choose which formant model to use
		int chooseModel(int new_model) {
						 
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

				case 3:
					N_bandpass = N_bandpass_3;
					N_table = N_table_3;
					table_F1 = table_F1_3; 
					table_F2 = table_F2_3; 
					table_F3 = table_F3_3;
					table_gain_F1 = table_gain_F1_3; 
					table_gain_F2 = table_gain_F2_3; 
					table_gain_F3 = table_gain_F3_3;
					break;

				case 4:
					N_bandpass = N_bandpass_4;
					N_table = N_table_4;
					table_F1 = table_F1_4; 
					table_F2 = table_F2_4; 
					table_F3 = table_F3_4;
					table_gain_F1 = table_gain_F1_4; 
					table_gain_F2 = table_gain_F2_4; 
					table_gain_F3 = table_gain_F3_4;
					break;
			}
			
			return model;
		}
		

  
  private:
		float low[3], band[3];
		float f[3], gain[3];
		float q;
		float overall_gain;
		int model;
		const float lfo_speed_scale = (1.0f/44100.0f)*2.0*10.0;  //fastest
		float lfo_increment = (1.0f/44100.0f)*0.5f;
		float lfo_val = 0.0f;
		float lfo_sign = 1.0; //switches between +1 and -1
		float vowel1, vowel2;

		  
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

		//https://engineering.purdue.edu/~ee649/notes/figures/formant_chart.gif
		//no gain adjustment
		const int N_bandpass_3 = 3; //how many bandpass filters to use
		const int N_table_3 = 10; //how long are the arrays below
		float table_F1_3[MAX_TABLE] = {270.0,	390.0,	530.0,	660.0,	660.0,	670.0,	440.0,	300.0,	640.0,	490.0};
		float table_F2_3[MAX_TABLE] = {2290.0,	1990.0,	1840.0,	1720.0,	1090.0,	840.0,	1020.0,	870.0,	1190.0,	1350.0};
		float table_F3_3[MAX_TABLE] = {3010.0,	2250.0, 2480.0,	2410.0,	2440.0,	2410.0,	2240.0,	2240.0,	2390.0,	1590.0};
		float table_gain_F1_3[MAX_TABLE] = {1.0,	1.0,	1.0,	1.0,	1.0,	1.0,	1.0,	1.0,	1.0,	1.0};	  
		float table_gain_F2_3[MAX_TABLE] = {1.0,	1.0,	1.0,	1.0,	1.0,	1.0,	1.0,	1.0,	1.0,	1.0};
		float table_gain_F3_3[MAX_TABLE] = {0.5,	0.5,	0.5,	0.5,	0.5,	0.5,	0.5,	0.5,	0.5,	0.5};		 

		//https://engineering.purdue.edu/~ee649/notes/figures/formant_chart.gif
		//with gain adjustment
		const int N_bandpass_4 = 3; //how many bandpass filters to use
		const int N_table_4 = 10; //how long are the arrays below
		float table_F1_4[MAX_TABLE] = {270.0,	390.0,	530.0,	660.0,	660.0,	670.0,	440.0,	300.0,	640.0,	490.0};
		float table_F2_4[MAX_TABLE] = {2290.0,	1990.0,	1840.0,	1720.0,	1090.0,	840.0,	1020.0,	870.0,	1190.0,	1350.0};
		float table_F3_4[MAX_TABLE] = {3010.0,	2250.0, 2480.0,	2410.0,	2440.0,	2410.0,	2240.0,	2240.0,	2390.0,	1590.0};
		float table_gain_F1_4[MAX_TABLE] = {1.0, 	1.0,	1.0,	1.0,	1.0,	1.0,	1.0,	1.0,	1.0,	1.0,	1.0};	  
		float table_gain_F2_4[MAX_TABLE] = {0.100,	0.100,	0.158,	0.282,	0.631,	0.447,	0.282,	0.158,	1.0,	0.316,	0.316};
		float table_gain_F3_4[MAX_TABLE] = {0.079,	0.063,	0.079,	0.100,	0.045,	0.020,	0.022,	0.010,	0.050,	0.178,	0.178};		 

	
};

#endif /* __StateVariableFilterPatch_hpp__ */



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
