
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
/* Adapted by Chip Audette to "VowelFilterWithTraj" 2021 */

////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef __VowelFilterWithTraj_hpp__
#define __VowelFilterWithTraj_hpp__



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

class VowelFilterWithTraj : public SampleBasedPatch {
	public:
		VowelFilterWithTraj(void) {
			registerParameter(PARAMETER_A, "Vowel"); //will be 0.0 to 1.0
			registerParameter(PARAMETER_B, "Trigger"); //will be 0.0 to 1.0
			registerParameter(PARAMETER_C, "Speed"); //will be 0.0 to 1.0
			registerParameter(PARAMETER_D, "Gain");

			//initialize states
			overall_gain = 1.0;
			for (int i=0; i<3; i++) {
				low[i] = 0.0;
				band[i]=0.0;
				
				f[i]= 0.1*((float)i);
				gain[i]=0.0;  //init to no gain (fully attenuated)
			}
			for (int i=0; i<n_ave; i++) {
				ave_buff[i] = 0.0f;
			}
			
			//choose baseline formant model
			chooseModel(2);  //this code has four models to choose from?

		}
		

		void updateFilters(float vowel_float, float time_float, float *_f, float *_gain) {
			float fc[3] = { 300., 1000., 3000.}; //dummy initial values
			
			//vowel_float is 0.0 to 1.0
			//time_float is 0.0 to 1.0
			
			// fc = cutoff freq in Hz 
			// fs = sampling frequency //(e.g. 44100Hz)
			// q = resonance/bandwidth [0 < q <= 1]  most res: q=1, less: q=0
			
			//prepare for interpolation for this vowel at this moment in time
			vowel_float = max(0.0f,min(1.0f, vowel_float));  //limit the value
			int vowel_int = (int)(((N_table-1)*vowel_float)+0.4999f); //get index of vowel that we want
			time_float = max(0.0f, min(1.0f, time_float));
			float interp_frac = time_float * (float)(N_time-1);
			int ind_low = (int)(interp_frac);
			int ind_high = (int)ceil(interp_frac);
			interp_frac = interp_frac - (float)ind_low;			
			

			//interpolate to get each formant's frequency for this vowel for this moment in time			
			fc[0] = interp_frac*(table_F1[vowel_int][ind_high]-table_F1[vowel_int][ind_low]) + table_F1[vowel_int][ind_low];
			fc[1] = interp_frac*(table_F2[vowel_int][ind_high]-table_F2[vowel_int][ind_low]) + table_F2[vowel_int][ind_low];
			fc[2] = interp_frac*(table_F3[vowel_int][ind_high]-table_F3[vowel_int][ind_low]) + table_F3[vowel_int][ind_low];
			
			//interpolate to get each formant's gain for this vowel for this moment in time
			/*
			_gain[0] = interp_frac*(table_gain_F1[vowel][ind_high]-table_gain_F1[vowel][ind_low]) + table_gain_F1[vowel][ind_low];
			_gain[1] = interp_frac*(table_gain_F2[vowel][ind_high]-table_gain_F2[vowel][ind_low]) + table_gain_F2[vowel][ind_low];
			_gain[2] = interp_frac*(table_gain_F3[vowel][ind_high]-table_gain_F3[vowel][ind_low]) + table_gain_F3[vowel][ind_low];
			*/
			_gain[0] = 1.0; //full gain, no attenuation
			_gain[1] = 1.0; //full gain, no attenuation
			_gain[2] = 1.0; //full gain, no attenuation

			

			for (int i=0; i<3; i++) { //only do two formants (two bandpass filters
				fc[i] = fc[i] / (44100.f / 2.0f);  //normalize by the nyquist
				_f[i] = sinf(M_PI * fc[i]);
			}
			
		}
		
		void prepare(void){
			vowel = getParameterValue(PARAMETER_A);   			//should be a value of 0.0 to 1.0
			trigger = getParameterValue(PARAMETER_B); 			//should be a value of 0.0 to 1.0
			float speed_frac = getParameterValue(PARAMETER_C);	//should be a value of 0.0 to 1.0
			overall_gain = getParameterValue(PARAMETER_D);		//should be a value of 0.0 to 1.0
			
			//choose the formant model to use
			//chooseModel(3);  //this code has four models to choose from?

			//set q and get it into the format that the algorithm needs
			q = 0.75; q = 1.0f - q;
			
			//set the trigger level (which is a power number) relative to full scale (FS = 1.0)
			float range_dB = 50.0f;  //here is the range that we would like to set for the knob
			float trigger_dBFS = trigger * range_dB - range_dB;  //should be negative and span -range_dB to 0.0
			trigger = powf(10.0f, trigger_dBFS / 10.0f); 
			
			//convert the speed into an lfo increment
			if (speed_frac > (1.0f-0.025f)) {
				//turn off the changing vowel...just hold it steady
				time_increment = 0.0f;
				time_val = 0.5f; //hold it at the midpoint of the vowel
			} else {
				
				//speed_frac is 0.0 to 1.0, linear.  Square it to get access to slower vowels.
				speed_frac = speed_frac*speed_frac;
				
				//compute time_increment, which is the amount added to the clock for evrey sample.
				//a smaller number makes the vowel evolve slowly.  A large number makes it evolve fast.
				const float shortest_time_sec = 0.2;
				const float longest_time_sec = 2.0;
				const float target_time_sec = speed_frac*(longest_time_sec - shortest_time_sec) + shortest_time_sec;
				const float sample_rate_Hz = 44100.;
				const float time_period_samples = target_time_sec * sample_rate_Hz;
				time_increment = 1.0 / time_period_samples;
			}
			
			//convert overall gain into logarithmic
			float desired_mid_point = 0.7f;  //without scaling, neutral volume appears to be about 75% of the knob
			if (overall_gain < 0.5f) {
				overall_gain = overall_gain / 0.5f * desired_mid_point;
			} else {
				overall_gain = ((overall_gain - 0.5f) / 0.5f) * (1.0f-desired_mid_point) + desired_mid_point;
			}
			overall_gain = overall_gain * 3.0f;  //make the center of the dial be zero gain.  max will be G=2 => 6dB
			overall_gain = overall_gain * overall_gain;  //max gain will be 4 => 12 dB
		}
			
		float processSample(float sample){
			//sample value should be -1.0 to +1.0
			
			//update running average estimate of signal amplitude
			float cur_pow = sample*sample; //square the signal to get power
			ave_ind = ave_ind + 1;  if (ave_ind >= n_ave) ave_ind = 0; //choose where to put the new data sample...this sample is the oldest in the buffer
			ave_ind = max(0,min(n_ave-1,ave_ind)); //make sure that we don't go out of bounds
			ave_sum = 0.9999f*ave_sum - ave_buff[ave_ind] + cur_pow;  //update the sum by removing the oldest value and adding the newest
			ave_buff[ave_ind] = cur_pow;  //save the new data sample over the oldest sample
			ave_pow = ave_sum / ((float) n_ave); //finish the calculation of the average
			ave_pow = max(0.000001f,min(1.0f,ave_pow)); //limit the value of the average power to reasonable values
			
			//based on the average signal power and retrigger threshold, decide whether to retrigger
			if (ave_pow >= trigger) {
				if (was_above_thresh == false) {
					//retrigger!
					time_val = 0.0f;
				} else {
					//do not reset.  we only want to reset on upward going transitions
				}
				was_above_thresh = true; //we are above the threshold, so note that
			} else {
				was_above_thresh = false; //we are below the threshold, so note that
			}

			//update the time
			time_val += time_increment;
			
			//update the filter parameters
			updateFilters(vowel, time_val, f, gain);
						
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
					N_time = N_time_1;
					table_F1 = traj_F1_1; 
					table_F2 = traj_F2_1; 
					table_F3 = traj_F3_1;
					//table_gain_F1 = table_gain_F1_1; 
					//table_gain_F2 = table_gain_F2_1; 
					//table_gain_F3 = table_gain_F3_1;
					break;
				case 2:
					N_bandpass = N_bandpass_2;
					N_table = N_table_2;
					N_time = N_time_2;
					table_F1 = traj_F1_2; 
					table_F2 = traj_F2_2; 
					table_F3 = traj_F3_2;
					//table_gain_F1 = table_gain_F1_1; 
					//table_gain_F2 = table_gain_F2_1; 
					//table_gain_F3 = table_gain_F3_1;
					break;
			}
			
			return model;
		}
		

  
  private:
		float low[3], band[3];
		float f[3], gain[3];
		float vowel;
		float q;
		float overall_gain;
		int model;
		//const float time_speed_scale = (1.0f/44100.0f)*20.0f;  //fastest is 20 per second
		float time_increment = (1.0f/44100.0f); //this will get overwritten in the methods
		float time_val = 0.0f; //time since the last trigger
		#define N_AVE (882*2)
		const int n_ave = N_AVE;
		float ave_buff[N_AVE];  //set for 20 msec, which should be a 50 Hz cutoff.  at 44.1kHz, that's about 882 samples
		int ave_ind = 0;
		float ave_sum = 0.0f;
		float ave_pow = 0.0f;
		float trigger = 0.01;
		float was_above_thresh = false; //state for the threshold detector
		
		  
		#define MAX_TABLE 16
		int N_bandpass, N_table, N_time;
		float (*table_F1)[3], (*table_F2)[3], (*table_F3)[3]; //pointer for a 2D array with the last dimension is size 3
		float (*table_gain_F1)[3], (*table_gain_F2)[3], (*table_gain_F3)[3];


		//https://web.nmsu.edu/~spsandov/papers/AverageFormantTrajectories.pdf
		//table 2, average vowel trajectories for monophthongs (regular vowels) for adult female subjects
		//20%, 50%, 80%
		int N_bandpass_1 = 2; //how many formants are we modeling here
		int N_table_1 = 12; //how many vowels (rows) are in the tables below
		int N_time_1 = 3;   //how many times (columns) are in the tables below
		float traj_F1_1[MAX_TABLE][3] = { {750., 794., 763.},
			{767., 815., 792.},
			{697., 723., 713.},
			{653., 683., 670.},
			{642., 604., 540.},
			{545., 558., 546.},
			{544., 560., 557.},
			{496., 484., 478.},
			{662., 665., 621.},
			{606., 608., 592.},
			{701., 724., 689.},
			{493., 492., 489.} };
		float traj_F2_1[MAX_TABLE][3] = { {1971., 1940., 1875.},
			{1382., 1355., 1413.},
			{1143., 1144., 1225.},
			{1927., 1900., 1828.},
			{2032., 2231., 2325.},
			{1528., 1562., 1609.},
			{2036., 2039., 1988.},
			{2234., 2388., 2356.},
			{1473., 1319., 1271.},
			{1512., 1501., 1479.},
			{1559., 1566., 1575.},
			{1526., 1352., 1271.}};
		float traj_F3_1[MAX_TABLE][3] = { {1.0, 1.0, 1.0},
			{1.0, 1.0, 1.0},
			{1.0, 1.0, 1.0},
			{1.0, 1.0, 1.0},
			{1.0, 1.0, 1.0},
			{1.0, 1.0, 1.0},
			{1.0, 1.0, 1.0},
			{1.0, 1.0, 1.0},
			{1.0, 1.0, 1.0},
			{1.0, 1.0, 1.0},
			{1.0, 1.0, 1.0},
			{1.0, 1.0, 1.0}};

		//https://web.nmsu.edu/~spsandov/papers/AverageFormantTrajectories.pdf
		//table 5, average vowel trajectories for diphthongs (regular vowels) for adult female subjects
		//20%, 50%, 80%
		int N_bandpass_2 = 2; //how many formants are we modeling here
		int N_table_2 = 9; //how many vowels (rows) are in the tables below
		int N_time_2 = 3;   //how many times (columns) are in the tables below
		float traj_F1_2[MAX_TABLE][3] = { {819., 819., 696.},
			{812., 839., 763.},
			{591., 596., 583.},
			{571., 571., 564.},
			{472., 469., 462.},
			{548., 551., 539.},
			{789., 727., 665.},
			{649., 659., 623.},
			{503., 483., 517.}}; //add "m" from Table 17
		float traj_F2_2[MAX_TABLE][3] = { {1469., 1667., 1955.},
			{1720., 1548., 1354.},
			{1581., 1562., 1599.},
			{1605., 1570., 1590.},
			{2054., 1957., 1856.},
			{1929., 1945., 1935.},
			{2030., 1991., 1942.},
			{1125., 1299., 1726,},
			{1490., 1502., 1567}};
		float traj_F3_2[MAX_TABLE][3] = { {1.0, 1.0, 1.0},
			{1.0, 1.0, 1.0},
			{1.0, 1.0, 1.0},
			{1.0, 1.0, 1.0},
			{1.0, 1.0, 1.0},
			{1.0, 1.0, 1.0},
			{1.0, 1.0, 1.0},
			{1.0, 1.0, 1.0},
			{1.0, 1.0, 1.0},
			{1.0, 1.0, 1.0},
			{1.0, 1.0, 1.0},
			{1.0, 1.0, 1.0},
			{1.0, 1.0, 1.0}};


	
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

