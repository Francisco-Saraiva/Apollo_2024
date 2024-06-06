// ~~~~~~~ INCLUDES ~~~~~~~
// ~~~~~~~~~~~~~~~~~


// ~~~~~~~ DEFINES ~~~~~~~
#define LED_RED PB14 
#define LED_GREEN PB0
#define LED_ORANGE PE1
#define BUFFER_SIZE 924
#define SAMPLING_FREQ 30000
#define PI 3.141592653589793238
#define FFT_WINDOW_SIZE 50
#define MIN_BETA 13
#define MAX_BETA 30
#define THRESHOLD 2.5
// ~~~~~~~~~~~~~~~~~~~~~~~~


// ~~~~~~~ VARIABLE CREATION ~~~~~~
String buffer = "";
int counter = 0;    // count the samples received from buffer 
int value;   // variable to hold the last received value by the buffer
double signalBuffer[BUFFER_SIZE];
double* signalLow;
double signalNotch[BUFFER_SIZE + 100] = {0};  // BUFFER_SIZE + FILTER_SIZE (101) - 1

unsigned long timeout = 30;  // in milliseconds

int num_bins = 1000;
int fft_sz = 0;
float* frequencies;  // will be resized in setup
float* freq_amps;

int beta_freqs[2] = {0,0};
int beta_size = 0;
float* beta_freq_amps;
float max_freq_value = 0;

double scaling_factor = 0.001;  // used to prevent overflow in the FFT


// ~~~~~~~ FUNCTIONS ~~~~~~~~

int FIR_Lowpass (double* filteredSignalBuffer, double inputSignalBuffer[]) {
  const int num_taps = 101;   // size of low-pass filter

  double lowpass_filterCoefficients[num_taps] = { 
    1.08192282e-34, -1.35855822e-06, -9.01414109e-06, -2.08586736e-05, -2.36424739e-05, 7.16701129e-20, 5.71022038e-05, 1.30806087e-04, 1.78120639e-04, 1.45558863e-04, -3.13502513e-19, -2.38584181e-04,
    -4.82272663e-04, -5.94830668e-04, -4.48486722e-04, 7.90413026e-19, 6.49450779e-04, 1.24993221e-03, 1.47686284e-03, 1.07197491e-03, -1.56521920e-18, -1.45512303e-03, -2.72369571e-03, -3.13760298e-03,
    -2.22520467e-03, 2.65066568e-18, 2.90012038e-03, 5.33287792e-03, 6.04489108e-03, 4.22506250e-03, -3.97433903e-18, -5.37390359e-03, -9.78620207e-03, -1.10049630e-02, -7.64562023e-03, 5.37283125e-18, 
    9.67157734e-03, 1.76334636e-02, 1.99166624e-02, 1.39510183e-02, -6.62066011e-18, -1.82149362e-02, -3.40946088e-02, -3.99334515e-02, -2.94109503e-02, 7.48617855e-18, 4.55752887e-02, 9.94462098e-02, 
    1.50382960e-01, 1.86788809e-01, 1.99993122e-01, 1.86788809e-01, 1.50382960e-01, 9.94462098e-02, 4.55752887e-02, 7.48617855e-18, -2.94109503e-02, -3.99334515e-02, -3.40946088e-02, -1.82149362e-02, 
    -6.62066011e-18, 1.39510183e-02, 1.99166624e-02, 1.76334636e-02, 9.67157734e-03, 5.37283125e-18, -7.64562023e-03, -1.10049630e-02, -9.78620207e-03, -5.37390359e-03, -3.97433903e-18, 4.22506250e-03,
     6.04489108e-03, 5.33287792e-03, 2.90012038e-03, 2.65066568e-18, -2.22520467e-03, -3.13760298e-03, -2.72369571e-03, -1.45512303e-03, -1.56521920e-18, 1.07197491e-03, 1.47686284e-03, 1.24993221e-03, 6.49450779e-04,
    7.90413026e-19, -4.48486722e-04, -5.94830668e-04, -4.82272663e-04, -2.38584181e-04, -3.13502513e-19, 1.45558863e-04, 1.78120639e-04, 1.30806087e-04, 5.71022038e-05, 7.16701129e-20, -2.36424739e-05, -2.08586736e-05, -9.01414109e-06, -1.35855822e-06, 1.08192282e-34
  };

  int new_size = Convolve(filteredSignalBuffer, inputSignalBuffer, lowpass_filterCoefficients, BUFFER_SIZE, num_taps);

  return new_size;
  
}

// ---------------------------------------------------------------------------

// 2nd Order IIR Notch filter
void Notch_Filter(double filteredSignalBuffer[], double inputSignalBuffer[], int input_size, float fc)
{
	// Notch filter parameters & constants
  // it is assumed a0 = 1
	double r = 0.80;
	double b0 = 1;
	double b1 = -2*cos(2*PI*fc/SAMPLING_FREQ);
	double b2 = 1;
	double a1 = -2*r*cos(2*PI*fc/SAMPLING_FREQ); 
	double a2 = r*r;
	
	
	double previous_values[2] = {0,0}; //stores previous filter outputs    //previous_values[0] -> 2 outputs ago; previous_values[1] -> 1 output ago
	double notch_value = 0; // Variable to hold most recent digital notch filter value
	
	
	//calculate filter output
	notch_value = (b0 * inputSignalBuffer[0]) + (b1 * 0) + (b2 * 0) - (a1 * previous_values[1]) - (a2 * previous_values[0]);

	//update filter output value
	previous_values[1] = notch_value;
	//update output array
	filteredSignalBuffer[0] = notch_value;
	
	
	//calculate filter output
	notch_value = (b0 * inputSignalBuffer[1]) + (b1 * inputSignalBuffer[0]) + (b2 * 0) - (a1 * previous_values[1]) - (a2 * previous_values[0]);
	//update filter output values
	previous_values[0] = previous_values[1];
	previous_values[1] = notch_value;
	//update output array
	filteredSignalBuffer[1] = notch_value;
	
	
	for(int n = 2; n < input_size; n++)
	{
		//calculate filter output
		notch_value = (b0 * inputSignalBuffer[n]) + (b1 * inputSignalBuffer[n-1]) + (b2 * inputSignalBuffer[n-2]) - (a1 * previous_values[1]) - (a2 * previous_values[0]);
		//update filter output values
		previous_values[0] = previous_values[1];
		previous_values[1] = notch_value;
		//update output array
		filteredSignalBuffer[n] = notch_value;
	}
}

// -----------------------------------------------------------------------------

int Convolve(double* convolutedSignal, double inputSignal[], double filter[], int signal_size, int filter_size){
  uint32_t i,j;  //iterators

  // initialize output arrray
  
  for (i=0; i<(signal_size + filter_size); i++){
    convolutedSignal[i] = 0;
  }

  for (i=0; i<signal_size; i++){
    for(j=0; j<filter_size; j++){
      convolutedSignal[i+j] = convolutedSignal[i+j] + (inputSignal[i] * filter[j]);
    }
  }

  return (signal_size + filter_size -1);

}

// ----------------------------------------------------

int FFT(float* frequencies, float* freq_amps, double inputSignalBuffer[], int input_size, int N, float Frequency){
  /*
  Inputs:
  1. inputSignalBuffer[]     : Data array, 
  2. N        : Number of samples (should be a power of 2!)
  3. Frequency: sampling frequency required as input (Hz)
  */

  unsigned int data[12]={1,2,4,8,16,32,64,128,256,512,1024,2048};
  int a,c1,f,o,x;
  a=N;  
                                  
  for(int i=0; i<12; i++) { 
    if(data[i]<=a){
      o=i;
    } 
  }

  int fft_size = data[o];

        
  int in_ps[fft_size]={};     //input for sequencing
  float out_r[fft_size]={};   //real part of transform
  float out_im[fft_size]={};  //imaginary part of transform       

  // kill the function if fft_size < input_size
  if (fft_size < input_size){
    return 1;
  }

  // Might not be padded if input_size = fft_size, but it will still be called that way
  double paddedInputSignal[fft_size] = {0};
  for (int k = 0; k < input_size; k++){
    paddedInputSignal[k] = inputSignalBuffer[k];
  }

  // bit reversal
  x=0; 
  for(int b=0; b<o; b++){
      c1=data[b];
      f=fft_size/(c1+c1);

      for(int j=0; j<c1; j++){   
        x=x+1;
        in_ps[x]=in_ps[j]+f;
      }
  }
  
  for(int i=0; i<fft_size; i++){  // update input array as per bit reverse order
    if(in_ps[i]<a){
      out_r[i]=paddedInputSignal[in_ps[i]];
    }
    if(in_ps[i]>a){
      out_r[i]=paddedInputSignal[in_ps[i]-a];
    }      
  }

  // Iterators and value-holders
  int i10,i11,n1;
  float e,c,s,tr,ti;

  // FFT calculations 
  for(int i=0;i<o;i++) {                                  

    i10=data[i];              
    i11=fft_size/data[i+1];    
    e=(2*PI)/data[i+1];
    e=0-e;
    n1=0;

    for(int j=0; j<i10; j++){
      c=cos(e*j);
      s=sin(e*j);    
      n1=j;
    
      for(int k=0; k<i11; k++){
        tr=c*out_r[i10+n1]-s*out_im[i10+n1];
        ti=s*out_r[i10+n1]+c*out_im[i10+n1];

        out_r[n1+i10]=out_r[n1]-tr;
        out_r[n1]=out_r[n1]+tr;

        out_im[n1+i10]=out_im[n1]-ti;
        out_im[n1]=out_im[n1]+ti;          

        n1=n1+i10+i10;
      }       
    }
  }

  

  //---> here onward out_r contains amplitude and out_im conntains frequency (Hz)
  for(int i=0; i<data[o-1]; i++){               // getting amplitude from complex number
    out_r[i]=sqrt(out_r[i]*out_r[i]+out_im[i]*out_im[i]); // to increase the speed delete sqrt
    out_im[i]=i*Frequency/N;
    /*
    Serial.print(out_im[i]); Serial.print("Hz");
    Serial.print("\t");                            // uncomment to print frequency bin    
    Serial.println(out_r[i]);
    */
         
  }

  Copy_vector(out_im, frequencies, (fft_size/2));
  Copy_vector(out_r, freq_amps, (fft_size/2));

  return 0;   

}

// ---------------------------------------------------
void scale_signal(double signal[], int signal_size, double scaling_factor){
  for (int i=0; i< signal_size; i++){
    signal[i] *= scaling_factor;
  }
}

// ----------------------------------------------------

int highest_power2(int number) {
  if (number < 1) {
    return 0; // Return 0 for invalid input
  }

  int pwr = 1;

  while (pwr < number) {
    pwr <<= 1; // Equivalent to pwr *= 2
  }

  return pwr; 
}

// --------------------------------------------------

// Copy_vector(original, to_copy, size)
void Copy_vector(float v1[], float* v2, int size){
  for (int i=0; i<size; i++){
    v2[i] = v1[i];
  }
}

// --------------------------------------------------
float Max(float* vector, int size){
  float max = -9999.99;
  for(int i=0; i< size; i++){
    if (vector[i] > max)
      max = vector[i];
  }
  return max;
}

// ---------------------------------------------------

int Max_idx(float* vector, int size){
  int max_idx = -1;
  float max_value = -9999.99;
  for(int i=0; i< size; i++){
    if (vector[i] > max_value){
      max_idx = i;
      max_value = vector[i];
    }
  }
  return max_idx;
}

// ---------------------------------------------------

void Subvector(float* subvector, int size_sub, float* original, int pos1){
  for (int i=0; i <= size_sub ;i++){
    subvector[i] = original[i+pos1];
  }
}

// ---------------------------------------------------

void cut_and_glue(double vector[], int vector_length, int sliding_length) {
    double mywindow[vector_length - sliding_length] = {0};

    for (int i = sliding_length; i < vector_length; i++) {
      mywindow[i - sliding_length] = vector[i];
    }

    for (int k = 0; k < vector_length - sliding_length; k++) {
      vector[k] = mywindow[k];
    }   
}

//-----------------------------------------------------
void find_indices(float* vector, int output[], int size, float min_value, float max_value){
  bool flag_min = true;
  
  for (int i=0; i<size; i++){
    if (vector[i] >= min_value && flag_min){
      output[0] = i;
      flag_min = false;
    }

    if (vector[i] > max_value){
      output[1] = i-1;
      break;
    }
  }
}

//----------------------------------------------------

void Blink_red(unsigned long time){
  digitalWrite(LED_RED, HIGH);
  delay(time);
  digitalWrite(LED_RED, LOW);
}

// -----------------------------------------------------

void Blink_green(unsigned long time){
  digitalWrite(LED_GREEN, HIGH);
  delay(time);
  digitalWrite(LED_GREEN, LOW);
}

// -----------------------------------------------------
void Blink_orange(unsigned long time){
  digitalWrite(LED_ORANGE, HIGH);
  delay(time);
  digitalWrite(LED_ORANGE, LOW);
}

// -------------------------------------------------------
void LightUpGreen(){
  digitalWrite(LED_GREEN, HIGH);
  digitalWrite(LED_RED, LOW);
}

// ---------------------------------------------------------

void LightUpRed(){
  digitalWrite(LED_RED, HIGH);
  digitalWrite(LED_GREEN, LOW);
}

// ----------------------------------------------------


// ~~~~~~ SETUP (code to be ran once) ~~~~~
void setup() {
  Serial.begin(115200);

  pinMode(LED_RED, OUTPUT);
  pinMode(LED_GREEN, OUTPUT);
  pinMode(LED_ORANGE, OUTPUT);

  num_bins = highest_power2(num_bins);
  fft_sz = num_bins/2;
  frequencies = new float[fft_sz];
  freq_amps = new float[fft_sz];
  Blink_orange(500);
}
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


// ~~~~ LOOP (code to be ran forever) ~~~~~
void loop() {
  if (Serial.available() > 0){

    buffer = Serial.readStringUntil('\n');
    value = buffer.toDouble();
    signalBuffer[counter] = value;
    counter++;

    if (counter == BUFFER_SIZE){

      // FILTERING
      int new_size = FIR_Lowpass(signalLow, signalBuffer);
      Notch_Filter(signalNotch, signalLow, new_size, 50);
      scale_signal(signalNotch, new_size, scaling_factor);

      int code = FFT(frequencies, freq_amps, signalNotch, new_size, num_bins, SAMPLING_FREQ);
      cut_and_glue(signalBuffer, BUFFER_SIZE, FFT_WINDOW_SIZE);  // For a sliding FFT

      // FEATURE EXTRACTION
      find_indices(frequencies, beta_freqs, fft_sz, MIN_BETA, MAX_BETA);
      beta_size = beta_freqs[1] - beta_freqs[0] + 1;
      beta_freq_amps = new float[beta_size];

      Subvector(beta_freq_amps, beta_size, freq_amps, beta_freqs[0]);
      max_freq_value = Max(beta_freq_amps, beta_size);

      int max_idx = Max_idx(freq_amps, fft_sz);
      float max_freq = frequencies[max_idx];

      // PRINTS AND RESETS
      counter = counter - FFT_WINDOW_SIZE;

      //Serial.println(max_freq_value);

      if (max_freq_value > THRESHOLD){
        LightUpGreen();
      }
      else {
        LightUpRed();
      }
    }
  }
}
// ~~~~~~~~~~~~~~~~~~
