# STFT
Short-time Fourier transform for IDL

Syntax:  
Result = stft(array [,window overlap] [,window lengh] [,timeline]  
　　　　　　　　[,sampling] [,window_function = scalar] [,timeline = variable]  
　　　　　　　　　　　[,frequency = variable] [,/cross] [,/inverse] [,/center])  

Return Value:  
　STFT return a two-dimensional vector which is Short-time Fourier  
transform result of input array or return a inverse STFT of input array.  
  
Arguments:  
Array  
　An one-dimensional vector of original time domaim signal.  
Timeline:  
　The time axis corresponding to the original signal.  
　If not spectified or 0, the converted timeline will not be output.  
Window length:  
　The length of window function, or the sampling length of STFT.  
　If not spectified, it will be preset to 512.  
Window overlap:  
　The length of overlap part between two window function.  
　If not spectified, it will be preset to 256.  
Sampling:  
　Sampling of input signal, in seconds.  
　If not spectified, it will be automatically determination.  
  
Keywords:  
Window_function:  
　0 or not spectified:Hanning  
　1:Rectangular  
　2:Hamming  
　3:Nuttall  
Timeline:  
　Return timeline of the spectrogram, in julday.  
Frequency:  
　Return frequency of the spectrogram. 
Cross:  
　Set this keyword to do STFT both fellow and reverse the timeline  
of original data, the result will be the average of the two.  
Inverse:  
　Set this keyword to do inverse STFT, meanwhile, other keywords and  
arguments except Window overlap are no effect, and default value of  
overlap will be 0.  
Center:  
　Set this keyword to center original signal befor FFT, can be  
regarded as a special high-pass filter.