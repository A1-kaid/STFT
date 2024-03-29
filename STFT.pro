;============================================================================
;=                                                                          =
;=                      Short-time Fourier transform                        =
;=                                                                          =
;============================================================================
;Short-time Fourier transform Function
;
;Syntax:
;Result = stft(input [,window overlap] [,window lengh] [,input timeline]
;                [,sampling] [,window_function = scalar] [,timeline = variable]
;                    [,frequency = variable] [,/two_way] [,/inverse] [,/remove_dc])
;
;Return Value:
;  STFT return a two-dimensional vector which is Short-time Fourier
;transform result of input array or return a inverse STFT of input array.
;
;Arguments:
;Input:
;  An one-dimensional vector of original time domaim signal.
;Input timeline:
;  The time axis corresponding to the original signal.
;  If not spectified or 0, the converted timeline will not be output.
;Window length:
;  The length of window function, or the sampling length of STFT.
;  If not spectified, it will be preset to 512.
;Window overlap:
;  The length of overlap part between two window function.
;  If not spectified, it will be preset to 256.
;Sampling:
;  Sampling period of input signal.
;  If not spectified, it will be automatically determination.
;
;Keywords:
;Window_function:
;  0 or not spectified:Hanning
;  1:Rectangular
;  2:Hamming
;  3:Nuttall
;Timeline:
;  Return timeline of the spectrogram, in julday.
;Frequency:
;  Return frequency of the spectrogram.
;Two_way:
;  Set this keyword to do STFT both fellow and reverse the timeline
;of original data, the result will be the average of the two.
;Inverse:
;  Set this keyword to do inverse STFT, meanwhile, other keywords and
;arguments except Window overlap are no effect, and default value of
;overlap will be 0.
;Remove_dc:
;  Remove the DC component of the input waveform.
;============================================================================

function stft, signal0, overlap, win, time0, sampling, $
  window_function = w_number,  timeline = time1, frequency = freq, $
  two_way = two_way, inverse = inverse, remove_dc = remove_dc

  on_error, 2
  
  if ~keyword_set(inverse) then begin
    
    if ~keyword_set(win) then win = 512
    if ~keyword_set(overlap) then overlap = 256
    if ~keyword_set(remove_dc) then remove_dc = 0
    
    n = floor((n_elements(signal0) - overlap) / (win - overlap))
    window_f = findgen(win) + 1
    signal   = signal0[ 0: (win-long(overlap))*n+overlap-1]
    spectral = complexarr(n, win, /nozero)
    if keyword_set(two_way) then spectral2 = complexarr(n, win, /nozero)


    if keyword_set(w_number) then begin
      case w_number of
        ;hanning
        1:window_f = hanning(win)
        ;hamming
        2:window_f = 0.53836- 0.46164*cos(2*!pi*temporary(window_f) / (win-1))
        ;nuttall
        3:window_f = 0.355768- 0.487396*cos(2*!pi*window_f / (win-1)) $
          + 0.144232*cos(4*!pi*window_f / (win-1))- 0.012604*cos(6*!pi*window_f / (win-1))
      endcase
    endif else window_f = replicate(1, win);rectangular
    
    
    for i = 1, n, 1 do spectral[i-1, *] = $
      fft( (signal[ (win-long(overlap))*(i-1): (win-long(overlap))*i+overlap-1] - $
        mean(signal[ (win-long(overlap))*(i-1): (win-long(overlap))*i+overlap-1])*remove_dc)*window_f)

    if keyword_set(two_way) then begin
      signal = reverse( temporary( signal))
      for i = 1, n, 1 do spectral2[n-i, *] = $
        fft( (signal[ (win-long(overlap))*(i-1): (win-long(overlap))*i+overlap-1] - $
          mean(signal[ (win-long(overlap))*(i-1): (win-long(overlap))*i+overlap-1])*remove_dc)*window_f)
      spectral = ( temporary( spectral) + spectral2) / 2
    endif


    if keyword_set(time0) then begin

      if n_elements(signal0) eq n_elements(time0) then begin
        time1 = interpol(time0[ 0: (win-long(overlap))*n+overlap-1], n)
      endif else print, 'STFT error: Length of timeline data must be equal to the input signal.'

      if ~keyword_set(sampling) then begin
        caldat, (time0[-1]-time0[0])/(n_elements(time0)-1), 0, 0, 0, h, m, s
        sampling = (h-12)*3600 + m*60 + s
      endif

      freq = findgen(win/2)/(win*sampling)

    endif
    
    return, spectral
    
  endif else begin
    
    if ~keyword_set(overlap) then overlap = 0
    
    s_size = size( signal0, /dimensions)
    time_domain = fltarr( n_elements(signal0)/2 + overlap)
    
    for i = 0, s_size[0]-1, 1 do $
      time_domain[ i*s_size[1]/2: i*s_size[1]/2+s_size[1]-1] = fft(signal0[i,*], /inverse)
    
    return, time_domain
    
  endelse
  
end
