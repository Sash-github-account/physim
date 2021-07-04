#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 27 14:20:15 2021

@author: sashwathnalinkanth
"""
import pyaudio
import wave

filename = 'file_example_WAV_1MG.wav'

# Set chunk size of 1024 samples per data frame
chunk = 1024  
optn = int(input("1.Play or 2.record ?:"))

if(optn == 1):
    # Open the sound file 
    wf = wave.open(filename, 'rb')
    
    # Create an interface to PortAudio
    p = pyaudio.PyAudio()
    
    # Open a .Stream object to write the WAV file to
    # 'output = True' indicates that the sound will be played rather than recorded
    stream = p.open(format = p.get_format_from_width(wf.getsampwidth()),
                    channels = wf.getnchannels(),
                    rate = wf.getframerate(),
                    output = True)
    
    # Read data in chunks
    data = wf.readframes(chunk)
    print(data)
    # Play the sound by writing the audio data to the stream
    while data != '':
        stream.write(data)
        data = wf.readframes(chunk)
    
    # Close and terminate the stream
    stream.close()
    p.terminate()
else :
    sample_format = pyaudio.paInt16  # 16 bits per sample
    channels = 2
    fs = 44100  # Record at 44100 samples per second
    seconds = 10
    filename = "output.wav"
    
    p = pyaudio.PyAudio()  # Create an interface to PortAudio
    
    print('Recording')
    
    stream = p.open(format=sample_format,
                    channels=channels,
                    rate=fs,
                    frames_per_buffer=chunk,
                    input=True)
    
    frames = []  # Initialize array to store frames
    
    # Store data in chunks for 3 seconds
    # for i in range(0, int(fs / chunk * seconds)):
    #     data = stream.read(chunk)
    #     frames.append(data)
    print("Keyboard interrupt to stop rec.")
    # while(True):
    #     try:
    #         data = stream.read(chunk)
    #         frames.append(data)
    #     except KeyboardInterrupt:
    #         # Stop and close the stream 
    #         stream.stop_stream()
    #         stream.close()
    #         # Terminate the PortAudio interface
    #         p.terminate()
    #             # Save the recorded data as a WAV file
    #         wf = wave.open(filename, 'wb')
    #         wf.setnchannels(channels)
    #         wf.setsampwidth(p.get_sample_size(sample_format))
    #         wf.setframerate(fs)
    #         wf.writeframes(b''.join(frames))
    #         wf.close()
    #         print('Finished recording')
    #         raise
    while(True):
        if keyboard.is_pressed('q' or 's'):  # if key 'q' is pressed 
            # Stop and close the stream 
            stream.stop_stream()
            stream.close()
            # Terminate the PortAudio interface
            p.terminate()
                # Save the recorded data as a WAV file
            wf = wave.open(filename, 'wb')
            wf.setnchannels(channels)
            wf.setsampwidth(p.get_sample_size(sample_format))
            wf.setframerate(fs)
            wf.writeframes(b''.join(frames))
            wf.close()
            print('Finished recording')       
            break  # finishing the loop            
