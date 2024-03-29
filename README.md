# SRT_sched

Create a list of targets for observations for the Sardinia Radio Telescope. We use a catalogue of sources from the TESS TOIs which is also supplied here as a .csv file 

To create a list of targets 9 hours of observations starting at a given UTC do following

>> python bl_srt_scripter.py --time "2024-03-16 15:30:00"

This will create a list of sources to observe in a sequence with 3 ONs and 3 OFFs.  
