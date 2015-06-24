# noaa-apt

This is a simple decoder for the APT signals that are transmitted by the
NOAA satellites.  It is based upon a very simple decoder that I wrote
a long time ago as part of an aborted article I was writing for AMSAT.
That encoder was quite simple, but also took a lot of short cuts.  My goal
is to have a simple command line program that you can run on a recording
(any format, bit depth, or sample rate) and have it do something straight
forward and produce the best image that it knows how.

It's fairly common for computers these days to have gigabytes of memory,
so I'm not working at all to reduce memory consumption.  The program works
by loading the entire signal into memory at once, and processing it as a 
whole.  
