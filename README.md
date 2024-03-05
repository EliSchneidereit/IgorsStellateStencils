Welcome to Igor's tool :)
=========================

Please compile the source file, using make.

Modify the existing input file according to your needs.
It is structured as follows, ordered by row:
Resolution in mm
Height of the pads in mm
Maximum circumfrerene lenth of a star jag in mm
Minimum free radius at the top surface in mm
Maximum overlapping radius at the top surface in mm
How long to recover for compressed material on the top surface (0.0, 1.0]

The following file would create 2 stars with 18 and 45mm diameter,
1mm resolution, 3.5mm thickness, 6.0mm star jags, 3mm guaranteed free
diameter on the top surface, and a desired overlap of 5mm. Hence, the
actual free diameter on the top surface is way greater than 3mm on both
pads.
The material compression on the cylinder sides is recovered on the first
2.5mm of the 5mm overlap:

 1.0
 3.5
 6.0
 1.5
 5.0
 0.5

 9
 22.5

Run the program using
 ./igor.exe igor.in
where mysettings.in is the filename of the input file you just created.
You should be able to open the SVG files created by igor.exe using
inkscape or a similar program.

Good luck!
