Welcome to Igor's tool :)
=========================

Please compile the source file, using make.

Modify the existing input file `igor.in` according to your needs.

Run the program using
 ./igor.exe
It will search for the configuration file on its own.
If you want to specify parameters without changing the configuration file use e.g.
  ./igor.exe --Rfirst=9.0
Run
  ./igor.exe --help
to get a help message telling you about all options and their default values.

You should be able to open the SVG files created by igor.exe using
inkscape or a similar program if you set `output_format=svg` else,
you can plot the results using gnuplot.
