/*
 * Sensirion NTCASCWE3 Test
 * 
 * Description: Implementation test of FEMTA's thermistor 
 * to determine functionality. NOTE - NOT A CALIBRATION TEST.
 *
 * Author: Vishal Ravi, Gouri Bellad, Mark Hartigan
 * Date: February 16 2021
 */
 
define enter start;

NTCASCWE3 { // need frequency from specific part number

[calibrate   | Temp, poly, raw, V | 1 0 ]
[conversions | Temp, raw, V, C    |     ] 


// options are Flow, Temperature, or Bubble
  [print | white, Temp | 1 ];
}

