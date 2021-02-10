/*
 * Sensirion SLF3S-1300F Test
 * 
 * Description: Implementation test of FEMTA's flow meter
 * to determine functionality. NOTE - NOT A CALIBRATION TEST.
 *
 * Author: Mark Hartigan
 * Date: November 5, 2020
 */

define enter start; 



Sensor slf3_left 100Hz {
  
  [calibrate   | Flow, poly, raw, ml/min | 1 0]
  [conversions | Flow, raw, ml/min]
  
  // options are Flow, Temperature, or Bubble
  [print | white, Flow | 1 ];
}



