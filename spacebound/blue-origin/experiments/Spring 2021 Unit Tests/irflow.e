/*
 * irflow Test //check sensor name
 * 
 * Description: Implementation test of FEMTA's IR flow meter
 * to determine functionality. NOTE - NOT A CALIBRATION TEST.
 *
 * Author: Vishal Ravi, Gouri Bellad, Mark Hartigan
 * Date: March 2 2021
 */
 
define enter start;

Sensor irflow 100Hz
{

  [calibrate   | Flow, poly, raw, ml/min | 1 0]
  [conversions | Flow, raw, ml/min]
  
  // options are Flow, Temperature, or Bubble
  [print | white, Flow | 1 ];
  
  
}
