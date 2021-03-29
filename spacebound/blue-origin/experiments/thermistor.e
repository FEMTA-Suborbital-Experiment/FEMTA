/*
 * Vishay NTCASCWE3 Test
 * 
 * Description: Implementation test of FEMTA's thermistor 
 * to determine functionality. NOTE - NOT A CALIBRATION TEST.
 *
 * Author: Vishal Ravi, Gouri Bellad, Mark Hartigan
 * Date: February 25 2021
 */
 
define enter start;

Sensor ad15_vcc 1 Hz { //vcc is placeholder adc. thermistors may be connected to different adc

 where A0 is thermistor1 //thermistors may be connected to different ports A0-A3 for this adc or they might be spread out among other adcs
 where A1 is thermistor2
 where A2 is thermistor3
 
 [calibrate   | thermistor1, poly, raw, V | 1 0 ]  //placeholder calibration curve, get from electronics, simbox electronics, or simbox software
 [conversions | thermistor1, raw, V, C    |     ]  
 
 [calibrate   | thermistor2, poly, raw, V | 1 0 ]
 [conversions | thermistor2, raw, V, C    |     ]  

 [calibrate   | thermistor3, poly, raw, V | 1 0 ]
 [conversions | thermistor3, raw, V, C    |     ]  
 
 // options are Flow, Temperature, or Bubble
   [print | white, thermistor1 | 1 ];
   [print | white, thermistor2 | 1 ];
   [print | white, thermistor3 | 1 ];
}

Sensor ad15_sda 1 Hz {

 where A0 is thermistor4 //thermistors may be connected to different ports A0-A3 for this adc or they might be spread out among other adcs
 where A1 is thermistor5
 
 [calibrate   | thermistor4, poly, raw, V | 1 0 ]
 [conversions | thermistor4, raw, V, C    |     ]  

 [calibrate   | thermistor5, poly, raw, V | 1 0 ]
 [conversions | thermistor5, raw, V, C    |     ]  
 
  // options are Flow, Temperature, or Bubble
  [print | red, thermistor4 | 1 ];
  [print | red, thermistor5 | 1 ];
}

 ds32 1Hz { 
  [calibrate   | Time, poly, raw, s | 0.0009765625, 0.0]
  [conversions | Time, raw, s       |                  ]
}

