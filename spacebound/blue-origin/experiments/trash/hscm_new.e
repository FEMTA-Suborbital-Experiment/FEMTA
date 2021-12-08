/*
 * HSCM Calibration Testing 
 *
 * Description: HSCM Calibration testing
 *
 * Author: Joseph Kawiecki 
 * Date: September 22, 2021
 */

define enter start;

Sensor hscm 1Hz {

  [calibrate   | Temperature, poly, raw, C | 1, 0];
  [conversions | Temperature, raw, C       |     ];

  [print | salmon, Pressure  | 0];
  [print | mint, Temperature | 3];
}
