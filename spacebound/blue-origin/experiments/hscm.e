/*
 * Honeywell HSCMNNN015PA2A3 Test
 *
 * Description: Test of board mounted pressure sensor to confirm functionality
 *
 * Author: Mark Hartigan
 * Date: June 10, 2021
 */

define enter start;

Sensor hscm 1Hz {
  [print | salmon, Pressure | 3];
}
