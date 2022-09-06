/*
** ds32 Clock Test
** 
** description: 
**      Read from the ds32 clock at a smaller timestep.
**
*/

define enter start


Sensor ds32 10Hz {

     [calibrate  | Time, poly, raw, s | 0.0009765]
     [conversion | Time, raw, s       | ]

}

