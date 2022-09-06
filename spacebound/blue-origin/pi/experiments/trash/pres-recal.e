/*
**  Pressure Recalibration Test
**  
**  Just recording raw sensor values from the pressure sensors
**  and the onboard HSCM. This is run twice; first, with two 
**  sensors to be tested on A1 and A2; the chamber sensor on A0,
**  and the HSCM.
**
**  Author: Evan Rittner
**  Date: 9/16/21
*/

define enter start;

Sensor hscm 10Hz {
	[print       | brown, Pressure      | 3       ];
}

Sensor ad15_vdd 10Hz {
	[conversions | A0, raw, V           |         ];
	[conversions | A1, raw, V           |         ];
	[conversions | A2, raw, V           |         ];
	[conversions | A3, raw, V           |         ];

	[calibrate   | A0, poly, raw, V     | 0.00041234314, -0.00198277794];
	[calibrate   | A1, poly, raw, V     | 0.00041234314, -0.00198277794];
	[calibrate   | A2, poly, raw, V     | 0.00041234314, -0.00198277794];
	[calibrate   | A3, poly, raw, V     | 0.00041234314, -0.00198277794];

	[print       | red, A0              | 3       ];
	[print       | blue, A1             | 3       ];
	[print       | green, A2            | 3       ];
	[print       | yellow, A3           | 3       ];
}

