/*
** Main Flight Experiment Script
**      author: Evan Rittner
**        date: July 20, 2021
** description:
**     Script to be run during the Blue Origin suborbital flight. Controls all
**     aspects of the FEMTA and propellant management experiments, and logs 
**     data. Sensing transitions between states must be done resiliently.
*/


/*
** Description of flight states:
**
** preflight:
**     The period from when the script begins to liftoff. All solenoids will
**     begin closed, and remain closed during this state.
**
** launch:
**     The transition to this state will likely come from UDP data sent by 
**     Blue Origin over an Ethernet connection. This should be sanity checked 
**     by acceleration and atmospheric pressure. If possible, the time of launch
**     should be recorded, to ensure that other transitions occur at reasonable 
**     times (i.e. the ascent phase will last longer than ten seconds, so if 
**     we get indication of MECO ten seconds after liftoff, something is wrong).
**     The script will remain in this state only briefly; just long enough to 
**     open the collection chamber solenoid. It will then transition to the 
**     ascent state while waiting for MECO.
**
** ascent:
**     The script will wait in this state for the duration of the boost phase.
**     
** meco:
**     When the New Shepard main engine cuts off (MECO), we need to set up to
**     begin the experiment. The transition to this phase will likely come from
**     UDP data, sanity-checked by acceleration, atmospheric pressure, and
**     the elapsed time since launch. In this state, the collection chamber
**     solenoid will be opened, then after a short delay, the run/flow valves
**     will be opened. It is to be determined whether the FEMTA experiment will
**     be started here, or later in the flight. (The propellant management 
**     experiment is the priority for this flight, and there is a small risk
**     that the high voltage from the FEMTA experiment will interfere with it.)
**     Once all outputs are made, this script will transition to the 
**     microgravity state.
**
** microgravity:
**     The script will wait in this state until the microgravity phase of the 
**     flight ends. If the FEMTA experiment is started with a delay, this phase
**     must be modified.
**
** descent:
**     Once the New Shepard vehicle leaves microgravity, the experiment must
**     finish. The transition to this state will likely be signalled by UDP
**     data, again sanity-checked by data we collect. At this point, all
**     solenoids will close, and the FEMTA experiment will end. Once both
**     experiments are secured, the script will transition to the postflight
**     state.
**
** postflight:
**     The script will remain in this state indefinitely after the experiments
**     have concluded (even if the New Shepard booster has not yet landed). 
**     Data will be recorded in this state until the flight computer loses its
**     power supply from the booster.
*/

define enter preflight;
define leave launch;
define leave ascent;
define leave meco;
define leave microgravity;
define leave descent;
define leave postflight;
