# SixNodeThermalAnalysis
Preliminary Six Node thermal analysis for CubeSats on orbit

To develop a baseline understanding of modeling the thermally dynamic environment of low earth orbit (LEO), a six-node analysis can be performed as introduced by Versteeg et al. [1]. This models each face of the satellite as one single node in a circular 600 km-1000 km orbit around Earth.  Internal temperature gradients are ignored. During the satellite’s orbit not all faces are simultaneously exposed to the sun. This is dependent on the orientation of the satellite and therefore the viewing factor of each face.

We first define the date and time of hypothetical launch, JD.
“The Julian day is the continuous count of days since the beginning of the Julian period, and is used primarily by astronomers, and in software for easily calculating elapsed days between two events.”
Then we construct all parameters having to do with the orbit (inclination, velocity, period, etc.)
We choose as base material that of Aluminium 6061-T6 and its properties (Cp, absorptivity, emissivity etc).
We initialize temperatures for each face of the CubeSat where we are going to solve the equations.
Define some counters to introduce periodicity for viewing factors and eclipse fraction calculations.
Create some files to write in.
We convert the JD parameter as the number of days since Greenwich noon, Terrestrial time, on 1 January 2000 (common time reference).
Using this new parameter, n, we can now define the mean longitude of the Sun and the mean anomaly of the Sun in the Ecliptic Coordinate System. We can then compute the ecliptic longitude of the Sun and the ecliptic latitude which is almost equal to zero since it never exceeds 0.00033o.
Then, the distance of the Sun from the Earth in astronomical units is defined as RR, and the obliquity of the ecliptic can be approximated.
Finally, we need to define the sun position in the Equatorial Coordinate System which is the one most frequently used and therefore we compute the well-known Sun’s ascension and Sun’s declination.
In the next step we start by adding a while loop to designate the time passing. Since we need to see how the eclipse fractions change in the day, I use a timestep of 1 sec, otherwise we would lose changes.
I compute the beta angle having first calculated the RAAN which changes due to the Earth’s oblateness.
By comparing the calculated beta angle with the critical beta angle which is a function only of the orbit’s height we can compute the eclipse fraction.
Both the albedo factor aa and the IR heating rate are affected by the beta angle. Providing 3.3-σ accuracy values, those are as seen. Important to note is that the IR flux has already been corrected for Earth’s emissivity (as done by Versteeg et al. [1]).
Then, basically, we need to know which faces are illuminated and when, by specifying the viewing factors.
This if statement contributes to the periodic calculation of the viewing factors because it basically changes the timeframe.
Finally, we solve one equation for each face using a backward Euler scheme. This means that we are using the calculated temperature of the previous timestep to specify the present one.
After that we only have to save the data in the created files. 

References
[1] “Preliminary Thermal Analysis of Small Satellites”, Casper Versteeg, David L. Cotton, PhD, Small Satellite Research Laboratory, The University of Georgia, Athens, Georgia, 30602

