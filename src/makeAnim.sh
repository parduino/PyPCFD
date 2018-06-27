#!/bin/bash

convert -delay 20 Pressure*png AnimPressure.gif
convert -delay 20 Stream*png   AnimStream.gif
convert -delay 20 Forces*png   AnimForces.gif
convert -delay 20 Velocity*png AnimVelocity.gif
convert -delay 20 ParticleVelocity*png AnimParticles.gif
