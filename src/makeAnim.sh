#!/bin/bash

CONVERT=/opt/local/bin/convert

${CONVERT} -delay 20 Pressure*png AnimPressure.gif
${CONVERT} -delay 20 Stream*png   AnimStream.gif
${CONVERT} -delay 20 Forces*png   AnimForces.gif
${CONVERT} -delay 20 Velocity*png AnimVelocity.gif
${CONVERT} -delay 20 ParticleVelocity*png AnimParticles.gif
