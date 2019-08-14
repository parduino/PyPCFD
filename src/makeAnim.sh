#!/bin/bash

if [ -x /opt/local/bin/convert ]
then
    CONVERT=/opt/local/bin/convert
elif [ -x /usr/local/bin/convert ]
then
    CONVERT=/usr/local/bin/convert
else
    exit 2
fi

if [ -x ${CONVERT} ]
then
	${CONVERT} -delay 20 Pressure*png AnimPressure.gif
	${CONVERT} -delay 20 Stream*png   AnimStream.gif
	${CONVERT} -delay 20 Forces*png   AnimForces.gif
	${CONVERT} -delay 20 Velocity*png AnimVelocity.gif
	${CONVERT} -delay 20 ParticleVelocity*png AnimParticles.gif
else
	echo "Install the ImageMagick package and adapt this script to locate convert"
	exit 1
fi

