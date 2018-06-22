#!/bin/bash

octave --no-init-file \
	-p ./mlib/ \
	--persist \
	--eval "pkg load geometry"
