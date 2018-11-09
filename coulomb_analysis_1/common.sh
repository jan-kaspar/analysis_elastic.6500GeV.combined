#!/bin/bash

function SetStandardExtraConfig()
{
	# defaults
	extra_config="-t-tr1 0.17 -t-tr2 0.25"

	if [ $t_max == "0.07" ]
	then
		extra_config="-t-tr1 0.09 -t-tr2 0.25"
	fi

	if [ $t_max == "0.15" ]
	then
		extra_config="-t-tr1 0.17 -t-tr2 0.25"
	fi
}
