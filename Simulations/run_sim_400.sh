#!/bin/bash
cd "${0%/*}"
more=true
while [ "$more" = true ]
do
    more = $(Rscript ./Run_simulation_400.R)
    echo $more
done
