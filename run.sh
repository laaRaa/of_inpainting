#!/bin/bash
set -e

output_method=$1
method=$2
lambda_lb=$3
g_lb=$4
lambda_amle=$5
g_amle=$6
epsilon=$7
S=$8
nt=$9
rr=${10}
virtualenv=${11}
BIN=${12}
MAT=${13}


if [ ! -d $virtualenv ]; then
  echo "Virtualenv not found" > demo_failure.txt
  exit 0
fi

source $virtualenv/bin/activate

if [ -f mask_0.png ]; then
	cp mask_0.png mask_d.png
	python $BIN/merge_mask.py
	cp mask.png input_2.png
fi

if [ -f mask_1.png ]; then
	cp mask_1.png mask_d.png
	python $BIN/merge_mask.py
	cp mask.png input_2.png
fi

if [ -f mask_2.png ]; then
	cp mask_2.png mask_d.png
	python $BIN/merge_mask.py
	cp mask.png input_2.png
fi

if [ ! -f input_3.flo ]; then
	if [ $method == "1" ]; then
		tvl1flow input_0.png input_1.png input.flo
	else
		main input_0.png input_1.png input.flo
	fi
else
	cp input_3.flo input.flo
fi

if [ $output_method == "1" ]; then

	run_inpaint_OF_with_laplace_beltrami_flow_interpolation.sh $MAT input.flo input_0.png input_2.png $lambda_lb $g_lb flow_out_lb.flo
	run_flo_to_png.sh $MAT input.flo input.flo input_flow.png
	run_flo_to_png.sh $MAT input.flo flow_out_lb.flo flow_out_lb.png

elif [ $output_method == "2" ]; then

	amle_recsep $S $lambda_amle $epsilon $g_amle $nt $rr input.flo input_2.png flow_out_amle.flo input_0.png
	run_flo_to_png.sh $MAT input.flo input.flo input_flow.png
	run_flo_to_png.sh $MAT input.flo flow_out_amle.flo flow_out_amle.png

elif [ $output_method == "3" ]; then

	amle_recsep $S $lambda_amle $epsilon $g_amle $nt $rr input.flo input_2.png flow_out_amle.flo input_0.png
	run_inpaint_OF_with_laplace_beltrami_flow_interpolation.sh $MAT input.flo input_0.png input_2.png $lambda_lb $g_lb flow_out_lb.flo
	run_flo_to_png.sh $MAT input.flo input.flo input_flow.png
	run_flo_to_png.sh $MAT input.flo flow_out_lb.flo flow_out_lb.png
	run_flo_to_png.sh $MAT input.flo flow_out_amle.flo flow_out_amle.png

fi
