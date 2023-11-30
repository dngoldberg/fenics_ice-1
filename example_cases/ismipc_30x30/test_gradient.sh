#!/bin/bash
set -e

#Generate the input data (100x100 grid, and use this to generate 'obs_vel')
#python $FENICS_ICE_BASE_DIR/aux/gen_rect_mesh.py -o ./input/momsolve_mesh.xml -xmax 40000 -ymax 40000 -nx 100 -ny 100
#python $FENICS_ICE_BASE_DIR/aux/gen_rect_mesh.py -o ./input/ismip_mesh.xml -xmax 40000 -ymax 40000 -nx 30 -ny 30

#python $FENICS_ICE_BASE_DIR/aux/gen_ismipC_domain.py -o ./input/ismipc_input.h5 -L 40000 -nx 100 -ny 100 --reflect
#python $FENICS_ICE_BASE_DIR/runs/run_momsolve.py momsolve.toml
#python $FENICS_ICE_BASE_DIR/aux/Uobs_from_momsolve.py -i "U.h5" -o "ismipc_U_obs.h5" -L 40000 -d ./output_momsolve

#cp output_momsolve/ismipc_U_obs.h5 input/

#Run each phase of the model in turn

#perturb_locx: float = np.nan
#    perturb_locy: float = np.nan
#    perturb_mag: float = np.nan
#    perturb_component: str = 'none'

xposition=(10000 10000)
yposition=(9000 30000)
perturb=(1 .5 .25 .125)
RUN_DIR=$FENICS_ICE_BASE_DIR/runs/
rm -f out_grad

toml set --to-float --toml-path ismipc_30x30.toml obs_sens.perturb_radius 5000

for ((i=1; i<=2; i++)); do
    xpos=${xposition[$((i-1))]}
    ypos=${yposition[$((i-1))]}

    toml set --toml-path ismipc_30x30.toml obs_sens.perturb_component 'v'
    toml set --to-float --toml-path ismipc_30x30.toml obs_sens.perturb_locx $xpos
    toml set --to-float --toml-path ismipc_30x30.toml obs_sens.perturb_locy $ypos

    toml set --to-float --toml-path ismipc_30x30.toml obs_sens.perturb_mag 0
    python $RUN_DIR/run_inv.py ismipc_30x30.toml
    python $RUN_DIR/run_forward_perturb.py ismipc_30x30.toml > ret
    echo $xpos, $ypos, 0, $(tail -n 1 ret) >> out_grad

    for pert in "${perturb[@]}"; do


	toml set --to-float --toml-path ismipc_30x30.toml obs_sens.perturb_mag $pert

	toml set --toml-path ismipc_30x30.toml obs_sens.perturb_sgn 'p'
	python $RUN_DIR/run_inv.py ismipc_30x30.toml 
        python $RUN_DIR/run_forward_perturb.py ismipc_30x30.toml > ret
	echo $xpos, $ypos, $pert, $(tail -n 1 ret) >> out_grad

	toml set --toml-path ismipc_30x30.toml obs_sens.perturb_sgn 'n'
	python $RUN_DIR/run_inv.py ismipc_30x30.toml 
        python $RUN_DIR/run_forward_perturb.py ismipc_30x30.toml > ret
	echo $xpos, $ypos, -$pert, $(tail -n 1 ret) >> out_grad

    done	    

done

exit

RUN_DIR=$FENICS_ICE_BASE_DIR/runs/
python $RUN_DIR/run_inv.py ismipc_30x30.toml
python $RUN_DIR/run_forward_perturb.py ismipc_30x30.toml
#python $RUN_DIR/run_eigendec.py ismipc_30x30.toml
#python $RUN_DIR/run_errorprop.py ismipc_30x30.toml
#python $RUN_DIR/run_obs_sens_prop.py ismipc_30x30.toml

# for verification, need to:

# 1. perturb single values of velocity data
# 2. run inversion
# 3. run forward, WITHOUT adjoint
# 4. record the QoI value
