# For fenics_ice copyright information see ACKNOWLEDGEMENTS in the fenics_ice
# root directory

# This file is part of fenics_ice.
#
# fenics_ice is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, version 3 of the License.
#
# fenics_ice is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with tlm_adjoint.  If not, see <https://www.gnu.org/licenses/>.

from fenics_ice.backend import compute_gradient, project

import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"

import sys
from pathlib import Path
import argparse

from fenics_ice import model, solver, inout
from fenics_ice import mesh as fice_mesh
from fenics_ice.config import ConfigParser
import fenics_ice.fenics_util as fu

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import datetime
import pickle
from IPython import embed


def run_forward_perturb(config_file):

    # Read run config file
    params = ConfigParser(config_file)
    log = inout.setup_logging(params)
    inout.log_preamble("forward", params)
    phase_sens = params.obs_sens.phase_name
    phase_suffix_sens = params.obs_sens.phase_suffix
    diag_dir = Path(params.io.diagnostics_dir)/phase_sens/phase_suffix_sens
    outdir = params.io.output_dir

    with open(os.path.join(diag_dir, 'vel_sens.pkl'), 'rb') as f:
       dict = pickle.load(f)

    gradient_val = np.nan
    if params.obs_sens.perturb_component != 'none':

       pert_mag = params.obs_sens.perturb_mag
       pert_locx = params.obs_sens.perturb_locx
       pert_locy = params.obs_sens.perturb_locy
       pert_rad = params.obs_sens.perturb_radius
       sign = params.obs_sens.perturb_sgn

       xpts = dict['uv_obs_pts'][:,0]
       ypts = dict['uv_obs_pts'][:,1]
       dObsU = dict['dObsU'][-1]
       dObsV = dict['dObsV'][-1]

       pert_dist = np.sqrt((xpts-pert_locx)**2 + (ypts-pert_locy)**2)
       Idx_update = (pert_dist<pert_rad)

       if (sign == 'p'):
           pert_mag = pert_mag
       elif (sign == 'n'):
           pert_mag = -1 * pert_mag
       else:
           pert_mag = 0.

       if (params.obs_sens.perturb_component == 'u'):
           dot = dObsU*pert_mag*np.exp(-pert_dist**2/pert_rad**2)
       elif (params.obs_sens.perturb_component == 'v'):
           dot = dObsV*pert_mag*np.exp(-pert_dist**2/pert_rad**2)
       dot[~Idx_update] = 0.0
       gradient_val = np.sum(dot)

    # Load the static model data (geometry, smb, etc)
    input_data = inout.InputData(params)

    # Get model mesh
    mesh = fice_mesh.get_mesh(params)

    # Define the model
    mdl = model.model(mesh, input_data, params)

    mdl.alpha_from_inversion()
    mdl.beta_from_inversion()

    # Solve
    slvr = solver.ssa_solver(mdl, mixed_space=params.inversion.dual)
    slvr.save_ts_zero()

    cntrl = slvr.get_control()

    qoi_func = slvr.get_qoi_func()

    # TODO here - cntrl now returns a list - so compute_gradient returns a list of tuples

    # Run the forward model
    Q = slvr.timestep(adjoint_flag=1, qoi_func=qoi_func)[-1].value()



    print(str(Q) + ',' + str(gradient_val))


    return mdl



if __name__ == "__main__":
    assert len(sys.argv) == 2, "Expected a configuration file (*.toml)"
    run_forward_perturb(sys.argv[1])
