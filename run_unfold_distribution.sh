#!/bin/bash

# The variable jobid is given by the environment

echo jobid is "$jobid"

# Suffix is an editable string used to put the results of different jobs in different directories
suffix="_27Jul2015"
# This array stores the variables for which the unfilding needs to be done
variable_array=('thadpt')
# Array of switches to turn on or off the use of the covariant matrix...
cov_matrix_array=('--cov_matrix=full')
# ...and corresponding strings to name the directories where the results will go
cov_matrix_name_array=('fullcov')
# Files in fit_file_array are "plots/${jobid}/ttxsec/${variable}/${variable}.${fit_file_array[$fit_file_index]}"...
# fit_file_array=('harvested.massdiscriminant.root' 'harvested.toy.massdiscriminant.root')
fit_file_array=('toy.harvested.root')
#  ...and corresponding strings to name the directories where the results will go
fit_file_name_array=('toy')
# Array of switches to turn on or off the use of the reco "truth" (i.e. the reco from the transport matrix) as the distribution to be unfolded
# Corresponding name is set automatically by the script
reco_truth_array=(' ')
# Array of transport (or "migration") matrices to use (this is also the directory name)...
scaled_matrix_array=('notscaled')
# ...and array of files in which each matrix is stored
scaled_matrix_file_array=("plots/${jobid}/ttxsec/migration_matrices.root")
regmodes=('Curvature')
# # Array of biases to apply to the true distribution
# bias_id_array=('_2015Jun26_skewed' '')
# #...and corresponding strings to name the directories where the results will go
# bias_id_name_array=('2015Jun26skewed' '')

for cov_matrix_index in "${!cov_matrix_array[@]}"
do
  for fit_file_index in "${!fit_file_array[@]}"
  do
    for variable in "${variable_array[@]}"
    do
      for reco_truth_index in "${!reco_truth_array[@]}"
      do
        for scaled_matrix_index in "${!scaled_matrix_array[@]}"
        do
          if [ "${reco_truth_array[$reco_truth_index]}" = "--use_reco_truth" ] && [ "${fit_file_name_array[$fit_file_index]}" = "toy" ]
          then
            continue
          fi
          if [ "${reco_truth_array[$reco_truth_index]}" = "--use_reco_truth" ]
          then
            fit_file_name="truth"
          else
            fit_file_name=${fit_file_name_array[$fit_file_index]}
          fi
					for regmode in "${regmodes[@]}"
					do
							target_dir=plots/${jobid}/ttxsec/${variable}_${cov_matrix_name_array[${cov_matrix_index}]}_reco${fit_file_name}_${scaled_matrix_array[$scaled_matrix_index]}_${regmode}${suffix}
							echo mkdir $target_dir
							mkdir plots/${jobid}/ttxsec/${variable}_${cov_matrix_name_array[${cov_matrix_index}]}_reco${fit_file_name}_${scaled_matrix_array[$scaled_matrix_index]}_${regmode}${suffix}
							cmd="unfold_distribution.py ${variable} plots/$jobid/ttxsec/${variable}/${variable}.${fit_file_array[$fit_file_index]} ${scaled_matrix_file_array[$scaled_matrix_index]} ${cov_matrix_array[$cov_matrix_index]} ${reco_truth_array[$reco_truth_index]} -d ${target_dir} --reg_mode=$regmode --runHandmade"
							echo python $cmd
							python $cmd
					done
        done
      done
    done
  done
done
