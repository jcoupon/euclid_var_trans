#!/bin/bash 
#$ -N compute_pz 
#$ -V 
#$ -q test 
#$ -cwd 
#$ -o /home/isdc/coupon/data/tmp/SGE_logs 
#$ -e /home/isdc/coupon/data/tmp/SGE_logs 
#$ -l h_rt=48:00:00 
#$ -l h_vmem=12G 
#$ -l h=!isdc-cn16.astro.unige.ch


echo $HOSTNAME
date


nnpz nnpz_base.conf         --reference_catalog /home/isdc/coupon/data/euclid/varTrans/fluxes_COSMOS15_East_training_vis24.5_corrected.fits         --reference_catalog_filters "[('u_obs', 'u_obs_err'),('g_obs', 'g_obs_err'),('r_tilt8_corr_obs', 'r_tilt8_corr_obs_err'),('i_obs', 'i_obs_err'),('z_obs', 'z_obs_err'),('vis_obs', 'vis_obs_err'),('Y_obs', 'Y_obs_err'),('J_obs', 'J_obs_err'),('H_obs', 'H_obs_err')]"         --target_catalog /home/isdc/coupon/data/euclid/varTrans/fluxes_COSMOS15_West_test_vis24.5/$( expr $SGE_TASK_ID - 1).fits         --target_catalog_filters "[('u_obs', 'u_obs_err'),('g_obs', 'g_obs_err'),('r_tilt8_obs', 'r_tilt8_obs_err'),('i_obs', 'i_obs_err'),('z_obs', 'z_obs_err'),('vis_obs', 'vis_obs_err'),('Y_obs', 'Y_obs_err'),('J_obs', 'J_obs_err'),('H_obs', 'H_obs_err')]"         --output_file /home/isdc/coupon/data/euclid/varTrans/fluxes_COSMOS15_West_test_vis24.5/$( expr $SGE_TASK_ID - 1)_pz_r_tilt8_corr.fits
photoz_metrics.py scatter -i /home/isdc/coupon/data/euclid/varTrans/fluxes_COSMOS15_West_test_vis24.5/$( expr $SGE_TASK_ID - 1)_pz_r_tilt8_corr.fits -o /home/isdc/coupon/data/euclid/varTrans/fluxes_COSMOS15_West_test_vis24.5/$( expr $SGE_TASK_ID - 1)_pz_r_tilt8_corr_density.pdf         -density 
photoz_metrics.py PDF -i /home/isdc/coupon/data/euclid/varTrans/fluxes_COSMOS15_West_test_vis24.5/$( expr $SGE_TASK_ID - 1)_pz_r_tilt8_corr.fits -o /home/isdc/coupon/data/euclid/varTrans/fluxes_COSMOS15_West_test_vis24.5/$( expr $SGE_TASK_ID - 1)_pz_r_tilt8_corr_PDF_tomo.pdf         -stats_output /home/isdc/coupon/data/euclid/varTrans/fluxes_COSMOS15_West_test_vis24.5/$( expr $SGE_TASK_ID - 1)_pz_r_tilt8_corr_PDF_tomo.csv
photoz_metrics.py PDF -i /home/isdc/coupon/data/euclid/varTrans/fluxes_COSMOS15_West_test_vis24.5/$( expr $SGE_TASK_ID - 1)_pz_r_tilt8_corr.fits -o /home/isdc/coupon/data/euclid/varTrans/fluxes_COSMOS15_West_test_vis24.5/$( expr $SGE_TASK_ID - 1)_pz_r_tilt8_corr_PDF_z0.63.pdf         -z_bins 0.55,0.70 -stats_output /home/isdc/coupon/data/euclid/varTrans/fluxes_COSMOS15_West_test_vis24.5/$( expr $SGE_TASK_ID - 1)_pz_r_tilt8_corr_PDF_z0.63.csv
ls /home/isdc/coupon/data/euclid/varTrans/fluxes_COSMOS15_West_test_vis24.5/$( expr $SGE_TASK_ID - 1)_pz_r_tilt8_corr_PDF_tomo.csv; ls /home/isdc/coupon/data/euclid/varTrans/fluxes_COSMOS15_West_test_vis24.5/$( expr $SGE_TASK_ID - 1)_pz_r_tilt8_corr_PDF_z0.63.csv

