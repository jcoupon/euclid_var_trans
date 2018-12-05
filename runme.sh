#! /bin/bash
# bash script to study the impact of variable transmisson
# on photo-z's

# configuration (machine dependent)
source scripts/config.sh

# set -o xtrace

#
# options
#

main() {

  # create filters and compute fluxes
  # -> VariableTransmissions.ipynb

  # to run in scripts
  # > rm -rf VariableTransmissions_files ; rm -f VariableTransmissions.{tex,pdf,out,log,aux} ; jupyter  nbconvert --to=latex VariableTransmissions.ipynb --template printviewlatex.tplx ; pdflatex VariableTransmissions.tex
  # or
  # > rm -rf VariableTransmissions_files ; rm -f VariableTransmissions.{tex,pdf,out,log,aux} ; jupyter  nbconvert --to=latex VariableTransmissions.ipynb --template template.tplx ; pdflatex VariableTransmissions.tex

  # test_nnpz
  test_plots

  return
}

# 130000 gals in 1.40 deg2
# functions
#

function test_plots {

  INPUT=$HOME/data/euclid/varTrans/fluxes_COSMOS15_West_test_vis24.5_pz_r.fits
  SELECT="(vis_obs_mag < 24.5) & (vis_obs_mag_err > 0.0) & (np.isfinite(TrueRedshiftPDZ_50)) & (0.2 < TrueRedshiftPDZ_50) & (TrueRedshiftPDZ_50 < 2.0)"

  # scripts/photoz_metrics.py scatter -i $INPUT -select "$SELECT" -title '$\mathrm{vis} < 24.5$' -zmin 0.0 -zmax 4.0 -stats_output stats.csv -density
  # photoz_metrics.py PDF -i $INPUT -select "$SELECT" -title '$\mathrm{vis} < 24.5$, ' -stats_output stats.csv  -z_bins 0.2,2.0  #  -z_bins 0.55,0.70
  # photoz_metrics.py bias -i $INPUT -select "$SELECT" -title '$\mathrm{vis} < 24.5$' -stats_output bias.csv #  -z_bins 0.55,0.70
  # photoz_metrics.py nz_bins -i $INPUT -select "$SELECT" -title '$\mathrm{vis} < 24.5$' # -z_bins 0.55,0.70

  photoz_metrics.py all -i $INPUT -select "$SELECT" -title '$\mathrm{vis} < 24.5$'  -zmin 0.0 -zmax 4.0 # -z_bins 0.55,0.70

  return
}

function test_nnpz
{

  TRAINING=$DATADIR/fluxes_COSMOS15_East_training_vis24.5.fits
  TEST=$DATADIR/fluxes_COSMOS15_West_test_vis24.5.fits

  # RESULTDIR=$DATADIR/photoz_training_no_noise
  # mkdir -p $RESULTDIR

  target_catalog_filters="[\
  ('u_obs', 'u_obs_err'),\
  ('g_obs', 'g_obs_err'),\
  ('r_obs', 'r_obs_err'),\
  ('i_obs', 'i_obs_err'),\
  ('z_obs', 'z_obs_err'),\
  ('vis_obs', 'vis_obs_err'),\
  ('Y_obs', 'Y_obs_err'),\
  ('J_obs', 'J_obs_err'),\
  ('H_obs', 'H_obs_err')\
  ]"

  reference_catalog_filters="[\
  ('u', None),\
  ('g', None),\
  ('r', None),\
  ('i', None),\
  ('z', None),\
  ('vis', None),\
  ('Y', None),\
  ('J', None),\
  ('H', None)\
  ]"

  if true; then

    nnpz config/nnpz_base.conf \
      --target_catalog $TEST \
      --target_catalog_filters "$target_catalog_filters" \
      --reference_catalog $TRAINING \
      --reference_catalog_filters "$reference_catalog_filters" \
      --output_file no_noise.fits --log_level 'DEBUG'

    #scripts/photoz_metrics.py scatter -i no_noise.fits -o no_noise.pdf \
    #  -select "(0.2<TrueRedshiftPDZ_50) & (TrueRedshiftPDZ_50<2.0)"

  fi



  #scripts/photoz_metrics.py PDF -i no_noise.fits -o no_noise_PDF.pdf \
  #    -select "(vis_obs_mag<24.5) & (0.2<TrueRedshiftPDZ_50) & (TrueRedshiftPDZ_50<2.0)"


  return
}

# ----------------------------------------------------------------------------------- #
# main
# ----------------------------------------------------------------------------------- #

main $@
