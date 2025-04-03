group_home="__________"
egfr_microsim_path=$group_home/egfr_microsim
cd $egfr_microsim_path

exp_name="full_100k"
R=100000

python egfr_microsim/analysis/plot_cover.py -e $exp_name
python egfr_microsim/analysis/plot_life_exp.py -e $exp_name

python egfr_microsim/analysis/plot_diag_time.py -e $exp_name --recalculate
python egfr_microsim/analysis/plot_param_prior.py -e $exp_name -r $R