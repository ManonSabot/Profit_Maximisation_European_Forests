#!/usr/bin/env bash

################################################################################
#                                                                              #
# This file is part of the TractLSM project.                                   #
#                                                                              #
# Copyright (c) 2019 Manon E. B. Sabot                                         #
#                                                                              #
# Please refer to the terms of the MIT License, which you should have received #
# along with the TractLSM.                                                     #
#                                                                              #
# N.B.: if running in windows, applying the "sed -i 's/\r$//' filename"        #
#       command on the file might be necessary before it can be made an        #
#       executable.                                                            #
#                                                                              #
# __author__ = "Manon E. B. Sabot"                                             #
# __version__ = "1.0 (12.07.2018)"                                             #
# __email__ = "m.e.b.sabot@gmail.com"                                          #
#                                                                              #
################################################################################


#### USER'S CHOICES ####
p1="var_kmax_adjust_average"
p2="var_kmax_adjust_extreme"
p3="var_kmax_sample"
p4="calib_g1"
p5="calib_fw"


#######################################
# Main: calls all the plotting scripts
#       used in this paper
# Globals:
#   lots
# Arguments:
#   None
# Returns:
#   figures
#######################################

main(){

# plot the performance scores of the experiments, for the data analysis period
${this_dir}${sl}plotting_scripts${sl}compare_best_performance_scores.py

# plot the relationships of the modelled kmax to MAP, for all experiments
${this_dir}${sl}plotting_scripts${sl}kmax_vs_climate.py ${p1} ${p2} -c ${p3} > \
                            ${data_dir}${sl}output${sl}kmax_vs_climate_stats.txt

# plot the relationships of the fitted D-sensitivities to g1
${this_dir}${sl}plotting_scripts${sl}supply_demand.py > \
                                 ${data_dir}${sl}output${sl}g1_sensitivities.txt

# plot the best calibrated kmax fluxes across sites and years
${this_dir}${sl}plotting_scripts${sl}flux_vs_years.py ${p3}

# plot the best calibrated kmax gs - soil moisture response functions
${this_dir}${sl}plotting_scripts${sl}gs_vs_soil_moisture.py ${p3}

# plot the best calibrated climate kmax experiments
${this_dir}${sl}plotting_scripts${sl}flux_vs_years.py ${p1} -c ${p2}
${this_dir}${sl}plotting_scripts${sl}gs_vs_soil_moisture.py ${p1} -c ${p2}

# compare the Control TractLSM with CABLE 
${this_dir}${sl}plotting_scripts${sl}Tract_LSM_vs_CABLE.py ${p1}

# plot the best calibrated Control TractLSM vs the best calibrated kmax 
${this_dir}${sl}plotting_scripts${sl}calib_control.py ${p4} ${p5} ${p3}

# conceptual plots that help to understand the instantaneous optimisation 
${this_dir}${sl}plotting_scripts${sl}how_the_model_works.py

}


#######################################
# Puts the figures of the main into a
# pdf document with legends
# Globals:
#   lots
# Arguments:
#   None
# Returns:
#   figure file
#######################################

fig_file(){

touch ${fdir}${sl}Figures_main.md

cat > ${fdir}${sl}Figures_main.md << EOF

![An example of the instantaneous profit maximisation algorithm. The carbon gain
 (green), hydraulic cost (purple), and net profit (blue) are shown as functions
 of the transpiration stream, which ranges between the soil water potential at 
saturation ($\Psi$<sub>sat</sub>) and the critical water potential 
($\Psi$<sub>crit</sub>). Maximum hydraulic conductance (k<sub>max</sub>) was 
calculated for each of the three behavioural solutions (i.e. 
k<sub>max,high</sub>, k<sub>max,opt</sub>, and k<sub>max,low</sub>; cf. Section 
2 of the Materials and Methods), before being used as an input the model. The 
dashed and dotted lines illustrate the impacts of alternative strategies for 
k<sub>max</sub> on the maximum profit; the optimal leaf water potentials 
($\Psi$<sub>leaf,opt</sub>) at which profit is maximised span a range of 0.4 MPa
 between the instantaneous model run using k<sub>max,high</sub> and the one 
using k<sub>max,low</sub>. The species used in this example is <i>Juniperus 
virginiana</i> (P<sub>50</sub>=-6.6 MPa and P<sub>88</sub>=-10.5 MPa; Choat 
<i>et al.</i>, 2012), with V<sub>cmax,25</sub>=100 µmol m<sup>-2</sup> 
s<sup>-1</sup>, T<sub>air</sub>=25 °C and D=1 kPa, $\Psi$<sub>s</sub>=-0.8 kPa, 
and LAI=2 m<sup>2</sup> m<sup>-2</sup>.](${fdir}${sl}profitmax_behaviours.png)

&nbsp;

![The sensitivity of the modelled optimal coordination between the maximum 
hydraulic conductance (k<sub>max,opt</sub>) and <b>(a)</b> the maximum 
carboxylation rate at 25°C (V<sub>cmax,25</sub>) and <b>(b)</b> air temperature 
(T<sub>air</sub>), both depending on vapour pressure deficit (D). In panel a, 
k<sub>max,opt</sub> increases near proportionally with V<sub>cmax,25</sub> and 
in a logarithmic fashion with D; T<sub>air</sub> is fixed to 25 °C. In panel b, 
k<sub>max,opt</sub> increases with T<sub>air</sub>, before decreasing (sharply 
at the two highest D) starting between 20 – 25 °C; V<sub>cmax,25</sub> is set to
 100 µmol m<sup>-2</sup> s<sup>-1</sup>. The valid range for k<sub>max,opt</sub>
 is constrained by physically plausible co-occurring values of T<sub>air</sub> 
and D under a relative humidity spanning 5 – 95%. The species used in this 
example is <i>Juniperus virginiana</i> (P<sub>50</sub>=-6.6 MPa and 
P<sub>88</sub>=-10.5 MPa; Choat <i>et al.</i>, 2012), with 
$\Psi$<sub>s</sub>=-0.8 kPa, and LAI=2 m<sup>2</sup> 
m<sup>-2</sup>.](${fdir}${sl}profitmax_2_forcing.png)

&nbsp;

![Quantile ranks of the best configurations of the Profit<sub>max</sub> model 
compared to the Control model, across drought (panels a and b) and non-drought 
years (panels c and d), and for Gross Primary Productivity (GPP; panels a and c)
 and Evapotranspiration (ET; panels b and d). The vertical lines – blue (Average
 scenario) or red (Extreme scenario) – correspond to the best climate 
configuration range of ranks across the three behavioural solutions for the 
maximum hydraulic conductance. The box and whisker plots to the right of the 
figure show summaries of the ranks across sites, but do not account for the 
behavioural range shown by the vertical lines. In the box and whiskers plots, 
the horizontal yellow line shows the average overall quantile rank, and the box 
shows the interquartile range; the whiskers extend to the 10<sup>th</sup> and 
90<sup>th</sup> percentile values of the ranks, with dots outside of the 
whiskers showing outliers.](${fdir}${sl}performance_scores.png)

&nbsp;

![A 14-day running average of the carbon and water fluxes predicted by the best 
selected calibration configuration from the Profit<sub>max</sub> model (green 
line) at the five northernmost eddy-covariance sites during the 2003 (panels a, 
b, e, f, i, j, m, n, o, p) and 2006 (panels c, d, g, h, k, l) European drought 
events, compared to the Control model (purple line), and to the observations 
(black line). Grey lines show the prescribed phenologies (LAI, m<sup>2</sup> 
m<sup>-2</sup>) and blue bars the daily precipitation (PPT, mm d<sup>-1</sup>). 
The Gross Primary Productivity (GPP) is in g C m<sup>-2</sup> d<sup>-1</sup> and
 the Evapotranspiration (ET) in mm 
d<sup>-1</sup>.](${fdir}${sl}calib_north_drought.png)

&nbsp;

![A 14-day running average of the carbon and water fluxes predicted by the best 
selected calibration configuration from the Profit<sub>max</sub> model (green 
line) at the five southernmost eddy-covariance sites during the 2003 (panels a, 
b, e, f, i, j, m, n, o, p) and 2006 (panels c, d, g, h, k, l) European drought 
events, compared to the Control model (purple line), and to the observations 
(black line). Grey lines show the prescribed phenologies (LAI, m<sup>2</sup> 
m<sup>-2</sup>) and blue bars the daily precipitation (PPT, mm d<sup>-1</sup>). 
The Gross Primary Productivity (GPP) is in g C m<sup>-2</sup> d<sup>-1</sup> and
 the Evapotranspiration (ET) in mm 
d<sup>-1</sup>.](${fdir}${sl}calib_south_drought.png)

&nbsp;

![A 14-day running average of the carbon and water fluxes predicted by the best 
selected calibration configuration from the Profit<sub>max</sub> model (green 
line) at the five northernmost eddy-covariance sites during the 2003 (panels a, 
b, e, f, i, j, m, n, o, p) and 2006 (panels c, d, g, h, k, l) European drought 
events, compared to the Control model (purple line), and to the observations 
(black line). The green line is the k<sub>max,opt</sub> strategy and the green 
shadings encompass the instantaneous range of fluxes predicted by the three 
behavioural solutions for k<sub>max</sub>. Grey lines show the prescribed 
phenologies (LAI, m<sup>2</sup> m<sup>-2</sup>) and blue bars the daily 
precipitation (PPT, mm d<sup>-1</sup>). The Gross Primary Productivity (GPP) is 
in g C m<sup>-2</sup> d<sup>-1</sup> and the Evapotranspiration (ET) in mm 
d<sup>-1</sup>.](${fdir}${sl}best_north_drought.png)

&nbsp;

![A 14-day running average of the carbon and water fluxes predicted by the best 
selected calibration configuration from the Profit<sub>max</sub> model (green 
line) at the five southernmost eddy-covariance sites during the 2003 (panels a, 
b, e, f, i, j, m, n, o, p) and 2006 (panels c, d, g, h, k, l) European drought 
events, compared to the Control model (purple line), and to the observations 
(black line). The green line is the k<sub>max,opt</sub> strategy and the green 
shadings encompass the instantaneous range of fluxes predicted by the three 
behavioural solutions for k<sub>max</sub>. Grey lines show the prescribed 
phenologies (LAI, m<sup>2</sup> m<sup>-2</sup>) and blue bars the daily 
precipitation (PPT, mm d<sup>-1</sup>). The Gross Primary Productivity (GPP) is 
in g C m<sup>-2</sup> d<sup>-1</sup> and the Evapotranspiration (ET) in mm 
d<sup>-1</sup>.](${fdir}${sl}best_south_drought.png)

&nbsp;

![A comparison of the sensitivity of the Control and calibrated 
Profit<sub>max</sub> models’ stomatal conductance (g<sub>s</sub>) to vapour 
pressure deficit (D) across the 10 eddy-covariance sites. Panel <b>(a)</b> shows
 the relationship between the implied water use efficiency (g<sub>1</sub>, 
kPa<sup>0.5</sup>) of the calibrated Profit<sub>max</sub> model and its 
sensitivity to D ($\sigma$, unitless). The values of g<sub>1</sub> obtained for 
the Profit<sub>max</sub> model were converted from unit kPa<sup>$\sigma$</sup> 
to unit kPa<sup>0.5</sup> for comparison with the values of g<sub>1</sub> used 
in the Control model. g<sub>1</sub> and $\sigma$ were obtained by least-square 
fitting of the g<sub>s</sub> simulated by the calibrated Profit<sub>max</sub> 
model to the Medlyn <i>et al.</i> (2011) stomatal conductance model. The 
estimates were produced using the site hydraulic and photosynthetic parameters, 
for temperatures ranging 10 – 40 °C and D ranging 0.05 – 3 kPa. The values of 
g<sub>1</sub> used in the Control model are plotted against the respective 
sites’ $\sigma$ for visual comparison with those of the Profit<sub>max</sub> 
model only, as they correspond to $\sigma$ = 0.5. Panel <b>(b)</b> shows the 
effect of the various $\sigma$ on the g<sub>s</sub> given by the Control model 
at 25 °C. The input parameter g<sub>1</sub> was set to 2 kPa<sup>0.5</sup> for 
all the generated curves, but it was transformed to kPa<sup>$\sigma$</sup> upon 
running the Control model with the site-specific values of $\sigma$ shown in 
panel a. The reference g<sub>s</sub> of the Control model ($\sigma$ = 0.5) is 
plotted for comparison. For both panels a and b, the models were run assuming 
well-watered conditions.](${fdir}${sl}sensitivities_to_VPD.png)

&nbsp;

![A comparison of the decline in stomatal conductance (g<sub>s</sub>) with 
predawn soil water potential ($\Psi$<sub>s</sub>), for the Control (plain lines)
 and the calibrated Profit<sub>max</sub> (dotted lines) models at a 
sub-selection of sites. The functional forms emerge from the soil parameters and
 the $\beta$ functions in the Control, and from the plant hydraulic traits and 
the profit maximisation algorithm in the best selected calibration. The inset 
zooms on the functional forms of the g<sub>s</sub> - $\Psi$<sub>s</sub> curves 
from the Control model for $\Psi$<sub>s</sub> > -0.3 MPa. The functional forms
 are made comparable by normalising g<sub>s</sub> by its maximum at a given 
site. Note that seemingly slow decreases in g<sub>s</sub> with 
$\Psi$<sub>s</sub> can be attributed to the non-linear relationship between 
$\Psi$<sub>s</sub> and volumetric soil moisture, whereby small variations in the
 latter can lead to large variations in the former. To avoid rainfall effects, 
the data up to 48 hours after rain were excluded. To avoid low solar radiation 
and low temperature effects, the g<sub>s</sub> data were restricted between 
9:00 h – 15:00 h from April – November across all years. The curves were fitted 
with a linear generalised additive model and the shadings show the 95% 
confidence interval of the fit.](${fdir}${sl}calib_all_water_stress.png)

&nbsp;

![The estimated site maximum hydraulic conductance (k<sub>max</sub>), for each 
climate configuration of the Profit<sub>max</sub> model and for the best 
selected calibration, shown as a function of Mean Annual Precipitation (MAP; as 
listed in Table 2). Note, the MAP was not used in the estimation of 
k<sub>max</sub>, however k<sub>max</sub> was multiplied by the sites’ weighted 
composite LAI, which normalises it to ground area and makes it comparable across
 sites. Linear regressions are used to show the positive relationship between 
k<sub>max</sub> and MAP, with a <i>r<sup>2</sup></i> of 0.53 and a 
<i>p-value</i> of 0.02 for the best selected calibration, a 
<i>r<sup>2</sup></i> of 0.21 and a <i>p-value</i> of 0.01 for the Average 
climate scenario, and a <i>r<sup>2</sup></i> of 0.30 and a <i>p-value</i> of 
0.002 for the Extreme climate scenario.](${fdir}${sl}kmax_per_area_2_MAP.png)

EOF

pandoc ${fdir}${sl}Figures_main.md -V geometry:margin=0.8in -s -o \
       ${fdir}${sl}Figures_main.pdf
wait # make sure the job is finished before moving on

rm ${fdir}${sl}Figures_main.md

}


#### execute functions #########################################################

# standard sanity checks and file locations
this_dir="$(cd "${BASH_SOURCE%/*}"; pwd)"

# deal with linux & windows compatibility
if [[ ${this_dir} == *\/* ]]; then
    sl="/"
elif [[ ${this_dir} == *\\* ]]; then
    sl="\\"
else
    echo "Unknown dir path format, now exiting."
    exit 1
fi

# make sure this is located above the TractLSM dir...
if [[ ${this_dir} == *"TractLSM"* ]]; then
    this_dir=$(echo ${this_dir} |awk -F ${sl}"TractLSM" '{print $1}' |xargs)
fi

if [[ ${this_dir} == *"src"* ]]; then
    data_dir=$(echo ${this_dir} |awk -F ${sl}"src" '{print $1}' |xargs)

else
    data_dir=${this_dir}
fi

fdir=${data_dir}${sl}output${sl}figures${sl}final_4_paper

if [[ ! -d ${fdir} ]]; then
	mkdir ${fdir}
fi

if [[ ! -d ${data_dir}${sl}output${sl}figures${sl}not_shown_in_paper ]]; then
	mkdir ${data_dir}${sl}output${sl}figures${sl}not_shown_in_paper
fi

# execute main
main
wait # make sure the job is finished before moving on

# some figures cannot be directly saved in eps due to latex
for i in ${fdir}${sl}*".pdf"; do

    OFNAME=$(echo ${i} | cut -d'.' -f1)
    gs -q -dNOCACHE -dNOPAUSE -dBATCH -dSAFER -sDEVICE=eps2write -r1800 \
       -sOutputFile=$OFNAME.eps ${i}
    rm ${i}

done

wait # make sure the job is finished before moving on

fig_file
