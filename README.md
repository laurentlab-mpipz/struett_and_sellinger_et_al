# struett_and_sellinger_et_al


This git provides access to all original workflows, scripts and data that are used in the publication; including creating the figures. To use tsABC, we recommend to use the provided workflow in https://github.com/sstruett/tsABC.git, to use teSMC please check https://github.com/TPPSellinger/eSMC2



Folder structure
==

+ figures_main
Scripots to reshape the analysed data for creating the main figures.

+ figures_suppl
Scripots to reshape the analysed data for creating the supplemental figures.

+ scripts
Workflows for the simulations and analyses from the publication

  + at_treeseq
  Workflow to generate the tree sequence from the already published A. thaliana data (Figure 5) which serves as input file for estimating the transition to selfing of A. thaliana.

  + bgs
  Workflow to analyse the robustness of tsABC to negative linked selection (Figure 4B)

  + simulations_for_teSMC
  Workflow to simulate analyse via teSMC (Figure 2)

  + tsABC_without_masking
  Workflow for the tsABC performance analyses (Figure 3)

  + TL_distr_TMwin_TMtrue
  Workflow to run simulations and visualise summarizing statistics (Figure 1)

  + tsabc_at_estimates
  Workflow to estimate parameters on A. thaliana using tsABC (Figure 5)

  + tsabc_with_masking
  Workflow to estimate parameters under negative linked selection if exonic regions are masked.

