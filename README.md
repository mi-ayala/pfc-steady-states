# pfc-rigorous-steady-states-2d

> The rigorous computation of Phase-Field-Crystal steady states in periodic rectangular domains

This repository contains code used to generate the data and results presented in https://arxiv.org/abs/2102.02338
You will need a working version of INTLAB and MATLAB to use this package. The code was tested with INTLAB 9.

The various states shown in the paper can be found in `CoefficientFiles/Table_X/state_rowY.mat` as coefficient arrays `A` along with the relevant PFC settings `pfc_g`.

**Tools**

The main tool is `TestAnsatz.m`, which runs the Newton solver and the radii polynomial approach on the four basic ansatz in three different PFC regimes. This used to verify the content of Table 1.

Next, `TestCoefficientFile.m` runs the radii polynomial approach for a given state and the given PFC settings. To use this file, make the file point to one of the states contained in `CoefficientFiles`. Energy and stabiblity can also be verified by uncommenting the relevant lines. Used to verify the contents of Table 2, 3, 4, 5.

Other general tools used in the paper are in `PaperTools/`
- `Explore_SS.m`: Used to explore the energy landscape by trying random initial states. This produces candidate steady states using the Newton solver and then verifies them.
- `Connect_SS.m`: Used to explore connections.
- `Continue_m_SS.m`: Used to continue states.
- The files in `PhaseDiagram/` can be used to generate the phase diagram shown in Figure 1.

Candidate steady states can also be produced through PFC simulations, for example with the tools in https://github.com/gmartinemath/pfc-gem. This can be useful to demonstrate that states that appear to be numerically stable are indeed (close to) steady states - or on the contrary only _metastable_. 
