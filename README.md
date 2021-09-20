# pfc-rigorous-steady-states-2d

This repository contains code used to generate shown in https://arxiv.org/abs/2102.02338
You will need a working version of INTLAB (tested with INTLAB 9)

The various states shown in the paper can be found in CoefficientFiles/Table_X/state_rowY.mat as coefficient arrays along with the relevant PFC settings.

TestAnsatz.m: Runs the Newton solver and the radii polynomial approach on the four basic ansatz in three different PFC regimes. Used to verify Table 1.

ExtractCoefficientFromAccumulate.m: The coefficients are saved in common files, so this file splits one of the selected states for easier labeling.
TestCoefficientFile.m: Runs the radii polynomial approach for a given state and the given PFC settings. Used to verify the contents of Table 2, 3, 4, 5. Can also extract energy and stability results if uncommented.

Other tools used in the paper are in PaperTools/
Explore_SS.m: Used to explore the energy landscape by trying random initial states.
Connect_SS.m: Used to explore connections.
Continue_m_SS.m: Used to continue states.

The files in PhaseDiagram/ can be used to generate the phase diagram shown in Figure 1.
