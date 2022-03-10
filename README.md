# Lepton Flavor Universality Violation at <img src="https://latex.codecogs.com/svg.latex?\Large&space;Z" title="Z" />-pole

Main reconstruction code, using CERN ROOT: `allinone.C`

Arguments: 

`type`: The input `.root` type, e.g. `'s1'` for muon mode data

`noise`: The detector location noise, default 10 micron 

'save': Boolean, save the `.root` features file or not.

`num_test`: For debug. If `num_test=0`, then run through all the events.


Example: `root -b -q 'allinone.C("s1")'`
