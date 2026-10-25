# run_all.py
#
# Re-run the attenuation-length fits and make the appendix figures, in order.
# Usage: python3 run_all.py
#
# STS + Claude, Oct 2026

import os
import subprocess
import sys

# all paths are relative to the folder containing this script
os.chdir(os.path.dirname(os.path.abspath(__file__)))

steps = [
    ['omega_extract_time.py'],          # OMEGA Us(t) and T(t) from the original records
    ['shot50_visible_model.py'],        # fused silica shot 50, isobaric release (appendix)
    ['shot50_visible_model.py', 'pyko'],  # fused silica shot 50, pressure behind the front from pyKO
    ['omega_quartz_attenuation.py'],    # alpha-quartz OMEGA shots, isobaric release (appendix)
]
for step in steps:
    print('\n==== python3 ' + ' '.join(step))
    subprocess.run([sys.executable] + step, check=True)
