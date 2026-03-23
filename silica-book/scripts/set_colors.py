# select colors
#print(plt.style.available)

# tableau-10
# https://www.tableau.com/blog/colors-upgrade-tableau-10-56782
tab10_colors=['#4d79a7', '#f28e2b', '#e15759', '#76b7b2', '#59a14f', '#edc948', '#b07aa1', '#ff9da7', '#9c755f', '#bab9ac']

# tableau matplotlib default
#tab:blue : #1f77b4
#tab:orange : #ff7f0e
#tab:green : #2ca02c
#tab:red : #d62728
#tab:purple : #9467bd
#tab:brown : #8c564b
#tab:pink : #e377c2
#tab:gray : #7f7f7f
#tab:olive : #bcbd22
#tab:cyan : #17becf

tab_colors= ['#1f77b4','#ff7f0e','#2ca02c','#d62728','#9467bd','#8c564b','#e377c2','#7f7f7f','#bcbd22','#17becf']


# blue, orange, red, grey, purple, brown, red-orange, light olive, dark grey, pale plue
# acc8 10

# https://github.com/mpetroff/accessible-color-cycles 8
# blue, orange, red, purple, light olive, pale blue, light blue, dark grey
acc8_colors=["#1845fb", "#ff5e02", "#c91f16", "#c849a9", "#adad7d", "#86c8dd", "#578dff", "#656364"]

# https://github.com/mpetroff/accessible-color-cycles 8 - swap red for green (3rd)
# blue, orange, green, purple, light olive, pale blue, light blue, dark grey
acc8g_colors=["#1845fb", "#ff5e02", '#2ca02c', "#c849a9", "#adad7d", "#86c8dd", "#578dff", "#656364"]

# Plotting items
fs_col = "#1845fb" # acc8 blue
cst_col = "#ff5e02" # acc8 orange
qtz_col = '#2ca02c' # tab:green
stv_col = "#adad7d" # acc8 light olive
bad_col = "#656364" # acc8 dark grey

mc_col = "#c849a9" # acc8 purple
sol_col = '#808080' # solidus color light grey alt "#94a4a2"
colsolidus = sol_col

# use for decaying shock data
custom_cycler = (cycler(color=['tab:cyan', 'tab:olive', 'tab:gray', 'tab:pink']) + 
                 cycler(lw=[1, 1, 1, 1])) # linewidth cycle
#fig, ax = plt.subplots(figsize=fsize)
#ax.set_prop_cycle(custom_cycler)

## plotting colors

# decaying gas gun shots fused silica
coldecayvishigh='deepskyblue'
coldecayvislow='gold'
coldecay48low = 'cornflowerblue'

coldecay=coldecayvishigh
coldecayt=coldecayvislow
col48 = coldecay48low
