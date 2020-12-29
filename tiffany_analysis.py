# Runs OpenSim Analyze for exotendon and natural running from corresponding setup xml files.
def analyzeSolution(state, subject, freq=""):
    traj = osim.MocoTrajectory(f'Subject_{subject}/running_{state}{freq}1/{state}{freq}1_MocoInverse_solution.sto')
    # traj = osim.MocoTrajectory('Jon/muscle_stateprescribe_grfprescribe_solution.sto')
    print("Creating controls file...")
    createControlsFile(traj, state, subject, freq)
    print("Analyzing...")
    analyze = osim.AnalyzeTool(f'Subject_{subject}/running_{state}{freq}1/setup_{state}{freq}.xml')
    # analyze = osim.AnalyzeTool('Jon/setup.xml')
    analyze.run()
    print("done!")


# Create a controls file for the AnalyzeTool's setup file. 
def createControlsFile(traj, state, subject=1, freq=None):
    controls = traj.exportToControlsTable()
    # cols = controls.getColumnLabels()
    # new_cols = list()
    # for i in range(len(cols)):
    #     new_cols.append(os.path.basename(cols[i:i+1][0]))
    # controls.setColumnLabels(new_cols)
    # osim.STOFileAdapter.write(controls, f'Subject_{subject}/running_{state}{freq}1/{state}{freq}_controls.sto')
    osim.STOFileAdapter.write(controls, f'Jon/{state}{freq}_controls.sto')

# Plots side-by-side bars plot for natural vs. exotendon costs, split swing vs. stance. Some code taken from matplotlib.org.
# @param type is a str that specifies whether to calculate the average metabolic power for the right leg, left leg, or averaging both. 
#        Accepts either ('r', 'l', or 'avg both legs')


# Calculates average metabolic power for either natural or exotendon running and separates the cost into stance vs. swing.
def calcAvgMetabolicPower(state, leg, subject, freq=""):
    muscles = ['glmed1', 'glmed2', 'glmed3', 'glmin1', 'glmin2', 'glmin3',
            'semimem', 'semiten', 'bflh', 'bfsh', 'sart', 'addlong', 'addbrev',
            'addmagDist', 'addmagIsch', 'addmagMid', 'addmagProx', 'tfl', 'grac',
            'glmax1', 'glmax2', 'glmax3', 'iliacus', 'psoas', 'piri', 'recfem',
            'vasmed', 'vasint', 'vaslat', 'gasmed', 'gaslat', 'soleus', 'tibpost',
            'fdl', 'fhl', 'tibant', 'perbrev', 'perlong', 'edl', 'ehl']
    total_cost = 0
    stance_cost = 0
    swing_cost = 0
    table_metabolics = osim.TimeSeriesTable(f"Subject_{subject}/running_{state}{freq}1/AnalyzeToolResults_{state}{freq}1/running_{state}{freq}1_analyze_MetabolicsReporter_probes.sto")
    time = table_metabolics.getIndependentColumn()
    if leg == 'r':
        swing = time[0:36]
        stance = time[36:len(time) - 1]
        stance_interval = time[len(time) - 1] - time[36]
        swing_interval = time[36] - time[0]
    else:
        stance = time[0:32]
        swing = time[32:len(time) - 1]
        stance_interval = time[32] - time[0]
        swing_interval = time[len(time) - 1] - time[32]
    interval = time[len(time) - 1]
    for muscle in muscles:
        metabolics_total = table_metabolics.getDependentColumn(
            f"metabolics_total_{muscle}_{leg}").to_numpy()
        if leg == 'r':
            stance_metabolics = metabolics_total[36:len(time)-1]
            swing_metabolics = metabolics_total[0:36]
        else:
            stance_metabolics = metabolics_total[0:32]
            swing_metabolics = metabolics_total[32:len(time) - 1]
        stance_cost += np.trapz(stance_metabolics, stance)
        swing_cost += np.trapz(swing_metabolics, swing)
        total_cost += np.trapz(metabolics_total, time)
    return stance_cost / stance_interval, swing_cost / swing_interval, total_cost / interval


def plotCostBars(type, subject):
    sum_natural_stance = sum_natural_swing = sum_natural_total = sum_exotendon_stance = sum_exotendon_swing = sum_exotendon_total = 0
    if type == 'avg both legs':
        legs = ['r', 'l']
        for leg in legs:
            natural_stance, natural_swing, natural_total = calcAvgMetabolicPower("natural", leg, subject)
            exotendon_stance, exotendon_swing, exotendon_total = calcAvgMetabolicPower("exotendon", leg, subject)
            sum_natural_stance += natural_stance / 2
            sum_natural_swing += natural_swing / 2
            sum_natural_total += natural_total / 2
            sum_exotendon_stance += exotendon_stance / 2
            sum_exotendon_swing += exotendon_swing / 2
            sum_exotendon_total += exotendon_total / 2
    else:
        sum_natural_stance, sum_natural_swing, sum_natural_total = calcAvgMetabolicPower("natural", type, subject)
        sum_exotendon_stance, sum_exotendon_swing, sum_exotendon_total = calcAvgMetabolicPower("exotendon", type, subject)
    print("Total percent diff:", percentDiff(sum_natural_total, sum_exotendon_total))
    print("Stance percent diff:", percentDiff(sum_natural_stance, sum_exotendon_stance))
    print("Swing percent diff:", percentDiff(sum_natural_swing, sum_exotendon_swing))
    fig, ax = pl.subplots()
    labels = ['stance', 'swing']
    x = np.arange(len(labels))
    width = 0.35
    dark_red = (0.7, 0.1, 0.2, 0.8)
    light_pink = (0.8, 0.3, 0.3, 0.3)
    # sum_natural_stance = 1286
    # sum_natural_swing = 966
    # sum_exotendon_stance = 882
    # sum_exotendon_swing = 624
    natural = ax.bar(x - width/2, [sum_natural_stance, sum_natural_swing], width, color=[(0.1, 0.2, 0.7, 0.8)], linewidth=2, label='natural')
    exotendon = ax.bar(x + width/2, [sum_exotendon_stance, sum_exotendon_swing], width, color=[(0.1, 0.2, 0.8, 0.3)], linewidth=2, label='exotendon')
    ax.set_title(f'Avg. Metabolic Power: Stance vs. Swing ({type})')
    ax.set_ylabel('average metabolic power (W)')
    ax.set_xticks(x)
    ax.set_xticklabels(['stance', 'swing'])
    ax.legend()
    def autolabel(rects):
        """Attach a text label above each bar in *rects*, displaying its height."""
        for rect in rects:
            height = int(rect.get_height())
            ax.annotate('{}'.format(height),
                        xy=(rect.get_x() + rect.get_width() / 2, height),
                        xytext=(0, 3),  # 3 points vertical offset
                        textcoords="offset points",
                        ha='center', va='bottom')
    autolabel(natural)
    autolabel(exotendon)
    fig.set_size_inches(6, 4.85)
    fig.tight_layout()
    pl.show()



if __name__ == '__main__':
    import os
    import copy
    import numpy as np
    import pylab as pl
    from scipy.signal import butter, filtfilt
    import opensim as osm
    
    subject = osm.Model('post_simple_model_all_the_probes.osim')
    # analyzeSolution("exotendon", 2)
    plotCostBars('l',1)