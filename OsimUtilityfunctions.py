import os
os.add_dll_directory('C:/Users/jonstingel/opensim-core-4.5.1-2024-08-23-cf3ef35/bin')
import opensim as osim
import matplotlib.pyplot as plt

def probeActivate(model):
    # get the probeset
    probeset = model.getProbeSet()
    numProbes = probeset.getSize()
    # need to loop through and set them all to be enabled hopefully
    for p in range(numProbes):
        probe = probeset.get(p)
        probe.setEnabled(True)
    # update the model to be returned. 
    model.updProbeSet()
    return model

# create an ID helper function that runs the analysis and plots the moments. 
def IDplotter(solution, tag):
    # get the column labels
    labels = solution.getColumnLabels()
    # get the time
    time = solution.getIndependentColumn()
    # get the number of columns
    numCols = solution.getNumColumns()
    # create a figure plotting all of the moments. 
    fig, axs = plt.subplots(1, 3, figsize=(10, 6))
    for i in range(numCols):
        temp = labels[i]
        if 'knee_angle' in temp and 'beta' not in temp:
            if '_l_' in temp:
                axs[0].plot(time, solution.getDependentColumn(temp).to_numpy(), label=temp, linestyle=':', color='red')
            else:
                axs[0].plot(time, solution.getDependentColumn(temp).to_numpy(), label=temp, linestyle='--', color='red')
        elif 'ankle_angle' in temp:
            if '_l_' in temp:
                axs[0].plot(time, solution.getDependentColumn(temp).to_numpy(), label=temp, linestyle=':', color='blue')
            else:
                axs[0].plot(time, solution.getDependentColumn(temp).to_numpy(), label=temp, linestyle='--', color='blue')
        else:
            axs[1].plot(time, solution.getDependentColumn(temp).to_numpy(), label=temp)
    
    axs[0].set_title('Knee/Ankle')
    axs[0].legend()
    axs[1].set_title('Other coordinates')
    # axs[1].legend()
    axs[2].axis('off')
    handles, heads = axs[1].get_legend_handles_labels()
    axs[2].legend(handles, heads, loc='center')

    plt.show()

    return
