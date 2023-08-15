import matplotlib.pyplot as plt

def show_plot_qs(df):
    fig, ax = plt.subplots()
    ax_2nd = ax.twinx()
    ax_3rd = ax.twinx()

    ax.plot(df['time'], df['qs'], label='qs', color='crimson')
    ax.plot(df['time'], df['qs_predicted'], label='qs calculated', color='violet')
    ax_2nd.plot(df['time'], df['biomass'], label='Biomass', color='blue')
    ax_3rd.plot(df['time'], df['glucose'], label='Substrate', color='orange')
    ax_3rd.spines['right'].set_position(('outward', 60))

    ax.set_xlabel('time [h]')
    ax.set_ylabel('qs [g/(Lh)]')
    ax_2nd.set_ylabel('Biomass [g/L]')
    ax_3rd.set_ylabel('Substrate [g/L]')

    handles, labels = ax.get_legend_handles_labels()
    handles_2nd, labels_2nd = ax_2nd.get_legend_handles_labels()
    handles_3rd, labels_3rd = ax_3rd.get_legend_handles_labels()
    all_handles = handles + handles_2nd + handles_3rd
    all_labels = labels + labels_2nd + labels_3rd

    # Create a single legend using the combined handles and labels
    ax.legend(all_handles, all_labels, loc='upper center', bbox_to_anchor=(0.5, 1.15), ncols=4)
    plt.show()