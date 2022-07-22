

def plotall(smarts,bincounts,func,min_E,rel_E,TEU_energies,degstatistics,xtaldeg,bestnofilter_deg,bestfilter_deg,cate,imgpath):
    # fig, axs = plt.subplots(1, 4, figsize=(60, 10), sharex=True, sharey=True)
    fig, axs = plt.subplots(1, 4, figsize=(60, 10), sharex=True)
    for num in range(4):
        if num == 0:
            # 1. plot_TEU
            ax1 = axs[num]
            ax1.set_ylim(0,10000)
            x = list(range(-175,180,10))
            ax1.bar(x,bincounts, width=10, color="b",align="center",edgecolor="black",label="CSD histogram")
            # 2. plot_QMoptE
            angles = list(range(-180,181,15))
            probe_angs = list(range(-180,181,1))
            ax2 = ax1.twinx()
            lns1 = ax2.plot(probe_angs, np.array([ func(ang) for ang in probe_angs ]) - min_E, color = "black",label="QM fragment", linewidth=5)
            ax2.plot(angles, rel_E, color="red", marker="^", markersize=10,linestyle="")
            ax2.set_ylim(0,20)

            # 3. plot CSD-TEU
            lns2 = ax2.plot(probe_angs, TEU_energies, color="blue", label="CSD TEU",linewidth=5)
            lns = lns1+lns2
            labs = [l.get_label() for l in lns]

            # plt.close()
            plt.xlim(-185,185)
            x_major_locator = MultipleLocator(30) # The sep is set to 30
            ax = plt.gca() # the instance of 2 axes
            ax.xaxis.set_major_locator(x_major_locator)
            # plt.figure()
            font1 = {"family": "Helvetica", "weight": "normal", "size": 30}
            font2 = {"family": "Helvetica", "weight": "normal", "size": 23}
            ax1.set_title(smarts, font1)
            ax1.set_xlabel("Angle(degree)",font2)
            ax1.set_ylabel("CSD histogram count",font2)
            ax2.set_ylabel("QM relE(kcal/mol) or CSD TEU",font2)
            # ax2.plot(probe_angs, TEU_energies, color="black", marker="o", markersize=5,linestyle="")
            ax1.tick_params(labelsize=23)
            ax2.tick_params(labelsize=23)
            ax2.legend(lns, labs, loc="upper right",prop={"size":23})
            
        else:
            ax1 = axs[num]
            degs = degstatistics[num-1]
            ax1.hist(degs, bins=int(360/10), range=[-180,180],edgecolor="black")
            ax1.set_ylim(0,1000)
            ax2 = ax1.twinx()
            lns1 = ax2.plot(probe_angs, np.array([ func(ang) for ang in probe_angs ]) - min_E, color = "black",label="QM fragment", linewidth=5)
            ax2.plot(angles, rel_E, color="red", marker="^", markersize=10,linestyle="")
            ax2.set_ylim(0,20)
            lns2 = ax2.plot(probe_angs, TEU_energies, color="blue", label="CSD TEU",linewidth=5)
            ax1.tick_params(labelsize=23)
            ax2.tick_params(labelsize=23)
        plt.axvline(x=xtaldeg,color="red",linestyle="--")
        plt.axhline(y=1.7,color="blue", linestyle="--")
    axs[2].axvline(x=bestnofilter_deg,color="green",linestyle="--")
    axs[3].axvline(x=bestfilter_deg,color="green",linestyle="--")
    axs[1].set_ylim(0,500)
    axs[1].set_title(cate,font1)
    axs[2].set_title("nofilter",font1)
    axs[3].set_title("filter",font1)
    
    plt.savefig(imgpath / f"{cate}.png",facecolor = "white", transparent = False)
    plt.show()
    plt.close()
