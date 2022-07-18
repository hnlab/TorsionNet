# %%
from pathlib import Path
from scipy import stats
import numpy as np
from matplotlib import pyplot as plt
import plotly
import plotly.express as px
import plotly.graph_objects as go
import os
import pandas as pd

def pxfigure(dataframe, xlabel, ylabel, title, color, hover_data, filename): # color to classify, hover_data to show extra label
    fig1 = px.scatter(
        dataframe, x=xlabel, y=ylabel, color=color,
        width=500, height=500,
        hover_data=hover_data,
        color_discrete_sequence=px.colors.qualitative.Alphabet,
        color_discrete_map={
            '1': 'red',
            # '0': 'rgba(255,0,0.4,0)'
            '0': 'blue'

        },
    )
    fig2 = px.line(
        dataframe, x=xlabel, y=xlabel,
        width=500, height=500,
    )
    # fig2['data'][0]['line']['color']="#00ff00"
    fig2.data[0].line.color = "#000000"

    fig = go.Figure(data=fig1.data + fig2.data)
    fig.update_layout(
        autosize=False,
        width=500,
        height=500,
        title=title
    )
    fig.show()

    if not os.path.exists("images"):
        os.mkdir("images")
    fig.to_image(format="png", engine="kaleido")
    fig.write_image("images/" + filename + ".png")

    return

def calc_hitr(df): # Given a df, calculate the hit ratio within
    try:
        hitr = len(df[df["Active"] == '1']) / len(df)
    except ZeroDivisionError:
        hitr = float(0)
    return hitr

def calc_allhitr(df,alldf): # Given a df, calculate the hit ratio within
    try:
        hitr = len(df[df["Active"] == '1']) / len(alldf)
    except ZeroDivisionError:
        hitr = float(0)
    return hitr

def QCX_plot(df, strain_col_name,scan_range,figname):
    # plot total active_ratio and molratio related with filter threshold
    scan_range = scan_range
    ratios = []
    nums = []
    allratios = []
    strained_ratios = []
    for strainf in scan_range:
        df_TEU = df[df[strain_col_name] < strainf]
        df_strained = df[df[strain_col_name] >= strainf]
        nums.append(len(df_TEU)/len(df))
        ratios.append(calc_hitr(df_TEU))
        strained_ratios.append(calc_hitr(df_strained))
        allratios.append(calc_allhitr(df_TEU,df))
        print(strainf)
        for df in [df_TEU, df_strained, df]:
            print(len(df[df["Active"] == '1']))
        print(len(df_TEU)/ len(df))
        print(len(df_TEU))

        
    plt.title(figname)
    plt.xlabel("filter threshold")
    plt.ylabel("active_ratio & molratio")
    plt.plot(scan_range, ratios,marker="*")
    plt.plot(scan_range, nums,marker="o")
    plt.plot(scan_range, allratios,marker="x")
    plt.plot(scan_range, strained_ratios,marker=".")
    plt.savefig(figname + ".png")
    plt.show()
    plt.close()
    return

def ShuoGu_plot(df,strain_col_name,scan_range,ylimit,figname):
    from itertools import cycle, islice
    data = []
    for strainf in scan_range:
        df_TEU_strain = df[df[strain_col_name] > strainf]
        # df_TEU_strain = df[df["Total strain(CSD)"] > strainf]
        hitr_strain = calc_hitr(df_TEU_strain)
        df_TEU_unstrain = df[df[strain_col_name] < strainf]
        hitr_unstrain = calc_hitr(df_TEU_unstrain)
        data.append((f"Total{strainf}", hitr_strain, hitr_unstrain))
        
    header = ["filterthresh","strained","unstrained"]
    tmpdf = pd.DataFrame(data, columns=header)
    my_colors = list(islice(cycle(['b', 'g', 'r', 'y', 'k']), None, len(df)))
    ax = tmpdf.plot.bar(x="filterthresh", figsize=(20,15), grid=True, title=figname, ylim=(0,ylimit), color=my_colors)
    fig = ax.get_figure()
    fig.savefig(figname +'.png')
    plt.show()
    plt.close()

    return

# %%
# if __name__ == "__main__":
# """
# for ShuoGu 44 AmpC, 5 binders + 39 nonbinders
if True:
    csvfile = "/pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/outputs/ShuoGu_AmpC_44/Shuo_Gu_AmpC.csv.sort"
    df = pd.read_csv(csvfile)

    df["Active"] = list([ str(i) for i in list(df["Active"]) ])
    pxfigure(df, xlabel="Total strain(CSD)", ylabel="Total strain(QM)",title="Total Strain (kcal/mol)",color="Active", hover_data=['ZINC ID', 'Docking score'], filename="Total_strain")
    pxfigure(df, xlabel="Maximum single strain(CSD)", ylabel="Maximum single strain(QM)",title="Maximum Single Strain (kcal/mol)",color="Active", hover_data=['ZINC ID', 'Docking score'], filename="Maximum_strain")

    # calculate hit rate

    # 1) Filtering
    ## 1. ShuoGu TEU
    scan_range = [float(i/2) for i in range(2,20,1)]
    ShuoGu_plot(df, strain_col_name="Total strain(CSD)",scan_range=scan_range,ylimit=0.2,figname="split_TotalS_CSD@ShuoGu_AmpC44")
    QCX_plot(df, strain_col_name="Total strain(CSD)",scan_range=scan_range,figname="all_TotalS_CSD@ShuoGu_AmpC44")
    scan_range = [float(i/10) for i in range(10,25,1)]
    ShuoGu_plot(df, strain_col_name="Maximum single strain(CSD)",scan_range=scan_range,ylimit=0.2,figname="split_MaxS_CSD@ShuoGu_AmpC44")
    QCX_plot(df, strain_col_name="Maximum single strain(CSD)",scan_range=scan_range,figname="all_MaxS_CSD@ShuoGu_AmpC44")

    scan_range = [float(i) for i in range(1,40,1)]
    ShuoGu_plot(df, strain_col_name="Total strain(QM)",scan_range=scan_range,ylimit=0.2,figname="split_TotalS_QM@ShuoGu_AmpC44")
    QCX_plot(df, strain_col_name="Total strain(QM)",scan_range=scan_range,figname="all_TotalS_QM@ShuoGu_AmpC44")
    scan_range = [float(i) for i in range(1,40,1)]
    ShuoGu_plot(df, strain_col_name="Maximum single strain(QM)",scan_range=scan_range,ylimit=0.2,figname="split_MaxS_QM@ShuoGu_AmpC44")
    QCX_plot(df, strain_col_name="Maximum single strain(QM)",scan_range=scan_range,figname="all_MaxS_QM@ShuoGu_AmpC44")




# %%

# for dude_AmpC_259(250 left) 59 binders + 200 decoys
if True:
    # csvfile = "/pubhome/qcxia02/git-repo/DL4molcst/scripts/dock_models/3-conf-docks/results/AmpC/studies/total/docking/outdock.csv"
    # csvfile = "/pubhome/qcxia02/git-repo/DL4molcst/scripts/dock_models/3-conf-docks/results/AmpC/studies/total/docking/outdock.csv.rmredund"
    # csvfile = "/pubhome/qcxia02/git-repo/DL4molcst/scripts/dock_models/3-conf-docks/results/AmpC/studies/total_strain/docking/outdock_strain.csv"
    csvfile = "/pubhome/qcxia02/git-repo/DL4molcst/scripts/dock_models/3-conf-docks/results/AmpC/studies/total_strain/outdock_strain.csv.rmredund"
    df = pd.read_csv(csvfile)
    
    df["Active"] = list([ str(i) for i in list(df["Active"]) ])
    
    pxfigure(df, xlabel="Docking score", ylabel="Total strain(QM)",title="Docking Score compared with Total Strain",color="Active", hover_data=['ZINC_ID', "MOL_NAME", 'Docking score', "Total strain(QM)", "Max strain(QM)"], filename="total_strain")


    scan_range = [float(i) for i in range(1,40,1)]
    ShuoGu_plot(df, strain_col_name="Total strain(QM)",scan_range=scan_range,ylimit=0.2,figname="split_TotalS_QM@dude_AmpC259")
    QCX_plot(df, strain_col_name="Total strain(QM)",scan_range=scan_range,figname="all_TotalS_QM@dude_AmpC259")
    scan_range = [float(i) for i in range(1,40,1)]
    ShuoGu_plot(df, strain_col_name="Max strain(QM)",scan_range=scan_range,ylimit=0.2,figname="split_MaxS_QM@dude_AmpC259")
    QCX_plot(df, strain_col_name="Max strain(QM)",scan_range=scan_range,figname="all_MaxS_QM@dude_AmpC259")


    DockingE = np.array(df["Docking score"])
    TotalStrain_QM = np.array(df["Total strain(QM)"])
    DockingE_p_TotalStrainQM = DockingE + TotalStrain_QM
    df["Docking score + Total Strain"] = DockingE_p_TotalStrainQM
    pxfigure(df, xlabel="Docking score", ylabel="Docking score + Total Strain",title="Docking Score w/ Total Strain",color="Active", hover_data=['ZINC_ID', "MOL_NAME", 'Docking score', "Total strain(QM)", "Max strain(QM)"], filename="Docking_Score")
    pxfigure(df, xlabel="Docking score", ylabel="Docking score with total strain",title="Total Strain in Docking Score compared with Docking Score without Torsion Strain",color="Active", hover_data=['ZINC_ID', "MOL_NAME", 'Docking score', "Total strain(QM)", "Max strain(QM)", "Total strain in docking(QM)", "Docking score with total strain"], filename="Strain_in_Docking_Score")
    pxfigure(df, xlabel="Docking score + Total Strain", ylabel="Docking score with total strain",title="Total Strain in Docking Score compared with Docking Score + Torsion Strain",color="Active", hover_data=['ZINC_ID', "MOL_NAME", 'Docking score', "Total strain(QM)", "Max strain(QM)", "Total strain in docking(QM)", "Docking score with total strain"], filename="Strain_in_Docking_Score_compare")


# """


# %%
import pandas as pd
# for ShuoGu_D4_256
if True:
    csvfile = "/pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/outputs/ShuoGu_D4_256/ShuoGu_D4_256.csv.sort"
    df = pd.read_csv(csvfile)

    df["Active"] = list([ str(i) for i in list(df["Active"]) ])
    pxfigure(df, xlabel="Total strain(CSD)", ylabel="Total strain(QM)",title="Total Strain (kcal/mol)",color="Active", hover_data=['ZINC ID', 'Docking score'], filename="Total_strain")
    pxfigure(df, xlabel="Maximum single strain(CSD)", ylabel="Maximum single strain(QM)",title="Maximum Single Strain (kcal/mol)",color="Active", hover_data=['ZINC ID', 'Docking score'], filename="Maximum_strain")

    # calculate hit rate

    # 1) Filtering
    ## 1. ShuoGu TEU
    scan_range = [float(i/2) for i in range(2,20,1)]
    ShuoGu_plot(df, strain_col_name="Total strain(CSD)",scan_range=scan_range,ylimit=0.35,figname="split_TotalS_CSD@ShuoGu_D4_256")
    QCX_plot(df, strain_col_name="Total strain(CSD)",scan_range=scan_range,figname="all_TotalS_CSD@ShuoGu_D4_256")
    scan_range = [float(i/10) for i in range(10,25,1)]
    ShuoGu_plot(df, strain_col_name="Maximum single strain(CSD)",scan_range=scan_range,ylimit=0.35,figname="split_MaxS_CSD@ShuoGu_D4_256")
    QCX_plot(df, strain_col_name="Maximum single strain(CSD)",scan_range=scan_range,figname="all_MaxS_CSD@ShuoGu_D4_256")

## NO QM DATA NOW 2022.04.20
# WTIH QM DATA NOW 2022.05.26
    scan_range = [float(i) for i in range(1,40,1)]
    ShuoGu_plot(df, strain_col_name="Total strain(QM)",scan_range=scan_range,ylimit=0.5,figname="split_TotalS_QM@ShuoGu_D4_256")
    QCX_plot(df, strain_col_name="Total strain(QM)",scan_range=scan_range,figname="all_TotalS_QM@ShuoGu_D4_256")
    scan_range = [float(i) for i in range(1,40,1)]
    ShuoGu_plot(df, strain_col_name="Maximum single strain(QM)",scan_range=scan_range,ylimit=0.5,figname="split_MaxS_QM@ShuoGu_D4_256")
    QCX_plot(df, strain_col_name="Maximum single strain(QM)",scan_range=scan_range,figname="all_MaxS_QM@ShuoGu_D4_256")