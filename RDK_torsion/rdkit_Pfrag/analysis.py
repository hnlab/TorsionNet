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
            '1': 'rgba(255,0,0,0.4)',
            # '0': 'rgba(255,0,0.4,0)'
            '0': 'rgba(255,0.4.0,0)'

        },
    )
    fig2 = px.line(
        dataframe, x=xlabel, y=xlabel,
        width=500, height=500,
    )
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
# %%

# if __name__ == "__main__":
if True:
    csvfile = "/pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/outputs/ShuoGu_AmpC_44/Shuo_Gu_AmpC.csv.sort"
    df = pd.read_csv(csvfile)

    df["Active"] = list([ str(i) for i in list(df["Active"]) ])
    pxfigure(df, xlabel="Total strain(CSD)", ylabel="Total strain(QM)",title="Total Strain (kcal/mol)",color="Active", hover_data=['ZINC ID', 'Docking score'], filename="Total_strain")
    pxfigure(df, xlabel="Maximum single strain(CSD)", ylabel="Maximum single strain(QM)",title="Maximum Single Strain (kcal/mol)",color="Active", hover_data=['ZINC ID', 'Docking score'], filename="Maximum_strain")

    # calculate hit rate

    # 1) Filtering
    ## 1. ShuoGu TEU6
    df_TEU6 = df[df["Total strain(CSD)"] < 6]
    print(len(df_TEU6[df_TEU6["Active"] == '1']) / len(df_TEU6[df_TEU6["Active"] == '0']))
    ## 2. qcxia QM TEU18.5
    df_TEU185 = df[df["Total strain(QM)"] < 18.5]
    print(len(df_TEU185[df_TEU185["Active"] == '1']) / len(df_TEU185[df_TEU185["Active"] == '0']))

    # 2) Ranking
    ## 3. Shuo Gu Docking Energy

    ## 4. qcxia Docking Energy + Total Strain (kcal/mol) QM
    DockingE = np.array(df["Docking score"])
    TotalStrain_QM = np.array(df["Total strain(QM)"])
    DockingE_p_TotalStrainQM = DockingE + TotalStrain_QM
    df["Docking score + Total Strain"] = DockingE_p_TotalStrainQM
    pxfigure(df, xlabel="Docking score", ylabel="Docking score + Total Strain",title="Docking Score w/ Total Strain",color="Active", hover_data=['ZINC ID', 'Docking score', "Total strain(QM)", "Total strain(CSD)"], filename="Docking_Score")

    # 'Docking Energy + Total Strain (kcal/mol)'


# %%
