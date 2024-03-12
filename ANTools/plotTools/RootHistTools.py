import os
import uproot3
import matplotlib.pyplot as plt
# plt.rcParams["font.family"] = "Times New Roman"
class RootPlotDrawer():
    def __init__(self) -> None:
        # self.cates = ["22", "1122", "11122", "1522", "151522", "2211122", "-99"]
        # self.cate_names = ["prompt \u0263", "FSR \u0263", "jet faked", "\u03C4 jet faked", "FSR from \u03C4", "\u0263 inside jet", "others"]
        self.cates = ["1122", "11122", "1522", "-99"]
        self.cate_names = ["FSR \u0263", "jet faked", "\u03C4 jet faked", "others"]
        self.colors = ["#66b3ff", "#bc711e", "#bc504e", "#368c4e", "#9336B4", "#000000", ]
        
    def set_plot_conf(self, path_to_conf):
        import json
        with open(path_to_conf, "r") as conf_file:
            plot_config = json.load(conf_file)
        self.conf = plot_config
    def set_plot_path(self, path_to_plots):
        self.plot_path = path_to_plots

    def set_out_path(self, path_to_save):
        self.out_path = path_to_save
    def set_corr_source(self, path_to_corr_table):
        self.corr_path = path_to_corr_table
        import numpy as np
        table_ = np.loadtxt(path_to_corr_table, dtype=str)
        table_dict = dict()
        # initial tables of correlation coefficient
        for erow in table_: 
            var1, var2, corr = erow.split("---")
            table_dict[var1+"---"+var2] = float(corr)
        self.corr_table = table_dict
        
    def draw_dupli_dpt(self):
        file = uproot3.open(self.plot_path)
        os.system(f"mkdir -p {self.out_path}/dupli")
        for region in ["EB", "EE"]:
            os.system(f"mkdir -p {self.out_path}/dupli/{region}")
            fig, ax = plt.subplots(1, 1, figsize=(8, 8))
            for cate, color, cate_name in zip(self.cates, self.colors, self.cate_names):
                histname = f"dupli_{region}_{cate}_delta_pt"
                h = file[histname]
                bins = h.edges
                vals = h.values
                vals = vals/sum(vals)
                ax.step((bins[:-1] + bins[1:])/2, vals,
                        # width=(bins[1] - bins[0]), 
                        ls="-", linewidth=3,
                        label=cate_name, color=color, where="post")
            ax.set_xlabel("$(p_{T}^{\u03B3_{1}}-p_{T}^{\u03B3_{2}})/p_{T}^{\u03B3_{1}}$", fontsize=20)
            ax.set_ylabel("a.u.", fontsize=20)
            ax.tick_params(axis = "x", direction="in", top=True, bottom = True, length=5, width=1.5, labelsize=15)
            ax.tick_params(axis = "y", direction="in", left=True, right=True, length=5, width=1.5, labelsize=15)
            ax.legend(fontsize=16, ncol = 2, frameon=False)
            fig.text(0.84, 0.9, region, fontsize=24)
            fig.text(0.12, 0.9, "CMS", fontsize=24, fontweight="bold")
            fig.text(0.23, 0.9, "Work-in-progress", fontsize=20, fontdict={"style": "italic"})
            fig.savefig(f"{self.out_path}/dupli/{region}/dpt_dist.png")
            fig.clf()
            plt.close("all")
    def draw_features(self):
        plot_config = self.conf["features"]
        # compare with all the cate hist
        file = uproot3.open(self.plot_path)
        for header in ["single", "dupli01", "dupli02"]:
            os.system(f"mkdir -p {self.out_path}/{header}")
            for region in ["EB", "EE"]:
                os.system(f"mkdir -p {self.out_path}/{header}/{region}")
                for var, pltset in plot_config.items():
                    # print(var)
                    xlabel, isylog = pltset
                    fig, ax = plt.subplots(1, 1, figsize=(8, 8))
                    for cate, color, cate_name in zip(self.cates, self.colors, self.cate_names):
                        histname = f"{header}_{region}_{cate}_{var}"
                        h = file[histname]
                        bins = h.edges
                        vals = h.values
                        vals = vals/sum(vals)
                        ax.step((bins[:-1] + bins[1:])/2, vals,
                                # width=(bins[1] - bins[0]), 
                                ls="-", linewidth=3,
                                label=cate_name, color=color, where="post")
                    ax.set_xlabel(xlabel, fontsize=20)
                    ax.set_ylabel("a.u.", fontsize=20)
                    ax.tick_params(axis = "x", direction="in", top=True, bottom = True, length=5, width=1.5, labelsize=15)
                    ax.tick_params(axis = "y", direction="in", left=True, right=True, length=5, width=1.5, labelsize=15)
                    if isylog:
                        ax.set_yscale("log")
                    if ("phi" in var or "Phi" in var )and "Width" not in var:
                        ax.set_ylim(0, 0.032)
                    # ax.set_xlim(1, 2.5) if var == "ETInCn" else ax.set_xlim(0, 2.5)
                    ax.legend(fontsize=16, ncol = 2, frameon=False)
                    fig.text(0.84, 0.9, region, fontsize=24)
                    fig.text(0.12, 0.9, "CMS", fontsize=24, fontweight="bold")
                    fig.text(0.23, 0.9, "Work-in-progress", fontsize=20, fontdict={"style": "italic"})
                    fig.savefig(f"{self.out_path}/{header}/{region}/{var}.png")
                    fig.clf()
                    plt.close("all")
    def draw_frixi_matrix(self):
        import pandas as pd
        import seaborn as sns
        file = uproot3.open(self.plot_path)
        for region in ["EB", "EE"]:
            for cate in self.cates:
                df = pd.DataFrame()
                fig, ax = plt.subplots(1, 1, figsize=(8, 8))
                r01, p01 = file[f"dupli_{region}_{cate}_IsoFrixioneMatrix"].values
                rr, rp = r01
                pr, pp = p01
                df.loc["pass", "reject"]   = pr
                df.loc["pass", "pass"]     = pp
                df.loc["reject", "pass"]   = rp
                df.loc["reject", "reject"] = rr
                sns.heatmap(df, annot=True, cmap=plt.cm.Blues,ax=ax,annot_kws={"size":20, "color":"k"}, fmt="7.0f", cbar=False)
                ax.collections[0].set_clim(0, 7*10**5) 
                # ax.collections[0].colorbar.ax.tick_params(labelsize=20)
                ax.tick_params(axis='x', labelsize=20)
                ax.tick_params(axis='y', labelsize=20)
                ax.tick_params(axis = "x", bottom=False, left=False, labelrotation = 85)
                ax.tick_params(axis = "y", bottom=False, left=False, labelrotation = 10)
                plt.margins(2.)
                plt.subplots_adjust(bottom=0.22, left=0.22)
                fig.text(0.8, 0.9, region, fontsize=24)
                fig.text(0.22, 0.9, "CMS", fontsize=24, fontweight="bold")
                fig.text(0.33, 0.9, "Work-in-progress", fontsize=20, fontdict={"style": "italic"})   
                fig.savefig(f"{self.out_path}/corr/dupli/{region}/{cate}_Frixione_matrix.png")
                plt.close("all")  
                  
    def draw_correlations_of(self, region):
        import pandas as pd
        import seaborn as sns
        # initial configuration
        conf_corr = self.conf["corr_plot"]
        conf_labels = self.conf["features"]
        
        os.system(f"mkdir -p {self.out_path}/corr")

        for header in ["single", "dupli01", "dupli02"]:
            os.system(f"mkdir -p {self.out_path}/corr/{header}")
            os.system(f"mkdir -p {self.out_path}/corr/{header}/{region}")
            for cate, cate_name in zip(self.cates, self.cate_names):
                # initial a new dataframe
                df = pd.DataFrame()
                size_fig = len(conf_corr["vars"][region])
                fig, ax = plt.subplots(1, 1, figsize=(size_fig/1.4, size_fig/1.4))
                for var1 in conf_corr["vars"][region]:
                    # label1, islog = conf_labels[var1]
                    label1 = str()
                    for s in var1.split("_")[1:]:
                        label1 += s + "_"
                    for var2 in conf_corr["vars"][region]:
                        # label2, islog = conf_labels[var2]
                        label2 = str()
                        for s in var2.split("_")[1:]:
                            label2 += s + "_"
                        corr_header = f"{header}_{region}_{cate}"
                        corr_val = self.corr_table[f"{corr_header}_{var1}---{corr_header}_{var2}"]
                        df.loc[label1, label2] = corr_val
                sns.heatmap(df, annot=True, cmap=plt.cm.bwr,ax=ax,annot_kws={"size":size_fig/2.1, "color":"k"}, fmt=".2f",)
                ax.collections[0].set_clim(-0.6, 0.6) 
                ax.collections[0].colorbar.ax.tick_params(labelsize=20)
                ax.tick_params(axis='x', labelsize=size_fig/1.2)
                ax.tick_params(axis='y', labelsize=size_fig/1.2)
                ax.tick_params(axis = "x", bottom=False, left=False, labelrotation = 85)
                ax.tick_params(axis = "y", bottom=False, left=False, labelrotation = 10)
                plt.margins(2.)
                plt.subplots_adjust(bottom=0.22, left=0.22)
                fig.text(0.74, 0.9, region, fontsize=35)
                fig.text(0.22, 0.9, "CMS", fontsize=35, fontweight="bold")
                fig.text(0.3, 0.9, "Work-in-progress", fontsize=25, fontdict={"style": "italic"})   
                fig.savefig(f"{self.out_path}/corr/{header}/{region}/{cate}_features_corr.png", dpi=1100)
                print(f"{self.out_path}/corr/{header}/{region}/{cate}_features_corr.png has been saved!")
                plt.close("all")                    
    def draw_correlation_duplicate(self, region):
        import pandas as pd
        import seaborn as sns
        conf_corr = self.conf["corr_plot"]
        conf_labels = self.conf["features"]
        header1 = "dupli01"
        header2 = "dupli02"
        os.system(f"mkdir -p {self.out_path}/corr/dupli")
        os.system(f"mkdir -p {self.out_path}/corr/dupli/{region}")
        for cate in self.cates:
            df = pd.DataFrame()
            size_fig = len(conf_corr["vars"][region])
            fig, ax = plt.subplots(1, 1, figsize=(size_fig/1.4, size_fig/1.4))
            for var1 in conf_corr["vars"][region]:
                label1 = str()
                for s in var1.split("_")[1:]:
                    label1 += s + "_"
                for var2 in conf_corr["vars"][region]:
                    label2 = str()
                    for s in var2.split("_")[1:]:
                        label2 += s + "_"
                    corr_header1 = f"{header1}_{region}_{cate}"
                    corr_header2 = f"{header2}_{region}_{cate}"
                    corr_val = self.corr_table[f"{corr_header1}_{var1}---{corr_header2}_{var2}"]
                    df.loc[label1, label2] = corr_val
            sns.heatmap(df, annot=True, cmap=plt.cm.bwr,ax=ax,annot_kws={"size":size_fig/2.1, "color":"k"}, fmt=".2f",)
            ax.collections[0].set_clim(-0.6, 0.6) 
            ax.collections[0].colorbar.ax.tick_params(labelsize=20)
            ax.tick_params(axis='x', labelsize=size_fig/1.2)
            ax.tick_params(axis='y', labelsize=size_fig/1.2)
            ax.tick_params(axis = "x", bottom=False, left=False, labelrotation = 85)
            ax.tick_params(axis = "y", bottom=False, left=False, labelrotation = 10)
            plt.margins(2.)
            plt.subplots_adjust(bottom=0.22, left=0.22)
            fig.text(0.74, 0.9, region, fontsize=35)
            fig.text(0.22, 0.9, "CMS", fontsize=35, fontweight="bold")
            fig.text(0.3, 0.9, "Work-in-progress", fontsize=25, fontdict={"style": "italic"})   
            fig.savefig(f"{self.out_path}/corr/dupli/{region}/{cate}_features_corr.png", dpi=1100)
            print(f"{self.out_path}/corr/dupli/{region}/{cate}_features_corr.png has been saved!")
            plt.close("all")                    
            
    def plot_number_of_single_categories(self,):
        import numpy as np
        file = uproot3.open(self.plot_path)
        # plot the number of single match photon with frixione isolation 
        def autolabel(rects, vals):
            """Attach a text label above each bar in *rects*, displaying its height."""
            for rect, val in zip(rects, vals):
                height = rect.get_height()
                ax.annotate('{:.1f}%'.format(val*100.),
                            xy=(rect.get_x() + rect.get_width() / 2, height),
                            xytext=(0, 3),  # 3 points vertical offset
                            fontsize=14,
                            textcoords="offset points",
                            ha='center', va='bottom')
        for header in ["single", "dupli01", "dupli02"]:
            for region in ["EB", "EE"]:
                cates, labels = [], []
                counts_pass, counts_rej = [], []
                ratio_pass, ratio_rej = [], []
                fig, ax = plt.subplots(1, 1, figsize=(8, 8))
                for cate, cate_name in zip(self.cates, self.cate_names):
                    histname = f"{header}_{region}_{cate}_gen_IsoFrixione"
                    rejev, passev = file[histname].values
                    cates.append(cate)
                    labels.append(cate_name)
                    counts_rej.append(rejev)
                    counts_pass.append(passev)
                    ratio_pass.append(passev/(rejev+passev))
                    ratio_rej.append(rejev/(rejev+passev))
                n = len(cates)
                r = np.arange(n)
                width = 0.4
                rects_pass = ax.bar(r-width/2, counts_pass, width=width, edgecolor="black", color = "#66B3FF", label="pass")
                rects_rej = ax.bar(r+width/2, counts_rej, width=width, edgecolor="black", color="#FFC78E", label="reject")
                autolabel(rects_pass, ratio_pass)
                autolabel(rects_rej, ratio_rej)
                ax.set_xticks(r, labels, fontsize=14)
                ax.set_ylabel("number of photons", fontsize=20)
                ax.tick_params(axis = "x", direction="in", top=True, bottom = True, length=5, width=1.5, labelsize=15)
                ax.tick_params(axis = "y", direction="in", left=True, right=True, length=5, width=1.5, labelsize=15)
                ax.legend(title="Frixione Isolation:", fontsize=16, title_fontsize=14, frameon=False)
                fig.text(0.12, 0.9, "CMS", fontsize=24, fontweight="bold")
                fig.text(0.23, 0.9, "Work-in-progress", fontsize=20, fontdict={"style": "italic"})    
                fig.savefig(f"{self.out_path}/{header}/{region}_photons.png", dpi=800)
                plt.close("all")    
                
        

  