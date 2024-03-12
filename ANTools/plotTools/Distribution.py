import json
import matplotlib.pyplot as plt
import numpy as np

categories = [22, 1522, 11122, 1122, -99]

class artist():
    def __init__(self, pltcfg_dir, save_dir) -> None:
        with open(pltcfg_dir+"plt_cfg.json", "r") as plot_config:
            plt_conf = json.load(plot_config)
        self.config = plt_conf
        self.dir_to_picture = save_dir
        self.colors = ["#66b3ff", "#368c4e", "#bc504e", "#bc711e", "#000000", ]
        
    def plot_number_of_each_categories(self, df):
        def autolabel(rects, vals):
            """Attach a text label above each bar in *rects*, displaying its height."""
            for rect, val in zip(rects, vals):
                height = rect.get_height()
                ax.annotate('{:.1f}%'.format(val*100.),
                            xy=(rect.get_x() + rect.get_width() / 2, height),
                            xytext=(0, 3),  # 3 points vertical offset
                            textcoords="offset points",
                            ha='center', va='bottom')

        categories_label = self.config["Categories"]
        cates, labels = [], []
        counts_pass, counts_rej = [], []
        ratio_pass, ratio_rej = [], []
        fig, ax = plt.subplots(1, 1, figsize=(8, 8))
        for cate, label in categories_label.items():
            num_rej = len(df[(df.genpho_Category == int(cate)) & (df.genpho_IsoFrixione == False)])
            num_pass = len(df[(df.genpho_Category == int(cate)) & (df.genpho_IsoFrixione == True)])
            cates.append(cate)
            counts_rej.append(num_rej)
            counts_pass.append(num_pass)
            ratio_pass.append(num_pass/(num_rej+num_pass))
            ratio_rej.append(num_rej/(num_rej+num_pass))
            labels.append(label)
            print(f"{label}: {num_pass}/{num_rej} and {num_pass + num_rej}")

        n = len(cates)
        r = np.arange(n)
        width = 0.4
        rects_pass = ax.bar(r-width/2, counts_pass, width=width, edgecolor="black", color = "#66B3FF", label="pass")
        rects_rej = ax.bar(r+width/2, counts_rej, width=width, edgecolor="black", color="#FFC78E", label="reject")
        autolabel(rects_pass, ratio_pass)
        autolabel(rects_rej, ratio_rej)
        ax.set_xticks(r, labels, fontsize=14)
        ax.set_ylabel("number of photons", fontsize=20)
        ax.tick_params(axis = "x", direction="in", top=True, bottom = True, length=5, width=1.5)
        ax.tick_params(axis = "y", direction="in", left=True, right=True, length=5, width=1.5)
        ax.legend(title="Frixione Isolation:", fontsize=14, title_fontsize=14)
        fig.text(0.12, 0.9, "CMS", fontsize=24, fontweight="bold")
        fig.text(0.23, 0.9, "Work-in-progress", fontsize=20)    
        fig.savefig(self.dir_to_picture+"counts_category.png", dpi=800)
        plt.clf()
        
        
    def compare_kinematics(self, df, region):
        # this aim to compare the fundamental kinematics among reco and gen level
        headers = ["genpho_", "Photon_"]
        kinematics = ["pt", "eta", "phi"]
        labels = ["$p_{T}$", "\u03B7", "\u03D5"]
        con_region = df.Photon_isScEtaEB == True if region == "EB" else df.Photon_isScEtaEE == True

        for kin, xlabel in zip(kinematics, labels):            
            fig, ax = plt.subplots(1, 1, figsize=(8, 8))
            xlowlim, xuplim, setlog = self.config["axis_database"][kin]
            for i, hat in enumerate(headers):
                df.loc[con_region, f"{hat}{kin}"].hist(
                    histtype = "step",
                    bins=40, range=(xlowlim, xuplim),
                    ax=ax, ls="-", color=self.colors[i],
                    linewidth = 2,
                    grid = False,
                    label = "$\u0263^{gen}$" if hat == "genpho_" else "$\u0263^{reco}$"
                )
            ax.set_xlabel(xlabel, fontsize=20)
            ax.set_ylabel("number of photons", fontsize=20)
            ax.tick_params(axis = "x", direction="in", top=True, bottom = True, length=5, width=1.5)
            ax.tick_params(axis = "y", direction="in", left=True, right=True, length=5, width=1.5)
            ax.legend(fontsize=14,)
            fig.text(0.12, 0.9, "CMS", fontsize=24, fontweight="bold")
            fig.text(0.23, 0.9, "Work-in-progress", fontsize=20)    
            fig.savefig(self.dir_to_picture+f"kin_genreco/compare_{kin}_{region}.png", dpi=800)
            plt.close("all")

    def plot_features_distribution(self, df, region):
        # this function aim to show the distribtion among all categories
        categories_label = self.config["Categories"]
        features = self.config["features_to_plot"]["common"]
        features.update(self.config["features_to_plot"][region])
        axis_choices = self.config["axis_database"]
        categories_label = self.config["Categories"]
        
        con_region = df.Photon_isScEtaEB == True if region == "EB" else df.Photon_isScEtaEE == True
        
        for ft, axis in features.items():
            axis, xlabel = axis
            xlowlim, xuplim, setlog = axis_choices[axis]
            fig, ax = plt.subplots(1, 1, figsize=(8, 8))
            for i, item in enumerate(categories_label.items()):
                cate, label = item
                con_cate = df.genpho_Category == int(cate)
                df.loc[(con_region) & (con_cate), ft].hist(
                    histtype = "step",
                    bins=40, range=(xlowlim, xuplim),
                    ax=ax, ls="-", color=self.colors[i],
                    linewidth = 2,
                    grid = False, density = True,
                    label = label
                )
            if axis[:3] == "iso" or axis[:12] == "shower_shape":
                ax.set_xlabel(ft.split("_")[-1]+xlabel, fontsize=20)
            else:
                ax.set_xlabel(xlabel, fontsize=20)
            ax.set_ylabel("a.u.", fontsize=20)
            ax.tick_params(axis = "x", direction="in", top=True, bottom = True, length=5, width=1.5)
            ax.tick_params(axis = "y", direction="in", left=True, right=True, length=5, width=1.5)
            ax.legend(fontsize=14,)
            if setlog:
                ax.set_yscale('log')
            fig.text(0.84, 0.9, region, fontsize=24)
            fig.text(0.12, 0.9, "CMS", fontsize=24, fontweight="bold")
            fig.text(0.23, 0.9, "Work-in-progress", fontsize=20)    
            fig.savefig(self.dir_to_picture+f"features/{region}/{ft}.png", dpi=800)
            print(self.dir_to_picture+f"features/{region}/{ft}.png has been saved!")
            plt.close("all")
    def plot_correlations(self, df):
        import seaborn as sns
        features_conf = self.config["features_corr"]
            
        for region in ["EB", "EE"]:
            con_region = df.Photon_isScEtaEB == True if region == "EB" else df.Photon_isScEtaEE == True
            features = features_conf["common"] + features_conf[region]
            fig, axes = plt.subplots(1, 1, figsize=(len(features)/1.4, len(features)/1.4))
            cor = df.loc[con_region, features].corr()
            new_name = {name: name.split("_")[1] for name in cor.index}
            cor = cor.rename(columns=new_name, index=new_name)
            sns.heatmap(cor, annot=True, cmap=plt.cm.bwr,ax=axes,annot_kws={"size":len(features)/2.1, "color":"k"}, fmt=".2f",)
            axes.collections[0].set_clim(-0.6, 0.6) 
            axes.collections[0].colorbar.ax.tick_params(labelsize=20)
            axes.tick_params(axis='x', labelsize=len(features)/1.2)
            axes.tick_params(axis='y', labelsize=len(features)/1.2)
            axes.tick_params(axis = "x", bottom=False, left=False, labelrotation = 85)
            axes.tick_params(axis = "y", bottom=False, left=False, labelrotation = 10)
            plt.margins(2.)
            plt.subplots_adjust(bottom=0.22, left=0.22)
            fig.text(0.74, 0.9, region, fontsize=35)
            fig.text(0.22, 0.9, "CMS", fontsize=35, fontweight="bold")
            fig.text(0.3, 0.9, "Work-in-progress", fontsize=25)   
            fig.savefig(f"{self.dir_to_picture}{region}_features_corr.png", dpi=1100)
            print(f"{self.dir_to_picture}{region}_features_corr.png has been saved!")
            plt.close("all")