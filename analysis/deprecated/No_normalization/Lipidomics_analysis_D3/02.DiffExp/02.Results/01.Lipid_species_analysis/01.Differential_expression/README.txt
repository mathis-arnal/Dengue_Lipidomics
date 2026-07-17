In lipid species analysis section, differentially expressed analysis is performed to find significant lipid species. In short, samples will be divided into two/multiple groups (independent) based on the Group Information of input. Two statistical methods are provided for different data types: the t-test and the Wilcoxon test (Wilcoxon rank-sum test)for two-group data, along with one-way ANOVA and the Kruskal-Wallis test for multi-group data. Additionally, the p-value will be adjusted using the Benjamini-Hochberg procedure. The condition and cut-offs for significant lipid species are also users selected.

The lollipop chart presents lipid species that meet the predefined cut-off criteria.
The x-axis indicates the log2 fold change for two-group data and the -log10 p-value for multi-group data. Lipid species are listed along the y-axis. The color of each point on the chart corresponds to the -log10 adjusted p-value or p-value.

For multiple-group data, the scatter plot of significant expressed lipid species in each class. The dot color is corresponding to the lipid class.

One-way ANOVA:
Post-hoc method: Tukey's HSD
Multiple Testing correction: Benjamini & Hochberg
p-value: 0.05

Kruskal-Wallis:
Post-hoc method: Dunn's test
Multiple Testing correction: Benjamini & Hochberg
p-value: 0.05