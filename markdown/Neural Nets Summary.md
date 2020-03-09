## Neural Nets Summary

*Data version: 10/04/2019 numeric only 

*307 compounds after cleaning

#### Goal 

Predict 11 chemistry categories by first 16 normalized PC's (more than 95% variance captured) using neural networks. 

#### Steps

- Table below provides the number of compounds of 11 chemistry categories. More than 90% compounds are coming from the first 5 categories, and only one compound comes from each of the last two groups. Thus instead of classifying all the compounds into 11 categories, we classify them into the first $k$ categories. For now, let $k=5$.

| Chem Gategory |  5   |  3   |  16  |  6   |  2   |  7   |  14  |  4   |  11  |  8   |  10  |
| :-----------: | :--: | :--: | :--: | :--: | :--: | :--: | :--: | :--: | :--: | :--: | :--: |
|     Count     | 109  | 102  |  36  |  27  |  12  |  7   |  5   |  4   |  3   |  1   |  1   |

- Use neural nets with four hidden layers, which have the number of nodes 20, 10, 10, 5 respectively. 

#### Results

- Below is a schematic of the neural nets. Correspond weights between each pair of nodes have been saved as *nn weights.csv* in Dropbox.

<img src="/Users/hanyin/Dropbox (CSU Fullerton)/Research Project/Solid State Chemistry/ChemStatsResearch/Chem Stats R code/result/screenshot/nn schematic.png" alt="nn schematic" style="zoom:21.5%;" />

- Histograms below plot the distributions of 200 (cross-validation) classification rates based on different $k$ between 3 and 11. The best $k$  renders classification rate varies from 0.55 to 0.60. 

<img src="/Users/hanyin/Dropbox (CSU Fullerton)/Research Project/Solid State Chemistry/ChemStatsResearch/Chem Stats R code/result/screenshot/dists of classification rate.png" alt="dists of classification rate" style="zoom:35%;" />

