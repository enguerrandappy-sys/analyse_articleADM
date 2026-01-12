# analyse_articleADM
Ce travail a été réalisé à partir d’une étude française s’intitulant: Microbial Ecology of French Dry Fermented Sausages and Mycotoxin Risk Evaluation During Storage réalisé par une équipe de 11 chercheurs.
Le saucisson sec français est un écosystème complexe où interagissent des bactéries lactiques, des staphylocoques à coagulase négative, des levures et des moisissures. Si cette microflore est essentielle au développement des qualités organoleptiques (goût, texture, odeur), elle joue aussi un rôle crucial de barrière protectrice.

Le problème central de l'article est l'évaluation de l'équilibre de l'écosystème microbien face au risque de contamination fongique toxique.

Plus précisément, l'étude cherche à comprendre :
- Identification des communautées fongiques et bactéries.
- Comparaison des profils microbiens avec/sans inoculation.

L'utilisation de DADA2 permettra de passer des séquences brutes (FastQ) à des ASV (Amplicon Sequence Variants), offrant une résolution taxonomique plus fine que les anciens OTU pour distinguer des espèces proches de Penicillium ou de Staphylococcus.

Une fois les séquences traitées par l’algorithme DADA2 et les tables d’ASV (Amplicon Sequence Variants) générées, l’intégration des données sous l’environnement Phyloseq constitue l'étape charnière de l’analyse. Cet outil permet de croiser les données taxonomiques avec les métadonnées de l'étude (type de saucisson, durée de stockage, présence de mycotoxines) pour répondre à trois objectifs fondamentaux :

### 1. Caractérisation de la Diversité Alpha

L’analyse de la diversité alpha permet de mesurer la richesse et l’équitabilité des espèces au sein de chaque échantillon. Dans le contexte du saucisson sec, il s’agit de déterminer si les processus de fermentation ou les conditions de stockage influencent le nombre d’espèces présentes. En utilisant des indices tels que Shannon ou Simpson, nous pouvons vérifier si une perte de diversité microbienne est corrélée à une vulnérabilité accrue face à la colonisation par des moisissures indésirables.

### 2. Comparaison de la Diversité Bêta

La diversité bêta est essentielle pour comparer la structure des communautés entre différents groupes d'échantillons. Par le biais de méthodes d'ordination comme le PCoA (Principal Coordinates Analysis) ou le NMDS basés sur les distances de Bray-Curtis ou Unifrac, Phyloseq permet de visualiser si les échantillons se regroupent en fonction de leur profil de sécurité sanitaire. Cette étape permet de distinguer clairement la signature microbienne des saucissons "sains" de ceux présentant un risque de contamination par les mycotoxines, mettant ainsi en évidence l'impact des variations de la flore sur la stabilité du produit.



## Mise en place de dada2

```{r, echo=FALSE}
library(dada2)
```
```{r, echo=FALSE}
csv <- read.csv2("~/dossier_1/analyse_article_ADM/data2/SraRunTable.csv", header = TRUE, sep = ",")
```

```{r}
path<- "~/dossier_1/analyse_article_ADM/data2/donnees"
list.files(path)
```
```{r}
 [1] "SRR15634651_1.fastq" "SRR15634651_2.fastq" "SRR15634652_1.fastq" "SRR15634652_2.fastq"
  [5] "SRR15634653_1.fastq" "SRR15634653_2.fastq" "SRR15634654_1.fastq" "SRR15634654_2.fastq"
  [9] "SRR15634655_1.fastq" "SRR15634655_2.fastq" "SRR15634656_1.fastq" "SRR15634656_2.fastq"
 [13] "SRR15634657_1.fastq" "SRR15634657_2.fastq" "SRR15634658_1.fastq" "SRR15634658_2.fastq"
 [17] "SRR15634659_1.fastq" "SRR15634659_2.fastq" "SRR15634660_1.fastq" "SRR15634660_2.fastq"
 [21] "SRR15634661_1.fastq" "SRR15634661_2.fastq" "SRR15634662_1.fastq" "SRR15634662_2.fastq"
 [25] "SRR15634663_1.fastq" "SRR15634663_2.fastq" "SRR15634664_1.fastq" "SRR15634664_2.fastq"
 [29] "SRR15634665_1.fastq" "SRR15634665_2.fastq" "SRR15634666_1.fastq" "SRR15634666_2.fastq"
 [33] "SRR15634667_1.fastq" "SRR15634667_2.fastq" "SRR15634668_1.fastq" "SRR15634668_2.fastq"
 [37] "SRR15634669_1.fastq" "SRR15634669_2.fastq" "SRR15634670_1.fastq" "SRR15634670_2.fastq"
 [41] "SRR15634671_1.fastq" "SRR15634671_2.fastq" "SRR15634672_1.fastq" "SRR15634672_2.fastq"
 [45] "SRR15634673_1.fastq" "SRR15634673_2.fastq" "SRR15634674_1.fastq" "SRR15634674_2.fastq"
 [49] "SRR15634675_1.fastq" "SRR15634675_2.fastq" "SRR15634676_1.fastq" "SRR15634676_2.fastq"
 [53] "SRR15634677_1.fastq" "SRR15634677_2.fastq" "SRR15634678_1.fastq" "SRR15634678_2.fastq"
 [57] "SRR15634679_1.fastq" "SRR15634679_2.fastq" "SRR15634680_1.fastq" "SRR15634680_2.fastq"
 [61] "SRR15634681_1.fastq" "SRR15634681_2.fastq" "SRR15634682_1.fastq" "SRR15634682_2.fastq"
 [65] "SRR15634683_1.fastq" "SRR15634683_2.fastq" "SRR15634684_1.fastq" "SRR15634684_2.fastq"
 [69] "SRR15634685_1.fastq" "SRR15634685_2.fastq" "SRR15634686_1.fastq" "SRR15634686_2.fastq"
 [73] "SRR15634687_1.fastq" "SRR15634687_2.fastq" "SRR15634688_1.fastq" "SRR15634688_2.fastq"
 [77] "SRR15634689_1.fastq" "SRR15634689_2.fastq" "SRR15634690_1.fastq" "SRR15634690_2.fastq"
 [81] "SRR15635525_1.fastq" "SRR15635525_2.fastq" "SRR15635526_1.fastq" "SRR15635526_2.fastq"
 [85] "SRR15635527_1.fastq" "SRR15635527_2.fastq" "SRR15635528_1.fastq" "SRR15635528_2.fastq"
 [89] "SRR15635529_1.fastq" "SRR15635529_2.fastq" "SRR15635530_1.fastq" "SRR15635530_2.fastq"
 [93] "SRR15635531_1.fastq" "SRR15635531_2.fastq" "SRR15635532_1.fastq" "SRR15635532_2.fastq"
 [97] "SRR15635533_1.fastq" "SRR15635533_2.fastq" "SRR15635534_1.fastq" "SRR15635534_2.fastq"
[101] "SRR15635535_1.fastq" "SRR15635535_2.fastq" "SRR15635536_1.fastq" "SRR15635536_2.fastq"
[105] "SRR15635537_1.fastq" "SRR15635537_2.fastq" "SRR15635538_1.fastq" "SRR15635538_2.fastq"
[109] "SRR15635539_1.fastq" "SRR15635539_2.fastq" "SRR15635540_1.fastq" "SRR15635540_2.fastq"
[113] "SRR15635541_1.fastq" "SRR15635541_2.fastq" "SRR15635542_1.fastq" "SRR15635542_2.fastq"
[117] "SRR15635543_1.fastq" "SRR15635543_2.fastq" "SRR15635544_1.fastq" "SRR15635544_2.fastq"
[121] "SRR15635545_1.fastq" "SRR15635545_2.fastq" "SRR15635546_1.fastq" "SRR15635546_2.fastq"
[125] "SRR15635547_1.fastq" "SRR15635547_2.fastq" "SRR15635548_1.fastq" "SRR15635548_2.fastq"
[129] "SRR15635549_1.fastq" "SRR15635549_2.fastq" "SRR15635550_1.fastq" "SRR15635550_2.fastq"
[133] "SRR15635551_1.fastq" "SRR15635551_2.fastq" "SRR15635552_1.fastq" "SRR15635552_2.fastq"
[137] "SRR15635553_1.fastq" "SRR15635553_2.fastq" "SRR15635554_1.fastq" "SRR15635554_2.fastq"
[141] "SRR15635555_1.fastq" "SRR15635555_2.fastq" "SRR15635556_1.fastq" "SRR15635556_2.fastq"
[145] "SRR15635557_1.fastq" "SRR15635557_2.fastq" "SRR15635558_1.fastq" "SRR15635558_2.fastq"
[149] "SRR15635559_1.fastq" "SRR15635559_2.fastq" "SRR15635560_1.fastq" "SRR15635560_2.fastq"
[153] "SRR15635561_1.fastq" "SRR15635561_2.fastq" "SRR15635562_1.fastq" "SRR15635562_2.fastq"
[157] "SRR15635563_1.fastq" "SRR15635563_2.fastq" "SRR15635564_1.fastq" "SRR15635564_2.fastq"
```

```{r}
system("gunzip ~/dossier_1/analyse_article_ADM/data2/donnees/*.gz")
```

```{r}
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
```

```{r}
plotQualityProfile(fnFs[1:2])
```

