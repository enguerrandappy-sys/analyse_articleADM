# analyse_articleADM
Ce travail a été réalisé à partir d’une étude française s’intitulant: Microbial Ecology of French Dry Fermented Sausages and Mycotoxin Risk Evaluation During Storage réalisé par une équipe de 11 chercheurs.
Le saucisson sec français est un écosystème complexe où interagissent des bactéries lactiques, des staphylocoques à coagulase négative, des levures et des moisissures. Si cette microflore est essentielle au développement des qualités organoleptiques (goût, texture, odeur), elle joue aussi un rôle crucial de barrière protectrice.

Le problème central de l'article est l'évaluation de l'équilibre de l'écosystème microbien face au risque de contamination fongique toxique.
Plus précisément, l'étude cherche à comprendre :
- Identification des communautées fongiques et bactéries.
- Comparaison des profils microbiens avec/sans inoculation. 
L'utilisation de DADA2 permettra de passer des séquences brutes (FastQ) à des ASV (Amplicon Sequence Variants), offrant une résolution taxonomique plus fine que les anciens OTU pour distinguer des espèces proches de Penicillium ou de Staphylococcus.
Une fois les séquences traitées par l’algorithme DADA2 et les tables d’ASV (Amplicon Sequence Variants) générées, l’intégration des données sous l’environnement Phyloseq constitue l'étape charnière de l’analyse. Cet outil permet de croiser les données taxonomiques avec les métadonnées de l'étude (type de saucisson, durée de stockage, présence de mycotoxines) pour répondre à trois objectifs fondamentaux :
## 1. Caractérisation de la Diversité Alpha
L’analyse de la diversité alpha permet de mesurer la richesse et l’équitabilité des espèces au sein de chaque échantillon. Dans le contexte du saucisson sec, il s’agit de déterminer si les processus de fermentation ou les conditions de stockage influencent le nombre d’espèces présentes. En utilisant des indices tels que Shannon ou Simpson, nous pouvons vérifier si une perte de diversité microbienne est corrélée à une vulnérabilité accrue face à la colonisation par des moisissures indésirables.
## 2. Comparaison de la Diversité Bêta
La diversité bêta est essentielle pour comparer la structure des communautés entre différents groupes d'échantillons. Par le biais de méthodes d'ordination comme le PCoA (Principal Coordinates Analysis) ou le NMDS basés sur les distances de Bray-Curtis ou Unifrac, Phyloseq permet de visualiser si les échantillons se regroupent en fonction de leur profil de sécurité sanitaire. Cette étape permet de distinguer clairement la signature microbienne des saucissons "sains" de ceux présentant un risque de contamination par les mycotoxines, mettant ainsi en évidence l'impact des variations de la flore sur la stabilité du produit.
## 3. Identification de Biomarqueurs et Interactions
Enfin, l’outil facilite l’identification de taxons spécifiques associés à des conditions particulières. En corrélant l’abondance relative de certaines bactéries (ex: Lactobacillus, Staphylococcus) avec l'absence de champignons toxinogènes (ex: Penicillium nordicum), il devient possible d'identifier des biomarqueurs de sécurité. Cette approche permet de mettre en lumière des phénomènes d'exclusion compétitive ou de synergie, offrant des pistes pour l'utilisation de cultures protectrices capables de limiter naturellement le risque de mycotoxines durant le stockage.

## Mise en place de dada2
```{r, echo=FALSE}
library(dada2)
```
