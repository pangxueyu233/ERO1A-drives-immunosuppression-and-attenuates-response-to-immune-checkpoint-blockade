# ***ERO1A* drives immunosuppression and attenuates response to immune checkpoint blockade**


This page recorded the codes and data used and mentioned in [*Cell Report Medicine*](https://www.sciencedirect.com/science/article/pii/S2666379123003737?via%3Dihub). And you could downloaded this paper by clicking [here](pdf/paper.pdf)

![Graphical abstract](README.assets/Graphical%20abstract.png)

Immunophenotyping of tumor microenvironment (TME) is crucial for improving immunotherapy efficacy in patients with solid tumours. However, strategies to characterize TME are largely limited with huge heterogeneity. Here, we show that endoplasmic reticular oxidoreductase-1α (ERO1A), which is ectopically overexpressed in tumours, induces an immune-suppressive TME and resistance to PD-1 blockade by promoting endoplasmic reticulum (ER) stress response. Single-cell RNA sequencing analyses confirm that ERO1A is correlated with immunosuppression and dysfunction of CD8+ T cell along anti-PD-1 treatment. Ablation of *Ero1a* in tumours promotes the infiltration of lymphocytes as well as cytotoxicity of CD8+ T cells and enhances response to anti-PD-1 treatment in mouse models. In human lung cancer, high expression level of ERO1A is negatively correlated with abundance of infiltrating effector T cells and is associated with higher recurrence risk after neoadjuvant immunotherapy. Mechanistically, disruption of ERO1A triggers lethal ER stress response in tumour cells and promotes host anti-tumour immunity via immunogenic cell death. Together, our study identifies a TME-based mechanism for immunosuppression and resistance to PD-1 blockade, suggesting an immunotherapeutic target for turning ‘cold’ tumours (immune-desert/exclusive) to ‘hot’ ones (immune-inflamed).

# **Codes of analyzing and visualization**

## Codes of Single cell RNA-seq analysis

To show how we analysis scRNA-seq steps by steps, we collected detail processed information and stored in Markdown files, recorded the quality controling, batch effect reducing, dimension reducing, cells clustering, pseudo time constructing, dynamics expression genes identifying and pathways enriching.

[Chapter1](scRNA-seq_Chapter1.md) is the code of data pre-analysis. And you could access this code by click [here](scRNA-seq_Chapter1.md).

[Chapter2](scRNA-seq_Chapter2.md) is the code of figure making. And you could access this code by click [here](scRNA-seq_Chapter2.md).

# **Citation**

Our paper has been published on [*Cell Report Medicine*](https://www.sciencedirect.com/science/article/pii/S2666379123003737?via%3Dihub)

You could downloaded raw data from [GEO Database GSE224525](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE224525)